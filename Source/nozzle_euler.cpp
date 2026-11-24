// Axisymmetric, cell-centred finite-volume Euler solver for the supplied air passage.
// Compile with: c++ -O3 -std=c++17 nozzle_euler.cpp -o nozzle_euler
// Run: nozzle_euler Nx Nr {subsonic|onset|shock|supersonic} output_directory [max_steps]
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using Vec = std::array<double,4>;
constexpr double gamma_air = 1.4, gm = gamma_air-1.0;
constexpr double pi = 3.14159265358979323846;
struct Face { double nx, nr, weight; };
struct Cell { double volume, area, x, r, dt; };
struct Balance { Vec l1{}; double inlet=0, outlet=0; };

// Primitive state: density, axial velocity, radial velocity, pressure.
// Reference units: L, p0, sqrt(R T0), p0/(R T0).
double radius(double x) { double c=std::cos(pi*x); return 0.05*std::sqrt(1+c*c); }
double area_mach(double m) { return std::pow((5+m*m)/6,3)/m; }
double exit_mach(bool supersonic) {
    double lo=supersonic?1.0:0.01, hi=supersonic?5.0:1.0;
    for(int i=0;i<60;++i) {
        double m=(lo+hi)/2;
        if((area_mach(m)<2)==supersonic) lo=m; else hi=m;
    }
    return (lo+hi)/2;
}
double pressure_ratio(double m) {return std::pow(1+0.2*m*m,-3.5);}
Vec conserved(const Vec& q) {
    return {q[0], q[0]*q[1], q[0]*q[2], q[3]/gm+0.5*q[0]*(q[1]*q[1]+q[2]*q[2])};
}
Vec primitive(const Vec& u) {
    double vx=u[1]/u[0], vr=u[2]/u[0];
    double p=gm*(u[3]-0.5*(u[1]*vx+u[2]*vr));
    if(!(p>0 && u[0]>0)) throw std::runtime_error("Nonphysical state");
    return {u[0],vx,vr,p};
}
Vec flux(const Vec& q, const Face& f) {
    double un=q[1]*f.nx+q[2]*f.nr, mf=q[0]*un;
    return {mf, mf*q[1]+q[3]*f.nx, mf*q[2]+q[3]*f.nr,
        un*(gamma_air*q[3]/gm+0.5*q[0]*(q[1]*q[1]+q[2]*q[2]))};
}
Vec hll(const Vec& l,const Vec& r,const Face& f) {
    double ul=l[1]*f.nx+l[2]*f.nr, ur=r[1]*f.nx+r[2]*f.nr;
    double al=std::sqrt(gamma_air*l[3]/l[0]), ar=std::sqrt(gamma_air*r[3]/r[0]);
    double sl=std::min(ul-al,ur-ar), sr=std::max(ul+al,ur+ar);
    if(sl>=0) return flux(l,f);
    if(sr<=0) return flux(r,f);
    Vec fl=flux(l,f), fr=flux(r,f), ql=conserved(l), qr=conserved(r);
    for(int k=0;k<4;++k) fl[k]=(sr*fl[k]-sl*fr[k]+sl*sr*(qr[k]-ql[k]))/(sr-sl);
    return fl;
}
Vec reflect(Vec q,const Face& f) {
    double v=q[1]*f.nx+q[2]*f.nr;
    q[1]-=2*v*f.nx; q[2]-=2*v*f.nr;
    return q;
}
// Limited primitive-variable increments in logical coordinates.
Vec slope(const Vec& l,const Vec& c,const Vec& r,bool smooth) {
    Vec ds{};
    for(int k=0;k<4;++k) {
        double a=c[k]-l[k], b=r[k]-c[k];
        if(smooth) ds[k]=0.5*(a+b);
        else if(a*b>0) ds[k]=std::copysign(std::min({1.25*std::abs(a),1.25*std::abs(b),0.5*std::abs(a+b)}),a);
    }
    return ds;
}
Vec extrapolate(Vec q,const Vec& ds,double factor) {
    for(int k=0;k<4;++k) q[k]+=factor*ds[k];
    return q;
}

class Nozzle {
    int nx,nr,n;
    double dx,back_pressure;
    bool smooth;
    std::string regime;
    std::vector<double> wall;
    std::vector<Face> axial,radial;
    std::vector<Cell> cell;
    std::vector<Vec> u,q,sx,sy,res,start;
    int id(int i,int j) const {return i*nr+j;}
    int rf(int i,int j) const {return i*(nr+1)+j;}
    Vec inlet(const Vec& inner) const {
        // Axial inflow. Outgoing acoustic invariant plus fixed upstream total T and p.
        double invariant=inner[1]-2*std::sqrt(gamma_air*inner[3]/inner[0])/gm;
        double a=gm*invariant;
        double disc=a*a-(gamma_air+1)*(gm*invariant*invariant-4*gamma_air/gm);
        double vel=(a+std::sqrt(std::max(0.0,disc)))/(gamma_air+1);
        double temp=1-gm*vel*vel/(2*gamma_air), p=std::pow(temp,gamma_air/gm);
        return {p/temp,vel,0,p};
    }
    Vec outlet(const Vec& inner) const {
        double a=std::sqrt(gamma_air*inner[3]/inner[0]);
        if(inner[1]>=a) return inner;
        double rho=inner[0]*std::pow(back_pressure/inner[3],1/gamma_air);
        double ab=std::sqrt(gamma_air*back_pressure/rho);
        return {rho,inner[1]+2*(a-ab)/gm,inner[2],back_pressure};
    }
    void add(int h,const Vec& f,double weight) {
        for(int k=0;k<4;++k) res[h][k]+=weight*f[k];
    }
    Balance residual(bool record,bool set_dt) {
        Balance bal;
        for(int h=0;h<n;++h) q[h]=primitive(u[h]);
        for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) {
            int h=id(i,j);
            // End cells use first-order normal extrapolation at open boundaries.
            sx[h]=(i>0 && i+1<nx)?slope(q[h-nr],q[h],q[h+nr],smooth):Vec{};
            const Face& bottom=radial[rf(i,j)], top=radial[rf(i,j+1)];
            Vec below=j?q[h-1]:reflect(q[h],bottom);
            Vec above=j+1<nr?q[h+1]:reflect(q[h],top);
            sy[h]=slope(below,q[h],above,smooth);
            res[h]={0,0,q[h][3]*cell[h].area,0};
            if(set_dt) {
                const double sound=std::sqrt(gamma_air*q[h][3]/q[h][0]);
                auto speed=[&](const Face& face) {
                    return face.weight*(sound+std::abs(q[h][1]*face.nx+q[h][2]*face.nr));
                };
                double sigma=speed(axial[h])+speed(axial[h+nr])+speed(top)+speed(bottom);
                cell[h].dt=0.55*2.0*cell[h].volume/sigma;
            }
        }
        for(int j=0;j<nr;++j) {
            int left=id(0,j), right=id(nx-1,j);
            Vec fi=flux(inlet(q[left]),axial[left]);
            Vec fo=flux(outlet(q[right]),axial[right+nr]);
            add(left,fi,axial[left].weight); add(right,fo,-axial[right+nr].weight);
            bal.inlet+=fi[0]*axial[left].weight; bal.outlet+=fo[0]*axial[right+nr].weight;
        }
        for(int i=1;i<nx;++i) for(int j=0;j<nr;++j) {
            int r=id(i,j),l=r-nr;
            const Face& f=axial[r];
            Vec ff=hll(extrapolate(q[l],sx[l],0.5),extrapolate(q[r],sx[r],-0.5),f);
            add(l,ff,-f.weight); add(r,ff,f.weight);
        }
        for(int i=0;i<nx;++i) {
            for(int j=1;j<nr;++j) {
                int hi=id(i,j),lo=hi-1;
                const Face& f=radial[rf(i,j)];
                Vec ff=hll(extrapolate(q[lo],sy[lo],0.5),extrapolate(q[hi],sy[hi],-0.5),f);
                add(lo,ff,-f.weight); add(hi,ff,f.weight);
            }
            int hi=id(i,nr-1);
            const Face& top=radial[rf(i,nr)];
            Vec inner=extrapolate(q[hi],sy[hi],0.5);
            add(hi,hll(inner,reflect(inner,top),top),-top.weight);
        }
        if(record) {
            for(const Vec& r:res) for(int k=0;k<4;++k) bal.l1[k]+=std::abs(r[k]);
            for(double& v:bal.l1) v/=std::abs(bal.inlet);
        }
        return bal;
    }

public:
    Nozzle(int nx_,int nr_,const std::string& name):
        nx(nx_),nr(nr_),n(nx*nr),dx(1.0/nx),regime(name),wall(nx+1),
        axial((nx+1)*nr),radial(nx*(nr+1)),cell(n),u(n),q(n),sx(n),sy(n),res(n),start(n) {
        smooth=(name=="subsonic");
        if(name=="subsonic") back_pressure=0.98;
        else if(name=="onset") back_pressure=pressure_ratio(exit_mach(false));
        else if(name=="shock") back_pressure=0.75;
        else if(name=="supersonic") back_pressure=pressure_ratio(exit_mach(true));
        else throw std::runtime_error("Unknown case");
        for(int i=0;i<=nx;++i) {
            wall[i]=radius(double(i)/nx);
            for(int j=0;j<nr;++j) axial[id(i,j)]={1,0,(j+0.5)*std::pow(wall[i]/nr,2)};
        }
        for(int i=0;i<nx;++i) {
            double r0=wall[i],r1=wall[i+1],avg=0.5*(r0+r1),dr=r1-r0;
            for(int j=0;j<=nr;++j) {
                double eta=double(j)/nr, len=std::hypot(dx,eta*dr);
                radial[rf(i,j)]={-eta*dr/len,dx/len,eta*avg*len};
            }
            for(int j=0;j<nr;++j) {
                int h=id(i,j);
                double x=(i+0.5)*dx, eta=(j+0.5)/nr;
                double vol=dx*(j+0.5)*(r0*r0+r0*r1+r1*r1)/(3.0*nr*nr);
                cell[h]={vol,dx*avg/nr,x,eta*avg,0};
                // Deliberately approximate start. The shock trial starts at 0.62 L.
                double mach=0.14+0.11*std::sin(pi*x);
                if(regime=="onset") mach=0.30+0.50*std::pow(std::sin(pi*x),2);
                if(regime=="supersonic" || regime=="shock") mach=0.31+1.80*x;
                bool downstream=regime=="shock" && x>0.62;
                if(downstream) mach=0.62-0.2*(x-0.62)/0.38;
                double t=1/(1+0.2*mach*mach), p=std::pow(t,3.5)*(downstream?0.80:1);
                double speed=mach*std::sqrt(gamma_air*t);
                u[h]=conserved({p/t,speed,0,p});
            }
        }
    }
    void run(const std::string& dir,int max_steps) {
        std::filesystem::create_directories(dir);
        std::ofstream hist(dir+"/history.csv");
        hist<<std::setprecision(12)<<"iteration,L1_mass,L1_axial_momentum,L1_radial_momentum,L1_energy,inlet_flux,outlet_flux\n";
        int step=0; Balance b;
        for(;;++step) {
            bool log=step%100==0 || step==max_steps;
            b=residual(log,true);
            if(log) {
                hist<<step;
                for(double val:b.l1) hist<<','<<val;
                hist<<','<<b.inlet<<','<<b.outlet<<'\n';
                if(step%1000==0) hist.flush();
                if(*std::max_element(b.l1.begin(),b.l1.end())<2e-5 || step==max_steps) break;
            }
            start=u;
            for(int h=0;h<n;++h) for(int k=0;k<4;++k)
                u[h][k]+=cell[h].dt/cell[h].volume*res[h][k];
            residual(false,false);
            for(int h=0;h<n;++h) for(int k=0;k<4;++k)
                u[h][k]=0.5*(start[h][k]+u[h][k]+cell[h].dt/cell[h].volume*res[h][k]);
        }
        std::ofstream field(dir+"/field.csv");
        field<<std::setprecision(12)<<"i,j,x_over_L,r_over_L,rho_over_rho0,u_over_uref,v_over_uref,p_over_p0,Mach\n";
        for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) {
            int h=id(i,j); auto z=primitive(u[h]);
            double m=std::hypot(z[1],z[2])/std::sqrt(gamma_air*z[3]/z[0]);
            field<<i<<','<<j<<','<<cell[h].x<<','<<cell[h].r;
            for(double v:z) field<<','<<v;
            field<<','<<m<<'\n';
        }
        std::cout<<regime<<' '<<nx<<'x'<<nr<<"  iterations "<<step
                 <<"  residual "<<*std::max_element(b.l1.begin(),b.l1.end())
                 <<"  signed mass mismatch "<<(b.outlet-b.inlet)/b.inlet<<'\n';
    }
};
int main(int argc,char** argv) {
    try {
        if(argc<5 || argc>6) throw std::runtime_error("Usage: nozzle_euler Nx Nr case output_dir [max_steps]");
        int nx=std::stoi(argv[1]), nr=std::stoi(argv[2]), steps=argc==6?std::stoi(argv[5]):40000;
        if(nx<8 || nr<2 || steps<1) throw std::runtime_error("Invalid grid or step count");
        Nozzle(nx,nr,argv[3]).run(argv[4],steps);
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
