"""Create dimensional ParaView fields, axial profiles and the verification tables."""
import csv
from math import pi
from xml.etree.ElementTree import Element, SubElement, ElementTree
import numpy as np
from model import ROOT, BACK_PRESSURE, RHO0, UREF, T0, L, P0, solution, SHOCK
from write_mesh import grid

OUT=ROOT/"Results"
CASES=[("subsonic",128,8),("onset",128,8),("shock",64,4),("shock",128,8),("shock",256,16),("supersonic",128,8)]

def save_table(path,columns,rows):
    with path.open("w",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(columns)
        writer.writerows(rows)

def export_vtu(path,xy,fields,nx,nr):
    root=Element("VTKFile",type="UnstructuredGrid",version="0.1",byte_order="LittleEndian")
    block=SubElement(SubElement(root,"UnstructuredGrid"),"Piece",
                     NumberOfPoints=str(xy.shape[0]*xy.shape[1]),NumberOfCells=str(nx*nr))
    def array(parent,name,data,kind="Float64",components=1):
        element=SubElement(parent,"DataArray",type=kind,Name=name,format="ascii",NumberOfComponents=str(components))
        vals=np.asarray(data).ravel()
        element.text=" ".join((f"{v:.11g}" if kind=="Float64" else str(int(v))) for v in vals)
    points=xy.reshape(-1,2)
    array(SubElement(block,"Points"),"Coordinates",np.column_stack((points,np.zeros(len(points)))),components=3)
    cells=SubElement(block,"Cells")
    a=np.arange((nx+1)*(nr+1)).reshape(nx+1,nr+1)
    indices=np.stack([a[:-1,:-1],a[1:,:-1],a[1:,1:],a[:-1,1:]],axis=2)
    array(cells,"connectivity",indices,"Int32")
    array(cells,"offsets",4*np.arange(1,nx*nr+1),"Int32")
    array(cells,"types",np.full(nx*nr,9),"UInt8")
    values=SubElement(block,"CellData",Scalars="Mach",Vectors="Velocity_m_s")
    for name,data in fields.items():
        array(values,name,data,components=(3 if name=="Velocity_m_s" else 1))
    ElementTree(root).write(path,encoding="UTF-8",xml_declaration=True)

def get_fields(case,nx,nr):
    directory=OUT/f"{case}_{nx}x{nr}"
    f=np.genfromtxt(directory/"field.csv",delimiter=",",names=True).reshape(nx,nr)
    xy=grid(nx,nr)
    centers=(xy[1:,1:]+xy[1:,:-1]+xy[:-1,1:]+xy[:-1,:-1])/4
    assert np.max(np.abs(centers[:,:,0]/L-f["x_over_L"]))<1e-10
    assert np.max(np.abs(centers[:,:,1]/L-f["r_over_L"]))<1e-10
    rho,u,v,p,m=(f[key] for key in f.dtype.names[4:])
    assert np.all(np.isfinite(m)) and np.all(p>0) and np.all(rho>0)
    assert np.max(np.abs(m-np.hypot(u,v)/np.sqrt(1.4*p/rho))) < 2e-10
    return directory,xy,f,dict(Mach=m,StaticPressure_Pa=P0*p,p_over_p0=p,
           Temperature_K=T0*p/rho,Density_kg_m3=RHO0*rho,
           Velocity_m_s=UREF*np.stack([u,v,np.zeros_like(v)],axis=-1),
           StagnationPressure_Pa=P0*p*(1+0.2*m*m)**3.5)

def main():
    records=[]
    for case,nx,nr in CASES:
        directory,xy,f,fields=get_fields(case,nx,nr)
        export_vtu(directory/"flow.vtu",xy,fields,nx,nr)
        x=f["x_over_L"][:,0]; w=(2*np.arange(nr)+1)/nr**2
        m=fields["Mach"]; p=fields["p_over_p0"]
        ma,pa=m@w,p@w; po=fields["StagnationPressure_Pa"]@w/P0
        ref=np.array([solution(q,case) for q in x])
        save_table(directory/"axial_profile.csv",
            ["x_over_L","x_mm","M_area_mean","M_near_axis","r_near_axis_mm","M_outermost_cell",
             "p_area_mean_over_p0","p_near_axis_over_p0","pt_area_mean_over_p0","M_quasi1D","p_quasi1D_over_p0"],
            zip(x,x*L*1000,ma,m[:,0],f["r_over_L"][:,0]*L*1000,m[:,-1],pa,p[:,0],po,*ref.T))
        hist=np.genfromtxt(directory/"history.csv",names=True,delimiter=",")[-1]
        l1=max(hist[key] for key in hist.dtype.names[1:5])
        imbalance=abs(hist["inlet_flux"]-hist["outlet_flux"])/hist["inlet_flux"]
        assert l1<2e-5 and imbalance<2e-5
        ishock=nx//2+np.argmax(np.diff(pa)[nx//2:])+1
        shock=ishock/nx if case=="shock" else ""
        mask=np.ones(nx,dtype=bool)
        if case=="shock": mask=np.abs(x-SHOCK[0])>0.025
        mt=np.interp(0.5,x,ma)
        rec=dict(case=case,Nx=nx,Nr=nr,cells=nx*nr,pb_over_p0=BACK_PRESSURE[case],
                 M_throat_CFD=mt,M_throat_1D=solution(0.5,case)[0],M_exit_cell_CFD=ma[-1],
                 M_exit_1D=solution(1,case)[0],p_exit_cell_over_p0=pa[-1],
                 shock_x_over_L=shock,p_at_080L_over_p0=np.interp(0.8,x,pa),
                 mean_abs_M_error_outside_shock=np.mean(np.abs(ma[mask]-ref[mask,0])),
                 mdot_in_kg_s=hist["inlet_flux"]*2*pi*RHO0*UREF*L**2,
                 mass_imbalance_pct=100*imbalance,max_L1_residual=l1,iterations=int(hist["iteration"]))
        records.append(rec)
    save_table(OUT/"case_summary.csv",list(records[0]),[r.values() for r in records])
    rows=[]
    for case,pb in BACK_PRESSURE.items():
        xs=list(np.linspace(0,1,501))
        if case=="shock": xs+= [SHOCK[0],SHOCK[0]]
        side=False
        for x in sorted(xs):
            m,p=solution(x,case,side)
            rows.append((case,pb,x,x*1000*L,m,p))
            if x==SHOCK[0]: side=True
    save_table(OUT/"quasi1D_reference.csv",["case","pb_over_p0","x_over_L","x_mm","Mach","p_over_p0"],rows)
    print("Wrote six verified CFD solutions and the independent quasi-1D curves.")
if __name__=="__main__": main()
