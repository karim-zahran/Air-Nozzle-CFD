"""Write and audit 2-D ASCII Fluent meshes; coordinates are in metres."""
from collections import Counter
import numpy as np
from model import ROOT,L,radius_m

def grid(nx,nr):
    xx=np.linspace(0,L,nx+1)
    rr=np.array([radius_m(x/L) for x in xx])
    return np.array([[(x,r*j/nr) for j in range(nr+1)] for x,r in zip(xx,rr)])

def write_mesh(nx,nr):
    xy=grid(nx,nr)
    node=lambda i,j: i*(nr+1)+j+1
    cell=lambda i,j: i*nr+j+1
    zones={name:[] for name in ("interior","inlet","outlet","wall","axis")}
    # Tangent x +z points to c1; c0 is left of each directed edge.
    for i in range(nx+1):
        for j in range(nr):
            lo,hi=node(i,j),node(i,j+1)
            if i==0: zones["inlet"].append((hi,lo,cell(i,j),0))
            elif i==nx: zones["outlet"].append((lo,hi,cell(i-1,j),0))
            else: zones["interior"].append((lo,hi,cell(i-1,j),cell(i,j)))
    for i in range(nx):
        for j in range(nr+1):
            left,right=node(i,j),node(i+1,j)
            if j==0: zones["axis"].append((left,right,cell(i,0),0))
            elif j==nr: zones["wall"].append((right,left,cell(i,j-1),0))
            else: zones["interior"].append((left,right,cell(i,j),cell(i,j-1)))
    uv=xy.reshape(-1,2)
    centers=((xy[:-1,:-1]+xy[1:,:-1]+xy[:-1,1:]+xy[1:,1:])/4).reshape(-1,2)
    degree=Counter()
    edges=set()
    for rows in zones.values():
        for a,b,c0,c1 in rows:
            edge=tuple(sorted((a,b)))
            assert edge not in edges
            edges.add(edge)
            vx,vy=uv[b-1]-uv[a-1]
            for c,sign in [(c0,1),(c1,-1)]:
                if c:
                    dx,dy=centers[c-1]-uv[a-1]
                    assert sign*(vx*dy-vy*dx)>0
                    degree[c]+=1
    assert len(degree)==nx*nr and set(degree.values())=={4}
    assert len(edges)==2*nx*nr+nx+nr
    codes={"interior":(2,2,"interior"),"inlet":(3,4,"pressure-inlet"),"outlet":(4,5,"pressure-outlet"),
           "wall":(5,3,"wall"),"axis":(6,37,"axis")}
    out=ROOT/"Mesh"/f"Air_Nozzle_{nx}x{nr}.msh"
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w",newline="\n") as f:
        f.write('(0 "Axisymmetric air passage | 2D quadrilateral mesh | coordinates in metres")\n(2 2)\n')
        f.write(f"(10 (0 1 {len(uv):x} 0 2))\n(10 (7 1 {len(uv):x} 1 2)(\n")
        for x,r in uv: f.write(f"{x:.13g} {r:.13g}\n")
        f.write(f"))\n(12 (0 1 {nx*nr:x} 0))\n(12 (1 1 {nx*nr:x} 1 3))\n")
        f.write(f"(13 (0 1 {len(edges):x} 0))\n")
        index=1
        for name,rows in zones.items():
            zid,bc,_=codes[name]
            f.write(f"(13 ({zid:x} {index:x} {index+len(rows)-1:x} {bc:x} 2)(\n")
            for row in rows: f.write(" ".join(f"{a:x}" for a in row)+"\n")
            f.write("))\n")
            index+=len(rows)
        f.write("(45 (1 fluid air)())\n")
        for name,(zid,_,ztype) in codes.items(): f.write(f"(45 ({zid} {ztype} {name})())\n")
    return out

if __name__=="__main__":
    for nx,nr in ((64,4),(128,8),(256,16)):
        print(write_mesh(nx,nr))
