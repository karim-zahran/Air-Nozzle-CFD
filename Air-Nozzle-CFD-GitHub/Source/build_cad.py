"""Build three STEP models in mm and verify their topology and dimensions on reimport."""
from math import pi, sqrt
from pathlib import Path
import json
import cadquery as cq
import numpy as np
from model import ROOT,L,RT,DISPLAY_WALL,radius_m

OUT=ROOT/"CAD"
OUT.mkdir(exist_ok=True)
mm=1000
length=L*mm
wall_thickness=DISPLAY_WALL*mm
samples=np.linspace(0,length,81)
V=cq.Vector
def point(x,y): return V(float(x),float(y),0)
def line(v1,v2): return cq.Edge.makeLine(v1,v2)
def build_wall(offset):
    points=[point(x,radius_m(x/length)*mm+offset) for x in samples]
    edge=cq.Edge.makeSpline(points,tangents=[V(1,0,0),V(1,0,0)],parameters=samples.tolist(),scale=False)
    return points,edge
pts,inner=build_wall(0)
outer_pts,outer=build_wall(wall_thickness)
a,b=point(0,0),point(length,0)
section=cq.Wire.assembleEdges([line(a,b),line(b,pts[-1]),inner,line(pts[0],a)])
annulus=cq.Wire.assembleEdges([inner,line(pts[-1],outer_pts[-1]),outer,line(outer_pts[0],pts[0])])
face=cq.Face.makeFromWires(section)
volume=cq.Solid.revolve(section,[],360,a,b)
nozzle=cq.Solid.revolve(annulus,[],360,a,b)
shapes={"Air_Nozzle_3D":nozzle, "Air_Volume_3D":volume, "CFD_Fluid_Section_2D":face}
verified={}
for name,shape in shapes.items():
    path=OUT/(name+".step")
    assert shape.isValid()
    shape.exportStep(str(path))
    imported=cq.importers.importStep(str(path)).val()
    assert imported.isValid()
    box=imported.BoundingBox()
    assert abs(box.xlen-length)<1e-4
    assert ("AUTOMOTIVE_DESIGN" in path.read_text() or "CONFIG_CONTROL_DESIGN" in path.read_text())
    values=dict(faces=len(imported.Faces()), solids=len(imported.Solids()),length_mm=box.xlen)
    if name=="CFD_Fluid_Section_2D":
        assert values["faces"]==1 and len(imported.Edges())==4 and box.zlen<1e-4
        curve=next(e for e in imported.Edges() if e.geomType()!="LINE")
        points=[curve.positionAt(float(t),mode="parameter") for t in np.linspace(*curve.bounds(),1001)]
        error=max(abs(p.y-mm*radius_m(p.x/length)) for p in points)
        assert error<0.0001
        values["max_sampled_wall_error_mm"]=error
        values["area_mm2"]=imported.Area()
    else:
        assert values["solids"]==1 and len(imported.Shells())==1
        # Intersect the imported part at the throat and near both ends.
        for x in (0.1, length/2, length-0.1):
            wires=cq.Workplane("YZ").newObject([imported]).section(x).wires().vals()
            radii=sorted((w.BoundingBox().ylen/2 for w in wires))
            expected=mm*radius_m(x/length)
            assert abs(radii[0]-expected)<0.0001
            if name=="Air_Nozzle_3D":
                assert len(radii)==2 and abs(radii[1]-expected-wall_thickness)<0.0001
            else:
                assert len(radii)==1
            if x==length/2: values["throat_bore_radius_mm"]=radii[0]
        values["volume_mm3"]=imported.Volume()
        if name=="Air_Volume_3D":
            exact=1.5*pi*(mm*RT)**2*length
            assert abs(imported.Volume()/exact-1)<1e-6
    verified[name]=values
(OUT/"cad_checks.json").write_text(json.dumps(verified,indent=2)+"\n")
print(json.dumps(verified,indent=2))
