"""Orthographic cutaway rendered from the reimported STEP solid."""
import cadquery as cq
import numpy as np
from PIL import Image
from model import ROOT,L,RT,radius_m

part=cq.importers.importStep(str(ROOT/"CAD"/"Air_Nozzle_3D.step")).val()
# The quadrant is removed only for this illustration.
cut=part.cut(cq.Solid.makeBox(L*1000+2,40,40,cq.Vector(-1,-40,0)))
eye=np.array([.57,-.75,.34]);eye/=np.linalg.norm(eye)
right=np.cross([0,0,1],eye);right/=np.linalg.norm(right)
up=np.cross(eye,right)
key=np.array([-.18,-.48,.86]);key/=np.linalg.norm(key)
spec=key+eye;spec/=np.linalg.norm(spec)
triangles=[];colors=[]
for face in cut.Faces():
    verts,idx=face.tessellate(.028,.06)
    pts=np.array([p.toTuple() for p in verts])
    idx=np.asarray(idx)
    if face.geomType()=="PLANE":
        n=np.array(face.normalAt().toTuple())
        normals=np.tile(n,(len(pts),1)); base=np.array([.87,.65,.36])
    else:
        rad=np.hypot(pts[:,1],pts[:,2])
        t=pts[:,0]/(L*1000); c=np.cos(np.pi*t)
        slope=-(RT/L)*np.pi*c*np.sin(np.pi*t)/np.sqrt(1+c*c)
        normals=np.column_stack([-slope,pts[:,1]/rad,pts[:,2]/rad])
        normals/=np.linalg.norm(normals,axis=1)[:,None]
        outer=np.median(rad-np.array([1000*radius_m(x) for x in t]))>1.
        if not outer: normals=-normals
        base=np.array([.45,.61,.69]) if outer else np.array([.21,.61,.63])
    shade=.52+.53*np.maximum(0,normals@key)
    rgb=np.clip(base[None,:]*shade[:,None]+.21*np.maximum(0,normals@spec)[:,None]**22,0,1)
    triangles.append(np.column_stack((pts@right,-pts@up,pts@eye))[idx]);colors.append(rgb[idx])
tri=np.concatenate(triangles);rgb=np.concatenate(colors)
low=tri[:,:,:2].min(axis=(0,1))-1.8; high=tri[:,:,:2].max(axis=(0,1))+1.8
scale=3200/(high[0]-low[0]); w=3200;h=int(np.ceil(scale*(high[1]-low[1])))
tri[:,:,:2]=(tri[:,:,:2]-low)*scale
depth=np.full((h,w),-np.inf)
pixels=np.zeros((h,w,4),dtype=np.uint8)
for t,col in zip(tri,rgb):
    x=t[:,0];y=t[:,1];z=t[:,2]
    det=(y[1]-y[2])*(x[0]-x[2])+(x[2]-x[1])*(y[0]-y[2])
    if abs(det)<1e-7: continue
    x0=max(0,int(x.min()));x1=min(w,int(x.max())+2)
    y0=max(0,int(y.min()));y1=min(h,int(y.max())+2)
    yy,xx=np.ogrid[y0:y1,x0:x1]
    dx=xx+.5-x[2];dy=yy+.5-y[2]
    a=((y[1]-y[2])*dx+(x[2]-x[1])*dy)/det
    b=((y[2]-y[0])*dx+(x[0]-x[2])*dy)/det
    c=1-a-b
    zz=a*z[0]+b*z[1]+c*z[2]
    window=depth[y0:y1,x0:x1]
    visible=(a>=-1e-8)&(b>=-1e-8)&(c>=-1e-8)&(zz>window)
    if visible.any():
        window[visible]=zz[visible]
        image=pixels[y0:y1,x0:x1]
        image[visible,:3]=255*(a[visible,None]*col[0]+b[visible,None]*col[1]+c[visible,None]*col[2])
        image[visible,3]=255
image=Image.fromarray(pixels).resize((1800,round(1800*h/w)),Image.Resampling.LANCZOS)
out=ROOT/"Report"/"figures";out.mkdir(parents=True,exist_ok=True)
image.save(out/"cad_cutaway.png")
