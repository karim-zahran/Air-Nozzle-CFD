"""Recreate the report figures from the saved finite-volume fields and histories."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import ROOT,L,SHOCK
from write_mesh import grid

out=ROOT/"Report"/"figures";out.mkdir(parents=True,exist_ok=True)
results=ROOT/"Results"
ink="#282c31"; gray="#44484b"; teal="#226979"; orange="#a34420"
palette={"subsonic":"#426882","onset":"#008778","shock":orange,"supersonic":"#634988"}
names={"subsonic":"Subsonic","onset":"Choking limit","shock":"Internal shock","supersonic":"Supersonic exit"}
plt.rcParams.update({"font.family":"serif","font.serif":["STIXGeneral","DejaVu Serif"],"mathtext.fontset":"stix","pdf.fonttype":42,"font.size":9,"text.color":ink,
    "axes.labelcolor":ink,"axes.edgecolor":"#a6b9bb","xtick.color":gray,"ytick.color":gray,
    "axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False,"lines.linewidth":2,
    "figure.facecolor":"white","axes.facecolor":"white","axes.labelpad":5,
    "savefig.dpi":240})
def load(path): return np.genfromtxt(path,delimiter=",",names=True)
def save(fig,name):
    fig.savefig(out/name,bbox_inches="tight",pad_inches=.04,facecolor="white")
    vector_dir=ROOT/"Report"/"LaTeX"/"figures"
    vector_dir.mkdir(parents=True,exist_ok=True)
    fig.savefig(vector_dir/(name[:-4]+".pdf"),bbox_inches="tight",pad_inches=.04,facecolor="white")
    plt.close(fig)
def ticks(ax):
    ax.grid(axis="y",alpha=.25,color="#abc0bf",lw=.7)
    ax.set_axisbelow(True);ax.tick_params(length=3,width=.7)

# A scientific field plot from the fine-grid cells. Reflect for a full-bore view.
xy=grid(256,16)*1000
X=xy[:,:,0].T;Y=xy[:,:,1].T
m=load(results/"shock_256x16"/"field.csv")["Mach"].reshape(256,16).T
fig,ax=plt.subplots(figsize=(8.2,2.15))
for sign in (1,-1):
    coll=ax.pcolormesh(X,sign*Y,m,cmap="viridis",vmin=.25,vmax=1.80,shading="flat",rasterized=True)
    ax.plot(X[-1],sign*Y[-1],c=ink,lw=1.15)
ax.set_xlim(0,200);ax.set_ylim(-16,23);ax.set_yticks([-10,0,10]);ax.set_xticks([0,50,100,150,200])
ax.set_xlabel("Axial distance (mm)");ax.set_ylabel("Transverse position (mm)")
for xx,txt in ((100,"Throat"),(143.75,"Shock")):
    ax.annotate(txt,(xx,10.8 if xx==100 else 13.0),xytext=(xx,19.3),ha="center",fontsize=9,
                arrowprops={"arrowstyle":"-","color":ink,"lw":.8})
cbar=fig.colorbar(coll,ax=ax,fraction=.033,pad=.023,aspect=16,ticks=[.3,.6,.9,1.2,1.5,1.8])
cbar.set_label("Mach number");cbar.outline.set_linewidth(.6)
save(fig,"mach_field.png")

# Positive-r computational domain and all edges of the coarse quadrilateral grid.
xy=grid(64,4)*1000
fig,ax=plt.subplots(figsize=(7.9,2.0))
x=xy[:,0,0];r=xy[:,-1,1]
ax.fill_between(x,0,r,color="#e5f3f2")
for row in xy: ax.plot(row[:,0],row[:,1],c="#8eb8b9",lw=.49)
for col in xy.transpose(1,0,2): ax.plot(col[:,0],col[:,1],c="#75a8aa",lw=.68)
ax.plot(x,r,c=ink,lw=1.7);ax.plot(x,0*x,c=ink,lw=1.2)
ax.plot([0,0],[0,r[0]],c=teal,lw=2);ax.plot([200,200],[0,r[-1]],c=orange,lw=2)
ax.text(2,15.1,"Inlet",color=teal,ha="left");ax.text(198,15.1,"Outlet",color=orange,ha="right")
ax.text(100,11.5,"Throat",ha="center");ax.text(145,14.4,"Slip wall",ha="center")
ax.text(100,-2.4,"Axis of revolution",ha="center")
ax.set_xlim(-3,203);ax.set_ylim(-3.8,18);ax.set_yticks([0,5,10,15]);ax.set_xticks([0,50,100,150,200])
ax.set_ylabel("Radius (mm)");ax.set_xlabel("Axial distance (mm)")
save(fig,"mesh.png")

fig,axes=plt.subplots(2,1,figsize=(7.8,4.5),sharex=True,gridspec_kw={"hspace":.16})
styles={"subsonic":"-","onset":(0,(4,2)),"shock":"-","supersonic":(0,(2,2))}
for case in names:
    data=load(results/f"{case}_128x8"/"axial_profile.csv")
    for ax,col in zip(axes,["M_area_mean","p_area_mean_over_p0"]):
        ax.plot(data["x_mm"],data[col],label=names[case],color=palette[case],linestyle=styles[case],
                lw=2.3 if case=="shock" else 1.8)
for ax in axes:
    ax.axvline(100,color="#93a6a6",lw=.9,ls=":");ticks(ax);ax.set_xlim(0,200)
axes[0].set_ylim(0,2.38);axes[1].set_ylim(0,1.06)
axes[0].set_ylabel("Mach number");axes[1].set_ylabel(r"Static pressure, $p/p_0$")
axes[1].set_xlabel("Axial distance (mm)");axes[1].set_xticks(np.arange(0,201,25))
axes[0].legend(ncol=4,loc="lower center",bbox_to_anchor=(.5,1.02),columnspacing=1.35,handlelength=2.2,fontsize=9)
save(fig,"four_cases.png")

fig,axes=plt.subplots(1,2,figsize=(8.2,2.42),gridspec_kw={"wspace":.32})
ref=np.genfromtxt(results/"quasi1D_reference.csv",delimiter=",",names=True,encoding="utf8",dtype=None)
ref=ref[ref["case"]=="shock"]
fine=load(results/"shock_256x16"/"axial_profile.csv")
for ax,var,mn,near,ylab in zip(axes,["Mach","p_over_p0"],["M_area_mean","p_area_mean_over_p0"],
                           ["M_near_axis","p_near_axis_over_p0"],["Mach number",r"Static pressure, $p/p_0$"]):
    aw,=ax.plot(ref["x_mm"],ref[var],c=ink,ls="--",lw=1.5,label="Quasi-1D")
    cw,=ax.plot(fine["x_mm"],fine[mn],c=teal,lw=2,label="CFD area mean")
    near_pts=(fine["x_mm"]>130)&(fine["x_mm"]<156)
    nw,=ax.plot(fine["x_mm"][near_pts][::2],fine[near][near_pts][::2],linestyle="none",marker="o",ms=3.4,
                mfc="white",mec=orange,label="CFD nearest axis")
    ax.set_xlim(131,156);ax.set_xticks([135,140,145,150,155]);ax.set_xlabel("Axial distance (mm)")
    ax.set_ylabel(ylab);ticks(ax)
axes[0].set_ylim(.45,1.85);axes[1].set_ylim(.13,.78)
fig.legend([aw,cw,nw],["Quasi-1D","CFD area mean","CFD nearest axis"],ncol=3,
           loc="lower center",bbox_to_anchor=(.5,.96),columnspacing=1.4)
save(fig,"shock_detail.png")

fig,ax=plt.subplots(figsize=(7.75,1.95))
for nx,nr,col in [(64,4,"#788a96"),(128,8,teal),(256,16,orange)]:
    h=load(results/f"shock_{nx}x{nr}"/"history.csv")
    worst=np.maximum.reduce([h[k] for k in h.dtype.names[1:5]])
    ax.semilogy(h["iteration"],worst,c=col,label=f"{nx} × {nr}",lw=1.4)
ax.axhline(2e-5,c=ink,ls=":",lw=.9,label="Stopping threshold")
ax.set_ylim(5e-6,30);ax.set_xlim(-100,14900);ax.set_yticks([1e-5,1e-3,1e-1,10])
ax.set_xlabel("Pseudo-time iterations");ax.set_ylabel(r"Largest $L_1$ residual")
ax.legend(ncol=4,loc="lower center",bbox_to_anchor=(.5,1.0),fontsize=8.5);ticks(ax)
save(fig,"residuals.png")
