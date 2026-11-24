"""Geometry, ideal-air reference units and quasi-1D analytical solutions."""
from math import cos, pi, sqrt, acos
from pathlib import Path
import os

ROOT = Path(os.environ.get("AIR_NOZZLE_OUTPUT", Path(__file__).resolve().parents[1])).resolve()
L = 0.200                      # metres
RT = 0.010                     # throat radius, metres
DISPLAY_WALL = 0.002           # radial thickness for the 3-D CAD view
P0, T0, R, GAMMA = 100000.0, 300.0, 287.05, 1.4
RHO0, UREF = P0/(R*T0), sqrt(R*T0)

def area_ratio(x_over_l):
    return 1.0 + cos(pi*x_over_l)**2

def radius_m(x_over_l):
    return RT*sqrt(area_ratio(x_over_l))

def area_mach(m):
    return ((5.0 + m*m)/6.0)**3/m

def mach_for_area(ratio, supersonic=False):
    if ratio < 1-1e-10:
        raise ValueError("A/A* must be at least one")
    if abs(ratio-1)<1e-12:
        return 1.0
    lo, hi = (1., 8.) if supersonic else (1e-8, 1.)
    for _ in range(65):
        mid = (lo+hi)/2
        if (area_mach(mid)<ratio)==supersonic:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2

def pressure(m):
    return (1+0.2*m*m)**-3.5

def normal_shock(m1):
    m2=sqrt((1+0.2*m1*m1)/(1.4*m1*m1-0.2))
    p2p1=1+7/6*(m1*m1-1)
    p02p01=pressure(m1)*p2p1/pressure(m2)
    return m2,p2p1,p02p01

def shock_state(pb):
    def trial(m1):
        m2,p2p1,p02p01=normal_shock(m1)
        # A_s/A_t = A_s/A*_1 ; downstream A*_2 changes with total pressure.
        exit_m=mach_for_area(2*area_mach(m2)/area_mach(m1))
        return p02p01*pressure(exit_m)
    lo,hi=1.,mach_for_area(2,True)
    if not trial(hi)<=pb<=trial(lo):
        raise ValueError("Back pressure outside the internal-shock regime")
    for _ in range(65):
        mid=(lo+hi)/2
        if trial(mid)>pb:
            lo=mid
        else:
            hi=mid
    m1=(lo+hi)/2
    x=acos(-sqrt(area_mach(m1)-1))/pi
    m2,p2p1,p02p01=normal_shock(m1)
    return x,m1,m2,p2p1,p02p01

BACK_PRESSURE={
    "subsonic":0.98, "onset":pressure(mach_for_area(2)), "shock":0.75,
    "supersonic":pressure(mach_for_area(2,True))
}
SHOCK=shock_state(BACK_PRESSURE["shock"])

def solution(x,case,downstream_at_shock=True):
    """Return local Mach number and static pressure divided by inlet total pressure."""
    ar=area_ratio(x)
    pb=BACK_PRESSURE[case]
    ptotal=1.0
    if case=="subsonic":
        exit_m=sqrt(5*(pb**(-2/7)-1))
        mach=mach_for_area(ar*area_mach(exit_m)/2)
    elif case=="onset":
        mach=mach_for_area(ar)
    elif case=="supersonic":
        mach=mach_for_area(ar,x>0.5)
    else:
        xs,m1,m2,_,p02=SHOCK
        if x>xs or (x==xs and downstream_at_shock):
            ptotal=p02
            mach=mach_for_area(ar*area_mach(m2)/area_mach(m1))
        else:
            mach=mach_for_area(ar,x>0.5)
    return mach,ptotal*pressure(mach)
