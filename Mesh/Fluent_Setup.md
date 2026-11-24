# Repeating the air-nozzle study in Ansys Fluent

The shipped CFD fields were computed with the included axisymmetric Euler solver. The following setup reproduces its geometry and physical assumptions in Fluent for an independent comparison.

1. Start Fluent in **2D, double precision**. Read **Air_Nozzle_128x8.msh** using **File → Read → Mesh**. The mesh is already in metres. Check for positive cell areas and bounds of x = 0 to 0.200 m and y = 0 to 0.01414214 m. There are 1,024 quadrilaterals. For the fine shock case, read **Air_Nozzle_256x16.msh** instead (4,096 cells).
2. Select a **steady, density-based implicit** solver and **axisymmetric** 2D space. Set the viscous model to **inviscid** and enable the energy equation. Leave gravity off.
3. Choose air with **ideal-gas density** and constant properties giving γ = 1.4. Use R = 287.05 J/(kg·K) and cp = 1004.675 J/(kg·K). With universal R = 8314.46 J/(kmol·K), the corresponding molecular weight is about 28.965 kg/kmol. Apply the material to the zone named **air**.
4. Set the **operating pressure to 0 Pa** so that gauge pressure inputs equal the absolute values below.

| Zone | Type | Setting |
|---|---|---|
| inlet | pressure-inlet | Gauge total pressure 100,000 Pa. Total temperature 300 K. Flow normal to inlet, directed into the domain |
| outlet | pressure-outlet | Use the absolute static pressure for the selected case. Backflow total temperature 300 K |
| wall | wall | Adiabatic, slip condition of the inviscid model |
| axis | axis | Along y = 0, the x-axis of revolution |
| interior | interior | Connects neighboring control volumes |

| Case | Pressure outlet (Pa) | Quasi-1D exit Mach |
|---|---:|---:|
| Subsonic | 98,000.0 | 0.17013 |
| Onset of choking | 93,716.25 | 0.30590 |
| Internal shock | 75,000.0 | 0.38034 |
| Matched supersonic exit | 9,393.26 | 2.19720 |

5. Initialize the flow from the inlet or a converged nearby case. The inlet always remains subsonic in the reference solutions. Where Fluent requests an initial static inlet pressure, about 93,700 Pa is a starting estimate for the choked cases, and 98,000 Pa for the subsonic case. For the supersonic-exit case, Fluent extrapolates the local outlet state once the exit becomes supersonic.
6. Use second-order upwind flow discretization for the comparison. Monitor the inlet/outlet mass-flow rates, throat Mach and static-pressure profile. Confirm the reported shock remains within one fine axial cell on further iterations and that the boundary mass mismatch is small relative to inlet flow. Fluent's scaled residuals and the custom solver's residuals have different definitions.
7. Plot Mach and static pressure along the axis. For the closest numerical comparison with **axial_profile.csv**, create axial stations and use **area-weighted** reports across them. A line on the axis and a section mean can differ, particularly near the shock. The independent fine-grid target is a pressure jump near x = 0.144 m. The theoretical jumps occur at x = 0.14361 m.

The initial back-pressure sweep shares one geometry. The normalized area is 1 + cos²(πx/L) relative to the throat. For CAD-based remeshing, import **../CAD/CFD_Fluid_Section_2D.step** as a surface body. Keep the positive-y half of the XY plane and its four boundaries.

## Sources

- [Ansys, axisymmetric model orientation](https://innovationspace.ansys.com/knowledge/forums/topic/what-are-the-requirements-for-an-axisymmetric-analysis/)
- [Fluent User’s Guide, Compressible Flows](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_sec_compressible.html)
- [Fluent User’s Guide, pressure outlet boundary conditions](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_bcs_sec_bound_cond.html)
- [Fluent User’s Guide, legacy mesh format](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_app_form_case_data.html)
