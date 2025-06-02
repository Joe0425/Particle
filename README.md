# Particle

A particle tracking module designed for multi-block overset grids [1] is presented here. By integrating the current particle module with an appropriate fluid solver, we can perform direct numerical simulations of particle-laden flows with complex geometries such as turbomachinery flows, and an example of the numerical configuration is shown in the figure below.

<div align=center><img width="1000" src="./images/fig1.png"/></div>

The present module for particle tracking is developed by Mr. Taiyang Wang, who is funded and advised by Prof. Yaomin Zhao from Peking University.

Programming language: Fortran

External libraries: MPI


The main features of this module are summarized:
1. Based on the particle location, each particle that is located in the valid domain is assigned to either the background H-type grid or the O-type grid. For particles located in overlapping regions, they are stored in the O-type grid. This procedure is conducted before the simulation by the subroutine read_particles_data_init, and the implementation details are shown in the figure below.
   <div align=center><img width="500" src="./images/fig2.png"/></div>
   In this process, owing to the curved geometry of the O-type grid, both the ray-casting [2] and angle summation algorithms [3] are utilized, corresponding to the functions particle_in_polygon and particle_in_convex_quadrilateral, respectively.
