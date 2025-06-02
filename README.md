# Particle

A particle tracking module designed for multi-block overset grids [1] is presented here. By integrating the current particle module with an appropriate fluid solver, we can perform direct numerical simulations of particle-laden flows with complex geometries such as turbomachinery flows, and an example of the numerical configuration is shown in the figure below.

<div align=center><img width="800" src="./images/fig1.png"/></div>

# Information
The present module for particle tracking is developed by Mr. Taiyang Wang, who is funded and advised by Prof. Yaomin Zhao from Peking University.

Programming language: Fortran

External libraries: MPI

# Features

The main features of this module are summarized:

1. Based on the particle location, each particle that is located in the valid domain is assigned to either the background H-type grid or the O-type grid. For particles located in overlapping regions, they are stored in the O-type grid. This procedure is conducted before the simulation by the subroutine ***read_particles_data_init***, and the implementation details are shown in the figure below.
   <div align=center><img width="400" src="./images/fig2.png"/></div>
   
   In this process, owing to the curved geometry of the O-type grid, both the ray-casting [2] and angle summation algorithms [3] are utilized, corresponding to the functions ***particle_in_polygon*** and ***particle_in_convex_quadrilateral***, respectively.
2. Neglecting particle rotation and temperature, each point particle is modeled as a spherical rigid body. Moreover, particle-wall interactions are described by the hard-sphere model, which is achieved by the subroutine ***bounce_particles***.
   <div align=center><img width="400" src="./images/fig3.png"/></div>

3. Affected by the drag force calculated in the subroutine ***compute_particles_force***, both the particle location and velocity need to be updated at each time step. In the present study, the second-order Velocity-Verlet algorithm is used, which is implemented by the subroutine ***timeadvance_redistribute_particles***. Besides, a trilinear scheme is applied for interpolation, where the weights are obtained by the function ***interpolation_weights***.

4. Utilizing the MPI parallelization strategy, particles in each block are handled separately on each processor unless they pass through the interface between blocks, namely the outer boundary of the O-type grid. Consequently, different situations may occur, which should be handled carefully. This tricky problem is solved in the subroutine ***timeadvance_redistribute_particles***, and the details are shown by the flow chart below.
   <div align=center><img width="400" src="./images/fig4.png"/, img width="400" src="./images/fig4.png"/>
