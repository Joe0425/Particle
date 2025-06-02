
      module particle

      use stdtypes

      implicit none

cty   Define the type arrays contain the information of the particle.
cty   Notice that the total global id is defined as the float rather than integer number,
cty   thus more particle numbers can be marked in the whole simulation. 
      type particle
        real(kind=wp) :: idpar
        real(kind=wp) :: dpar, rhopar
        real(kind=wp) :: transpar
        real(kind=wp) :: xpar, ypar
        real(kind=wp) :: zpar
        real(kind=wp) :: upar, vpar, wpar
        real(kind=wp) :: axpar, aypar, azpar
        integer :: procpar
        integer :: ipar, jpar, kpar 
cty     The integer is_monitor is used to monitor the particles!
cty     is_monitor .eq. 999 -> don't monitor the particle
cty     is_monitor .eq. 1 -> monitor the particle
        integer :: is_monitor
      end type
      type(particle), dimension(:), pointer :: Pointer_par

      integer, dimension(:), pointer :: idpar_proc_del
      integer, allocatable :: mark_block_id(:)
     
      integer, save :: Npar_local_collision
      real(kind=wp), allocatable, save :: pardata_real_collision(:,:)
      integer, allocatable, save :: pardata_integer_collision(:,:)

cty   Npar is the particles number at each time step 
      integer, parameter :: Npar_Max = 10000000, 
     &                      Npar_del_Max = 10000000, 
     &                      Npar_send_Max = 10000000
      integer :: Npar_local, Npar_global
      integer :: Npar_del, Npar_send

      integer, parameter :: Npar_read_onetime = 10000000
 
cty   1: idpar 
cty   2: dpar
cty   3: rhopar, 
cty   4: transpar
cty   5: xpar
cty   6: ypar
cty   7: zpar
cty   8: upar
cty   9: vpar
cty   10: wpar
cty   11: axpar
cty   12: aypar
cty   13: azpar 
      integer, parameter :: Nvar_par_real = 13
cty   procpar, ipar, jpar, kpar, is_monitor 
      integer, parameter :: Nvar_par_integer = 5

cty   Control dimensional parameters 
      type Para_dimensional
        real(kind=wp) :: rho_ref
        real(kind=wp) :: U_ref
        real(kind=wp) :: L_ref
      end type
      type(Para_dimensional) :: Para_dim

cty   Control non-dimensional parameters
      type Para_non_dimensional 
        real(kind=wp) :: diameter_par
        real(kind=wp) :: density_par
      end type
      type(Para_non_dimensional) :: Para_non_dim

cty   Control parameters of particle computation 
      type Para_particle 
        integer :: imode
        integer :: Initial_velocity_particle
        integer :: Monitor_particle_steps
        integer :: Capture_particle_steps
        integer :: Store_particle_steps
        integer :: Average_particle_steps
        integer :: Collision_particle_steps
        integer :: Add_particle_steps
        integer :: Add_particle_number
        real(kind=wp) :: diameter_par
        real(kind=wp) :: density_par
      end type
      type(Para_particle) :: Para_par

cty   Id of the monitor particles 
      integer, save :: Npar_monitor
      real(kind=wp), allocatable, save :: id_monitor(:)
      integer :: Npar_monitor_local, Npar_monitor_global

cty   Number of the processors corresponding to background and O-type grids
      integer, save :: B_nprocs
      integer, save :: O_nprocs

cty   Boundary of the background and O-type grids
      real(kind=wp), save :: B_xmin, B_xmax
     &                     , B_ymin, B_ymax
      real(kind=wp), save :: B_zmin, B_zmax
      real(kind=wp), save :: O_outer_xmin, O_outer_xmax, 
     &                       O_inner_xmin, O_inner_xmax
      real(kind=wp), allocatable, save :: O_inner_px(:) 
     &                                  , O_inner_py(:) 
     &                                  , O_inner_py_trans(:)
     &                                  , O_outer_px(:)
     &                                  , O_outer_py(:) 
     &                                  , O_outer_py_trans(:)

cty   Boundary of the subblock of the background and the O-type grids
      type subblock_coordinates
cty     size of the subblock_coordinates
        integer :: xsize, ysize
        real(kind=wp), allocatable :: coords(:,:,:)
cty     xmin, xmax, ymin, ymax is specially for the background grid
        real(kind=wp) :: xmin, xmax
        real(kind=wp) :: ymin, ymax
cty     px(:), py(:), py_trans(:) is specially for the O-type grid
        real(kind=wp), allocatable :: px(:)
        real(kind=wp), allocatable :: py(:)
     &                              , py_trans(:)
      end type
      type(subblock_coordinates), allocatable, save :: subblock_proc(:)
      real(kind=wp), allocatable, save :: B_z_x(:)
      real(kind=wp), allocatable, save :: B_r_y(:)
      real(kind=wp), allocatable, save :: B_theta(:)

cty   Mapping relationships between the H- type and O-type blocks
      type mapping_particles
        integer, allocatable :: nearest_proc(:,:)
        integer, allocatable :: nearest_i(:,:)
        integer, allocatable :: nearest_j(:,:)
      end type
      type(mapping_particles), allocatable, save :: mapping(:)

cty   Process-related parameters 
      integer, allocatable, save :: block_proc(:)
      integer, allocatable, save :: ind_i_proc(:), ind_j_proc(:)
      integer, allocatable, save :: ibeg_proc(:), iend_proc(:),
     &                              jbeg_proc(:), jend_proc(:)
      integer, allocatable, save :: xm_proc(:), xp_proc(:),
     &                              ym_proc(:), yp_proc(:)
      integer, allocatable, save :: xmym_proc(:), xpym_proc(:),
     &                              xpyp_proc(:), xmyp_proc(:)
cty   Assemble 3 nearby proccessors in the list
      type nearby_proc_list
        integer :: id(3)
        integer :: length
      end type
      type(nearby_proc_list), allocatable, save :: xmym_proc_list(:),
     &                                             xpym_proc_list(:),
     &                                             xmyp_proc_list(:),
     &                                             xpyp_proc_list(:)
      type(nearby_proc_list), allocatable, save :: ymxm_proc_list(:),
     &                                             ypxm_proc_list(:),
     &                                             ymxp_proc_list(:),
     &                                             ypxp_proc_list(:)

cty   Blade surface vectors
      real(kind=wp), allocatable, save :: 
     &               normal_x_global(:), normal_y_global(:), 
     &               tangen_x_global(:), tangen_y_global(:)

cty   Arrays contain the added particles data 
      integer, save :: Npar_added_local, Npar_added_global
      real(kind=wp), allocatable, save :: added_pardata_real(:,:)
      integer, allocatable, save :: added_pardata_integer(:,:)

cty   Define the statistics for particles
      real(kind=wp), allocatable, save :: qstat_par(:,:,:,:)
      real(kind=wp), allocatable, save :: B_qstat_par_gather(:,:,:,:)
      real(kind=wp), allocatable, save :: O_qstat_par_gather(:,:,:,:)
      integer, parameter :: n_qstat_par_vel = 10
      integer, save :: n_qstat_par
      integer, parameter :: nsprint_par = 1
cty   (1) n -> n
cty   (2) u -> u
cty   (3) v -> v
cty   (4) w -> v
cty   (5)  u*u -> u'u'
cty   (6)  u*v -> v'v'
cty   (7)  w*w -> w'w'
cty   (8)  u*v -> u'v'
cty   (9)  u*w -> u'w'
cty   (10) v*w -> v'w'

      real(kind=wp), parameter :: PI = 3.141592653589792_wp

      contains



cty========================================================================
cty     subroutine init_particles
cty========================================================================
cty
cty     Description:
cty     Initialize the particle framework.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine init_particles

        use mppvars, only : ioproc
        use corevars, only : qph, phys_deriv, nzp, nvar, xhalo, yhalo
        use deriv, only : calc_grad_nsrhs      

        implicit none

        call read_dimensional_parameters
        call read_particle_parameters
        call set_non_dimensional_parameters
        call mkdir_particles_subfolders_PGI 
        call allocate_particles_memory
        call set_block_coordinates
        call init_glob_index
        call init_grids_nprocs
        call set_subblock_coordinates
        call set_B_subblock_proc_coordinates
        call set_mapping_indexes
        call init_surface_vector
        call read_particles_data_init  
        call read_monitor_particles
        call init_monitor_particles
        call init_adjacent_proc
        call init_corner_proc
        call init_nearby_proc_list
        call fill_corner_halopoints_qph
        call set_particles_velocity   
        call read_added_particles_init
        call init_particles_statistics

      end subroutine ! init_particles



cty========================================================================
cty     subroutine read_dimensional_parameters
cty========================================================================
cty
cty     Description:
cty     Read the dimensional parameters.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty       1. rho_ref:
cty          dimensional (reference) density 
cty       2. U_ref:
cty          dimensional (reference) velocity
cty       3. L_ref:
cty          dimensional (reference) length
cty
cty========================================================================

      subroutine read_dimensional_parameters
        
        use mppvars, only : real_mp_type, MPI_comm_fluid, ioid, ioproc
        
        implicit none

        real(kind=wp) :: rho_ref
        real(kind=wp) :: U_ref
        real(kind=wp) :: L_ref
        integer :: ierr

        namelist /control_dimensional/ rho_ref, U_ref, L_ref

        real(kind=wp) :: rparameters(5) 

        logical :: lex


        inquire(file='dimensional_parameter.in',exist=lex)
        if (lex) then        
          if (ioproc) then

            open(unit=99,file='dimensional_parameter.in')
            read(99,nml=control_dimensional)
            close(99)
            rparameters(1) = rho_ref
            rparameters(2) = U_ref
            rparameters(3) = L_ref

cty         rho_ref
            write(*,'(A,F20.16)') 
     &              'rho_ref = ', rparameters(1)
            if (rparameters(1) .gt. 0) then
              write(*,'(A,F20.16,A7)') 
     &                ' The dimensional density is', rparameters(1),
     &                ' kg/m^3'
            else
              call stopnow('The rho_ref must be positive')
            end if ! (rparameters(1) .gt. 0)
 
cty         U_ref
            write(*,'(A,F20.16)') 
     &              'U_ref = ', rparameters(2)
            if (rparameters(2) .gt. 0) then
              write(*,'(A,F20.16,A4)') 
     &                ' The dimensional velocity is', rparameters(2),
     &                ' m/s'
            else
              call stopnow('The U_ref must be positive')
            end if ! (rparameters(2) .gt. 0)

cty         L_ref
            write(*,'(A,F20.16)') 
     &              'L_ref = ', rparameters(3)
            if (rparameters(3) .gt. 0) then
              write(*,'(A,F20.16,A2)') 
     &                ' The dimensional length is', rparameters(3),
     &                ' m'
            else
              call stopnow('The L_ref must be positive')
            end if ! (rparameters(3) .gt. 0)

          end if ! (ioproc)

        else
          call stopnow('file dimensional_parameter.in does not exist')
        end if ! (lex)
        
        call MPI_bcast(rparameters,5,real_mp_type,
     &                 ioid,MPI_comm_fluid,ierr) 
       
        Para_dim%rho_ref = rparameters(1)
        Para_dim%U_ref = rparameters(2)
        Para_dim%L_ref = rparameters(3)

      end subroutine ! read_dimensional_parameters



cty========================================================================
cty     subroutine read_particle_parameters
cty========================================================================
cty
cty     Description:
cty     Read the parameters of particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty       1. nparameters(1)
cty          imode: new or subsequent simulation for particles?
cty            0: new particles simulation
cty            1: subsequent particles simulation
cty       2. nparameters(2)
cty          Initial_velocity_particle: how to set the initial velocity
cty            imode = 0
cty              0: zero initial velocity
cty              1: follow the fluid
cty            imode = 1
cty              given by the last simulation
cty       3. nparameters(3)
cty          Initial_vorticity_particle: how to set the initial vorticity
cty            imode = 0
cty              0: zero initial vorticity
cty              1: follow the fluid
cty            imode = 1
cty              given by the last simulation
cty       4. nparameters(4)
cty          Initial_temperature_particle: how to set the initial temperature
cty            imode = 0
cty              0: one initial temperature
cty              1: follow the fluid
cty            imode = 1
cty              given by the last simulation
cty       5. nparameters(5)
cty          Monitor_particle_steps: 
cty            the interval of time steps to monitor particles data
cty            particles data will be saved in the single file
cty       6. nparameters(6)
cty          Capture_particle_steps: 
cty            the interval of time steps to capture particles data
cty            particles data will be saved in one file
cty       7. nparameters(7)
cty          Store_particle_steps: 
cty            the interval of time steps to store particles data 
cty            for subsequent particles simulation
cty       8. nparameters(8)
cty          Average_particle_steps: 
cty            the interval of time steps to gather particles statistics data 
cty       9. nparameters(9)
cty          Collision_particle_steps: 
cty            the interval of time steps to collect data of particles which
cty            collide with the wall.
cty       10. nparameters(10)
cty           Impaction_particle_steps: 
cty             the interval of time steps to collect data of particles which
cty             impact the wall.
cty       11. nparameters(11)
cty           Deposition_particle_steps: 
cty             the interval of time steps to collect data of particles which
cty             deposite on the wall.
cty       12. nparameters(12)
cty           Add_particle_steps: 
cty             the interval of time steps that particles were added into 
cty             the flow field from the inlet 
cty       13. nparameters(13)
cty           Add_particle_number: 
cty             the number that particles were added into the flow field 
cty             from the inlet
cty       14. rparameters(1)
cty           diameter_par:
cty             the dimensional diameter of the material of the single type
cty             of particle  
cty       15. rparameters(2)
cty           density_par:
cty             the dimensional density of the material of the single type 
cty             of particle
cty
cty========================================================================

      subroutine read_particle_parameters
        
        use corevars, only : nsteps, ncapt, nnstat
        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_comm_fluid, 
     &                      ioid, ioproc
        
        implicit none

        integer :: imode
        integer :: Initial_velocity_particle
        integer :: Monitor_particle_steps
        integer :: Capture_particle_steps
        integer :: Store_particle_steps
        integer :: Average_particle_steps
        integer :: Collision_particle_steps
        integer :: Add_particle_steps
        integer :: Add_particle_number
        real(kind=wp) :: diameter_par
        real(kind=wp) :: density_par

        integer :: ierr

        namelist /control_particle/ imode
     &                             ,Initial_velocity_particle
     &                             ,Monitor_particle_steps
     &                             ,Capture_particle_steps
     &                             ,Store_particle_steps
     &                             ,Average_particle_steps
     &                             ,Collision_particle_steps
     &                             ,Add_particle_steps
     &                             ,Add_particle_number
     &                             ,diameter_par
     &                             ,density_par

        integer :: nparameters(13)
        real(kind=wp) :: rparameters(12) 
        logical :: lex

        imode = 0
        Initial_velocity_particle = 1
        Monitor_particle_steps = 0
        Capture_particle_steps = 0
        Store_particle_steps = 0
        Average_particle_steps = 0
        Collision_particle_steps = 0
        Add_particle_steps = 0
        Add_particle_number = 0

        inquire(file='particle_parameter.in',exist=lex)
        if (lex) then        
          if (ioproc) then

            open(unit=99,file='particle_parameter.in')
            read(99,nml=control_particle)
            close(99)
            nparameters(1) = imode
            nparameters(2) = Initial_velocity_particle
            nparameters(5) = Monitor_particle_steps
            nparameters(6) = Capture_particle_steps
            nparameters(7) = Store_particle_steps
            nparameters(8) = Average_particle_steps
            nparameters(9) = Collision_particle_steps
            nparameters(12) = Add_particle_steps
            nparameters(13) = Add_particle_number
            rparameters(1) = diameter_par
            rparameters(2) = density_par

cty         imode
            write(*,'(A,I1)') 'imode = ', nparameters(1)
            if (nparameters(1) .eq. 0) then
              write(*,*) 'The particles simulation is new'
            else if (nparameters(1) .eq. 1) then
              write(*,*) "The particles simulation continues, "//
     &                   "and we don't change the processors"
            else 
              call stopnow
     &        ('The present imode is not defined')
            end if ! (nparameters(1) .eq. 0)

cty         Initial_velocity_particle 
            write(*,'(A,I1)') 
     &              'Initial_velocity_particle = ', nparameters(2)
            if (nparameters(1) .eq. 0) then
              if (nparameters(2) .eq. 0) then
                write(*,*) 
     &          'The initial velocity of added particles '//
     &          'are set to zero'
                write(*,*) 
     &          'The initial velocity is set to zero'
              else if (nparameters(2) .eq. 1) then
                write(*,*) 
     &          'The initial velocity of added particles '//
     &          'are interpolated from the flow field'
                write(*,*) 
     &          'The initial velocity is interpolated from '//
     &          'the flow field'
              else 
                call stopnow
     &          ('The present Initial_velocity_particle is not defined')
              end if ! (nparameters(2) .eq. 0)
            else if (nparameters(1) .eq. 1) then
              write(*,*) 
     &        'The initial velocity is given by '//
     &        'the last simulation'
            end if ! (nparameters(1) .eq. 0) 

            write(*,*) 
     &      'The collision between particles and the '//
     &      'blade surface is handled with the '//
     &      'hard-sphere model'

cty         Monitor_particle_steps
            write(*,'(A,I7)') 
     &              'Monitor_particle_steps = ', nparameters(5)
            if (nparameters(5) .le. 0) then
              write(*,*) 
     &        'The Monitor_particle_steps is invalid, '//
     &        'so the subroutine monitor_particles_data is invalid'
            else 
              if (mod(nsteps,nparameters(5)) .ne. 0) then
                call stopnow
     &          ('The Monitor_particle_steps must be '//
     &           'a divider of nsteps')
              else 
                write(*,'(A,I7,A6)') 
     &          ' The particle data will be saved in single file'//
     &          ' at every', nparameters(5), 'steps'
              end if ! (mod(nsteps,nparameters(5)) .ne. 0)
            end if ! (nparameters(5) .le. 0)

cty         Capture_particle_steps
            write(*,'(A,I7)') 
     &              'Capture_particle_steps = ', nparameters(6)
            if (nparameters(6) .le. 0) then
              write(*,*) 
     &        'The Capture_particle_steps is invilid, '//
     &        'so the subroutine capture_particles_data is invalid'
            else 
              if (mod(nsteps,nparameters(6)) .ne. 0) then
                call stopnow
     &          ('The Capture_particle_steps must be '//
     &           'a divider of nsteps')
              else 
                write(*,'(A,I7,A6)') 
     &          ' The particle data will be saved in one file'//
     &          ' at every', nparameters(6), 'steps'
              end if ! (mod(nsteps,nparameters(6)) .ne. 0)
            end if ! (nparameters(6) .le. 0)

cty         Store_particle_steps
            write(*,'(A,I7)') 
     &              'Store_particle_steps = ', nparameters(7)
            if (nparameters(7) .ne. ncapt) then 
              call stopnow
     &       ('The Store_particle_steps must be equal to ncapt')
            else 
              write(*,'(A,I7,A6)') 
     &        ' The particle data is stored for each processor'//
     &        ' at every', nparameters(7), 'steps'
            end if

cty         Average_particle_steps
            write(*,'(A,I7)') 
     &              'Average_particle_steps = ', nparameters(8)
            if (nparameters(8) .eq. 0) then 
              write(*,*) 
     &        'The Average_particle_steps is invilid, '//
     &        'so the statistics of particles will not be computed'
            else 
              if (nparameters(8) .ne. nnstat) then 
                call stopnow
     &         ('The Average_particle_steps must be equal to nnstat')
              else
                write(*,'(A,I7,A6)') 
     &          ' The particle statistics will be saved in one file'//
     &          ' at every', nparameters(8), 'steps'
              end if ! (nparameters(8) .ne. nnstat)
            end if ! (nparameters(8) .eq. 0)

cty         Collision_particle_steps
            write(*,'(A,I7)') 
     &              'Collision_particle_steps = ', nparameters(9)
            if (nparameters(9) .eq. 0) then 
              write(*,*) 
     &        'The Collision_particle_steps is invilid, so the '//
     &        'subroutine collect_collision_particles_data is invalid'
            else 
              if (mod(nnstat,nparameters(9)) .ne. 0) then 
                call stopnow
     &         ('The Collision_particle_steps must be '//
     &          'a divider of nnstat')
              else
                write(*,'(A,I7,A6)') 
     &          ' The Collision particle data will be saved in '//
     &          ' one file at every', nparameters(9), 'steps'
              end if ! (mod(nnstat,nparameters(9)) .ne. 0)
            end if ! (nparameters(9) .eq. 0)

            write(*,*) 'Continuous particles will be added'
            write(*,*) 'Added particles are given by the input file'//
     &                 ' at the fixed position'
     &                ,'but the velocities are interpolated'

cty         Add_particle_steps
            write(*,'(A,I7)') 
     &              'Add_particle_steps = ', nparameters(12)
            if (nparameters(12) .le. 0) then
              call stopnow
     &        ('The Add_particle_steps must be positive')
            else 
              if (mod(nsteps,nparameters(12)) .ne. 0) then
                call stopnow
     &          ('The Add_particle_steps must be '//
     &           'a divider of nsteps')
              else 
                write(*,'(A,I7,A6)') 
     &          ' Particles are added at every',
     &           nparameters(12), 'steps'
              end if ! (mod(nsteps,nparameters(12)) .ne. 0)
            end if ! (nparameters(12) .le. 0)
cty         Add_particle_number
            write(*,'(A,I7)') 
     &              'Add_particle_number = ', nparameters(13)
            if (nparameters(13) .le. 0) then
              call stopnow
     &        ('The Add_particle_number must be positive')
            else 
              write(*,'(A,I7,A20)')
     &        ' Each time,', nparameters(13), 'particles are added'
            end if ! (nparameters(13) .le. 0)

cty         diameter_par
            write(*,'(A,F20.16)') 
     &              'diameter_par = ', rparameters(1)
            if (rparameters(1) .gt. 0) then
              write(*,'(A,F20.16,A2)') 
     &                ' The dimensional diameter of the material of'//
     &                ' the single type of particle is',
     &                rparameters(1), ' m'
            else
              call stopnow('The diameter_par must be positive')
            end if ! (rparameters(1) .gt. 0)

cty         density_par
            write(*,'(A,F20.12)') 
     &              'density_par = ', rparameters(2)
            if (rparameters(2) .gt. 0) then
              write(*,'(A,F20.12,A7)') 
     &                ' The dimensional density of the material of'//
     &                ' the single type of particle is',
     &                rparameters(2), ' kg/m^3'
            else
              call stopnow('The density_par must be positive')
            end if ! (rparameters(2) .gt. 0)           

          end if ! (ioproc)

        else
          call stopnow('file particle_parameter.in does not exist')
        end if ! (lex)

        call MPI_bcast(nparameters,13,MPI_INTEGER,
     &                 ioid,MPI_comm_fluid,ierr) 
        call MPI_bcast(rparameters,12,real_mp_type,
     &                 ioid,MPI_comm_fluid,ierr) 

        Para_par%imode = nparameters(1)
        Para_par%Initial_velocity_particle = nparameters(2)     
        Para_par%Monitor_particle_steps = nparameters(5)
        Para_par%Capture_particle_steps = nparameters(6)
        Para_par%Store_particle_steps = nparameters(7)
        Para_par%Average_particle_steps = nparameters(8)
        Para_par%Collision_particle_steps = nparameters(9)  
        Para_par%Add_particle_steps = nparameters(12)
        Para_par%Add_particle_number = nparameters(13)
        Para_par%diameter_par = rparameters(1)
        Para_par%density_par = rparameters(2)

      end subroutine ! read_particle_parameters



cty========================================================================
cty     subroutine set_non_dimensional_parameters 
cty========================================================================
cty
cty     Description:
cty     Set non-dimensional parameters used in the simulation.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine set_non_dimensional_parameters

        use corevars, only : nzp, thl
        use mppvars, only : ioproc

        implicit none

        real(kind=wp) :: temp1, temp2

        Para_non_dim%diameter_par = Para_par%diameter_par / 
     &                              Para_dim%L_ref
        Para_non_dim%density_par = Para_par%density_par / 
     &                             Para_dim%rho_ref
        if (ioproc) then
          write(*,'(A,F20.16)') 
     &            ' The non-dimensional diameter of the material of'//
     &            ' the single type of particle is', 
     &            Para_non_dim%diameter_par
          write(*,'(A,F25.16)') 
     &            ' The non-dimensional density of the material of'//
     &            ' the single type of particle is',
     &            Para_non_dim%density_par
        end if

      end subroutine ! set_non_dimensional_parameters



cty========================================================================
cty     subroutine mkdir_particles_subfolders_PGI
cty========================================================================
cty
cty     Description:
cty
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine mkdir_particles_subfolders_PGI

        use mppvars, only : MPI_comm_fluid, procid_global, ioproc

        implicit none

        integer :: ierr

        if (procid_global .eq. 0) then

          if (Para_par%Monitor_particle_steps .gt. 0) then
            call mkdir_implementation('Monitor_particles',
     &                                 ierr,.false.)
            if (ierr .ne. 0) then 
              call stopnow('Error mkdir Monitor_particles')
            end if
            call mkdir_implementation('Monitor_particles/out',
     &                                 ierr,.false.)
            if (ierr .ne. 0) then 
              call stopnow('Error mkdir Monitor_particles/out')
            end if
          end if
          
          if (Para_par%Capture_particle_steps .gt. 0) then
            call mkdir_implementation('Capture_particles',
     &                                 ierr,.false.)
            if (ierr .ne. 0) then 
              call stopnow('Error mkdir Capture_particles')
            end if
          end if

          if (Para_par%Store_particle_steps .gt. 0) then
            call mkdir_implementation('Continue_particles',
     &                                 ierr,.false.)
            if (ierr .ne. 0) then 
              call stopnow('Error mkdir Continue_particles')
            end if
          end if

          if (Para_par%Average_particle_steps .gt. 0) then
            call mkdir_implementation('Average_particles',
     &                                 ierr,.false.)
            if (ierr .ne. 0) then 
              call stopnow('Error mkdir Average_particles')
            end if
          end if

          if (Para_par%Collision_particle_steps .gt. 0) then
            call mkdir_implementation('Collision_particles',
     &                                 ierr,.false.)
            if (ierr .ne. 0) then 
              call stopnow('Error mkdir Collision_particles')
            end if
          end if

        end if

        if (ioproc) then
          call mkdir_implementation('Mapping_particles',
     &                               ierr,.false.)
          if (ierr .ne. 0) then 
            call stopnow('Error mkdir Mapping_particles')
          end if
        end if
        call MPI_barrier(MPI_comm_fluid,ierr) 

      end subroutine ! mkdir_particles_subfolders_PGI



cty========================================================================
cty     subroutine allocate_particles_memory
cty========================================================================
cty
cty     Description:
cty     Allocate the memories of particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================
 
      subroutine allocate_particles_memory

        use corevars, only : nxp, nyp, nzp, kmodw, nvar, ngrad_pre, 
     &                       xhalo, yhalo, dt
        use mppvars, only : nprocs

        implicit none
 
        allocate(Pointer_par(Npar_Max))
cty     Arrays used for monitor particles
        allocate(idpar_proc_del(Npar_del_Max))
        idpar_proc_del = 0
cty     Arrays used for exchange particles
        allocate(mark_block_id(Npar_del_Max))
        mark_block_id = 0

        Npar_local_collision = 0
        allocate(pardata_real_collision(Nvar_par_real,
     &                                  Npar_del_Max))
        pardata_real_collision = 0._wp
        allocate(pardata_integer_collision(Nvar_par_integer-1,
     &                                     Npar_del_Max))
        pardata_integer_collision = 0

      end subroutine ! allocate_particles_memory



cty========================================================================
cty     subroutine set_block_coordinates
cty========================================================================
cty
cty     Description:
cty     Set the boundary coordinates of block grids.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty       B_xmin, B_xmax, B_ymin, B_ymax: Boundary of the background grid
cty       B_zmin, B_zmax
cty       O_outer_xmin, O_outer_xmax: boundary of the outer surface of the O-type grid 
cty       O_inner_xmin, O_inner_xmax: boundary of the inner surface of the O-type grid 
cty
cty     Notice: 
cty       1. This subroutine is used to identify whether the particle locates
cty          in the flow field coarsely.
cty
cty========================================================================

      subroutine set_block_coordinates

        use corevars, only : thl
        use mbvars, only : num_blocks
     &                   , ov_bgrid_block, ov_ogrid_block
        use mppvars, only : ioproc
        use global_grid, only : grids
        use overset, only : ov_mesh

        implicit none

        integer :: block_id
        integer :: nx, ny

        integer :: i, j

        do block_id = 1, num_blocks 
          if (ov_bgrid_block(block_id)) then
            nx = grids(block_id)%nxp_block
            ny = grids(block_id)%nyp_block
            if (ioproc) then
              write(*,'(A45,A8,I8,A8,I8)') 
     &                'For the background grid, the dimension is',
     &                'nx = ', nx, 'ny = ', ny
            end if
            B_xmin = minval(grids(block_id)%coords(1:nx,1:ny,1))
            B_xmax = maxval(grids(block_id)%coords(1:nx,1:ny,1))
            B_ymin = minval(grids(block_id)%coords(1:nx,1:ny,2))
            B_ymax = maxval(grids(block_id)%coords(1:nx,1:ny,2))
            B_ymax = B_ymax + (B_ymax-B_ymin)/dble(ny-1)
            B_zmin = 0._wp
            B_zmax = thl
          end if

          if (ov_ogrid_block(block_id)) then
            nx = grids(block_id)%nxp_block
            ny = grids(block_id)%nyp_block
            if (ioproc) then
              write(*,'(A45,A8,I8,A8,I8)') 
     &                'For the O-type grid, the dimension is',
     &                'nx = ', nx, 'ny = ', ny
            end if
            O_outer_xmin = minval(grids(block_id)%coords(1:nx,ny,1))
            O_outer_xmax = maxval(grids(block_id)%coords(1:nx,ny,1))
            O_inner_xmin = minval(grids(block_id)%coords(1:nx,1,1))
            O_inner_xmax = maxval(grids(block_id)%coords(1:nx,1,1))
          end if
        end do
        
        if (ioproc) then
          write(*,'(2(A10,F20.16))') 'B_xmin = ', B_xmin, 
     &                               'B_xmax = ', B_xmax
          write(*,'(2(A10,F20.16))') 'B_ymin = ', B_ymin,
     &                               'B_ymax = ', B_ymax
          write(*,'(2(A10,F20.16))') 'B_zmin = ', B_zmin,
     &                               'B_zmax = ', B_zmax
          write(*,'(2(A16,F20.16))') 'O_outer_xmin = ', O_outer_xmin,
     &                               'O_outer_xmax = ', O_outer_xmax
          write(*,'(2(A16,F20.16))') 'O_inner_xmin = ', O_inner_xmin,
     &                               'O_inner_xmax = ', O_inner_xmax
        end if 

        do block_id = 1, num_blocks 
          if (ov_ogrid_block(block_id)) then 

            nx = grids(block_id)%nxp_block
            ny = grids(block_id)%nyp_block

            allocate(O_inner_px(nx))
            allocate(O_inner_py(nx))
            allocate(O_inner_py_trans(nx))
            ! ym
            j = 1
            do i = 1, nx
              O_inner_px(i) = grids(block_id)%coords(i,j,1)
              O_inner_py(i) = grids(block_id)%coords(i,j,2)
              O_inner_py_trans(i) = ov_mesh(block_id)%y_trans(i,j)
            end do

            allocate(O_outer_px(nx))
            allocate(O_outer_py(nx))
            allocate(O_outer_py_trans(nx))
            ! yp
            j = ny
            do i = 1, nx
              O_outer_px(i) = grids(block_id)%coords(i,j,1)
              O_outer_py(i) = grids(block_id)%coords(i,j,2)
              O_outer_py_trans(i) = ov_mesh(block_id)%y_trans(i,j)
            end do

          end if
        end do

      end subroutine ! set_block_coordinates



cty========================================================================
cty     subroutine init_glob_index
cty========================================================================
cty
cty     Description:
cty     Initialize arrays used for the mapping between blocks and processors.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty       1. For each processor:
cty          (1) block_proc(:): block_id
cty       2. For each processor, each block:
cty          (1) ind_i_proc(:): the i index in Cartesian processor network 
cty              ind_j_proc(:): the j index in Cartesian processor network 
cty          (2) ibeg_proc(:): the start i index of the global grid 
cty              iend_proc(:): the end i index of the global grid 
cty              jbeg_proc(:): the start j index of the global grid 
cty              jend_proc(:): the end j index of the global grid 
cty
cty     Notice: 
cty
cty========================================================================

      subroutine init_glob_index

        use mppvars, only : procid_fluid, nprocs, ioproc
        use mbvars, only : num_blocks
        use corevars, only : nxp, nyp, nzp
        use global_grid, only : grids

        implicit none

        integer :: block_id, p
        integer :: i, j

        allocate(ind_i_proc(0:nprocs-1))
        allocate(ibeg_proc(0:nprocs-1))
        allocate(iend_proc(0:nprocs-1))
        allocate(ind_j_proc(0:nprocs-1))
        allocate(jbeg_proc(0:nprocs-1))
        allocate(jend_proc(0:nprocs-1))
        allocate(block_proc(0:nprocs-1))

        do block_id = 1, num_blocks 
          do j = 1, grids(block_id)%nproc_y
            do i = 1, grids(block_id)%nproc_x
              
              p = grids(block_id)%proc_map(i,j)
              block_proc(p) = block_id
              ind_i_proc(p) = i
              ind_j_proc(p) = j

              ibeg_proc(p) = grids(block_id)%ibeg(i)
              if (i .eq. grids(block_id)%nproc_x) then
                iend_proc(p) = grids(block_id)%nxp_block
              else
                iend_proc(p) = grids(block_id)%ibeg(i+1)-1
              end if
    
              jbeg_proc(p) = grids(block_id)%jbeg(j)
              if (j .eq. grids(block_id)%nproc_y) then
                jend_proc(p) = grids(block_id)%nyp_block
              else
                jend_proc(p) = grids(block_id)%jbeg(j+1)-1
              end if
  
            end do ! i
          end do ! j       
        end do ! block_id

        if (ioproc) then
          do p = 0, nprocs-1
            write(*,'(A,I4)') 'For the processor', p
            write(*,'(A15,I3,2(A12,I3),4(A11,I5))') 
     &              'block id = ', block_proc(p),
     &              'ind_i = ', ind_i_proc(p), 
     &              'ind_j = ', ind_j_proc(p), 
     &              'ibeg = ', ibeg_proc(p), 
     &              'iend = ', iend_proc(p), 
     &              'jbeg = ', jbeg_proc(p), 
     &              'jend = ', jend_proc(p)
          end do
        end if

      end subroutine ! init_glob_index



cty========================================================================
cty     subroutine init_grids_nprocs
cty========================================================================
cty
cty     Description:
cty     Initialize the number of processors corresponding to the background
cty     and O-type grids.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty       B_nprocs: the number of processors corresponding to background grid
cty       O_nprocs: the number of processors corresponding to O-type grid
cty 
cty     Notice: 
cty       Global processors id: 
cty       1. background grid: 0:B_nprocs-1
cty       2. O-type grid: B_nprocs:nprocs-1
cty
cty========================================================================

      subroutine init_grids_nprocs

        use mppvars, only : nprocs, ioproc
        use mbvars, only : num_blocks
     &                   , ov_bgrid_block, ov_ogrid_block
        use global_grid, only : grids

        implicit none

        integer :: block_id

        do block_id = 1, num_blocks
          if (ov_bgrid_block(block_id)) then 
            B_nprocs = grids(block_id)%nproc_x * grids(block_id)%nproc_y
          end if

          if (ov_ogrid_block(block_id)) then
            O_nprocs = grids(block_id)%nproc_x * grids(block_id)%nproc_y
          end if
        end do ! block_id

        if (ioproc) then
          write(*,'(A12,I5)') 'B_nprocs = ', B_nprocs
          write(*,'(A12,I5)') 'O_nprocs = ', O_nprocs 
          write(*,'(A12,I5)') 'nprocs = ', nprocs 
        end if

      end subroutine ! init_grids_nprocs



cty========================================================================
cty     subroutine set_subblock_coordinates
cty========================================================================
cty
cty     Description:
cty     Set the boundary coordinates of the subblock grids.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine set_subblock_coordinates

        use mbvars, only : ov_ogrid_block, ov_bgrid_block
        use mppvars, only : procid_fluid, nprocs, ioproc
        use global_grid, only : grids
        use overset, only : ov_mesh

        implicit none

        integer :: p
        integer :: block_id, ind_i, ind_j
        integer :: sx, sy
        integer :: nx, ny
        integer :: i, j
        integer :: l
        character(3) :: idproc_name 
        character(300) :: particle_file  

        allocate(subblock_proc(0:nprocs-1))

        do p = 0, nprocs-1

cty       Set parameters associated with the present processor: block_id, ind_i, ind_j
          block_id = block_proc(p)
          ind_i = ind_i_proc(p)
          ind_j = ind_j_proc(p)

cty       Set sx, sy 
          sx = grids(block_id)%ibeg(ind_i)
          sy = grids(block_id)%jbeg(ind_j)

cty       Here nx and ny save the number of points in each direction

cty       Set nx
          if (ind_i .lt. grids(block_id)%nproc_x) then
            nx = grids(block_id)%ibeg(ind_i+1) -sx +1
          else if (ind_i .eq. grids(block_id)%nproc_x) then
cty         This layer need to be considered spectially
            if (ov_bgrid_block(block_id)) then 
              nx = grids(block_id)%nxp_block -sx +1
            end if

            if (ov_ogrid_block(block_id)) then 
              nx = (grids(block_id)%nxp_block -sx +1) + 1
            end if ! (ov_bgrid_block(block_id))
          end if ! (ind_i .lt. grids(block_id)%nproc_x)

cty       Set ny
          if (ind_j .lt. grids(block_id)%nproc_y) then
            ny = grids(block_id)%jbeg(ind_j+1) -sy +1
          else if (ind_j .eq. grids(block_id)%nproc_y) then
cty         This layer need to be considered spectially
            if (ov_bgrid_block(block_id)) then 
              ny = (grids(block_id)%nyp_block -sy +1) + 1
            end if

            if (ov_ogrid_block(block_id)) then 
              ny = grids(block_id)%nyp_block -sy +1
            end if ! (ov_bgrid_block(block_id))
          end if ! (ind_j .lt. grids(block_id)%nproc_y)

          subblock_proc(p)%xsize = nx
          subblock_proc(p)%ysize = ny

          allocate(subblock_proc(p)%coords(nx,ny,2))

          do j = 1, ny
            do i = 1, nx
              subblock_proc(p)%coords(i,j,1) = 
     &                         grids(block_id)%coords(sx+i-1,sy+j-1,1)
              subblock_proc(p)%coords(i,j,2) = 
     &                         grids(block_id)%coords(sx+i-1,sy+j-1,2)
            end do ! i = 1, nx
          end do ! j = 1, ny 

          if (ov_bgrid_block(block_id)) then

            if (ind_j .eq. grids(block_id)%nproc_y) then
              subblock_proc(p)%coords(:,ny,1) = 
     &        subblock_proc(p)%coords(:,ny-1,1)
              subblock_proc(p)%coords(:,ny,2) = 
     &        subblock_proc(p)%coords(:,ny-1,2) +
     &       (B_ymax-B_ymin)/dble(grids(block_id)%nyp_block)
            end if ! (ind_j .eq. grids(block_id)%nproc_y)

          end if !(ov_bgrid_block(block_id))

          if (ov_ogrid_block(block_id)) then

            if (ind_i .eq. grids(block_id)%nproc_x) then
              do j = 1, ny
                subblock_proc(p)%coords(nx,j,1) = 
     &                           grids(block_id)%coords(1,sy+j-1,1)
                subblock_proc(p)%coords(nx,j,2) =
     &                           grids(block_id)%coords(1,sy+j-1,2)
              end do ! j = 1, ny 
            end if ! (ind_i .eq. grids(block_id)%nproc_x)

          end if !(ov_ogrid_block(block_id))

          if (ov_bgrid_block(block_id)) then

            subblock_proc(p)%xmin = subblock_proc(p)%coords(1,ny,1) 
            subblock_proc(p)%xmax = subblock_proc(p)%coords(nx,1,1)
            subblock_proc(p)%ymin = subblock_proc(p)%coords(nx,1,2)
            subblock_proc(p)%ymax = subblock_proc(p)%coords(1,ny,2)

          end if !(ov_bgrid_block(block_id))

          if (ov_ogrid_block(block_id)) then

            allocate(subblock_proc(p)%px(2*(nx+ny)-4))
            allocate(subblock_proc(p)%py(2*(nx+ny)-4))
            allocate(subblock_proc(p)%py_trans(2*(nx+ny)-4))
        
            l = 0

            ! ym
            j = sy
            do i = sx, sx+nx-2
              l = l + 1
              subblock_proc(p)%px(l) = grids(block_id)%coords(i,j,1)
              subblock_proc(p)%py(l) = grids(block_id)%coords(i,j,2)
              subblock_proc(p)%py_trans(l) = 
     &                                 ov_mesh(block_id)%y_trans(i,j)
            end do

            ! xp
            i = sx+nx-1
            i = mod(i,grids(block_id)%nxp_block)
            do j = sy, sy+ny-1
              l = l + 1
              subblock_proc(p)%px(l) = grids(block_id)%coords(i,j,1)
              subblock_proc(p)%py(l) = grids(block_id)%coords(i,j,2)
              subblock_proc(p)%py_trans(l) = 
     &                                 ov_mesh(block_id)%y_trans(i,j)
            end do

            ! yp
            j = sy+ny-1
            do i = sx+nx-2, sx, -1
              l = l + 1
              subblock_proc(p)%px(l) = grids(block_id)%coords(i,j,1)
              subblock_proc(p)%py(l) = grids(block_id)%coords(i,j,2)
              subblock_proc(p)%py_trans(l) = 
     &                                 ov_mesh(block_id)%y_trans(i,j)
            end do

            ! xm
            i = sx
            do j = sy+ny-2, sy+1, -1
              l = l + 1
              subblock_proc(p)%px(l) = grids(block_id)%coords(i,j,1)
              subblock_proc(p)%py(l) = grids(block_id)%coords(i,j,2)
              subblock_proc(p)%py_trans(l) = 
     &                                 ov_mesh(block_id)%y_trans(i,j)
            end do
 
          end if !(ov_ogrid_block(block_id))

        end do ! p

        nx = subblock_proc(procid_fluid)%xsize
        ny = subblock_proc(procid_fluid)%ysize

cty     Set the file name of the particle_file
        write(idproc_name,'(I2)') procid_fluid
cty     Here we use 2 places to mark the id of the processors, 
cty     and fill the blanks with zeros.
        do i = 1, len_trim(idproc_name)
          if (idproc_name(i:i) .eq. '') then
            idproc_name(i:i) = '0' 
          else 
            exit
          end if
        end do 
        particle_file = 'Mapping_particles/'//
     &                  'Subblock_coordinates'//
     &                  trim(idproc_name)//'.dat'
        particle_file = trim(adjustl(particle_file))

        open(unit=procid_fluid+1000,file=particle_file,
     &       form='formatted',status='replace',
     &       action='write')

        write(procid_fluid+1000,*) nx, ny

        do j = 1, ny 
          do i = 1, nx 
            write(procid_fluid+1000,*) 
     &            subblock_proc(procid_fluid)%coords(i,j,1),
     &            subblock_proc(procid_fluid)%coords(i,j,2)
          end do
        end do

        close(procid_fluid+1000) 

      end subroutine ! set_subblock_coordinates



cty========================================================================
cty     subroutine set_B_subblock_proc_coordinates
cty========================================================================
cty
cty     Description:
cty     Set the boundary coordinates of the subblock grids of the processor.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine set_B_subblock_proc_coordinates

        use corevars, only : nzp, theta, thl
        use mppvars, only : procid_fluid, ioproc

        implicit none

        integer :: nx, ny
        integer :: i, j, k

        nx = subblock_proc(procid_fluid)%xsize
        ny = subblock_proc(procid_fluid)%ysize

        allocate(B_z_x(nx))
        allocate(B_r_y(ny))

        allocate(B_theta(nzp+1))

        if (block_proc(procid_fluid) .eq. 1) then
          do i = 1, nx
            B_z_x(i) = subblock_proc(procid_fluid)%coords(i,1,1)
          end do ! i = 1, nx
          do j = 1, ny
            B_r_y(j) = subblock_proc(procid_fluid)%coords(1,j,2)
          end do ! j = 1, ny
        end if ! (block_proc(procid_fluid) .eq. 1)

        do k = 1, nzp
          B_theta(k) = theta(k)
        end do
        B_theta(nzp+1) = thl

      end subroutine ! set_B_subblock_proc_coordinates



cty========================================================================
cty     subroutine set_mapping_indexes
cty========================================================================
cty
cty     Description:
cty     Set the mapping indexes for the H- and O-type blocks to accelerate
cty     the procedure of locating the particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty 
cty     Notice: 
cty       1. The present mapping relationships are only valid for the outer 
cty          surface of the O-type block.
cty
cty========================================================================

      subroutine set_mapping_indexes

        use mppvars, only : nprocs, ioproc

        implicit none

        integer :: p
        integer :: nx, ny
        integer :: i
        character(3) :: idproc_name 
        character(300) :: particle_file        

        logical, allocatable :: lex(:)


        allocate(mapping(0:nprocs-1))
        allocate(lex(0:nprocs-1))

        do p = 0, nprocs-1

          nx = subblock_proc(p)%xsize 
          ny = subblock_proc(p)%ysize 

          allocate(mapping(p)%nearest_proc(nx,ny))
          allocate(mapping(p)%nearest_i(nx,ny))
          allocate(mapping(p)%nearest_j(nx,ny))

cty       Set the file name of the particle_file
          write(idproc_name,'(I2)') p
cty       Here we use 2 places to mark the id of the processors, 
cty       and fill the blanks with zeros.
          do i = 1, len_trim(idproc_name)
            if (idproc_name(i:i) .eq. '') then
              idproc_name(i:i) = '0' 
            else 
              exit
            end if
          end do 
          particle_file = 'Mapping_particles/'//
     &                    'Mapping_indexes'//
     &                    trim(idproc_name)//'.dat'
          particle_file = trim(adjustl(particle_file))

          inquire(file=particle_file,exist=lex(p))

        end do ! p

        if (all(lex(0:nprocs-1) .eq. .true.)) then
          call read_mapping_indexes
        else
          call init_mapping_indexes
          call write_mapping_indexes
        end if
        deallocate(lex)

cty     Check the mapping relationships
        call check_mapping_indexes

      end subroutine ! set_mapping_indexes



cty========================================================================
cty     subroutine read_mapping_indexes
cty========================================================================
cty
cty     Description:
cty     Read the mapping indexes.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty 
cty     Notice: 
cty
cty========================================================================

      subroutine read_mapping_indexes

        use mppvars, only : MPI_INTEGER, MPI_comm_fluid, ioid, nprocs, 
     &                      ioproc

        implicit none

        integer :: p 
        integer :: nx, ny
        integer :: i, j
        character(3) :: idproc_name 
        character(300) :: particle_file        

        integer :: ierr

        if (ioproc) then

          do p = 0, nprocs-1

cty         Set the file name of the particle_file
            write(idproc_name,'(I2)') p
cty         Here we use 2 places to mark the id of the processors, 
cty         and fill the blanks with zeros.
            do i = 1, len_trim(idproc_name)
              if (idproc_name(i:i) .eq. '') then
                idproc_name(i:i) = '0' 
              else 
                exit
              end if
            end do 
            particle_file = 'Mapping_particles/'//
     &                      'Mapping_indexes'//
     &                      trim(idproc_name)//'.dat'
            particle_file = trim(adjustl(particle_file))

            open(unit=99,file=particle_file,
     &           form='formatted',status='old',action='read')

            read(99,*) nx, ny
            if ( (nx .ne. subblock_proc(p)%xsize) .or.
     &           (ny .ne. subblock_proc(p)%ysize) ) then
              call stopnow('Dimensions of mapping relationships are '//
     &                     'not matching for subblock')
            end if

            do j = 1, ny 
              do i = 1, nx 
                read(99,*) mapping(p)%nearest_proc(i,j),
     &                     mapping(p)%nearest_i(i,j),
     &                     mapping(p)%nearest_j(i,j)
              end do
            end do

            close(99) 

          end do ! p

          do p = 0, nprocs-1
            nx = subblock_proc(p)%xsize
            ny = subblock_proc(p)%ysize
            call MPI_bcast(mapping(p)%nearest_proc,nx*ny,
     &                     MPI_INTEGER,ioid,MPI_comm_fluid,ierr)
            call MPI_bcast(mapping(p)%nearest_i,nx*ny,
     &                     MPI_INTEGER,ioid,MPI_comm_fluid,ierr)
            call MPI_bcast(mapping(p)%nearest_j,nx*ny,
     &                     MPI_INTEGER,ioid,MPI_comm_fluid,ierr)
          end do ! p

        end if

      end subroutine ! read_mapping_indexes



cty========================================================================
cty     subroutine init_mapping_indexes
cty========================================================================
cty
cty     Description:
cty     Initialize the mapping indexes for the H- and O-type blocks to accelerate
cty     the procedure of locating the particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty 
cty     Notice: 
cty       1. The present mapping relationships are only valid for the outer 
cty          surface of the O-type block.
cty
cty========================================================================

      subroutine init_mapping_indexes

        use mppvars, only : nprocs, ioproc, nullproc
        use global_grid, only : grids
        use overset, only : ov_mesh, ov_pitch, trans_mode
 
        implicit none

        real(kind=wp) :: dis_threshold
        integer :: p 
        integer :: nx_B, ny_B
        integer :: nx_O, ny_O
        integer :: nx, ny
        integer :: i, j
        integer :: nearest_i, nearest_j
        integer :: k
        integer :: l, m
        real(kind=wp) :: xc, yc, yc_trans
        real(kind=wp) :: xc_B, yc_B 
        real(kind=wp) :: xc_O, yc_O, yc_O_trans
        real(kind=wp) :: dis, dis_trans, dis_min
        real(kind=wp), allocatable :: dis_local(:)

        dis_threshold = 3._wp*(B_ymax-B_ymin)/dble(grids(1)%nyp_block-1)

        nx_B = grids(1)%nxp_block
        ny_B = grids(1)%nyp_block
        nx_O = grids(2)%nxp_block
        ny_O = grids(2)%nyp_block

        do p = 0, nprocs-1

          nx = subblock_proc(p)%xsize 
          ny = subblock_proc(p)%ysize 

cty       Inatialize the nearest_proc, nearest_i and nearest_j
          mapping(p)%nearest_proc = -999
          mapping(p)%nearest_i = 0
          mapping(p)%nearest_j = 0
cty       Only the values i=1,nx-1, j=1,ny-1 are valid
          mapping(p)%nearest_proc(nx,:) = nullproc
          mapping(p)%nearest_proc(:,ny) = nullproc
          mapping(p)%nearest_i(nx,:) = -999
          mapping(p)%nearest_i(:,ny) = -999
          mapping(p)%nearest_j(nx,:) = -999
          mapping(p)%nearest_j(:,ny) = -999

        end do ! p = 0, nprocs-1

cty     First find the non-mapping cells in the outer surface
        do p = 0, B_nprocs-1

          nx = subblock_proc(p)%xsize 
          ny = subblock_proc(p)%ysize 

          do j = 1, ny-1
            do i = 1, nx-1

cty           Four ponits of the cell are located in the outer surface
              if ( particle_in_outer_surface(
     &               subblock_proc(p)%coords(i  ,j  ,1), 
     &               subblock_proc(p)%coords(i  ,j  ,2)) .and.
     &             particle_in_outer_surface(
     &               subblock_proc(p)%coords(i+1,j  ,1),
     &               subblock_proc(p)%coords(i+1,j  ,2)) .and.
     &             particle_in_outer_surface(
     &               subblock_proc(p)%coords(i+1,j+1,1),
     &               subblock_proc(p)%coords(i+1,j+1,2)) .and.
     &             particle_in_outer_surface(
     &               subblock_proc(p)%coords(i  ,j+1,1),
     &               subblock_proc(p)%coords(i  ,j+1,2)) ) then
cty             The non-mapping cells are marked by the nearest_proc = nullproc
cty             nearest_i = -999, nearest_j = -999
                mapping(p)%nearest_proc(i,j) = nullproc
                mapping(p)%nearest_i(i,j) = -999
                mapping(p)%nearest_j(i,j) = -999
              end if

            end do ! i
          end do ! j

        end do ! p = 0, B_nprocs-1

cty     Then further identify the mapping and non-mapping cells, which are 
cty     marked by the nearest_proc = -999 now
        do p = 0, B_nprocs-1

          allocate(dis_local(nx_O))
          nx = subblock_proc(p)%xsize 
          ny = subblock_proc(p)%ysize 

          do j = 1, ny-1
            do i = 1, nx-1
              if (mapping(p)%nearest_proc(i,j) .eq. -999) then

cty             Firstly, we need to calculate the dis_local
                xc = 0.25_wp*(subblock_proc(p)%coords(i  ,j  ,1)+
     &                        subblock_proc(p)%coords(i+1,j  ,1)+
     &                        subblock_proc(p)%coords(i+1,j+1,1)+
     &                        subblock_proc(p)%coords(i  ,j+1,1))
                yc = 0.25_wp*(subblock_proc(p)%coords(i  ,j  ,2)+
     &                        subblock_proc(p)%coords(i+1,j  ,2)+
     &                        subblock_proc(p)%coords(i+1,j+1,2)+
     &                        subblock_proc(p)%coords(i  ,j+1,2))
                do k = 1, nx_O-1
                  xc_O = 0.25_wp*(grids(2)%coords(k  ,ny_O-1,1)+
     &                            grids(2)%coords(k+1,ny_O-1,1)+
     &                            grids(2)%coords(k+1,ny_O  ,1)+
     &                            grids(2)%coords(k  ,ny_O  ,1))
                  yc_O = 0.25_wp*(grids(2)%coords(k  ,ny_O-1,2)+
     &                            grids(2)%coords(k+1,ny_O-1,2)+
     &                            grids(2)%coords(k+1,ny_O  ,2)+
     &                            grids(2)%coords(k  ,ny_O  ,2))
                  yc_O_trans = yc_O - trans_mode*ov_pitch
                  dis = dsqrt( (xc-xc_O)**2._wp + 
     &                         (yc-yc_O)**2._wp )
                  dis_trans = dsqrt( (xc-xc_O)**2._wp + 
     &                               (yc-yc_O_trans)**2._wp )
                  dis_local(k) = min(dis,dis_trans)
                end do ! k 
                xc_O = 0.25_wp*(grids(2)%coords(nx_O,ny_O-1,1)+
     &                          grids(2)%coords(1   ,ny_O-1,1)+
     &                          grids(2)%coords(1   ,ny_O  ,1)+
     &                          grids(2)%coords(nx_O,ny_O  ,1))
                yc_O = 0.25_wp*(grids(2)%coords(nx_O,ny_O-1,2)+
     &                          grids(2)%coords(1   ,ny_O-1,2)+
     &                          grids(2)%coords(1   ,ny_O  ,2)+
     &                          grids(2)%coords(nx_O,ny_O  ,2))
                yc_O_trans = yc_O - trans_mode*ov_pitch
                dis = dsqrt( (xc-xc_O)**2._wp + 
     &                       (yc-yc_O)**2._wp )
                dis_trans = dsqrt( (xc-xc_O)**2._wp + 
     &                             (yc-yc_O_trans)**2._wp )
                dis_local(nx_O) = min(dis,dis_trans)
cty             Now, we have obtained the dis_local                

cty             Use dis_local to make a further classfication
                nearest_i = minloc(dis_local,dim=1) 
                if (dis_local(nearest_i) .le. dis_threshold) then

                  do k = B_nprocs, nprocs-1
                    if (ind_j_proc(k) .eq. grids(2)%nproc_y) then
                      if ( (nearest_i .ge. ibeg_proc(k)) .and.
     &                     (nearest_i .le. iend_proc(k)) ) then
                        mapping(p)%nearest_proc(i,j) = k
                        mapping(p)%nearest_i(i,j) = 
     &                                    nearest_i-ibeg_proc(k)+1 
                        mapping(p)%nearest_j(i,j) = jend_proc(k)-1
                        exit
                      end if
                    end if
                  end do ! k

                else

                  mapping(p)%nearest_proc(i,j) = nullproc 
                  mapping(p)%nearest_i(i,j) = -999
                  mapping(p)%nearest_j(i,j) = -999

                end if ! (dis_local(nearest_i) .le. dis_threshold)

              end if ! (mapping(p)%nearest_proc(i,j) .eq. -999)
            end do ! i
          end do ! j

          deallocate(dis_local)

        end do ! p = 0, B_nprocs-1

cty     For the O-type block, we only consider the outer surface
        do p = B_nprocs, nprocs-1

          nx = subblock_proc(p)%xsize 
          ny = subblock_proc(p)%ysize 

          if (ind_j_proc(p) .eq. grids(2)%nproc_y) then

            do i = 1, nx-1
              xc_O = 0.25_wp*(subblock_proc(p)%coords(i  ,ny-1,1)+
     &                        subblock_proc(p)%coords(i+1,ny-1,1)+
     &                        subblock_proc(p)%coords(i+1,ny  ,1)+
     &                        subblock_proc(p)%coords(i  ,ny  ,1))
              yc_O = 0.25_wp*(subblock_proc(p)%coords(i  ,ny-1,2)+
     &                        subblock_proc(p)%coords(i+1,ny-1,2)+
     &                        subblock_proc(p)%coords(i+1,ny  ,2)+
     &                        subblock_proc(p)%coords(i  ,ny  ,2))
cty           mapping%nearest_proc(i,j) must not be equal to nullproc for j = ny-1
cty           Thus we remain them equals to -999, which are handled in the next procedure
              do j = 1, ny-2
                xc = 0.25_wp*(subblock_proc(p)%coords(i  ,j  ,1)+
     &                        subblock_proc(p)%coords(i+1,j  ,1)+
     &                        subblock_proc(p)%coords(i+1,j+1,1)+
     &                        subblock_proc(p)%coords(i  ,j+1,1))
                yc = 0.25_wp*(subblock_proc(p)%coords(i  ,j  ,2)+
     &                        subblock_proc(p)%coords(i+1,j  ,2)+
     &                        subblock_proc(p)%coords(i+1,j+1,2)+
     &                        subblock_proc(p)%coords(i  ,j+1,2))
                dis = dsqrt( (xc-xc_O)**2._wp + 
     &                       (yc-yc_O)**2._wp )
                if (dis .gt. dis_threshold) then
                  mapping(p)%nearest_proc(i,j) = nullproc 
                  mapping(p)%nearest_i(i,j) = -999
                  mapping(p)%nearest_j(i,j) = -999
                end if
              end do ! j
            end do ! i

          else

            mapping(p)%nearest_proc = nullproc 
            mapping(p)%nearest_i = -999
            mapping(p)%nearest_j = -999

          end if

        end do ! p = B_nprocs, nprocs-1

        do p = B_nprocs, nprocs-1

          nx = subblock_proc(p)%xsize 
          ny = subblock_proc(p)%ysize 

          do j = 1, ny-1
            do i = 1, nx-1

cty           The mapping cells are marked by the nearest_proc = -999,
cty           and we need to find its mapping processor 
              if (mapping(p)%nearest_proc(i,j) .eq. -999) then

                xc = 0.25_wp*(subblock_proc(p)%coords(i  ,j  ,1)+
     &                        subblock_proc(p)%coords(i+1,j  ,1)+
     &                        subblock_proc(p)%coords(i+1,j+1,1)+
     &                        subblock_proc(p)%coords(i  ,j+1,1))
                yc = 0.25_wp*(subblock_proc(p)%coords(i  ,j  ,2)+
     &                        subblock_proc(p)%coords(i+1,j  ,2)+
     &                        subblock_proc(p)%coords(i+1,j+1,2)+
     &                        subblock_proc(p)%coords(i  ,j+1,2))
                yc_trans = yc - trans_mode*ov_pitch
cty             Set a large value for dis_min at first
                dis_min = B_ymax - B_ymin 

                do k = 0, B_nprocs-1

                  do m = 1, subblock_proc(k)%ysize-1
                    do l = 1, subblock_proc(k)%xsize-1
                      if (mapping(k)%nearest_proc(l,m) .ne. nullproc) 
     &                  then 

                        xc_B = 0.25_wp*(
     &                         subblock_proc(k)%coords(l  ,m  ,1)+
     &                         subblock_proc(k)%coords(l+1,m  ,1)+
     &                         subblock_proc(k)%coords(l+1,m+1,1)+
     &                         subblock_proc(k)%coords(l  ,m+1,1))
                        yc_B = 0.25_wp*(
     &                         subblock_proc(k)%coords(l  ,m  ,2)+
     &                         subblock_proc(k)%coords(l+1,m  ,2)+
     &                         subblock_proc(k)%coords(l+1,m+1,2)+
     &                         subblock_proc(k)%coords(l  ,m+1,2))
                        dis = dsqrt( (xc_B-xc)**2._wp + 
     &                               (yc_B-yc)**2._wp )
                        dis_trans = dsqrt( (xc_B-xc)**2._wp + 
     &                                     (yc_B-yc_trans)**2._wp )
                        if (min(dis,dis_trans) .le. dis_min) then
                          dis_min = min(dis,dis_trans)
                          mapping(p)%nearest_proc(i,j) = k 
                          mapping(p)%nearest_i(i,j) = l
                          mapping(p)%nearest_j(i,j) = m
                        end if
cty                     Theoretically, we can must find such a minimum?

                      end if ! (mapping(k)%nearest_proc(l,m) .ne. nullproc)
                    end do ! l
                  end do ! m

                end do ! k

              end if ! (mapping(p)%nearest_proc(i,j) .eq. -999)

            end do ! i
          end do ! j

        end do ! p = B_nprocs, nprocs-1
 
      end subroutine ! init_mapping_indexes



cty========================================================================
cty     subroutine write_mapping_indexes
cty========================================================================
cty
cty     Description:
cty     Write the mapping indexes.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty 
cty     Notice: 
cty
cty========================================================================

      subroutine write_mapping_indexes

        use mppvars, only : nprocs, ioproc

        implicit none

        integer :: p 
        integer :: nx, ny
        integer :: i, j
        character(3) :: idproc_name 
        character(300) :: particle_file        

        if (ioproc) then

          do p = 0, nprocs-1
            nx = subblock_proc(p)%xsize 
            ny = subblock_proc(p)%ysize 

cty         Set the file name of the particle_file
            write(idproc_name,'(I2)') p
cty         Here we use 2 places to mark the id of the processors, 
cty         and fill the blanks with zeros.
            do i = 1, len_trim(idproc_name)
              if (idproc_name(i:i) .eq. '') then
                idproc_name(i:i) = '0' 
              else 
                exit
              end if
            end do 
            particle_file = 'Mapping_particles/'//
     &                      'Mapping_indexes'//
     &                      trim(idproc_name)//'.dat'
            particle_file = trim(adjustl(particle_file))

            open(unit=99,file=particle_file,
     &           form='formatted',status='replace',action='write')

            write(99,*) nx, ny

            do j = 1, ny 
              do i = 1, nx 
                write(99,*) mapping(p)%nearest_proc(i,j),
     &                      mapping(p)%nearest_i(i,j),
     &                      mapping(p)%nearest_j(i,j)
              end do
            end do

            close(99) 
          end do ! p

        end if

      end subroutine ! write_mapping_indexes



cty========================================================================
cty     subroutine check_mapping_indexes
cty========================================================================
cty
cty     Description:
cty     Check the mapping indexes.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty 
cty     Notice: 
cty
cty========================================================================

      subroutine check_mapping_indexes

        use mppvars, only : nprocs

        implicit none

        integer :: p 
        integer :: nx, ny 
        integer :: i, j

        do p = 0, nprocs-1

          nx = subblock_proc(p)%xsize 
          ny = subblock_proc(p)%ysize 

          do j = 1, ny
            do i = 1, nx
              if ( (mapping(p)%nearest_proc(i,j) .eq. -999) .or.
     &             (mapping(p)%nearest_i(i,j) .eq. 0) .or.
     &             (mapping(p)%nearest_j(i,j) .eq. 0) ) then
                call stopnow('There are somthings wrong with the '//
     &                       'mapping relationships ')
              end if
            end do ! i
          end do ! j

        end do ! p

      end subroutine ! check_mapping_indexes



cty========================================================================
cty     subroutine init_surface_vector
cty========================================================================
cty
cty     Description:
cty     For the Hertz-Mindlin no-slip model, we need to calculate the blade 
cty     surface vector for the subsequent calculation.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty 
cty     Notice: 
cty
cty========================================================================

      subroutine init_surface_vector

        use corevars, only : nxp, z_x, r_y
        use mbvars, only : num_blocks
     &                   , ov_ogrid_block, ov_bgrid_block
        use mppvars, only : procid_fluid
        use global_grid, only : grids

        implicit none

        integer :: block_id
        integer :: i, nx
        real(kind=wp) :: length


        do block_id = 1, num_blocks
          if (ov_ogrid_block(block_id)) then

            nx = grids(block_id)%nxp_block
            allocate(tangen_x_global(nx)) 
            allocate(tangen_y_global(nx)) 
            allocate(normal_x_global(nx)) 
            allocate(normal_y_global(nx)) 

cty         Calculate the normal and tangential vectors             
            do i = 1, nx-1
              if (i .lt. nx) then
                tangen_x_global(i) = grids(block_id)%coords(i+1,1,1)
     &                             - grids(block_id)%coords(i,1,1)
                tangen_y_global(i) = grids(block_id)%coords(i+1,1,2)
     &                             - grids(block_id)%coords(i,1,2)
              else if (i .eq. nx) then
                tangen_x_global(i) = grids(block_id)%coords(1,1,1)
     &                             - grids(block_id)%coords(nx,1,1)
                tangen_y_global(i) = grids(block_id)%coords(1,1,2)
     &                             - grids(block_id)%coords(nx,1,2)
              end if
              length = dsqrt(tangen_x_global(i)**2._wp + 
     &                       tangen_y_global(i)**2._wp)
              tangen_x_global(i) = tangen_x_global(i) / length
              tangen_y_global(i) = tangen_y_global(i) / length
              normal_x_global(i) = -tangen_y_global(i)
              normal_y_global(i) = tangen_x_global(i)
            end do ! i

          end if
        end do

      end subroutine ! init_surface_vector



cty========================================================================
cty     subroutine read_particles_data_init
cty========================================================================
cty
cty     Description:
cty     Read the data of particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty       1. Parameters named as Npar*
cty          (1) Npar_local: 
cty              the number of particles in every processor at each time step
cty          (2) Npar_global:
cty              the sum of Npar_local at each time step
cty          (3) Npar_init:
cty              the initial total particle number given by PARTICLES_INPUT.dat
cty          (4) Npar_init_inside: 
cty              the initial total particle number inside the background grid
cty          (5) Npar_read: 
cty              the number of particles been read at each read times
cty          (6) Npar_read_onetime: 
cty              the max number of particles been read at each read times
cty          (7) Npar_local_proc(procid_fluid): 
cty              saves Npar_local of each processor in ioid (only used to print)
cty
cty========================================================================

      subroutine read_particles_data_init

        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_SUM,
     &                      MPI_comm_fluid, ioid, procid_fluid, 
     &                      nprocs, ioproc
        use mbvars, only : num_blocks
     &                   , ov_ogrid_block, ov_bgrid_block
        use global_grid, only : grids
        use overset, only : ov_mesh
        use corevars, only : z_x, r_y
 
        implicit none 

        integer :: Npar_init, Npar_init_inside 
        integer, allocatable, dimension(:) :: Npar_local_proc
        integer :: Npar_read 
        integer :: nk, read_times
        integer :: j, k, ierr

        real(kind=wp), allocatable, dimension(:,:) :: pardata
        real(kind=wp) :: transpar
        real(kind=wp) :: xpar, ypar
        real(kind=wp) :: zpar
        integer :: block_id
     &           , which_block
        logical :: is_in

        logical :: lex

        character(3) :: idproc_name
        character(300) :: particle_file

        integer, parameter :: nprocs_max = 80
        integer :: p
        real(kind=wp), allocatable, dimension(:,:) :: pardata_real
        integer, allocatable, dimension(:,:) :: pardata_integer

        if (Para_par%imode .eq. 0) then

          inquire(file='PARTICLES_INPUT.dat',exist=lex)
          if (lex) then

            if (ioproc) then
              open(unit=1000,file='PARTICLES_INPUT.dat',
     &                       status='old',action='read')

              read(1000,*) Npar_init
              write(*,*) 'Begin to read PARTICLES_INPUT.dat'
              write(*,'(A,I11)') 'Particles number = ', Npar_init
            end if 
            call MPI_bcast(Npar_init,1,MPI_INTEGER,
     &                     ioid,MPI_comm_fluid,ierr) 
          
cty         Allocate pardata arrays used for reading data, which is covered when do k = 1, Npar_read  
            allocate(pardata(4,Npar_read_onetime))

cty         Obtain the read_times
            read_times = int(Npar_init/Npar_read_onetime) + 1

cty         Initialize Npar_local and Npar_init_inside for each processor
            Npar_local = 0
            Npar_init_inside = Npar_init

            do nk = 1, read_times

cty           Obtain the Npar_read
              if (nk<read_times) then 
                Npar_read = Npar_read_onetime 
              else 
                Npar_read = Npar_init - (read_times-1)*Npar_read_onetime
              end if

cty           Notice the data is only be read in io processor, 
cty           the data is then be bcasted to other processors. 
              if (ioproc) then  
                do k = 1, Npar_read

                  read(1000,*) (pardata(j,k), j = 1, 4)  

                end do ! k 
              end if ! (ioproc)

              call MPI_bcast(pardata,Npar_read*4,
     &                       real_mp_type,ioid,MPI_comm_fluid,ierr) 

cty           Find which subblock is governed by the current processor
              do k = 1, Npar_read
                xpar = pardata(2,k)
                ypar = pardata(3,k)
                zpar = pardata(4,k)

cty             Delete input particles outside the background grid previously
                if (   
     &                  (xpar.lt.B_xmin .or. xpar.ge.B_xmax) 
     &             .or. (ypar.le.B_ymin .or. ypar.gt.B_ymax) 
     &             .or. (zpar.le.B_zmin .or. zpar.gt.B_zmax) 
     &             ) then

                  Npar_init_inside = Npar_init_inside - 1
                  cycle

                else  

cty               Notice that every point must be in background grid now,
cty               then we need to carefully identify which block does the particle locate

cty               Firstly, we need to delete the points located in the blade surface
                  block_id = block_proc(procid_fluid)

cty               Secondly, identify which block does the point locate, and O-type grid is prior
cty               For the sake of simplicity, we only consider the case with two blocks !!!
                  if (.not. particle_in_inner_surface(xpar,ypar)) then

                    do j = 1, num_blocks
                      if (ov_bgrid_block(j)) then
                        which_block = j
cty                     Here we can coarsely identify does the particle only locate in backgroun grid
                        if ( (xpar .lt. O_outer_xmin) .or.
     &                       (xpar .gt. O_outer_xmax) ) then
                          exit
                        end if
                      else if (ov_ogrid_block(j) .and. 
     &                         particle_in_outer_surface(xpar,ypar))
     &                  then
                        which_block = j
                      end if
                    end do

cty                 Subsequently, save particles in the corresponding processor
                    if (which_block .eq. block_id) then
                      call particle_in_subblock(xpar,ypar,
     &                                          procid_fluid,
     &                                          is_in,transpar)

                      if (is_in) then

                        Npar_local = Npar_local + 1
                        Pointer_par(Npar_local)%idpar = pardata(1,k) 
                        Pointer_par(Npar_local)%dpar = 
     &                                 Para_non_dim%diameter_par 
                        Pointer_par(Npar_local)%rhopar = 
     &                                 Para_non_dim%density_par
cty                     transpar is used to simplify the process of compute_ij_particle
                        Pointer_par(Npar_local)%transpar = transpar
                        Pointer_par(Npar_local)%xpar = xpar 
                        Pointer_par(Npar_local)%ypar = ypar
                        Pointer_par(Npar_local)%zpar = zpar
cty                     The initial velocity are setted to zero, but nextly these numerical values 
cty                     may be changed when call set_particles_velocity
                        Pointer_par(Npar_local)%upar = 0._wp
                        Pointer_par(Npar_local)%vpar = 0._wp
                        Pointer_par(Npar_local)%wpar = 0._wp 
cty                     Initially, we can set the acceleration as zero 
                        Pointer_par(Npar_local)%axpar = 0._wp  
                        Pointer_par(Npar_local)%aypar = 0._wp 
                        Pointer_par(Npar_local)%azpar = 0._wp
                        Pointer_par(Npar_local)%procpar = procid_fluid 
                        Pointer_par(Npar_local)%is_monitor = 999

cty                     Notice that ipar, jpar and kpar are not set up to now.
cty                     In the end, find which cell does the particle locate.
cty                     This part is divided into two procedure owing to the 
cty                     different homogeneity in the x-y plane and z direction.
                        call compute_ij_particle(
     &                       transpar,
     &                       xpar,ypar, 
     &                       Pointer_par(Npar_local)%ipar, 
     &                       Pointer_par(Npar_local)%jpar)
                        call compute_k_particle(zpar,
     &                       Pointer_par(Npar_local)%kpar)

                      end if ! (is_in)
                    end if ! (which_block .eq. block_id)
                   
                  else

                  end if ! (.not. particle_in_inner_surface(xpar,ypar))
                end if ! outside the background grid

              end do ! k 

            end do ! nk

            allocate(Npar_local_proc(0:nprocs-1))
cty         Collect Npar_local in ioid, which is save in Npar_local_proc(procid_fluid)
            call MPI_gather(Npar_local,1,MPI_INTEGER,
     &                      Npar_local_proc(procid_fluid),1,MPI_INTEGER,
     &                      ioid,MPI_comm_fluid,ierr)
cty         Check the particle number Npar_local in each processor 
            if (ioproc) then
              do k = 0, nprocs-1
                write(*,'(A,I4,A9,I4,A4,I11)') 
     &                  'The particle number in processor', k, 
     &                  'of block', block_proc(k), 
     &                  'is', Npar_local_proc(k)
              end do
            end if
            deallocate(Npar_local_proc)
           
cty         Sum Npar_local in ioid, which results in Npar_global
            call MPI_reduce(Npar_local,Npar_global,1,MPI_INTEGER,
     &                      MPI_SUM,ioid,MPI_comm_fluid,ierr)

            if (ioproc) then
              close(1000) 
cty           print Npar_init_inside to check the function particle_in_subblock 
              write(*,'(A,I11)') 
     &                'Intially, the total particle number' // 
     &                ' inside the background grid is',
     &                 Npar_init_inside
cty           print Npar_global to check the function particle_in_subblock 
              write(*,'(A,I11)') 
     &                'Intially, the total particle number' // 
     &                ' in the flow field is',
     &                 Npar_global 
              write(*,*) 'PARTICLES_INPUT.dat has been read'
            end if

            deallocate(pardata)

          else
            call stopnow('file PARTICLES_INPUT.dat does not exist')
          end if      

        else if (Para_par%imode .eq. 1) then

cty       Read the particles data from the file PARTICLES_INPUT_CONTINUE*.bin
cty       Set the file name of the particle_file
          write(idproc_name,'(I2)') procid_fluid
cty       Here we use 2 places to mark the id of the processors, 
cty       and fill the blanks with zeros.
          do j = 1, len_trim(idproc_name)
            if (idproc_name(j:j) .eq. '') then
              idproc_name(j:j) = '0' 
            else 
              exit
            end if
          end do 
          particle_file = 'Continue_particles/'//
     &                    'PARTICLES_INPUT_CONTINUE'//
     &                    trim(idproc_name)//'.bin'
          particle_file = trim(adjustl(particle_file))
          inquire(file=particle_file,exist=lex)
          if (lex) then
            open(unit=procid_fluid+1000,file=particle_file,
     &           form='unformatted',status='old',action='read') 

            read(procid_fluid+1000) Npar_local
            
            do k = 1, Npar_local
              read(procid_fluid+1000) Pointer_par(k)%idpar,
     &                                Pointer_par(k)%dpar,    
     &                                Pointer_par(k)%rhopar,  
     &                                Pointer_par(k)%transpar, 
     &                                Pointer_par(k)%xpar,    
     &                                Pointer_par(k)%ypar,    
     &                                Pointer_par(k)%zpar,    
     &                                Pointer_par(k)%upar,    
     &                                Pointer_par(k)%vpar,    
     &                                Pointer_par(k)%wpar,    
     &                                Pointer_par(k)%axpar,    
     &                                Pointer_par(k)%aypar,    
     &                                Pointer_par(k)%azpar    
              read(procid_fluid+1000) Pointer_par(k)%procpar,
     &                                Pointer_par(k)%ipar, 
     &                                Pointer_par(k)%jpar,
     &                                Pointer_par(k)%kpar, 
     &                                Pointer_par(k)%is_monitor

            end do ! k

            allocate(Npar_local_proc(0:nprocs-1))
cty         Collect Npar_local in ioid, which is save in Npar_local_proc(procid_fluid)
            call MPI_gather(Npar_local,1,MPI_INTEGER,
     &                      Npar_local_proc(procid_fluid),1,MPI_INTEGER,
     &                      ioid,MPI_comm_fluid,ierr)
cty         Check the particle number Npar_local in each processor 
            if (ioproc) then
              do k = 0, nprocs-1
                write(*,'(A,I4,A9,I4,A4,I11)') 
     &                  'The particle number in processor', k, 
     &                  'of block', block_proc(k), 
     &                  'is', Npar_local_proc(k)
              end do
            end if
            deallocate(Npar_local_proc)
           
cty         Sum Npar_local in ioid, which results in Npar_global
            call MPI_reduce(Npar_local,Npar_global,1,MPI_INTEGER,
     &                      MPI_SUM,ioid,MPI_comm_fluid,ierr)

            if (ioproc) then
cty           print Npar_global to check the function particle_in_subblock 
              write(*,'(A,I11)') 
     &                'Intially, the total particle number' // 
     &                ' in the flow field is',
     &                 Npar_global 
              write(*,*) 'PARTICLES_INPUT_CONTINUE*.dat has been read'
            end if

            close(procid_fluid+1000)

          else
            call stopnow('file PARTICLES_INPUT_CONTINUE*.dat'//
     &                   ' does not exist')
          end if      

        end if ! Para_par%imode

      end subroutine ! read_particles_data_init



cty========================================================================
cty     subroutine read_monitor_particles
cty========================================================================
cty
cty     Description:
cty     Read the monitor particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine read_monitor_particles

        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_comm_fluid,
     &                      ioid, ioproc
      
        implicit none

        integer :: k
        integer :: ierr
        logical :: lex

        if (Para_par%Monitor_particle_steps .gt. 0) then       

          inquire(file='PARTICLES_MONITOR.dat',exist=lex)
          if (lex) then

            if (ioproc) then
              open(unit=1000,file='PARTICLES_MONITOR.dat',
     &                       status='old',action='read')
              read(1000,*) Npar_monitor
              write(*,*) 'Begin to read PARTICLES_MONITOR.dat'
              write(*,'(A,I11)') 'Initial monitor particles number = ', 
     &                            Npar_monitor
            end if

            call MPI_bcast(Npar_monitor,1,MPI_INTEGER,
     &                     ioid,MPI_comm_fluid,ierr)
            allocate(id_monitor(Npar_monitor))

            if (ioproc) then
              do k = 1, Npar_monitor
                read(1000,*) id_monitor(k)
              end do
              close(1000)
              write(*,*) 'PARTICLES_MONITOR.dat has been read'
            end if

            call MPI_bcast(id_monitor,Npar_monitor,real_mp_type,
     &                     ioid,MPI_comm_fluid,ierr)

          else
            call stopnow('file PARTICLES_MONITOR.dat does not exist')
          end if ! (lex)

        end if ! (Para_par%Monitor_particle_steps .gt. 0)

      end subroutine ! read_monitor_particles



cty========================================================================
cty     subroutine init_monitor_particles
cty========================================================================
cty
cty     Description:
cty     Init the particle%is_monitor to identify monitor particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine init_monitor_particles
      
        use mppvars, only : MPI_INTEGER, MPI_SUM, MPI_comm_fluid, 
     &                      ioid, procid_fluid, nprocs, ioproc

        implicit none

        type(particle), pointer :: Pt_par
        integer :: k, j, p
        integer, allocatable, dimension(:) :: Npar_monitor_local_proc
        integer :: ierr

        if (Para_par%Monitor_particle_steps .gt. 0) then       

          Npar_monitor_local = 0

          do k = 1, Npar_local
            Pt_par => Pointer_par(k)
           
            if (Pt_par%is_monitor .eq. 1) then 
              Npar_monitor_local = Npar_monitor_local + 1
            else if (Pt_par%is_monitor .eq. 999) then
              do j = 1, Npar_monitor
                if (Pt_par%idpar .eq. id_monitor(j)) then
                  Pt_par%is_monitor = 1
                  Npar_monitor_local = Npar_monitor_local + 1
                end if ! (Pt_par%idpar .eq. id_monitor(j))      
              end do ! j
            end if ! (Pt_par%is_monitor .eq. 1)

          end do ! k

          allocate(Npar_monitor_local_proc(0:nprocs-1))
cty       Collect Npar_monitor_local in ioid, 
cty       which is save in Npar_monitor_local_proc(procid_fluid)
          call MPI_gather
     &        (Npar_monitor_local,1,MPI_INTEGER,
     &         Npar_monitor_local_proc(procid_fluid),1,MPI_INTEGER,
     &         ioid,MPI_comm_fluid,ierr)
cty       Check the monitor particles number Npar_monitor_local in each processor 
          if (ioproc) then
            do p = 0, nprocs-1
              write(*,'(A,I4,A9,I4,A4,I11)') 
     &                'The monitor particles number in processor', p, 
     &                'of block', block_proc(p), 
     &                'is', Npar_monitor_local_proc(p)
            end do
          end if
          deallocate(Npar_monitor_local_proc)

cty       Sum Npar_monitor_local in ioid, which results in Npar_monitor_global
          call MPI_reduce(Npar_monitor_local,Npar_monitor_global,1,
     &                    MPI_INTEGER,MPI_SUM,ioid,MPI_comm_fluid,ierr)
          if (ioproc) then
            write(*,'(A,I11)') 
     &              'Intially, the total monitor particles number' // 
     &              ' in the flow field is',
     &               Npar_monitor_global 
          end if
  
        end if ! (Para_par%Monitor_particle_steps .gt. 0)

      end subroutine ! init_monitor_particles



cty========================================================================
cty     function particle_in_block
cty========================================================================
cty
cty     Description:
cty     Determine whether a point lies inside a certain block.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Physical coordinate of the point
cty       block_id: block id to be searched
cty
cty     Output arguments:
cty       is_in: True if (x, y) lies inside this block
cty     
cty     Notice: 
cty       1. The point (x, y) may locate both in the background and O-type 
cty          (or the transiformed O-type) grid.
cty
cty========================================================================

      pure function particle_in_block(x,y,
     &                                block_id) result (is_in)

        use mbvars, only : ov_ogrid_block, ov_bgrid_block
        use overset, only : ov_mesh

        implicit none

        real(kind=wp), intent(in) :: x, y
        integer, intent(in) :: block_id
        logical :: is_in
        logical :: is_in_tmp

cty     For background grid          
        if (ov_bgrid_block(block_id)) then

          if ( (x .ge. B_xmin) .and. 
     &         (x .lt. 1.65_wp)
     &       ) then
            is_in = .true.
          else
            is_in = .false.
          end if ! (x.ge.B_xmin .and. x.lt.B_xmax)
        end if !(ov_bgrid_block(block_id))
       
cty     For O-type grid          
        if (ov_ogrid_block(block_id)) then

          if ( (x .lt. O_outer_xmin) .or. 
     &         (x .gt. O_outer_xmax) ) then
            is_in = .false.
            return
          else
            ! ym
            is_in = particle_in_polygon(x,y,
     &              O_inner_px,O_inner_py,size(O_inner_px))
            if (ov_mesh(block_id)%wrap_block) then
              is_in_tmp = particle_in_polygon(x,y,
     &                    O_inner_px,O_inner_py_trans,size(O_inner_px))
            end if

            ! yp
            is_in = is_in .neqv. particle_in_polygon(x,y,
     &                           O_outer_px,O_outer_py,size(O_outer_px))
            if (ov_mesh(block_id)%wrap_block) then
              is_in_tmp = is_in_tmp .neqv. 
     &                    particle_in_polygon(x,y,
     &                    O_outer_px,O_outer_py_trans,size(O_outer_px))
            end if

            if (is_in .or. is_in_tmp) then
              is_in = .true.
            else
              is_in = .false.
            end if

          end if ! (x.lt.O_outer_xmin .or. x.gt.O_outer_xmax)

        end if ! (ov_ogrid_block(block_id))

      end function ! particle_in_block



cty========================================================================
cty     function particle_in_inner_surface
cty========================================================================
cty
cty     Description:
cty     Determine whether a point lies inside the blade surface or not, if true, 
cty     then the particle will be omitted or be dealed with boundary condition.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Physical coordinate of the point
cty
cty     Output arguments:
cty       is_in: True if (x, y) lies inside the blade surface 
cty              (contains the blade surface) 
cty     
cty     Notice: 
cty
cty========================================================================

      pure function particle_in_inner_surface(x,y) result (is_in)

        use mbvars, only : num_blocks, ov_ogrid_block
        use overset, only : ov_mesh

        implicit none

        real(kind=wp), intent(in) :: x, y

        logical :: is_in
     &           , is_in_tmp

        integer :: block_id

        block_id = 2        

        if ( (x .lt. O_inner_xmin) .or. 
     &       (x .gt. O_inner_xmax) ) then
          is_in = .false.
          return
        else
          ! ym
          is_in = particle_in_polygon(x,y,
     &            O_inner_px,O_inner_py,size(O_inner_px))
          if (ov_mesh(block_id)%wrap_block) then
            is_in_tmp = particle_in_polygon(x,y,
     &                  O_inner_px,O_inner_py_trans,size(O_inner_px))
          end if
  
          if (is_in .or. is_in_tmp) then
            is_in = .true.
          else
            is_in = .false.
          end if

        end if ! (x.lt.O_inner_xmin .or. x.gt.O_inner_xmax)

      end function ! particle_in_inner_surface



cty========================================================================
cty     function particle_in_outer_surface
cty========================================================================
cty
cty     Description:
cty     Determine whether a point lies inside the outer boundary of the O-type
cty     grid or not.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Physical coordinate of the point
cty
cty     Output arguments:
cty       is_in: True if (x, y) lies inside the outer boundary of the O-type 
cty              grid. (contains the surface) 
cty     
cty     Notice: 
cty
cty========================================================================

      pure function particle_in_outer_surface(x,y) result (is_in)

        use mbvars, only : num_blocks, ov_ogrid_block
        use overset, only : ov_mesh

        implicit none

        real(kind=wp), intent(in) :: x, y

        logical :: is_in
     &           , is_in_tmp

        integer :: block_id

        block_id = 2

        if ( (x .lt. O_outer_xmin) .or. 
     &       (x .gt. O_outer_xmax) ) then
          is_in = .false.
          return
        else
          ! yp
          is_in = particle_in_polygon(x,y,
     &            O_outer_px,O_outer_py,size(O_outer_px))
          if (ov_mesh(block_id)%wrap_block) then
            is_in_tmp = particle_in_polygon(x,y,
     &                  O_outer_px,O_outer_py_trans,size(O_outer_px))
          end if
  
          if (is_in .or. is_in_tmp) then
            is_in = .true.
          else
            is_in = .false.
          end if

        end if ! (x.lt.O_outer_xmin .or. x.gt.O_outer_xmax)

      end function ! particle_in_outer_surface



cty========================================================================
cty     subroutine particle_in_subblock
cty========================================================================
cty
cty     Description:
cty     Determine if a point lies inside a subvolume of the certain block.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Physical coordinate of the point
cty       proc: Processor id
cty
cty     Output arguments:
cty       is_in: True if (x, y) lies inside the subblock governed by the
cty              present processor
cty       transpar: trans_mode if (x, y) lies inside the y transformed subblock of 
cty              the O-type grid.
cty     
cty     Notice: 
cty
cty========================================================================

      pure subroutine particle_in_subblock(x,y,
     &                                     proc,is_in,transpar)

        use mbvars, only : ov_ogrid_block, ov_bgrid_block
        use overset, only : ov_mesh, trans_mode

        implicit none

        real(kind=wp), intent(in) :: x, y
        integer, intent(in) :: proc
        logical, intent(out) :: is_in
        real(kind=wp), intent(out) :: transpar

        logical :: is_in_ori, is_in_trans
        logical :: on_y
     &           , on_y_trans

        integer :: block_id
        integer :: nx, ny
        integer :: i, j

        real(kind=wp), allocatable :: px(:)
        real(kind=wp), allocatable :: py(:)
        real(kind=wp), allocatable :: py_trans(:) 

cty     Set parameters associated with the present processor: block_id, ind_i, ind_j
        block_id = block_proc(proc)

cty     Here nx and ny save the number of points in each direction
        nx = subblock_proc(proc)%xsize 
        ny = subblock_proc(proc)%ysize 

        if (ov_bgrid_block(block_id)) then
cty     Background grid: xperi_block is false while yperi_block is ture

          if ( 
     &             (x .ge. subblock_proc(proc)%xmin)  
     &       .and. (y .gt. subblock_proc(proc)%ymin)  
     &       .and. (y .le. subblock_proc(proc)%ymax) 
     &       .and. (x .lt. subblock_proc(proc)%xmax) 
     &       ) then
            is_in = .true. 
          else 
            is_in = .false.
          end if
          transpar = 0._wp

        end if ! (ov_bgrid_block(block_id))

        if (ov_ogrid_block(block_id)) then
cty     O-type grid: xperi_block is true while yperi_block is false

          allocate(px(2*(nx+ny)-4))
          allocate(py(2*(nx+ny)-4))
          allocate(py_trans(2*(nx+ny)-4))

          px = subblock_proc(proc)%px
          py = subblock_proc(proc)%py
          py_trans = subblock_proc(proc)%py_trans
         
          is_in_ori = particle_in_inner_polygon(x,y,px,py,2*(nx+ny)-4)
          if(ov_mesh(block_id)%wrap_block) then
            is_in_trans = 
     &      particle_in_inner_polygon(x,y,px,py_trans,2*(nx+ny)-4)
          endif

          if (is_in_ori .or. is_in_trans) then
            is_in = .true. 
            if (is_in_trans) then
              transpar = trans_mode
            else
              transpar = 0._wp 
            end if
          else
            is_in = .false. 
          end if
         
          j = 2*(nx+ny)-4
          do i = 1, 2*(nx+ny)-4

            on_y =  particle_on_segment
     &              (x,y,px(i),py(i),px(j),py(j))
            on_y_trans = particle_on_segment
     &                   (x,y,px(i),py_trans(i),px(j),py_trans(j))

            if (on_y .or. on_y_trans) then

              if (i .eq. 1) then
                if ( (dsqrt((x-px(i))**2._wp
     &                     +(y-py(i))**2._wp) .le. 0._wp) 
     &          .or. (dsqrt((x-px(i))**2._wp
     &                     +(y-py_trans(i))**2._wp) .le. 0._wp) 
     &              )
     &             then
                  is_in = .false.
                  transpar = 0._wp
                else
                  is_in = .true.
                  if (on_y_trans) then
                    transpar = trans_mode
                  else
                    transpar = 0._wp 
                  end if
                end if
                return

              else if ( (i .ge. 2) .and. 
     &                  (i .le. (nx+ny-1)) ) then
                is_in = .false.
                transpar = 0._wp
                return

              else if (i .eq. (nx+ny)) then
                if ( (dsqrt((x-px(j))**2._wp
     &                     +(y-py(j))**2._wp) .le. 0._wp) 
     &          .or. (dsqrt((x-px(j))**2._wp
     &                     +(y-py_trans(j))**2._wp) .le. 0._wp)
     &              )
     &             then
                  is_in = .false.
                  transpar = 0._wp
                else
                  is_in = .true.
                  if (on_y_trans) then
                    transpar = trans_mode
                  else 
                    transpar = 0._wp
                  end if
                end if
                return

              else if ( (i .gt. (nx+ny)) .and. 
     &                  (i .le. 2*(nx+ny)-4) ) then
                is_in = .true.
                if (on_y_trans) then
                  transpar = trans_mode
                else  
                  transpar = 0._wp
                end if
                return
              end if         
       
            end if ! (on_y .or. on_y_trans)

            j = i
          end do ! i

          deallocate(px)
          deallocate(py)
          deallocate(py_trans)

        end if ! (ov_ogrid_block(block_id))

      end subroutine ! particle_in_subblock



cty========================================================================
cty     function particle_in_inner_polygon
cty========================================================================
cty
cty     Description:
cty     Determine if a point lies inside the inner part of the given polygon.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Point to be searched
cty       px, py: Coordinates of the vertices of the polygon, order must be
cty               correct and extremities should not be repeated
cty
cty     Output arguments:
cty       is_in: True if (x,y) lies inside the inner part of the polygon
cty
cty     Notice: 
cty
cty========================================================================

      pure function particle_in_inner_polygon(x,y,px,py,n) result(is_in)

        implicit none

        real(kind=wp), intent(in) :: x, y
        real(kind=wp), intent(in) :: px(n), py(n)
        integer, intent(in) :: n

        logical :: is_in

        integer :: i, j

        is_in = .false.

        j = n

        do i = 1, n
          if ( ( (py(i) .lt. y) .and. 
     &           (py(j) .ge. y) ) .or.
     &         ( (py(j) .lt. y) .and. 
     &           (py(i) .ge. y) ) ) then
            if (px(i)+(y-py(i))/(py(j)-py(i))*(px(j)-px(i)) .lt. x) then
              is_in = .not. is_in
            end if
          end if
          j = i
        end do

      end function ! particle_in_inner_polygon



cty========================================================================
cty     function particle_in_polygon
cty========================================================================
cty
cty     Description:
cty     Determine if a point lies inside a given polygon.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Point to be searched
cty       px, py: Coordinates of the vertices of the polygon, order must be
cty               correct and extremities should not be repeated
cty
cty     Output arguments:
cty       is_in: True if (x,y) lies inside the polygon
cty
cty     Notice: 
cty
cty========================================================================

      pure function particle_in_polygon(x,y,px,py,n) result(is_in)

        implicit none

        real(kind=wp), intent(in) :: x, y
        real(kind=wp), intent(in) :: px(n), py(n)
        integer, intent(in) :: n

        logical :: is_in

        integer :: i, j

        is_in = .false.

        j = n

        do i = 1, n
          if (particle_on_segment(x,y,px(i),py(i),px(j),py(j))) then
            is_in = .true.
            return
          end if
          if ( ( (py(i) .lt. y) .and. 
     &           (py(j) .ge. y) ) .or.
     &         ( (py(j) .lt. y) .and. 
     &           (py(i) .ge. y) ) ) then
            if (px(i)+(y-py(i))/(py(j)-py(i))*(px(j)-px(i)) .lt. x) then
              is_in = .not. is_in
            end if
          end if
          j = i
        end do

      end function ! particle_in_polygon



cty========================================================================
cty     function particle_on_segment
cty========================================================================
cty
cty     Description:
cty     Determine if a point lies on a given segment.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Point to determine
cty       x1, y1: First extremity of the segment
cty       x2, y2: Second extremity of the segment
cty
cty     Output arguments:
cty       is_on: True if (x,y) lies on the segment (x1,y1)--(x2,y2)
cty
cty     Notice: 
cty
cty========================================================================

      pure function particle_on_segment(x,y,x1,y1,x2,y2) result(is_on)

        implicit none

        real(kind=wp), intent(in) :: x, y, x1, y1, x2, y2

        logical :: is_on

        real(kind=wp) :: d1, d2, d

        d1 = dsqrt((x-x1)**2._wp + (y-y1)**2._wp)
        d2 = dsqrt((x-x2)**2._wp + (y-y2)**2._wp)
        d = dsqrt((x1-x2)**2._wp + (y1-y2)**2._wp)

        is_on = (d1+d2) .le. d

      end function ! particle_on_segment



cty========================================================================
cty     subroutine compute_ij_particle
cty========================================================================
cty
cty     Description:
cty     Compute which the cell does the particle locate. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar: Coordinates of the particle
cty
cty     Output arguments:
cty       ipar, jpar: Local indices of the cell where the particle locates in
cty
cty     Notice:
cty
cty========================================================================

      pure subroutine compute_ij_particle(transpar,xpar,ypar,ipar,jpar)

        use mppvars, only : procid_fluid
        use mbvars, only : ov_ogrid_block, ov_bgrid_block

        implicit none 

        real(kind=wp), intent(in) :: transpar, xpar, ypar
        integer, intent(out) :: ipar, jpar

        integer :: i, j
        integer :: nx, ny
        integer :: ibpar
        real(kind=wp) :: y

        logical :: is_in

cty     Set the original or transformed coordinates of the particle
cty     for the subsequent calculation
        y = transformed_y(transpar,ypar)

cty     Set parameters associated with the present processor: ibpar, ind_i, ind_j
        ibpar = block_proc(procid_fluid)

cty     Here nx and ny save the number of points in each direction
        nx = subblock_proc(procid_fluid)%xsize
        ny = subblock_proc(procid_fluid)%ysize

        if (ov_bgrid_block(ibpar)) then

          do i = 1, nx-1
            if ( (xpar .ge. B_z_x(i)) .and. 
     &           (xpar .lt. B_z_x(i+1)) ) then
              ipar = i 
              exit
            end if
          end do ! i

          do j = 1, ny-1
            if ( (y .gt. B_r_y(j)) .and. 
     &           (y .le. B_r_y(j+1)) ) then
              jpar = j 
              exit
            end if
          end do ! j

        end if

        if (ov_ogrid_block(ibpar)) then  
cty     O-type grid: xperi_block is true while yperi_block is false
     
          do j = 1, ny-1
            do i = 1, nx-1

              call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                         ipar,jpar)
              if (is_in) return

            end do ! i = 1, nx-1
          end do ! j = 1, ny-1

        end if ! (ov_bgrid_block(ibpar))

      end subroutine ! compute_ij_particle



cty========================================================================
cty     subroutine update_B_grid_ij_particle
cty========================================================================
cty
cty     Description:
cty     For the background grid, compute which the cell does the particle 
cty     locate.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar: Coordinates of the particle
cty       upar, vpar: Velocities of the particle
cty       procpar_ori: Original processor id of the particle
cty       proc: Present processor id of the particle
cty       ipar_ori, jpar_ori: Local indices of the cell where the particle 
cty                           originally locates in
cty
cty     Output arguments:
cty       ipar, jpar: Local indices of the cell where the particle 
cty                   presently locates in
cty
cty     Notice:
cty
cty========================================================================

      subroutine update_B_grid_ij_particle
     &          (transpar,xpar,ypar,upar,vpar,procpar_ori,proc,
     &           ipar_ori,jpar_ori,ipar,jpar)

        use mppvars, only : procid_fluid

        implicit none 

        real(kind=wp), intent(in) :: transpar, xpar, ypar
        real(kind=wp), intent(in) :: upar, vpar
        integer, intent(inout) :: procpar_ori
        integer, intent(in) :: proc
        integer, intent(in) :: ipar_ori, jpar_ori
        integer, intent(out) :: ipar, jpar

        integer :: i, j
        integer :: nx, ny
        integer :: start_i, end_i
        integer :: start_j, end_j
        real(kind=wp) :: y

cty     Set the original or transformed coordinates of the particle
cty     for the subsequent calculation
        y = transformed_y(transpar,ypar)

cty     Here nx and ny save the number of points in each direction
        nx = subblock_proc(procid_fluid)%xsize
        ny = subblock_proc(procid_fluid)%ysize

cty     Set impossible initial values of ipar and jpar for debug
        ipar = -9999
        jpar = -9999

        if (upar .ge. 0._wp) then
          if (procpar_ori .eq. proc) then
            start_i = ipar_ori
          else
            start_i = 1
          end if ! (procpar_ori .eq. proc)
        else
          if (procpar_ori .eq. proc) then
            end_i = ipar_ori
          else
            end_i = nx-1
          end if ! (procpar_ori .eq. proc)
        end if ! (upar .ge. 0._wp)

        if (vpar .ge. 0._wp) then
          if (procpar_ori .eq. proc) then
            start_j = jpar_ori
          else
            start_j = 1
          end if ! (procpar_ori .eq. proc)
        else
          if (procpar_ori .eq. proc) then
            end_j = jpar_ori
          else
            end_j = ny-1
          end if ! (procpar_ori .eq. proc)
        end if ! (vpar .ge. 0._wp)

        if (upar .ge. 0._wp) then
          do i = start_i, nx-1
            if ( (xpar .ge. B_z_x(i)) .and.
     &           (xpar .lt. B_z_x(i+1)) ) then
              ipar = i
              exit
            end if
          end do ! i
        else
          do i = end_i, 1, -1
            if ( (xpar .ge. B_z_x(i)) .and.
     &           (xpar .lt. B_z_x(i+1)) ) then
              ipar = i
              exit
            end if
          end do ! i
        end if

        if (vpar .ge. 0._wp) then
          do j = start_j, ny-1
            if ( (y .gt. B_r_y(j)) .and.
     &           (y .le. B_r_y(j+1)) ) then
              jpar = j
              exit
            end if
          end do ! j
        else
          do j = end_j, 1, -1
            if ( (y .gt. B_r_y(j)) .and.
     &           (y .le. B_r_y(j+1)) ) then
              jpar = j
              exit
            end if
          end do ! j
        end if

        procpar_ori = procid_fluid 

        if ( (ipar .eq. (-9999)) .or. 
     &       (jpar .eq. (-9999)) ) then
          call compute_ij_particle(transpar,xpar,ypar,ipar,jpar)
        end if

      end subroutine ! update_B_grid_ij_particle



cty========================================================================
cty     subroutine identify_B_grid_ij_particle
cty========================================================================
cty
cty     Description:
cty     Indentify which the cell does the particle locate. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Coordinates
cty       i, j: Local indices
cty
cty     Output arguments:
cty       is_in: True if (x,y) lies inside the cell
cty       ipar, jpar: Local indices of the cell where the particle locates in
cty
cty     Notice:
cty
cty========================================================================

      pure subroutine identify_B_grid_ij_particle(x,y,i,j,is_in,
     &                                            ipar,jpar)

        implicit none 

        real(kind=wp), intent(in) :: x, y
        integer, intent(in) :: i, j
        integer, intent(out) :: ipar, jpar
        logical, intent(out) :: is_in

        if ( (x .ge. B_z_x(i)) .and.
     &       (x .lt. B_z_x(i+1)) .and.
     &       (y .gt. B_r_y(j)) .and.
     &       (y .le. B_r_y(j+1)) ) then
          is_in = .true.
          ipar = i
          jpar = j
          return
        else
          is_in = .false.
          return
        end if

      end subroutine ! identify_B_grid_ij_particle



cty========================================================================
cty     subroutine identify_O_grid_ij_particle
cty========================================================================
cty
cty     Description:
cty     Indentify which the cell does the particle locate. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Coordinates
cty       i, j: Local indices
cty
cty     Output arguments:
cty       is_in: True if (x,y) lies inside the cell
cty       ipar, jpar: Local indices of the cell where the particle locates in
cty
cty     Notice:
cty
cty========================================================================

      pure subroutine identify_O_grid_ij_particle(x,y,i,j,is_in,
     &                                            ipar,jpar)

        use corevars, only : z_x, r_y

        implicit none 

        real(kind=wp), intent(in) :: x, y
        integer, intent(in) :: i, j
        integer, intent(out) :: ipar, jpar
        logical, intent(out) :: is_in

        if (particle_in_convex_quadrilateral
     &      (           x,            y,
     &       z_x(i  ,j  ), r_y(i  ,j  ),
     &       z_x(i+1,j  ), r_y(i+1,j  ),
     &       z_x(i+1,j+1), r_y(i+1,j+1),
     &       z_x(i  ,j+1), r_y(i  ,j+1) )) then
          is_in = .true.
          ipar = i
          jpar = j 
          return
        else if ( (particle_on_segment
     &             (           x,            y,
     &              z_x(i  ,j+1), r_y(i  ,j+1),
     &              z_x(i  ,j  ), r_y(i  ,j  ) )) .and.
     &            (dsqrt((x-z_x(i,j))**2._wp + 
     &                   (y-r_y(i,j))**2._wp) .gt. 0._wp) 
     &          ) then
          is_in = .true. 
          ipar = i
          jpar = j
          return 
        else if ( (particle_on_segment
     &             (           x,            y,
     &              z_x(i+1,j+1), r_y(i+1,j+1),
     &              z_x(i  ,j+1), r_y(i  ,j+1) )) .and.
     &            (dsqrt((x-z_x(i+1,j+1))**2._wp + 
     &                   (y-r_y(i+1,j+1))**2._wp) .gt. 0._wp)
     &          ) then
          is_in = .true.
          ipar = i
          jpar = j 
          return
        else
          is_in = .false.
          return
        end if

      end subroutine ! identify_O_grid_ij_particle



cty========================================================================
cty     subroutine update_O_grid_ij_particle_prediction_first_layer
cty========================================================================
cty
cty     Description:
cty     For the O-type grid, compute which the cell does the particle locate
cty     based on the present velocity and the first layer. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar: Coordinates of the particle
cty       vel_tangen, vel_normal: Velocities of the particle for the prediction
cty       procpar: Original processor id of the particle
cty       proc: Present processor id of the particle
cty       ipar_ori, jpar_ori: Local indices of the cell where the particle 
cty                           originally locates in
cty
cty     Output arguments:
cty       ipar, jpar: Local indices of the cell where the particle 
cty                   presently locates in
cty
cty     Notice:
cty
cty========================================================================

      subroutine update_O_grid_ij_particle_prediction_first_layer
     &          (transpar,xpar,ypar,vel_tangen,vel_normal,
     &           procpar_ori,proc,ipar_ori,jpar_ori,ipar,jpar)

        use corevars, only : nxp, nyp, z_x, r_y
        use global_grid, only : grids

        implicit none 

        real(kind=wp), intent(in) :: transpar, xpar, ypar
        real(kind=wp), intent(in) :: vel_tangen, vel_normal
        integer, intent(in) :: procpar_ori
        integer, intent(in) :: proc
        integer, intent(in) :: ipar_ori, jpar_ori
        integer, intent(out) :: ipar, jpar

        integer :: i, j
        integer :: nx, ny
        real(kind=wp) :: y

        logical :: is_in

cty     Set the original or transformed coordinates of the particle
cty     for the subsequent calculation
        y = transformed_y(transpar,ypar)

cty     Here nx and ny save the number of points in each direction
        nx = subblock_proc(proc)%xsize 
        ny = subblock_proc(proc)%ysize 

cty     Set impossible initial values of ipar and jpar for debug
        ipar = -9999
        jpar = -9999

        if (procpar_ori .eq. proc) then
 
cty       O-type grid: xperi_block is true while yperi_block is false

          if (vel_tangen .ge. 0._wp) then

            if (vel_normal .ge. 0._wp) then

              do j = jpar_ori, min(jpar_ori+1,ny-1)
                do i = ipar_ori, min(ipar_ori+1,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j

            else

              do j = jpar_ori, max(jpar_ori-1,1), -1
                do i = ipar_ori, min(ipar_ori+1,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 

            end if ! (vel_normal .ge. 0._wp)

          else

            if (vel_normal .ge. 0._wp) then

              do j = jpar_ori, min(jpar_ori+1,ny-1)
                do i = ipar_ori, max(ipar_ori-1,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j

            else

              do j = jpar_ori, max(jpar_ori-1,1), -1
                do i = ipar_ori, max(ipar_ori-1,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 

            end if ! (vel_normal .ge. 0._wp)

          end if ! (vel_tangen .ge. 0._wp)

        else
 
cty       O-type grid: xperi_block is true while yperi_block is false
     
          if (vel_tangen .ge. 0._wp) then

            if (vel_normal .ge. 0._wp) then

              do j = jpar_ori, min(jpar_ori+1,ny-1)
                do i = 1, 1 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = 1, 1
                do i = ipar_ori, min(ipar_ori+1,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = 1, 1 
                do i = 1, 1 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j

            else

              do j = jpar_ori, max(jpar_ori-1,1), -1
                do i = 1, 1 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = ny-1, ny-1, -1 
                do i = ipar_ori, min(ipar_ori+1,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = ny-1, ny-1, -1 
                do i = 1, 1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 

            end if ! (vel_normal .ge. 0._wp)

          else

            if (vel_normal .ge. 0._wp) then

              do j = jpar_ori, min(jpar_ori+1,ny-1)
                do i = nx-1, nx-1, -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = 1, 1
                do i = ipar_ori, max(ipar_ori-1,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = 1, 1 
                do i = nx-1, nx-1, -1 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j

            else

              do j = jpar_ori, max(jpar_ori-1,1), -1
                do i = nx-1, nx-1, -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = ny-1, ny-1, -1
                do i = ipar_ori, max(ipar_ori-2,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = ny-1, ny-1, -1
                do i = nx-1, nx-1, -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 

            end if ! (vel_normal .ge. 0._wp)

          end if ! (vel_tangen .ge. 0._wp)

        end if ! (procpar_ori .eq. proc)

      end subroutine ! update_O_grid_ij_particle_prediction_first_layer



cty========================================================================
cty     subroutine update_O_grid_ij_particle_prediction_second_layer
cty========================================================================
cty
cty     Description:
cty     For the O-type grid, compute which the cell does the particle locate
cty     based on the present velocity and the second layer. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar: Coordinates of the particle
cty       vel_tangen, vel_normal: Velocities of the particle for the prediction
cty       procpar: Original processor id of the particle
cty       proc: Present processor id of the particle
cty       ipar_ori, jpar_ori: Local indices of the cell where the particle 
cty                           originally locates in
cty
cty     Output arguments:
cty       ipar, jpar: Local indices of the cell where the particle 
cty                   presently locates in
cty
cty     Notice:
cty
cty========================================================================

      subroutine update_O_grid_ij_particle_prediction_second_layer
     &          (transpar,xpar,ypar,vel_tangen,vel_normal,
     &           procpar_ori,proc,ipar_ori,jpar_ori,ipar,jpar)

        use corevars, only : nxp, nyp, z_x, r_y
        use global_grid, only : grids

        implicit none 

        real(kind=wp), intent(in) :: transpar, xpar, ypar
        real(kind=wp), intent(in) :: vel_tangen, vel_normal
        integer, intent(in) :: procpar_ori
        integer, intent(in) :: proc
        integer, intent(in) :: ipar_ori, jpar_ori
        integer, intent(out) :: ipar, jpar

        integer :: i, j
        integer :: nx, ny
        real(kind=wp) :: y

        logical :: is_in

cty     Set the original or transformed coordinates of the particle
cty     for the subsequent calculation
        y = transformed_y(transpar,ypar)

cty     Here nx and ny save the number of points in each direction
        nx = subblock_proc(proc)%xsize 
        ny = subblock_proc(proc)%ysize 

cty     Set impossible initial values of ipar and jpar for debug
        ipar = -9999
        jpar = -9999

        if (procpar_ori .eq. proc) then
 
cty       O-type grid: xperi_block is true while yperi_block is false

          if (vel_tangen .ge. 0._wp) then

            if (vel_normal .ge. 0._wp) then

              do j = jpar_ori, min(jpar_ori+1,ny-1)
                do i = min(ipar_ori+2,nx-1), min(ipar_ori+2,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = min(jpar_ori+2,ny-1), min(jpar_ori+2,ny-1)
                do i = ipar_ori, min(ipar_ori+1,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = min(jpar_ori+2,ny-1), min(jpar_ori+2,ny-1)
                do i = min(ipar_ori+2,nx-1), min(ipar_ori+2,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j

            else

              do j = jpar_ori, max(jpar_ori-1,1), -1
                do i = min(ipar_ori+2,nx-1), min(ipar_ori+2,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = max(jpar_ori-2,1), max(jpar_ori-2,1), -1
                do i = ipar_ori, min(ipar_ori+1,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = max(jpar_ori-2,1), max(jpar_ori-2,1), -1
                do i = min(ipar_ori+2,nx-1), min(ipar_ori+2,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 

            end if ! (vel_normal .ge. 0._wp)

          else

            if (vel_normal .ge. 0._wp) then

              do j = jpar_ori, min(jpar_ori+1,ny-1)
                do i = max(ipar_ori-2,1), max(ipar_ori-2,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = min(jpar_ori+2,ny-1), min(jpar_ori+2,ny-1)
                do i = ipar_ori, max(ipar_ori-1,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = min(jpar_ori+2,ny-1), min(jpar_ori+2,ny-1)
                do i = max(ipar_ori-2,1), max(ipar_ori-2,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j

            else

              do j = jpar_ori, max(jpar_ori-1,1), -1
                do i = max(ipar_ori-2,1), max(ipar_ori-2,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = max(jpar_ori-2,1), max(jpar_ori-2,1), -1
                do i = ipar_ori, max(ipar_ori-1,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = max(jpar_ori-2,1), max(jpar_ori-2,1), -1
                do i = max(ipar_ori-2,1), max(ipar_ori-2,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 

            end if ! (vel_normal .ge. 0._wp)

          end if ! (vel_tangen .ge. 0._wp)

        else
 
cty       O-type grid: xperi_block is true while yperi_block is false
     
          if (vel_tangen .ge. 0._wp) then

            if (vel_normal .ge. 0._wp) then

              do j = jpar_ori, min(jpar_ori+1,ny-1)
                do i = 2, 2 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = 2, 2
                do i = ipar_ori, min(ipar_ori+1,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = 2, 2 
                do i = 2, 2 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j

            else

              do j = jpar_ori, max(jpar_ori-1,1), -1
                do i = 2, 2 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = ny-2, ny-2, -1 
                do i = ipar_ori, min(ipar_ori+1,nx-1)
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = ny-2, ny-2, -1 
                do i = 2, 2 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 

            end if ! (vel_normal .ge. 0._wp)

          else

            if (vel_normal .ge. 0._wp) then

              do j = jpar_ori, min(jpar_ori+1,ny-1)
                do i = nx-2, nx-2, -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = 2, 2 
                do i = ipar_ori, max(ipar_ori-1,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j
              do j = 2, 2 
                do i = nx-2, nx-2, -1 
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i 
              end do ! j

            else

              do j = jpar_ori, max(jpar_ori-1,1), -1
                do i = nx-2, nx-2, -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = ny-2, ny-2, -1
                do i = ipar_ori, max(ipar_ori-1,1), -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 
              do j = ny-2, ny-2, -1
                do i = nx-2, nx-2, -1
                  
                  call identify_O_grid_ij_particle(xpar,y,i,j,is_in,
     &                                             ipar,jpar)
                  if (is_in) return

                end do ! i
              end do ! j 

            end if ! (vel_normal .ge. 0._wp)

          end if ! (vel_tangen .ge. 0._wp)

        end if ! (procpar_ori .eq. proc)

      end subroutine ! update_O_grid_ij_particle_prediction_second_layer



cty========================================================================
cty     subroutine update_O_grid_ij_particle_prediction
cty========================================================================
cty
cty     Description:
cty     For the O-type grid, compute which the cell does the particle locate
cty     based on the present velocity. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar: Coordinates of the particle
cty       vel_tangen, vel_normal: Velocities of the particle for the prediction
cty       procpar: Original processor id of the particle
cty       proc: Present processor id of the particle
cty       ipar_ori, jpar_ori: Local indices of the cell where the particle 
cty                           originally locates in
cty
cty     Output arguments:
cty       ipar, jpar: Local indices of the cell where the particle 
cty                   presently locates in
cty
cty     Notice:
cty
cty========================================================================

      subroutine update_O_grid_ij_particle_prediction
     &          (transpar,xpar,ypar,vel_tangen,vel_normal,
     &           procpar,proc,ipar_ori,jpar_ori,ipar,jpar)

        use mppvars, only : procid_fluid

        implicit none 

        real(kind=wp), intent(in) :: transpar, xpar, ypar
        real(kind=wp), intent(in) :: vel_tangen, vel_normal
        integer, intent(inout) :: procpar
        integer, intent(in) :: proc
        integer, intent(in) :: ipar_ori, jpar_ori
        integer, intent(out) :: ipar, jpar

        call update_O_grid_ij_particle_prediction_first_layer
     &      (transpar,xpar,ypar,vel_tangen,vel_normal,
     &       procpar,procid_fluid,ipar_ori,jpar_ori,ipar,jpar)

        if ( (ipar .eq. (-9999)) .or. 
     &       (jpar .eq. (-9999)) ) then
          call update_O_grid_ij_particle_prediction_second_layer
     &        (transpar,xpar,ypar,vel_tangen,vel_normal,
     &         procpar,procid_fluid,ipar_ori,jpar_ori,ipar,jpar)
        end if

        procpar = procid_fluid

        if ( (ipar .eq. (-9999)) .or. 
     &       (jpar .eq. (-9999)) ) then
          call compute_ij_particle(transpar,xpar,ypar,ipar,jpar)
        end if

      end subroutine ! update_O_grid_ij_particle_prediction



cty========================================================================
cty     subroutine update_B_O_grid_ij_particle
cty========================================================================
cty
cty     Description:
cty     For the particle from the background grid flow into the O-type grid, 
cty     compute which the cell does the particle locate.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar: Coordinates of the particle
cty       procpar_ori: Original processor id of the particle
cty       proc: Present processor id of the particle
cty       ipar_ori, jpar_ori: Local indices of the cell where the particle 
cty                           originally locates in
cty
cty     Output arguments:
cty       ipar, jpar: Local indices of the cell where the particle 
cty                   presently locates in
cty
cty     Notice:
cty
cty========================================================================

      subroutine update_B_O_grid_ij_particle
     &          (transpar,xpar,ypar,procpar_ori,proc,
     &           ipar_ori,jpar_ori,ipar,jpar)
       
        use mppvars, only : procid_fluid
        use corevars, only : nxp, nyp, z_x, r_y
        use global_grid, only : grids

        implicit none 

        real(kind=wp), intent(in) :: transpar, xpar, ypar
        integer, intent(inout) :: procpar_ori
        integer, intent(in) :: proc
        integer, intent(in) :: ipar_ori, jpar_ori
        integer, intent(out) :: ipar, jpar

        integer :: i, j
        integer :: nx, ny
        real(kind=wp) :: y

        logical :: is_in

cty     Set the original or transformed coordinates of the particle
cty     for the subsequent calculation
        y = transformed_y(transpar,ypar)

cty     Here nx and ny save the number of points in each direction
        nx = subblock_proc(proc)%xsize 
        ny = subblock_proc(proc)%ysize 

cty     Set impossible initial values of ipar and jpar for debug
        ipar = -9999
        jpar = -9999
 
cty     O-type grid: xperi_block is true while yperi_block is false
        do j = ny-1, 1, -1
          do i = 1, nx-1
            
            call identify_O_grid_ij_particle(xpar,y,i,j,is_in,ipar,jpar)
            if (is_in) return

          end do ! i = 1, nx-1
        end do ! j = ny-1, 1, -1

        procpar_ori = procid_fluid 

        if ( (ipar .eq. (-9999)) .or. 
     &       (jpar .eq. (-9999)) ) then
          call compute_ij_particle(transpar,xpar,ypar,ipar,jpar)
        end if

      end subroutine ! update_B_O_grid_ij_particle



cty========================================================================
cty     function transformed_y
cty========================================================================
cty
cty     Description:
cty     Determine the transformed y location of the particle.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       y: The original y location of the particle
cty
cty     Output arguments:
cty       trans_y: The transformed y location of the particle
cty
cty     Notice: 
cty       1. The transformation used here is the opposite of the transformation 
cty          used in the ov_mesh.
cty       2. trans_mode is a new varible defiend in overset.F
cty          (1) if the particle locates in the original O-type grid, 
cty              then trans_mode = 0._wp 
cty          (2) if (yo_min.lt.y_min+tol), 
cty              then trans_mode = -1._wp, y need to be moved down by one pitch
cty          (3) if (yo_max.gt.y_max-tol), 
cty              then trans_mode =  1._wp, y need to be moved up by one pitch
cty
cty========================================================================

      pure function transformed_y(transpar,y) result(trans_y)

        use overset, only : ov_pitch

        implicit none

        real(kind=wp), intent(in) :: transpar, y

        real(kind=wp) :: trans_y
   
        trans_y = y + transpar*ov_pitch

      end function ! transformed_y



cty========================================================================
cty     function particle_in_convex_quadrilateral
cty========================================================================
cty
cty     Description:
cty     Determine if a point lies inside a certain convex quadrilateral.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Point to be searched
cty       px, py: Coordinates of the vertices of the convex quadrilateral
cty
cty     Output arguments:
cty       is_in: True if (x,y) lies inside a certain convex quadrilateral
cty
cty     Notice: 
cty
cty========================================================================

      pure function particle_in_convex_quadrilateral
     &              (x,y,x1,y1,x2,y2,x3,y3,x4,y4) result(is_in)

        implicit none

        real(kind=wp), intent(in) :: x, y
        real(kind=wp), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4

        logical :: is_in 

        if ( (((x2-x1)*(y-y1)-(x-x1)*(y2-y1)) .gt. 0) .and.
     &       (((x3-x2)*(y-y2)-(x-x2)*(y3-y2)) .gt. 0) .and.
     &       (((x4-x3)*(y-y3)-(x-x3)*(y4-y3)) .gt. 0) .and.
     &       (((x1-x4)*(y-y4)-(x-x4)*(y1-y4)) .gt. 0) ) then
          is_in = .true.
        else
          is_in = .false.
        end if

      end function ! particle_in_convex_quadrilateral



cty========================================================================
cty     subroutine compute_k_particle
cty========================================================================
cty
cty     Description:
cty     Compute which the cell of the particle is located in. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       zpar: Coordinates of the particle
cty
cty     Output arguments:
cty       kpar: Local indices of the cell where the particle locates in
cty
cty     Notice:
cty
cty========================================================================

      pure subroutine compute_k_particle(zpar,kpar)

        use corevars, only : nzp

        implicit none 

        real(kind=wp), intent(in) :: zpar 
        integer, intent(out) :: kpar

        integer :: k

cty     Here we set kpar for both two block together,
cty     since they are the same in the spanwsie direction

        do k = 1, nzp
          if ( (zpar .gt. B_theta(k)) .and. 
     &         (zpar .le. B_theta(k+1)) ) then
            kpar = k
            return
          end if
        end do ! k

      end subroutine ! compute_k_particle



cty========================================================================
cty     subroutine update_k_particle
cty========================================================================
cty
cty     Description:
cty     Compute which the cell does the particle locate based on the second
cty     layer. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       zpar: Coordinates of the particle
cty       wpar: Velocities of the particle
cty       kpar_ori: Local indices of the cell where the particle 
cty                 originally locates in
cty
cty     Output arguments:
cty       kpar: Local indices of the cell where the particle 
cty             presently locates in
cty
cty     Notice:
cty
cty========================================================================

      pure subroutine update_k_particle(zpar,wpar,kpar_ori,kpar)

        use corevars, only : nzp

        implicit none 

        real(kind=wp), intent(in) :: zpar
        real(kind=wp), intent(in) :: wpar
        integer, intent(in) :: kpar_ori
        integer, intent(out) :: kpar

        integer :: k
        integer :: nz

cty     Set nz, which saves the number of points in each direction
cty     Periodical spanwise boundary condition
        nz = nzp + 1

cty     Set impossible initial values of kpar for debug
        kpar = -9999

cty     Here we set kpar for both two block together,
cty     since they are the same in the spanwsie direction

        if (wpar .gt. 0._wp) then

          do k = kpar_ori, nz-1
            if ( (zpar .gt. B_theta(k)) .and. 
     &           (zpar .le. B_theta(k+1)) ) then
              kpar = k
              return
            end if 
          end do 
          do k = 1, kpar_ori
            if ( (zpar .gt. B_theta(k)) .and. 
     &           (zpar .le. B_theta(k+1)) ) then
              kpar = k
              return
            end if 
          end do 

        else

          do k = kpar_ori, 1, -1
            if ( (zpar .gt. B_theta(k)) .and. 
     &           (zpar .le. B_theta(k+1)) ) then
              kpar = k
              return
            end if 
          end do
          do k = nz-1, kpar_ori, -1
            if ( (zpar .gt. B_theta(k)) .and. 
     &           (zpar .le. B_theta(k+1)) ) then
              kpar = k
              return
            end if 
          end do

        end if ! (wpar .gt. 0._wp)

        if (kpar .eq. (-9999)) then
          call compute_k_particle(zpar,kpar)
        end if

      end subroutine ! update_k_particle



cty========================================================================
cty     subroutine set_particles_velocity
cty========================================================================
cty
cty     Description:
cty     Set the initial velocity of particles. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty       refer to subroutine read_particle_parameters
cty       1. Initial_velocity_particle: how to set the initial velocity
cty          0: zero initial velocity
cty          1: follow the fluid
cty
cty========================================================================

      subroutine set_particles_velocity

        use mppvars, only : procid_fluid
        use, intrinsic :: ieee_arithmetic

        implicit none 
       
        integer :: num
        type(particle), pointer :: Pt_par
       
        if (Para_par%imode .eq. 0) then 
          if (Para_par%Initial_velocity_particle .eq. 1) then

            do num = 1, Npar_local
              Pt_par => Pointer_par(num)

              call interpolate_velocity_information(
     &             Pt_par%transpar,
     &             Pt_par%xpar,Pt_par%ypar,
     &             Pt_par%zpar,
     &             Pt_par%procpar,Pt_par%ipar,Pt_par%jpar,Pt_par%kpar,
     &             Pt_par%upar,Pt_par%vpar,Pt_par%wpar)


            end do ! num

          end if ! (Para_par%Initial_velocity_particle .eq. 1) 
        end if ! (Para_par%imode .eq. 0 )
        
        if (Npar_local .gt. 0) then
          write(*,'(A,I4,A9,I4,A4,I11,/,
     &              2(A9,F20.16),/,
     &              2(A9,F20.16),/,
     &              2(A9,F20.16))') 
     &            'The particle number in processor', procid_fluid, 
     &            'of block', block_proc(procid_fluid), 
     &            'is', Npar_local,
     &            'max u = ', 
     &             maxval(Pointer_par(1:Npar_local)%upar),
     &            'min u = ', 
     &             minval(Pointer_par(1:Npar_local)%upar),
     &            'max v = ', 
     &             maxval(Pointer_par(1:Npar_local)%vpar),
     &            'min v = ', 
     &             minval(Pointer_par(1:Npar_local)%vpar),
     &            'max w = ', 
     &             maxval(Pointer_par(1:Npar_local)%wpar),
     &            'min w = ', 
     &             minval(Pointer_par(1:Npar_local)%wpar)
        end if ! (Npar_local .gt. 0)

      end subroutine ! set_particles_velocity


 
cty========================================================================
cty     function weighted_average
cty========================================================================
cty
cty     Description:
cty     Weighted average flow varibles. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       var_number: Number of the varibles of the input arrays
cty       arrays: The source data we need
cty       i, j, k: i, j, k indices of the input arrays
cty       var_id: ID of the varibles of the input arrays
cty       c_ijk: weight coefficients includes
cty
cty     Output arguments:
cty       weighted_mean: flow varibles gotten by weighted average.
cty
cty     Notice:
cty
cty========================================================================

      pure function weighted_average(var_number,arrays,i,j,k,var_id,
     &                               c_ijk)
     &                       result (weighted_mean)

        use corevars, only : z_x, r_y, nxp, nyp, nzp, xhalo, yhalo

        implicit none

        integer, intent(in) :: var_number
        real(kind=wp), intent(in) :: arrays(1-xhalo:nxp+xhalo,
     &                                      1-yhalo:nyp+yhalo,
     &                                      nzp,var_number)
        integer, intent(in) :: i, j, k
        integer, intent(in) :: var_id
        real(kind=wp), intent(in) :: c_ijk(-2:2,-2:2,-2:2)

        integer :: l, m
        real(kind=wp) :: weighted_mean 
 
          weighted_mean = 
     &           c_ijk(0,0,0) * arrays(i  ,j  ,k           ,var_id)
     &         + c_ijk(1,0,0) * arrays(i+1,j  ,k           ,var_id)
     &         + c_ijk(0,1,0) * arrays(i  ,j+1,k           ,var_id) 
     &         + c_ijk(0,0,1) * arrays(i  ,j  ,mod(k,nzp)+1,var_id) 
     &         + c_ijk(0,1,1) * arrays(i  ,j+1,mod(k,nzp)+1,var_id) 
     &         + c_ijk(1,0,1) * arrays(i+1,j  ,mod(k,nzp)+1,var_id) 
     &         + c_ijk(1,1,0) * arrays(i+1,j+1,k           ,var_id) 
     &         + c_ijk(1,1,1) * arrays(i+1,j+1,mod(k,nzp)+1,var_id) 

      end function ! weighted_average



cty========================================================================
cty     subroutine interpolate_velocity_information
cty========================================================================
cty
cty     Description:
cty     Interpolate flow velocities at the location of particles. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar, zpar: Coordinates of the particle
cty       procpar: Processor id of the particle
cty       ipar, jpar, kpar: Local indices of the cell where the particle locates 
cty
cty     Output arguments:
cty       uflow, vflow, wflow: velocity varibles gotten by interpolation
cty
cty     Notice:
cty
cty========================================================================

      subroutine interpolate_velocity_information(
     &           transpar,xpar,ypar,zpar,
     &           procpar,ipar,jpar,kpar,
     &           uflow,vflow,wflow)

        use mppvars, only : procid_fluid, ioproc
        use corevars, only : qph, nzp

        implicit none

        real(kind=wp), intent(in) :: transpar, xpar, ypar, zpar
        integer, intent(in) :: procpar
        integer, intent(in) :: ipar, jpar, kpar
        real(kind=wp), intent(out) :: uflow, vflow, wflow
 
        integer :: i, j, k
        integer :: ibpar
        real(kind=wp) :: y
        real(kind=wp) :: c_ijk(-2:2,-2:2,-2:2)

cty     Set parameters associated with the present processor: ibpar
        ibpar = block_proc(procpar)

cty     Set indices of the cell where the particle locates
        i = ipar
        j = jpar
        k = kpar

cty     Set the original or transformed coordinates of the particle
        y = transformed_y(transpar,ypar)

cty     Calculate the interpolation weights
        c_ijk = interpolation_weights(xpar,y,zpar,ibpar,i,j,k)

        uflow = weighted_average(5,qph,i,j,k,2,c_ijk)
        vflow = weighted_average(5,qph,i,j,k,3,c_ijk)
        wflow = weighted_average(5,qph,i,j,k,4,c_ijk)

      end subroutine ! interpolate_velocity_information



cty========================================================================
cty     subroutine interpolate_flow_part_information
cty========================================================================
cty
cty     Description:
cty     Interpolate flow varibles at the location of particles. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar, zpar: Coordinates of the particle
cty       procpar: Processor id of the particle
cty       ipar, jpar, kpar: Local indices of the cell where the particle locates 
cty
cty     Output arguments:
cty       uflow, vflow, wflow, 
cty
cty     Notice:
cty
cty========================================================================

      subroutine interpolate_flow_part_information(transpar,
     &                                             xpar,ypar,zpar,
     &                                             procpar,
     &                                             ipar,jpar,kpar,
     &                                             uflow,vflow,wflow
     &                                             )

        use mppvars, only : procid_fluid, ioproc
        use corevars, only : qph, nzp

        implicit none

        real(kind=wp), intent(in) :: transpar, xpar, ypar, zpar
        integer, intent(in) :: procpar
        integer, intent(in) :: ipar, jpar, kpar
        real(kind=wp), intent(out) :: uflow, vflow, wflow
 
        integer :: i, j, k
        integer :: ibpar
        real(kind=wp) :: y
        real(kind=wp) :: c_ijk(-2:2,-2:2,-2:2)

cty     Set parameters associated with the present processor: ibpar
        ibpar = block_proc(procpar)

cty     Set indices of the cell where the particle locates
        i = ipar
        j = jpar
        k = kpar

cty     Set the original or transformed coordinates of the particle
        y = transformed_y(transpar,ypar)

cty     Calculate the interpolation weights
        c_ijk = interpolation_weights(xpar,y,zpar,ibpar,i,j,k)

        uflow = weighted_average(5,qph,i,j,k,2,c_ijk)
        vflow = weighted_average(5,qph,i,j,k,3,c_ijk)
        wflow = weighted_average(5,qph,i,j,k,4,c_ijk)

      end subroutine ! interpolate_flow_part_information



cty========================================================================
cty     subroutine interpolate_flow_information
cty========================================================================
cty
cty     Description:
cty     Interpolate flow varibles at the location of particles. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       xpar, ypar, zpar: Coordinates of the particle
cty       procpar: Processor id of the particle
cty       ipar, jpar, kpar: Local indices of the cell where the particle locates 
cty
cty     Output arguments:
cty       uflow, vflow, wflow, 
cty       rhoflow, muflow, 
cty       pflow_deriv_x, pflow_deriv_y, pflow_deriv_z, 
cty       vorflow_x, vorflow_y, vorflow_z: flow varibles gotten by interpolation
cty
cty     Notice:
cty
cty========================================================================

      subroutine interpolate_flow_information(transpar,xpar,ypar,zpar,
     &                                        procpar,ipar,jpar,kpar,
     &                                        c_ijk,
     &                                        uflow,vflow,wflow,
     &                                        rhoflow,muflow
     &                                        )

        use mppvars, only : procid_fluid, ioproc
        use corevars, only : qph, nzp
        use flow_properties, only : rrey, rsu

        implicit none

        real(kind=wp), intent(in) :: transpar, xpar, ypar, zpar
        integer, intent(in) :: procpar
        integer, intent(in) :: ipar, jpar, kpar
        real(kind=wp), intent(out) :: c_ijk(-2:2,-2:2,-2:2)
        real(kind=wp), intent(out) :: uflow, vflow, wflow
        real(kind=wp), intent(out) :: rhoflow, muflow
        real(kind=wp) :: Tflow

        integer :: i, j, k
        integer :: ibpar
        real(kind=wp) :: y

cty     Set parameters associated with the present processor: ibpar
        ibpar = block_proc(procpar)

cty     Set indices of the cell where the particle locates
        i = ipar
        j = jpar
        k = kpar

cty     Set the original or transformed coordinates of the particle
        y = transformed_y(transpar,ypar)

cty     Calculate the interpolation weights
        c_ijk = interpolation_weights(xpar,y,zpar,ibpar,i,j,k)

        uflow = weighted_average(5,qph,i,j,k,2,c_ijk)
        vflow = weighted_average(5,qph,i,j,k,3,c_ijk)
        wflow = weighted_average(5,qph,i,j,k,4,c_ijk)

        rhoflow = weighted_average(5,qph,i,j,k,1,c_ijk)
        Tflow = weighted_average(5,qph,i,j,k,5,c_ijk)

cty     After we obtain Tflow, we can calculate muflow with the Sutherlands law
        muflow = Tflow*dsqrt(Tflow)*(1._wp+rsu)/(Tflow+rsu)

      end subroutine ! interpolate_flow_information



cty========================================================================
cty     function interpolation_weights 
cty========================================================================
cty
cty     Description:
cty     Calculate the weights coefficient of the interpolation.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty       c_ijk: weight coefficients includes
cty              c_ijk(0,0,0) = c_ic0jc0kc0
cty              c_ijk(1,0,0) = c_ir1jc0kc0
cty              c_ijk(0,1,0) = c_ic0jr1kc0
cty              c_ijk(0,0,1) = c_ic0jc0kr1
cty              c_ijk(0,1,1) = c_ic0jr1kr1
cty              c_ijk(1,0,1) = c_ir1jc0kr1
cty              c_ijk(1,1,0) = c_ir1jr1kc0
cty              c_ijk(1,1,1) = c_ir1jr1kr1
cty
cty========================================================================

      pure function interpolation_weights(x,y,z,ibpar,ipar,jpar,kpar)
     &                            result (c_ijk)

        use mbvars, only : ov_ogrid_block, ov_bgrid_block
        use corevars, only : qph, nxp, nyp, nzp, z_x, r_y
     &                     , theta, thl
        use mppvars, only : procid_fluid
        use global_grid, only : grids
        use overset, only : ov_pitch
        use, intrinsic :: ieee_arithmetic

        implicit none

        real(kind=wp), intent(in) :: x, y, z
        integer, intent(in) :: ibpar, ipar, jpar, kpar      

        real(kind=wp) :: c_ijk(-2:2,-2:2,-2:2)
        real(kind=wp) :: c_i(-2:2)
        real(kind=wp) :: c_j(-2:2)
        real(kind=wp) :: c_k(-2:2)
        real(kind=wp) :: dz

        integer :: i, j, k
        integer :: l, m, n

        real(kind=wp) :: x_xm, y_xm,
     &                   x_yp, y_yp

cty     Set indices of the cell where the particle locates
        i = ipar
        j = jpar
        k = kpar
cty     Assume the particle is located in the box of
cty     i, i+1, j, j+1, k, k+1     

        c_i = 1._wp 
        c_j = 1._wp        
        c_k = 1._wp
cty     We need to avoid the situation that the 4th Langrange
cty     is invalid in the near-wall region!
        c_ijk = 0._wp


          if (ov_bgrid_block(ibpar)) then
            c_i(1) = (x-z_x(i,j))/(z_x(i+1,j+1)-z_x(i,j+1))
            c_i(0) = 1._wp - c_i(1)
cty         Notice that the y scale of each cell is equal for background grid
            c_j(1) = (y-r_y(i,j))/(r_y(i,j+1)-r_y(i,j))
            c_j(0) = 1._wp - c_j(1)

          end if

          if (ov_ogrid_block(ibpar)) then 
cty         Notice that each cell is not rectangle for O-type grid,
cty         but we assume that they are rectangle here,
cty         The grid point is actually locating at (z_x(i,j+1),r_y(i,j+1))
            call line_xmyp_intersection(x           ,y           ,
     &                                  z_x(i  ,j  ),r_y(i  ,j  ),
     &                                  z_x(i+1,j  ),r_y(i+1,j  ),
     &                                  z_x(i  ,j+1),r_y(i  ,j+1),
     &                                  z_x(i+1,j+1),r_y(i+1,j+1),
     &                                  x_xm        ,   y_xm     ,
     &                                  x_yp        ,   y_yp     )
            if ( abs(z_x(i+1,j+1) - z_x(i,j+1)) .gt.
     &           abs(r_y(i+1,j+1) - r_y(i,j+1)) ) then
              c_i(1) = (        x_yp - z_x(i  ,j+1)) /
     &                 (z_x(i+1,j+1) - z_x(i  ,j+1))
              c_i(0) = 1._wp - c_i(1)
              c_j(1) = (        y_xm - r_y(i  ,j  )) /
     &                 (r_y(i  ,j+1) - r_y(i  ,j  ))
              c_j(0) = 1._wp - c_j(1)
            else
              c_i(1) = (        y_yp - r_y(i  ,j+1)) /
     &                 (r_y(i+1,j+1) - r_y(i  ,j+1))
              c_i(0) = 1._wp - c_i(1)
              c_j(1) = (        x_xm - z_x(i  ,j  )) /
     &                 (z_x(i  ,j+1) - z_x(i  ,j  ))
              c_j(0) = 1._wp - c_j(1)
            end if

          end if

cty       Notice that the z scale of each cell is equal
          c_k(1) = (z-theta(k))/(thl/nzp)
          c_k(0) = 1._wp - c_k(1)

cty       Set coefficient of the interpolation formula
          c_ijk(0,0,0) = c_i(0) * c_j(0) * c_k(0) 
          c_ijk(1,0,0) = c_i(1) * c_j(0) * c_k(0)
          c_ijk(0,1,0) = c_i(0) * c_j(1) * c_k(0)
          c_ijk(0,0,1) = c_i(0) * c_j(0) * c_k(1)
          c_ijk(0,1,1) = c_i(0) * c_j(1) * c_k(1)
          c_ijk(1,0,1) = c_i(1) * c_j(0) * c_k(1)
          c_ijk(1,1,0) = c_i(1) * c_j(1) * c_k(0)
          c_ijk(1,1,1) = c_i(1) * c_j(1) * c_k(1)

      end function ! interpolation_weights



cty========================================================================
cty     subroutine fill_corner_halopoints_qph
cty========================================================================
cty
cty     Description:
cty     Fill the halo points of qph in the corner, which is needed for the 
cty     interpolation of the particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty       1. data_send from xmym_proc to data_recv in xpyp_proc:
cty                    qph(nxp+1:nxp+xhalo,nyp+1:nyp+yhalo,:,:)
cty          data_send from xpym_proc to data_recv in xmyp_proc:
cty                    qph(1-xhalo:0,nyp+1:nyp+yhalo,:,:)
cty          data_send from xpyp_proc to data_recv in xmym_proc:
cty                    qph(1-xhalo:0,1-yhalo:0,:,:)
cty          data_send from xmyp_proc to data_recv in xpym_proc:
cty                    qph(nxp+1:nxp+xhalo,1-yhalo:0,:,:)
cty
cty========================================================================

      subroutine fill_corner_halopoints_qph

        use mppvars, only : real_mp_type, MPI_comm_fluid,
     &                      MPI_STATUS_IGNORE, procid_fluid, 
     &                      ioproc, nullproc
        use corevars, only : qph, nxp, nyp, nzp, kmodw, nvar, 
     &                       xhalo, yhalo

        implicit none
        
        integer :: i, j, k, var
        integer :: l, data_number
        integer :: ierr 
        real(kind=wp), allocatable, dimension(:) :: data_send, data_recv
 
        integer, parameter :: tag = 8888

        call swap_mb(qph,nvar,kmodw,.false.)
!$acc update host(qph)

        data_number = xhalo*yhalo*nzp*nvar
        allocate(data_send(data_number))
        allocate(data_recv(data_number))

cty     data_send from xmym_proc: physical points 
        if (xmym_proc(procid_fluid) .ne. nullproc) then
          l = 0 
          do var = 1, nvar
            do k = 1, nzp
              do j = 1, yhalo 
                do i = 1, xhalo 
                  l = l + 1
                  data_send(l) = qph(i,j,k,var)
                end do
              end do
            end do
          end do
        end if
        call MPI_sendrecv(
     &       data_send,data_number,real_mp_type,
     &       xmym_proc(procid_fluid),tag,
     &       data_recv,data_number,real_mp_type,
     &       xpyp_proc(procid_fluid),tag,
     &       MPI_comm_fluid,MPI_STATUS_IGNORE,ierr)
cty     to data_recv in xpyp_proc: halo points 
        if (xpyp_proc(procid_fluid) .ne. nullproc) then 
          l = 0
          do var = 1, nvar
            do k = 1, nzp
              do j = nyp+1, nyp+yhalo
                do i = nxp+1, nxp+xhalo
                  l = l + 1
                  qph(i,j,k,var) = data_recv(l)
                end do
              end do
            end do
          end do  
        end if

cty     data_send from xpym_proc: physical points
        if (xpym_proc(procid_fluid) .ne. nullproc) then
          l = 0 
          do var = 1, nvar
            do k = 1, nzp
              do j = 1, yhalo 
                do i = nxp-xhalo+1, nxp 
                  l = l + 1
                  data_send(l) = qph(i,j,k,var)
                end do
              end do
            end do
          end do
        end if
        call MPI_sendrecv(
     &       data_send,data_number,real_mp_type,
     &       xpym_proc(procid_fluid),tag,
     &       data_recv,data_number,real_mp_type,
     &       xmyp_proc(procid_fluid),tag,
     &       MPI_comm_fluid,MPI_STATUS_IGNORE,ierr)
cty     to data_recv in xmyp_proc: halo points 
        if (xmyp_proc(procid_fluid) .ne. nullproc) then 
          l = 0
          do var = 1, nvar
            do k = 1, nzp
              do j = nyp+1, nyp+yhalo
                do i = 1-xhalo, 0
                  l = l + 1
                  qph(i,j,k,var) = data_recv(l)
                end do
              end do
            end do
          end do
        end if

cty     data_send from xpyp_proc: physical points 
        if (xpyp_proc(procid_fluid) .ne. nullproc) then
          l = 0 
          do var = 1, nvar
            do k = 1, nzp
              do j = nyp-yhalo+1, nyp
                do i = nxp-xhalo+1, nxp 
                  l = l + 1
                  data_send(l) = qph(i,j,k,var)
                end do
              end do
            end do
          end do
        end if
        call MPI_sendrecv(
     &       data_send,data_number,real_mp_type,
     &       xpyp_proc(procid_fluid),tag,
     &       data_recv,data_number,real_mp_type,
     &       xmym_proc(procid_fluid),tag,
     &       MPI_comm_fluid,MPI_STATUS_IGNORE,ierr)
cty     to data_recv in xmym_proc: halo points  
        if (xmym_proc(procid_fluid) .ne. nullproc) then 
          l = 0
          do var = 1, nvar
            do k = 1, nzp
              do j = 1-yhalo, 0
                do i = 1-xhalo, 0
                  l = l + 1
                  qph(i,j,k,var) = data_recv(l)
                end do
              end do
            end do
          end do
        end if

cty     data_send from xmyp_proc: physical points
        if (xmyp_proc(procid_fluid) .ne. nullproc) then
          l = 0 
          do var = 1, nvar
            do k = 1, nzp
              do j = nyp-yhalo+1, nyp
                do i = 1, xhalo
                  l = l + 1
                  data_send(l) = qph(i,j,k,var)
                end do
              end do
            end do
          end do
        end if
        call MPI_sendrecv(
     &       data_send,data_number,real_mp_type,
     &       xmyp_proc(procid_fluid),tag,
     &       data_recv,data_number,real_mp_type,
     &       xpym_proc(procid_fluid),tag,
     &       MPI_comm_fluid,MPI_STATUS_IGNORE,ierr)
cty     to data_recv in xpym_proc: halo points  
        if (xpym_proc(procid_fluid) .ne. nullproc) then 
          l = 0
          do var = 1, nvar
            do k = 1, nzp
              do j = 1-yhalo, 0
                do i = nxp+1, nxp+xhalo
                  l = l + 1
                  qph(i,j,k,var) = data_recv(l)
                end do
              end do
            end do
          end do
        end if

      end subroutine ! fill_corner_halopoints_qph



cty========================================================================
cty     subroutine init_adjacent_proc
cty========================================================================
cty
cty     Description:
cty     Initialize arrays used for the communication between processors and 
cty     processors in the adjacent region.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty       1. xm_proc: processor in xm direction of the present processor
cty          xp_proc: processor in xp direction of the present processor
cty          ym_proc: processor in ym direction of the present processor
cty          yp_proc: processor in yp direction of the present processor
cty
cty========================================================================

      subroutine init_adjacent_proc

        use mppvars, only : nprocs, ioproc, nullproc
        use mbvars, only : num_blocks, ov_ogrid_block, ov_bgrid_block
        use global_grid, only : grids

        implicit none

        integer :: p, ind_i, ind_j
        integer :: block_id, nx, ny

        allocate(xm_proc(0:nprocs-1))
        allocate(xp_proc(0:nprocs-1))
        allocate(ym_proc(0:nprocs-1))
        allocate(yp_proc(0:nprocs-1)) 
        xm_proc = nullproc 
        xp_proc = nullproc
        ym_proc = nullproc
        yp_proc = nullproc
        
        do p = 0, nprocs-1

          ind_i = ind_i_proc(p)
          ind_j = ind_j_proc(p)
          block_id = block_proc(p)
          nx = grids(block_id)%nproc_x
          ny = grids(block_id)%nproc_y

          if (ov_bgrid_block(block_id)) then

            if ( (nx .eq. 1) .and.
     &           (ny .eq. 1) ) then
            else if ( (nx .eq. 1) .and. 
     &                (ny .gt. 1) ) then

              if (ind_j .eq. 1) then
                ym_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ny     )
                yp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j+1) 
              else if (ind_j .eq. ny) then
                ym_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j-1)
                yp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,1      )
              else
                ym_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j-1)
                yp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j+1)
              end if ! (ind_j .eq. 1)

            else if ( (nx .gt. 1) .and. 
     &                (ny .eq. 1) ) then

              if (ind_i .eq. 1) then
                xp_proc(p) = grids(block_id)%
     &                       proc_map(ind_i+1,1      )
                ym_proc(p) = grids(block_id)%
     &                       proc_map(ind_i  ,1      )
                yp_proc(p) = grids(block_id)%
     &                       proc_map(ind_i  ,1      )
              else if (ind_i .eq. nx) then
                xm_proc(p) = grids(block_id)%
     &                       proc_map(ind_i-1,1      )
                ym_proc(p) = grids(block_id)%
     &                       proc_map(ind_i  ,1      )
                yp_proc(p) = grids(block_id)%
     &                       proc_map(ind_i  ,1      )
              else
                xm_proc(p) = grids(block_id)%
     &                       proc_map(ind_i-1,1      )
                xp_proc(p) = grids(block_id)%
     &                       proc_map(ind_i+1,1      )
                ym_proc(p) = grids(block_id)%
     &                       proc_map(ind_i  ,1      )
                yp_proc(p) = grids(block_id)%
     &                       proc_map(ind_i  ,1      )
              end if ! (ind_i .eq. 1)

            else if ( (nx .gt. 1) .and. 
     &                (ny .gt. 1) ) then

              if (ind_i .eq. 1) then 

                if (ind_j .eq. 1) then
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ny     )
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                else if (ind_j .eq. ny) then 
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,1      )
                else
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                end if ! (ind_j .eq. 1)

              else if (ind_i .eq. nx) then

                if (ind_j .eq. 1) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ny     )
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                else if (ind_j .eq. ny) then 
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,1      )
                else
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                end if ! (ind_j .eq. 1) 

              else

                if (ind_j .eq. 1) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ny     )
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                else if (ind_j .eq. ny) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,1      )
                else
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                end if ! (ind_j .eq. 1) 

              end if ! (ind_i .eq. 1)

            end if ! ( (nx .eq. 1) .and. (ny .eq. 1) )

          end if

          if (ov_ogrid_block(block_id)) then
cty       Notice xperi_block(block_id) .eq. .ture. here
cty       Notice yperi_block(block_id) .ne. .ture. here
            if ( (nx .eq. 1) .and. 
     &           (ny .eq. 1) ) then

              xm_proc(p) = grids(block_id)%
     &                     proc_map(1      ,1      )
              xp_proc(p) = grids(block_id)%
     &                     proc_map(1      ,1      )

            else if ( (nx .eq. 1) .and. 
     &                (ny .gt. 1) ) then

              if (ind_j .eq. 1) then
                xm_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j  )
                xp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j  )
                yp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j+1)
              else if (ind_j .eq. ny) then
                xm_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j  )
                xp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j  )
                ym_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j-1)
              else
                xm_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j  )
                xp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j  )
                ym_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j-1)
                yp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j+1)
              end if ! (ind_j .eq. 1)

            else if ( (nx .gt. 1) .and. 
     &                (ny .eq. 1) ) then

              if (ind_i .eq. 1) then
                xm_proc(p) = grids(block_id)%
     &                       proc_map(nx     ,ind_j  )
                xp_proc(p) = grids(block_id)%
     &                       proc_map(ind_i+1,ind_j  )
              else if (ind_i .eq. nx) then
                xm_proc(p) = grids(block_id)%
     &                       proc_map(ind_i-1,ind_j  )
                xp_proc(p) = grids(block_id)%
     &                       proc_map(1      ,ind_j  )
              else
                xm_proc(p) = grids(block_id)%
     &                       proc_map(ind_i-1,ind_j  )
                xp_proc(p) = grids(block_id)%
     &                       proc_map(ind_i+1,ind_j  )
              end if ! (ind_i .eq. 1)

            else if ( (nx .gt. 1) .and. 
     &                (ny .gt. 1) ) then
 
              if (ind_j .eq. 1) then

                if (ind_i .eq. 1) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(nx     ,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                else if (ind_i .eq. nx) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j  )
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                else
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                end if ! (ind_j .eq. 1)

              else if (ind_j .eq. ny) then

                if (ind_i .eq. 1) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(nx     ,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                else if (ind_i .eq. nx) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                else
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                end if ! (ind_j .eq. 1)

              else

                if (ind_i .eq. 1) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(nx     ,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                else if (ind_i .eq. nx) then
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                else
                  xm_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,ind_j  )
                  xp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,ind_j  )
                  ym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j-1)
                  yp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i  ,ind_j+1)
                end if ! (ind_j .eq. 1)

              end if ! (ind_j .eq. 1)
            
            end if ! ( (nx .gt. 1) .and. (ny .gt. 1) )

          end if ! (ov_bgrid_block(block_id)) 

        end do

cty     Check the meaning of xm_proc, xp_proc, ym_proc, yp_proc
        if (ioproc) then
          write(*,*) 'Adjacent processors:' 
          do p = 0, nprocs-1
            write(*,'(A,I4,A9,I4)') 
     &              'For the processor', p, 
     &              'of block', block_proc(p)
            write(*,'(4(A13,I5))') 
     &              'xm_proc = ', xm_proc(p),
     &              'xp_proc = ', xp_proc(p),
     &              'ym_proc = ', ym_proc(p),
     &              'yp_proc = ', yp_proc(p)
          end do
        end if

      end subroutine ! init_adjacent_proc



cty========================================================================
cty     subroutine init_corner_proc
cty========================================================================
cty
cty     Description:
cty     Initialize arrays used for the communication between processors and 
cty     processors in the corner region.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty       1. xmym_proc: processor in xm&ym direction of the present processor
cty          xpym_proc: processor in xp&ym direction of the present processor
cty          xpyp_proc: processor in xp&yp direction of the present processor
cty          xmyp_proc: processor in xm&yp direction of the present processor
cty
cty========================================================================

      subroutine init_corner_proc

        use mppvars, only : procid_fluid, nprocs, ioproc, nullproc
        use mbvars, only : num_blocks, ov_ogrid_block, ov_bgrid_block
        use global_grid, only : grids

        implicit none

        integer :: p, ind_i, ind_j
        integer :: block_id, nx, ny

        allocate(xmym_proc(0:nprocs-1))
        allocate(xpym_proc(0:nprocs-1))
        allocate(xpyp_proc(0:nprocs-1))
        allocate(xmyp_proc(0:nprocs-1)) 
        xmym_proc = nullproc 
        xpym_proc = nullproc
        xpyp_proc = nullproc
        xmyp_proc = nullproc
        
        do p = 0, nprocs-1

          ind_i = ind_i_proc(p)
          ind_j = ind_j_proc(p)
          block_id = block_proc(p)
          nx = grids(block_id)%nproc_x
          ny = grids(block_id)%nproc_y

          if (ov_bgrid_block(block_id)) then

            if ( (nx .eq. 1) .and. 
     &           (ny .eq. 1) ) then
            else if ( (nx .eq. 1) .and. 
     &                (ny .gt. 1) ) then

              if (ind_j .eq. 1) then
                xmym_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ny     )
                xpym_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ny     )
              else if (ind_j .eq. ny) then
                xmyp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,1      )
                xpyp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,1      )
              else
              end if

            else if ( (nx .gt. 1) .and. 
     &                (ny .eq. 1) ) then

              if (ind_i .eq. 1) then
                xpym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,1      )
                xpyp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,1      )
              else if (ind_i .eq. nx) then
                xmym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,1      )
                xmyp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,1      )
              else
                xmym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,1      )
                xmyp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i-1,1      )
                xpym_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,1      )
                xpyp_proc(p) = grids(block_id)%
     &                         proc_map(ind_i+1,1      )
              end if

            else if ( (nx .gt. 1) .and. 
     &                (ny .gt. 1) ) then

              if (ind_i .eq. 1) then 

                if (ind_j .eq. 1) then
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ny     )
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j+1)
                else if (ind_j .eq. ny) then 
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j-1)
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,1      )
                else
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j-1)
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j+1)
                end if ! (ind_j .eq. 1)

              else if (ind_i .eq. nx) then

                if (ind_j .eq. 1) then
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ny     )
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j+1)
                else if (ind_j .eq. ny) then 
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j-1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,1      )
                else
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j-1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j+1)
                end if ! (ind_j .eq. 1) 

              else

                if (ind_j .eq. 1) then
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ny     )
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ny     )
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j+1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j+1)
                else if (ind_j .eq. ny) then
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j-1)
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j-1)
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,1      )
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,1      )
                else
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j-1)
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j-1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j+1)
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j+1)
                end if ! (ind_j .eq. 1) 

              end if ! (ind_i .eq. 1)

            end if ! ( (nx .eq. 1) .and. (ny .eq. 1) )

          end if

          if (ov_ogrid_block(block_id)) then
cty       Notice xperi_block(block_id) .eq. .ture. here
cty       Notice yperi_block(block_id) .ne. .ture. here
            if ( (nx .eq. 1) .and. 
     &           (ny .eq. 1) ) then

            else if ( (nx .eq. 1) .and. 
     &                (ny .gt. 1) ) then

              if (ind_j .eq. 1) then
                xpyp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j+1)
                xmyp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j+1)
              else if (ind_j .eq. ny) then
                xpym_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j-1)
                xmym_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j-1)
              else
                xmym_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j-1)
                xpym_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j-1)
                xpyp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j+1)
                xmyp_proc(p) = grids(block_id)%
     &                         proc_map(1      ,ind_j+1)
              end if

            else if ( (nx .gt. 1) .and. 
     &                (ny .eq. 1) ) then

            else if ( (nx .gt. 1) .and. 
     &                (ny .gt. 1) ) then
 
              if (ind_j .eq. 1) then

                if (ind_i .eq. 1) then
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j+1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(nx     ,ind_j+1)
                else if (ind_i .eq. nx) then
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(1      ,ind_j+1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j+1)
                else
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j+1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j+1)
                end if ! (ind_j .eq. 1)

              else if (ind_j .eq. ny) then

                if (ind_i .eq. 1) then
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(nx     ,ind_j-1)
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j-1)
                else if (ind_i .eq. nx) then
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j-1)
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(1      ,ind_j-1)
                else
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j-1)
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j-1)
                end if ! (ind_j .eq. 1)

              else

                if (ind_i .eq. 1) then
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(nx     ,ind_j-1)
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j-1)
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j+1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(nx     ,ind_j+1)
                else if (ind_i .eq. nx) then
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j-1)
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(1      ,ind_j-1)
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(1      ,ind_j+1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j+1)
                else
                  xmym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j-1)
                  xpym_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j-1)
                  xpyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i+1,ind_j+1)
                  xmyp_proc(p) = grids(block_id)%
     &                           proc_map(ind_i-1,ind_j+1)
                end if ! (ind_j .eq. 1)

              end if ! (ind_j .eq. 1)
            
            end if ! ( (nx .gt. 1) .and. (ny .gt. 1) )

          end if ! (ov_bgrid_block(block_id)) 

        end do

cty     Check the meaning of xmym_proc, xpym_proc, xpyp_proc, xmyp_proc
        if (ioproc) then
          write(*,*) 'Corner processors:' 
          do p = 0, nprocs-1
            write(*,'(A,I4,A9,I4)') 
     &              'For the processor', p, 
     &              'of block', block_proc(p)
            write(*,'(4(A13,I5))') 
     &              'xmym_proc = ', xmym_proc(p),
     &              'xpym_proc = ', xpym_proc(p),
     &              'xpyp_proc = ', xpyp_proc(p),
     &              'xmyp_proc = ', xmyp_proc(p)
          end do
        end if


      end subroutine ! init_corner_proc



cty========================================================================
cty     subroutine init_nearby_proc_list 
cty========================================================================
cty
cty     Description:
cty     Initialize arrays used for the communication between processors and 
cty     processors in the nearby region.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty       xmym_proc_list(proc)%id = [xm_proc(proc), ym_proc(proc), xmym_proc(proc)]
cty       xpym_proc_list(proc)%id = [xp_proc(proc), ym_proc(proc), xpym_proc(proc)]
cty       xmyp_proc_list(proc)%id = [xm_proc(proc), yp_proc(proc), xmyp_proc(proc)]
cty       xpyp_proc_list(proc)%id = [xp_proc(proc), yp_proc(proc), xpyp_proc(proc)]
cty
cty       ymxm_proc_list(proc)%id = [ym_proc(proc), xm_proc(proc), xmym_proc(proc)]
cty       ypxm_proc_list(proc)%id = [yp_proc(proc), xm_proc(proc), xmyp_proc(proc)]
cty       ymxp_proc_list(proc)%id = [ym_proc(proc), xp_proc(proc), xpym_proc(proc)]
cty       ypxp_proc_list(proc)%id = [yp_proc(proc), xp_proc(proc), xpyp_proc(proc)]
cty
cty     Notice: 
cty
cty========================================================================

      subroutine init_nearby_proc_list

        use mppvars, only : nprocs, ioproc

        implicit none

        integer :: p

        allocate(xmym_proc_list(0:nprocs-1))
        allocate(xpym_proc_list(0:nprocs-1))
        allocate(xmyp_proc_list(0:nprocs-1)) 
        allocate(xpyp_proc_list(0:nprocs-1))
 
        allocate(ymxm_proc_list(0:nprocs-1))
        allocate(ypxm_proc_list(0:nprocs-1))
        allocate(ymxp_proc_list(0:nprocs-1)) 
        allocate(ypxp_proc_list(0:nprocs-1))

        do p = 0, nprocs-1

          xmym_proc_list(p)%id = [xm_proc(p), ym_proc(p), xmym_proc(p)]
          xpym_proc_list(p)%id = [xp_proc(p), ym_proc(p), xpym_proc(p)]
          xmyp_proc_list(p)%id = [xm_proc(p), yp_proc(p), xmyp_proc(p)]
          xpyp_proc_list(p)%id = [xp_proc(p), yp_proc(p), xpyp_proc(p)]

          ymxm_proc_list(p)%id = [ym_proc(p), xm_proc(p), xmym_proc(p)]
          ypxm_proc_list(p)%id = [yp_proc(p), xm_proc(p), xmyp_proc(p)]
          ymxp_proc_list(p)%id = [ym_proc(p), xp_proc(p), xpym_proc(p)]
          ypxp_proc_list(p)%id = [yp_proc(p), xp_proc(p), xpyp_proc(p)]

        end do


        do p = 0, nprocs-1

          call reduce_proc_list(xmym_proc_list(p)%id,
     &                          xmym_proc_list(p)%length)
          call reduce_proc_list(xpym_proc_list(p)%id,
     &                          xpym_proc_list(p)%length)
          call reduce_proc_list(xmyp_proc_list(p)%id,
     &                          xmyp_proc_list(p)%length)
          call reduce_proc_list(xpyp_proc_list(p)%id,
     &                          xpyp_proc_list(p)%length)

          call reduce_proc_list(ymxm_proc_list(p)%id,
     &                          ymxm_proc_list(p)%length)
          call reduce_proc_list(ypxm_proc_list(p)%id,
     &                          ypxm_proc_list(p)%length)
          call reduce_proc_list(ymxp_proc_list(p)%id,
     &                          ymxp_proc_list(p)%length)
          call reduce_proc_list(ypxp_proc_list(p)%id,
     &                          ypxp_proc_list(p)%length)

        end do


      end subroutine ! init_nearby_proc_list



cty========================================================================
cty     subroutine line_xmyp_intersection 
cty========================================================================
cty
cty     Description:
cty     Calculate the intersection of parallel line passing through the point 
cty     location (or the transformed location) of the particle with another line.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       x, y: Coordinates (or the transformed coordinates) of the particle
cty       x0, y0: 
cty       x1, y1: 
cty       x2, y2: 
cty       x3, y3: 
cty
cty     Output arguments:
cty       x_xm, y_xm: Coordinates of the point of intersection
cty       x_yp, y_yp: 
cty
cty     Notice: 
cty
cty========================================================================

      pure subroutine line_xmyp_intersection(x,y,x0,y0,
     &                                       x1,y1,x2,y2,x3,y3,
     &                                       x_xm,y_xm,x_yp,y_yp)

        implicit none

        real(kind=wp), intent(in) :: x, y
        real(kind=wp), intent(in) :: x0, y0
        real(kind=wp), intent(in) :: x1, y1
        real(kind=wp), intent(in) :: x2, y2
        real(kind=wp), intent(in) :: x3, y3
        real(kind=wp), intent(out) :: x_xm, y_xm
        real(kind=wp), intent(out) :: x_yp, y_yp
         
        real(kind=wp) :: l, m

        l = ((x -x0)*(y2-y0) - (y -y0)*(x2-x0)) /
     &      ((x2-x0)*(y2-y3) - (y2-y0)*(x2-x3))
        x_xm = x + l*(x2-x3)
        y_xm = y + l*(y2-y3)

        m = ((x -x2)*(y3-y2) - (y -y2)*(x3-x2)) /
     &      ((x3-x2)*(y2-y0) - (y3-y2)*(x2-x0))
        x_yp = x + m*(x2-x0)
        y_yp = y + m*(y2-y0)

      end subroutine ! line_xmyp_intersection 



cty========================================================================
cty     subroutine monitor_particles_data
cty========================================================================
cty
cty     Description:
cty     This routine monitor the collected particle data.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty       1. (1) nstart: The timestamp returned when the previous simulation
cty                      was completed
cty          (2) nsteps: The total number of calculation steps
cty          (3) nit: Counter of the time step when HiPSTAR is running
cty
cty========================================================================

      subroutine monitor_particles_data

        use corevars, only : nstart, nsteps, nit, time
        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_comm_fluid, 
     &                      ioid, procid_fluid, nprocs, ioproc
        use ompvars, only : start_master, end_master

        implicit none

        type(particle), pointer :: Pt_par         
        integer :: k, i, p
        character(21) :: idpar_name
        character(300) :: particle_file        
        logical :: lex


        if (Para_par%Monitor_particle_steps .gt. 0) then

          if ( (mod(nit-nstart-1,Para_par%Monitor_particle_steps).eq.0)
     &         .and. (nit .gt. nstart+1) ) then

            if (start_master()) then

cty           Write the particles data into the corresponding file
              do k = 1, Npar_local
                Pt_par => Pointer_par(k)

                if (Pt_par%is_monitor .eq. 1) then

cty               Set the file name of the particle_file
                  write(idpar_name,'(F9.0)') Pt_par%idpar
cty               Here we use 8 places to mark the total global id of the 
cty               particles, and fill the blanks with zeros.
                  do i = 1, len_trim(idpar_name)
                    if (idpar_name(i:i) .eq. '') then
                      idpar_name(i:i) = '0' 
                    else 
                      exit
                    end if
                  end do 
                  particle_file = 'Monitor_particles/'//
     &                            'Monitor_particle_'//
     &                             trim(idpar_name)//'bin'
                  particle_file = trim(adjustl(particle_file))

                  inquire(file=particle_file,exist=lex)

cty               Open the particle_file
                  if (lex) then
                    open(unit=procid_fluid+1000,file=particle_file,
     &                   form='unformatted',status='old',
     &                   action='write',position='append') 
                  else
                    open(unit=procid_fluid+1000,file=particle_file,
     &                   form='unformatted',status='new',
     &                   action='write',position='append') 
                  end if ! (lex)

cty               Write particle data
                  write(procid_fluid+1000) time,
     &                                     Pt_par%dpar,    
     &                                     Pt_par%rhopar,  
     &                                     Pt_par%transpar, 
     &                                     Pt_par%xpar,    
     &                                     Pt_par%ypar,    
     &                                     Pt_par%zpar,    
     &                                     Pt_par%upar,    
     &                                     Pt_par%vpar,    
     &                                     Pt_par%wpar,    
     &                                     Pt_par%axpar,    
     &                                     Pt_par%aypar,    
     &                                     Pt_par%azpar    
                  write(procid_fluid+1000) Pt_par%procpar,
     &                                     Pt_par%ipar, 
     &                                     Pt_par%jpar,
     &                                     Pt_par%kpar 

cty               Close the particle_file
                  close(procid_fluid+1000)

                end if ! (Pt_par%is_monitor .eq. 1)

              end do ! k

            end if ! start_master
            call end_master()

          end if ! (mod(nit-nstart-1,Para_par%Monitor_particle_steps).eq.0) 

        end if ! (Para_par%Monitor_particle_steps .gt. 0)         

      end subroutine ! monitor_particles_data



cty========================================================================
cty     subroutine collect_collision_particles_data
cty========================================================================
cty 
cty     Description:
cty     This routine save the data of particles which collide with the blade wall.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine collect_collision_particles_data

        use corevars, only : nstart, nsteps, nit, time
        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_status_size,
     &                      MPI_SUM, MPI_comm_fluid, ioid, procid_fluid,
     &                      nprocs, ioproc
        use ompvars, only : start_master, end_master

        implicit none

        integer :: k, i, p
        character(15) :: id_name
        character(300) :: particle_file, particle_file_temp 
        logical :: lex
        integer :: status(MPI_status_size)
        integer :: ierr
        character(700) :: exec_arg

        integer, allocatable, dimension(:) :: Npar_local_collision_proc
        integer :: Npar_collision

        if (start_master()) then
 
          allocate(Npar_local_collision_proc(0:nprocs-1))

cty       Collect Npar_local_collision in ioid, which is save in 
cty       Npar_local_collision_proc(procid_fluid)
          call MPI_gather
     &        (Npar_local_collision,1,MPI_INTEGER,
     &         Npar_local_collision_proc(procid_fluid),1,MPI_INTEGER,
     &         ioid,MPI_comm_fluid,ierr)

cty       Set the file name of the particle_file
          particle_file_temp = 'Collision_particles/'//
     &                         'Collision_particle_temp.bin'
          particle_file_temp = trim(adjustl(particle_file_temp))

          if (ioproc) then
            inquire(file=particle_file_temp,exist=lex)

cty         Open the particle_file
            if (lex) then
cty           It means that some particles collide with the blade wall before
              open(unit=99,file=particle_file_temp,
     &             form='unformatted',status='old',
     &             action='write',position='append') 
            else
cty           It means that no particles collide with the blade wall before
              open(unit=99,file=particle_file_temp,
     &             form='unformatted',status='new',
     &             action='write',position='append') 
            end if ! (lex)

cty         Write the particle data of io processor
            do k = 1, Npar_local_collision_proc(procid_fluid)
              write(99) time,
     &                  pardata_real_collision(1,k),
     &                  pardata_real_collision(2,k),
     &                  pardata_real_collision(3,k),
     &                  pardata_real_collision(4,k),
     &                  pardata_real_collision(5,k),
     &                  pardata_real_collision(6,k),
     &                  pardata_real_collision(7,k),
     &                  pardata_real_collision(8,k),
     &                  pardata_real_collision(9,k),
     &                  pardata_real_collision(10,k),
     &                  pardata_real_collision(11,k),
     &                  pardata_real_collision(12,k),
     &                  pardata_real_collision(13,k)
              write(99) pardata_integer_collision(1,k),
     &                  pardata_integer_collision(2,k),
     &                  pardata_integer_collision(3,k),
     &                  pardata_integer_collision(4,k)
            end do ! k

          end if ! (ioproc)

          if (.not. ioproc) then

            if (Npar_local_collision .ne. 0) then
cty           Send the particle data of other processors
              call MPI_send
     &            (pardata_real_collision,
     &             Nvar_par_real*Npar_local_collision,
     &             real_mp_type,ioid,procid_fluid,
     &             MPI_comm_fluid,ierr)
              call MPI_send
     &            (pardata_integer_collision,
     &            (Nvar_par_integer-1)*Npar_local_collision,
     &             MPI_INTEGER,ioid,procid_fluid,
     &             MPI_comm_fluid,ierr)
            end if ! (Npar_local_collision .ne. 0)
          else
            do p = 1, nprocs-1
              Npar_collision = Npar_local_collision_proc(p)   
              if (Npar_collision .ne. 0) then 
cty             Receive the particle data of other processors
                call MPI_recv
     &              (pardata_real_collision,
     &               Nvar_par_real*Npar_collision,
     &               real_mp_type,p,p,
     &               MPI_comm_fluid,status,ierr)
                call MPI_recv
     &              (pardata_integer_collision,
     &              (Nvar_par_integer-1)*Npar_collision,
     &               MPI_INTEGER,p,p,
     &               MPI_comm_fluid,status,ierr)
cty             Write the particle data of other processors
                do k = 1, Npar_collision
                  write(99) time,
     &                      pardata_real_collision(1,k),
     &                      pardata_real_collision(2,k),
     &                      pardata_real_collision(3,k),
     &                      pardata_real_collision(4,k),
     &                      pardata_real_collision(5,k),
     &                      pardata_real_collision(6,k),
     &                      pardata_real_collision(7,k),
     &                      pardata_real_collision(8,k),
     &                      pardata_real_collision(9,k),
     &                      pardata_real_collision(10,k),
     &                      pardata_real_collision(11,k),
     &                      pardata_real_collision(12,k),
     &                      pardata_real_collision(13,k)
                  write(99) pardata_integer_collision(1,k),
     &                      pardata_integer_collision(2,k),
     &                      pardata_integer_collision(3,k),
     &                      pardata_integer_collision(4,k)
                end do ! k
              end if ! (Npar_collision .ne. 0)
            end do ! p

          end if ! (.not. ioproc)

          if (ioproc) then
cty         Close the particle_file
            close(99)
          end if

          deallocate(Npar_local_collision_proc)

cty       Rename the temporary files, and them collect the particles
cty       that collide with the blade wall during the time until 
cty       (mod(nit-nstart-1,Para_par%Collision_particle_steps).eq.0)
          if ( (mod(nit-nstart-1,Para_par%Collision_particle_steps)
     &         .eq.0)
     &         .and. (nit .gt. nstart+1) ) then

            if (ioproc) then
              inquire(file=particle_file_temp,exist=lex)

              if (lex) then 
                write(id_name,'(I15)') nit-1

                particle_file = 'Collision_particles/'//
     &                          'Collision_particle_'//
     &                          trim(adjustl(id_name))//'.bin'
                particle_file = trim(adjustl(particle_file))
          
                write(exec_arg,'(A,A,A,A)') 'mv ', particle_file_temp, 
     &                                      ' ', particle_file
                call system(exec_arg)
              end if ! (lex)

            end if ! (ioproc)

          end if

cty       It is important to make them zero to conduct the collection
cty       in the next timestep
          Npar_local_collision = 0
          pardata_real_collision = 0._wp
          pardata_integer_collision = 0

        end if ! start_master
        call end_master()

      end subroutine ! collect_collision_particles_data



cty========================================================================
cty     subroutine capture_particles_number
cty========================================================================
cty
cty     Description:
cty     This routine capture the total particles number.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine capture_particles_number(time,Npar_global)

        use corevars, only : nstart, nit
        use mppvars, only : ioproc
        use ompvars, only : start_master, end_master

        implicit none

        real(kind=wp), intent(in) :: time
        integer, intent(in) :: Npar_global

        character(300) :: particle_file
        logical :: lex

        if (Para_par%Capture_particle_steps .gt. 0) then

          if ( (mod(nit-nstart,Para_par%Capture_particle_steps).eq.0)
     &         .and. (nit .gt. nstart) ) then

            if (start_master()) then

              if (ioproc) then
cty             Set the file name of the particle_file
                particle_file = 'Npar_global.csv'
                particle_file = trim(adjustl(particle_file))

                inquire(file=particle_file,exist=lex)

cty             Open the particle_file
                if (lex) then
                  open(unit=99,file=particle_file,
     &                 form='formatted',status='old',
     &                 action='write',position='append') 
                else
                  open(unit=99,file=particle_file,
     &                 form='formatted',status='new',
     &                 action='write',position='append') 
                end if ! (lex)

cty             Write the particle data of io processor
                write(99,'(F25.16,I11)') time, Npar_global 

cty             Close the particle_file
                close(99)

              end if

            end if ! start_master
            call end_master()

          end if ! (mod(nit-nstart-1,Para_par%Capture_particle).eq.0) 

        end if ! (Para_par%Capture_particle .gt. 0) 

      end subroutine ! capture_particles_number



cty========================================================================
cty     subroutine capture_particles_data
cty========================================================================
cty
cty     Description:
cty     This routine capture the collected particle data.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine capture_particles_data

        use corevars, only : nstart, nsteps, nit, time
        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_status_size,
     &                      MPI_comm_fluid, ioid, procid_fluid, nprocs,
     &                      ioproc
        use ompvars, only : start_master, end_master

        implicit none

        type(particle), pointer :: Pt_par         
        integer :: k, i, p
        character(15) :: id_name
        character(300) :: particle_file        
        logical :: lex
        integer :: status(MPI_status_size)
        integer :: ierr
       
        real(kind=wp), allocatable, dimension(:,:) :: pardata_real 
        integer, allocatable, dimension(:,:) :: pardata_integer
        integer, allocatable, dimension(:) :: Npar_local_proc
        integer :: Npar


        if (Para_par%Capture_particle_steps .gt. 0) then

          if ( (mod(nit-nstart-1,Para_par%Capture_particle_steps).eq.0)
     &         .and. (nit .gt. nstart+1) ) then

            if (start_master()) then

              allocate(Npar_local_proc(0:nprocs-1))
cty           Collect Npar_local in ioid, which is save in Npar_local_proc(procid_fluid)
              call MPI_gather
     &            (Npar_local,1,MPI_INTEGER,
     &             Npar_local_proc(procid_fluid),1,MPI_INTEGER,
     &             ioid,MPI_comm_fluid,ierr)

              allocate(pardata_real(Nvar_par_real,Npar_Max))
              allocate(pardata_integer(Nvar_par_integer-1,Npar_Max))
              do k = 1, Npar_local
                Pt_par => Pointer_par(k)
                call transfer_real_data(Pt_par,pardata_real(:,k))
                call transfer_integer_data(Pt_par,pardata_integer(:,k))
              end do

cty           Set the file name of the particle_file
              write(id_name,'(I15)') nit-1
cty           Here we use 15 places to mark and fill the blanks with zeros.
              particle_file = 'Capture_particles/'//
     &                        'Capture_particle_'//
     &                         trim(adjustl(id_name))//'.bin'
              particle_file = trim(adjustl(particle_file))

              if (ioproc) then
cty             Open the particle_file
                open(unit=99,file=particle_file,
     &               form='unformatted',status='replace',
     &               action='write') 
                write(99) time
                write(99) sum(Npar_local_proc)
cty             Write the particle data of io processor
                do k = 1, Npar_local_proc(procid_fluid)
                  write(99) pardata_real(1,k),
     &                      pardata_real(2,k),
     &                      pardata_real(3,k),
     &                      pardata_real(4,k),
     &                      pardata_real(5,k),
     &                      pardata_real(6,k),
     &                      pardata_real(7,k),
     &                      pardata_real(8,k),
     &                      pardata_real(9,k),
     &                      pardata_real(10,k),
     &                      pardata_real(11,k),
     &                      pardata_real(12,k),
     &                      pardata_real(13,k)
                  write(99) pardata_integer(1,k),
     &                      pardata_integer(2,k),
     &                      pardata_integer(3,k),
     &                      pardata_integer(4,k)
                end do ! k
              end if

              if (.not. ioproc) then
                if (Npar_local .ne. 0) then
cty               Send the particle data of other processors
                  call MPI_send
     &                (pardata_real,Nvar_par_real*Npar_local,
     &                 real_mp_type,ioid,procid_fluid,
     &                 MPI_comm_fluid,ierr)
                  call MPI_send
     &                (pardata_integer,(Nvar_par_integer-1)*Npar_local,
     &                 MPI_INTEGER,ioid,procid_fluid,
     &                 MPI_comm_fluid,ierr)
                end if ! (Npar_local .ne. 0)
              else
                do p = 1, nprocs-1
                  Npar = Npar_local_proc(p)   
                  if (Npar .ne. 0) then 
cty                 Receive the particle data of other processors
                    call MPI_recv
     &                  (pardata_real,Nvar_par_real*Npar,
     &                   real_mp_type,p,p,
     &                   MPI_comm_fluid,status,ierr)
                    call MPI_recv
     &                  (pardata_integer,(Nvar_par_integer-1)*Npar,
     &                   MPI_INTEGER,p,p,
     &                   MPI_comm_fluid,status,ierr)
cty                 Write the particle data of other processors
                    do k = 1, Npar
                      write(99) pardata_real(1,k),
     &                          pardata_real(2,k),
     &                          pardata_real(3,k),
     &                          pardata_real(4,k),
     &                          pardata_real(5,k),
     &                          pardata_real(6,k),
     &                          pardata_real(7,k),
     &                          pardata_real(8,k),
     &                          pardata_real(9,k),
     &                          pardata_real(10,k),
     &                          pardata_real(11,k),
     &                          pardata_real(12,k),
     &                          pardata_real(13,k)
                      write(99) pardata_integer(1,k),
     &                          pardata_integer(2,k),
     &                          pardata_integer(3,k),
     &                          pardata_integer(4,k)
                    end do ! k
                  end if ! (Npar .ne. 0)
                end do ! p
              end if

              if (ioproc) then
cty             Close the particle_file
                close(99)
              end if
              
              deallocate(Npar_local_proc)
              deallocate(pardata_real)
              deallocate(pardata_integer)

            end if ! start_master
            call end_master()

          end if ! (mod(nit-nstart-1,Para_par%Capture_particle).eq.0) 

        end if ! (Para_par%Capture_particle .gt. 0)         

      end subroutine ! capture_particles_data



cty========================================================================
cty     subroutine transfer_real_data
cty========================================================================
cty
cty     Description:
cty     Copy all the real(kind=wp) data of the particles into the arrays pardata.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       Pt_par: Pointer to the particle
cty
cty     Output arguments:
cty       pardata_real: Arrays contains the real(kind=wp) data of the particles
cty
cty     Notice: 
cty
cty========================================================================

      pure subroutine transfer_real_data(Pt_par,pardata_real)
       
        implicit none 

        type(particle), intent(in) :: Pt_par
        real(kind=wp), intent(out) :: pardata_real(Nvar_par_real)

        pardata_real(1) = Pt_par%idpar 
        pardata_real(2) = Pt_par%dpar
        pardata_real(3) = Pt_par%rhopar 
        pardata_real(4) = Pt_par%transpar
        pardata_real(5) = Pt_par%xpar 
        pardata_real(6) = Pt_par%ypar
        pardata_real(7) = Pt_par%zpar
        pardata_real(8) = Pt_par%upar  
        pardata_real(9) = Pt_par%vpar
        pardata_real(10) = Pt_par%wpar
        pardata_real(11) = Pt_par%axpar 
        pardata_real(12) = Pt_par%aypar
        pardata_real(13) = Pt_par%azpar

      end subroutine ! transfer_real_data



cty========================================================================
cty     subroutine transfer_integer_data
cty========================================================================
cty
cty     Description:
cty     Copy the integer data of the particles except the is_monitor  into 
cty     the arrays pardata.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       Pt_par: Pointer to the particle
cty
cty     Output arguments:
cty       pardata_integer: Arrays contains the integer data of the particles
cty
cty     Notice: 
cty
cty========================================================================

      pure subroutine transfer_integer_data(Pt_par,pardata_integer)
       
        implicit none 

        type(particle), intent(in) :: Pt_par
        integer, intent(out) :: pardata_integer(Nvar_par_integer-1)

        pardata_integer(1) = Pt_par%procpar
        pardata_integer(2) = Pt_par%ipar
        pardata_integer(3) = Pt_par%jpar
        pardata_integer(4) = Pt_par%kpar

      end subroutine ! transfer_integer_data



cty========================================================================
cty     subroutine store_particles_data
cty========================================================================
cty
cty     Description:
cty     This routine store the collected particle data for the subsequent 
cty     particles simulation.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine store_particles_data

        use corevars, only : nstart, nsteps, nit, time
        use mppvars, only : procid_fluid, ioproc
        use ompvars, only : start_master, end_master

        implicit none

        type(particle), pointer :: Pt_par         
        integer :: k, i
        character(3) :: idproc_name 
        character(300) :: particle_file        


        if (Para_par%Store_particle_steps .gt. 0) then

          if ( (mod(nit-nstart-1,Para_par%Store_particle_steps).eq.0)
     &         .and. (nit .gt. nstart+1) ) then

            if (start_master()) then

cty           Write the particles data into the same file for continued simulation
cty           Set the file name of the particle_file
              write(idproc_name,'(I2)') procid_fluid
cty           Here we use 2 places to mark the id of the processors, 
cty           and fill the blanks with zeros.
              do i = 1, len_trim(idproc_name)
                if (idproc_name(i:i) .eq. '') then
                  idproc_name(i:i) = '0' 
                else 
                  exit
                end if
              end do 
              particle_file = 'Continue_particles/'//
     &                        'PARTICLES_INPUT_CONTINUE'//
     &                        trim(idproc_name)//'.bin'
              particle_file = trim(adjustl(particle_file))

              open(unit=procid_fluid+1000,file=particle_file,
     &             form='unformatted',status='replace',
     &             action='write') 

cty           Write particle data
              write(procid_fluid+1000) Npar_local
              do k = 1, Npar_local
                Pt_par => Pointer_par(k)
                write(procid_fluid+1000) Pt_par%idpar,
     &                                   Pt_par%dpar,    
     &                                   Pt_par%rhopar,  
     &                                   Pt_par%transpar, 
     &                                   Pt_par%xpar,    
     &                                   Pt_par%ypar,    
     &                                   Pt_par%zpar,    
     &                                   Pt_par%upar,    
     &                                   Pt_par%vpar,    
     &                                   Pt_par%wpar,    
     &                                   Pt_par%axpar,    
     &                                   Pt_par%aypar,    
     &                                   Pt_par%azpar 
                write(procid_fluid+1000) Pt_par%procpar,
     &                                   Pt_par%ipar, 
     &                                   Pt_par%jpar,
     &                                   Pt_par%kpar, 
     &                                   Pt_par%is_monitor
              end do ! k
cty           Close the particle_file
              close(procid_fluid+1000)

            end if ! start_master
            call end_master()

          end if ! (mod(nit-nstart-1,Para_par%Store_particle_steps).eq.0) 

        end if ! (Para_par%Store_particle_steps .gt. 0)         

      end subroutine ! store_particles_data



cty========================================================================
cty     subroutine timeadvance_redistribute_particles 
cty========================================================================
cty
cty     Description:
cty     Time advances particles according to Newton's third law. Besides, 
cty     in this subroutine, we also need to redistribute particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine timeadvance_redistribute_particles

        use mbvars, only : num_blocks
     &                   , ov_bgrid_block, ov_ogrid_block
        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_SUM, 
     &                      MPI_comm_fluid, MPI_STATUS_IGNORE,
     &                      ioid, procid_fluid, nprocs, ioproc, nullproc
        use corevars, only : time, dt, qph, nit, nstart, nxp, nyp, nzp, 
     &                       nvar, xhalo, yhalo
        use global_grid, only : grids

        implicit none

        real(kind=wp) :: c_ijk(-2:2,-2:2,-2:2)
        real(kind=wp) :: Fpar(3)

        type(particle), pointer :: Pt_par         
        integer :: k, kdel
        logical :: is_in
        logical :: is_in_temp
        real(kind=wp) :: transpar
        
        real(kind=wp) :: uflow, vflow, wflow
        real(kind=wp) :: rhoflow, muflow

        integer :: i, j
        integer :: which_block

        integer :: p
        integer :: B_proc_list(3)
        integer :: B_proc_list_length
        integer :: O_proc_list(3)
        integer :: O_proc_list_length

        integer, allocatable :: Npar_out(:), Npar_in(:)
cty     Npar_out(0:nprocs-1)
cty     number of particles flow out of the present processor to processors 0:nprocs-1
cty     Npar_in(0:nprocs-1) 
cty     number of particles flow into the present processor from processors 0:nprocs-1
        integer :: Npar_out_total, Npar_in_total
cty     Npar_out_total: Total number of particles that flow out of the original processor
cty     Npar_in_total: Total number of particles that flow into the new processor
        real(kind=wp), allocatable :: pardata_send_real(:,:,:),
     &                                pardata_recv_real(:,:,:)
cty     pardata_send_real(Nvar_par_real,Npar_send_Max,0:nprocs-1)
cty     save the real(kind=wp) data corresponding to the particle in Npar_out(0:nprocs-1)
cty     pardata_recv_real(Nvar_par_real,Npar_send_Max,0:nprocs-1)
cty     save the real(kind=wp) data corresponding to the particle in Npar_in(0:nprocs-1)
        integer, allocatable :: pardata_send_integer(:,:,:),
     &                          pardata_recv_integer(:,:,:)
cty     pardata_send_integer(Nvar_par_integer,Npar_send_Max,0:nprocs-1)
cty     save the integer data corresponding to the particle in Npar_out(0:nprocs-1)
cty     pardata_recv_integer(Nvar_par_integer,Npar_send_Max,0:nprocs-1)
cty     save the integer data corresponding to the particle in Npar_in(0:nprocs-1)

        integer :: ibpar
        integer :: ind_i, ind_j
        integer, parameter :: tag = 8888
        integer :: ierr 
        integer :: procid_send, procid_recv
        integer :: Npar_local_settled, Npar_del_settled

        real(kind=wp) :: xpar_ori, ypar_ori
     &                 , zpar_ori        
        real(kind=wp) :: upar_ori, vpar_ori
     &                 , wpar_ori        
        integer :: ipar_ori, jpar_ori, kpar_ori        
cty     We need to make sure that the O-type grid is clockwise rotation
        real(kind=wp) :: vel_tangen, vel_normal
        integer :: ipar_predict
        integer :: procpar_mapping
        integer :: ipar_mapping, jpar_mapping
        logical :: is_in_mapping
        integer :: nx, ny
        real(kind=wp) :: y 
        integer :: start_i, end_i, step_i
        integer :: start_j, end_j, step_j

        kdel = 0
cty     mark_block_id = 0: represents the particle still exists in
cty                        the original processor
cty     mark_block_id = 999: represents the non-physical particle
cty                          the corresponding output file will be
cty                          moved to other folder
cty     mark_block_id) = -999: represents the particle impact with
cty                            the wall blade, because these data 
cty                            need to be saved before compress the 
cty                            data, we need to handle specially
cty     mark_block_id = 1: represents the particle that flow into 
cty                        other processors of the background grid
cty     mark_block_id = 2: represents the particle that flow into 
cty                        other processors of the O-type grid
        mark_block_id = 0

        do k = 1, Npar_local
          Pt_par => Pointer_par(k)

          call compute_particles_force(Pointer_par(k)
     &                                ,c_ijk
     &                                ,Fpar
     &                                 )


          Pt_par%axpar = Fpar(1)
          Pt_par%aypar = Fpar(2)
          Pt_par%azpar = Fpar(3)

cty       Time advances: update the coordinate and velocity of the particle
          Pt_par%xpar = Pt_par%xpar + Pt_par%upar*dt
     &                              + 0.5_wp*Pt_par%axpar*dt**2._wp
          Pt_par%ypar = Pt_par%ypar + Pt_par%vpar*dt
     &                              + 0.5_wp*Pt_par%aypar*dt**2._wp
          Pt_par%zpar = Pt_par%zpar + Pt_par%wpar*dt 
     &                              + 0.5_wp*Pt_par%azpar*dt**2._wp

cty       The coordinate of the particle need to be modified with the
cty       periodical boundary condition.
          call transform_particles_xyz_coordinates
     &         (Pt_par%xpar,Pt_par%ypar,Pt_par%zpar)


cty       Because spanwise direction is homogeneous, thus we can determine 
cty       Pt_par%kpar in the initial moment.
          kpar_ori = Pt_par%kpar
          call update_k_particle
     &        (Pt_par%zpar,Pt_par%wpar,kpar_ori,Pt_par%kpar)

cty       We need to carefully settle the particles that flow into other processors, 
cty       and the corresponding id are recorded in idpar_proc_del(:)

cty       Set parameters associated with the present processor: ibpar, ind_i, ind_j
          ibpar = block_proc(procid_fluid)
          ind_i = ind_i_proc(procid_fluid)
          ind_j = ind_j_proc(procid_fluid)

          which_block = 2

cty       1. If the particle originally locates in the background grid, 
          if (ov_bgrid_block(ibpar)) then

cty         (1): Flow out of the background grid
            if (.not. particle_in_block
     &                (Pt_par%xpar,Pt_par%ypar,ibpar)
     &         ) then 

              kdel = kdel + 1
              idpar_proc_del(kdel) = k
              mark_block_id(k) = 999

cty         (2): Flows into the O-type grid
            else if (particle_in_outer_surface
     &               (Pt_par%xpar,Pt_par%ypar)) then 

cty           The inner_surface is the blade surface
              kdel = kdel + 1
              idpar_proc_del(kdel) = k
              mark_block_id(k) = 2

cty         (3): If the particle also doesn't flow into the O-type grid

            else


              call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                  procid_fluid,is_in,transpar)

cty           (a): Still exists in the present grid
              if (is_in) then

cty             Update transpar, ipar, jpar and kpar
                Pt_par%transpar = transpar
                ipar_ori = Pt_par%ipar
                jpar_ori = Pt_par%jpar
                call update_B_grid_ij_particle
     &              (Pt_par%transpar,Pt_par%xpar,Pt_par%ypar,
     &               Pt_par%upar,Pt_par%vpar,
     &               Pt_par%procpar,procid_fluid,
     &               ipar_ori,jpar_ori,
     &               Pt_par%ipar,Pt_par%jpar)

cty           (b): Flow into other subblock of the background grid
              else

                kdel = kdel + 1
                idpar_proc_del(kdel) = k
                mark_block_id(k) = 1

              end if ! (is_in)

            end if ! (.not. particle_in_block(Pt_par%xpar,Pt_par%ypar,ibpar))

          end if ! (ov_bgrid_block(ibpar))

cty       2. If the particle originally locates in the O-type grid
          if (ov_ogrid_block(ibpar)) then

cty         (1): Flow into the blade surface
            if (particle_in_inner_surface
     &          (Pt_par%xpar,Pt_par%ypar)
     &         ) then    
 
              Pt_par%is_monitor = -Pt_par%is_monitor
cty           The inner_surface is the blade surface
              call bounce_particles(Pt_par%transpar,Pt_par%ipar
     &                             ,Pt_par%xpar,Pt_par%ypar
     &                             ,Pt_par%upar,Pt_par%vpar
     &                              )
              call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                  procid_fluid,is_in,transpar)

cty           (a): Still exists in the present grid
              if (is_in) then

cty             Update transpar, ipar, jpar and kpar
                Pt_par%transpar = transpar
                ipar_ori = Pt_par%ipar
                jpar_ori = Pt_par%jpar
                ipar_predict = ibeg_proc(Pt_par%procpar)+ipar_ori-1
                vel_tangen = Pt_par%upar*tangen_x_global(ipar_predict)
     &                     + Pt_par%vpar*tangen_y_global(ipar_predict)
                vel_normal = Pt_par%upar*normal_x_global(ipar_predict)
     &                     + Pt_par%vpar*normal_y_global(ipar_predict)
                call update_O_grid_ij_particle_prediction
     &              (Pt_par%transpar,
     &               Pt_par%xpar,Pt_par%ypar,
     &               vel_tangen,vel_normal,
     &               Pt_par%procpar,procid_fluid,
     &               ipar_ori,jpar_ori,
     &               Pt_par%ipar,Pt_par%jpar)
                Pt_par%dpar = -Pt_par%dpar

cty           (b): Flow into other subblock of the O-type grid
              else

cty             For the particle which bounces with the blade surface, 
cty             we should not to change its velocity when exchanges 
cty             the particles between different processors.
                Pt_par%dpar = -Pt_par%dpar
                kdel = kdel + 1
                idpar_proc_del(kdel) = k
                mark_block_id(k) = 2

              end if ! (is_in)

cty         (2): If the particle still locates in the O-type grid
            else if (particle_in_outer_surface
     &               (Pt_par%xpar,Pt_par%ypar)
     &              ) then

              call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                  procid_fluid,is_in,transpar)

cty           (a): Still exists in the present grid
              if (is_in) then

cty             Update transpar, ipar, jpar and kpar
                Pt_par%transpar = transpar
                ipar_ori = Pt_par%ipar
                jpar_ori = Pt_par%jpar
                ipar_predict = ibeg_proc(Pt_par%procpar)+ipar_ori-1
                vel_tangen = Pt_par%upar*tangen_x_global(ipar_predict)
     &                     + Pt_par%vpar*tangen_y_global(ipar_predict)
                vel_normal = Pt_par%upar*normal_x_global(ipar_predict)
     &                     + Pt_par%vpar*normal_y_global(ipar_predict)
                call update_O_grid_ij_particle_prediction
     &              (Pt_par%transpar,
     &               Pt_par%xpar,Pt_par%ypar,
     &               vel_tangen,vel_normal,
     &               Pt_par%procpar,procid_fluid,
     &               ipar_ori,jpar_ori,
     &               Pt_par%ipar,Pt_par%jpar)

cty           (b): Flow into other subblock of the O-type grid
              else

                kdel = kdel + 1
                idpar_proc_del(kdel) = k
                mark_block_id(k) = 2

              end if ! (is_in)

cty         (3): Flow out of the O-type grid and into the background grid
            else

              kdel = kdel + 1
              idpar_proc_del(kdel) = k
              mark_block_id(k) = 1

            end if ! (particle_in_inner_surface(Pt_par%xpar,Pt_par%ypar))

          end if ! (ov_ogrid_block(ibpar))

        end do ! k = 1, Npar_local

        do k = 1, Npar_local
          Pt_par => Pointer_par(k)

          if (Pt_par%is_monitor .lt. 0) then 

            Npar_local_collision = Npar_local_collision + 1
            call transfer_real_data(Pt_par,
     &           pardata_real_collision(:,Npar_local_collision))
            call transfer_integer_data(Pt_par,
     &           pardata_integer_collision(:,Npar_local_collision))

            Pt_par%is_monitor = -Pt_par%is_monitor
          end if 

        end do ! k = 1, Npar_local

cty     Npar_del is the number of particles which flow out of the present
cty     processor at each time step.
        Npar_del = kdel

cty     For the particles which flow out of the present processor, we need 
cty     to carefully identify which processor does the particle flow into, 
cty     thus we obtain Npar_out, pardata_send_real and pardata_send_integer.
cty     Furthermore, exchange the data with 1 communication.
        allocate(Npar_out(0:nprocs-1))
        allocate(Npar_in(0:nprocs-1))
        Npar_out = 0
        Npar_in = 0
       
        allocate(pardata_send_real
     &           (Nvar_par_real,Npar_send_Max,0:nprocs-1))
        allocate(pardata_recv_real
     &           (Nvar_par_real,Npar_send_Max,0:nprocs-1))
        allocate(pardata_send_integer
     &           (Nvar_par_integer,Npar_send_Max,0:nprocs-1))
        allocate(pardata_recv_integer
     &           (Nvar_par_integer,Npar_send_Max,0:nprocs-1))

        do k = 1, Npar_del
          Pt_par => Pointer_par(idpar_proc_del(k))

cty       1. If the particle originally locates in the background grid, 
          if (ov_bgrid_block(ibpar)) then

            if (mark_block_id(idpar_proc_del(k)) .eq. 999) then

              if (Pt_par%is_monitor .eq. 1) then
                call move_monitor_particles_file(Pt_par%idpar)
              end if
              Pt_par%idpar = -Pt_par%idpar

            else if (mark_block_id(idpar_proc_del(k)) .eq. 2) then

              procpar_mapping = 
     &        mapping(Pt_par%procpar)%nearest_proc(Pt_par%ipar,
     &                                             Pt_par%jpar)
              call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                  procpar_mapping,is_in_mapping,
     &                                  transpar)
              if (is_in_mapping) then
cty             Update transpar, ipar, jpar and kpar
                Pt_par%idpar = -Pt_par%idpar
                Pt_par%transpar = transpar
cty             Upadate Npar_out, pardata_send_real and pardata_send_integer
                Npar_out(procpar_mapping) = 
     &          Npar_out(procpar_mapping) + 1
                call copyin_real_data(Pt_par,
     &               pardata_send_real(:,Npar_out(procpar_mapping),
     &                                            procpar_mapping))
                call copyin_integer_data(Pt_par,
     &               pardata_send_integer(:,Npar_out(procpar_mapping),
     &                                               procpar_mapping))
              else

                do p = B_nprocs, nprocs-1

cty               The inner_surface is the blade surface
cty               or the inner_surface is the interpolation face
                  if ( (ind_j_proc(p) .eq. grids(which_block)%nproc_y) 
     &            .and. (p .ne. procpar_mapping)
     &               ) then

                    call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                        p,is_in,transpar) 

                    if (is_in) then
cty                   Update transpar, ipar, jpar and kpar
                      Pt_par%idpar = -Pt_par%idpar
                      Pt_par%transpar = transpar
cty                   Upadate Npar_out, pardata_send_real and pardata_send_integer
                      Npar_out(p) = Npar_out(p) + 1
                      call copyin_real_data(Pt_par,
     &                     pardata_send_real(:,Npar_out(p),p))
                      call copyin_integer_data(Pt_par,
     &                     pardata_send_integer(:,Npar_out(p),p))
                      exit
                    end if ! (is_in)
                  
                  end if ! (ind_j_proc(p) .eq. grids(which_block)%nproc_y)

                end do ! p

              end if ! (is_in_mapping)



            else if (mark_block_id(idpar_proc_del(k)) .eq. 1) then

              call obtain_proc_list(Pt_par%procpar,
     &                              Pt_par%ipar,Pt_par%jpar,
     &                              B_proc_list,B_proc_list_length)

              do i = 1, B_proc_list_length

                p = B_proc_list(i)
                call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                    p,is_in,transpar) 

                if (is_in) then
cty               Update transpar, ipar, jpar and kpar
                  Pt_par%idpar = -Pt_par%idpar
                  Pt_par%transpar = transpar
cty               Upadate Npar_out, pardata_send_real and pardata_send_integer
                  Npar_out(p) = Npar_out(p) + 1
                  call copyin_real_data(Pt_par,
     &                 pardata_send_real(:,Npar_out(p),p))
                  call copyin_integer_data(Pt_par,
     &                 pardata_send_integer(:,Npar_out(p),p))
                  exit
                end if ! (is_in)

              end do ! i

            end if ! (mark_block_id(idpar_proc_del(k)) .eq. 999) 

          end if


cty       2. If the particle originally locates in the O-type grid
          if (ov_ogrid_block(ibpar)) then

            if (mark_block_id(idpar_proc_del(k)) .eq. 999) then
              
              if (Pt_par%is_monitor .eq. 1) then
                call move_monitor_particles_file(Pt_par%idpar)
              end if
              Pt_par%idpar = -Pt_par%idpar

            else if (mark_block_id(idpar_proc_del(k)) .eq. 2) then

              call obtain_proc_list(Pt_par%procpar,
     &                              Pt_par%ipar,Pt_par%jpar,
     &                              O_proc_list,O_proc_list_length)

              do i = 1, O_proc_list_length

                p = O_proc_list(i)
                call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                    p,is_in,transpar) 

                if (is_in) then
cty               Update transpar, ipar, jpar and kpar
                  Pt_par%idpar = -Pt_par%idpar
                  Pt_par%transpar = transpar
cty               Upadate Npar_out, pardata_send_real and pardata_send_integer
                  Npar_out(p) = Npar_out(p) + 1
                  call copyin_real_data(Pt_par,
     &                 pardata_send_real(:,Npar_out(p),p))
                  call copyin_integer_data(Pt_par,
     &                 pardata_send_integer(:,Npar_out(p),p))
                  exit
                end if ! (is_in)

              end do ! i

            else if (mark_block_id(idpar_proc_del(k)) .eq. 1) then

              procpar_mapping = 
     &        mapping(Pt_par%procpar)%nearest_proc(Pt_par%ipar,
     &                                             Pt_par%jpar)
              call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                  procpar_mapping,is_in_mapping,
     &                                  transpar)
              if (is_in_mapping) then
cty             Update transpar, ipar, jpar and kpar
                Pt_par%idpar = -Pt_par%idpar
                Pt_par%transpar = transpar
cty             Upadate Npar_out, pardata_send_real and pardata_send_integer
                Npar_out(procpar_mapping) = 
     &          Npar_out(procpar_mapping) + 1
                call copyin_real_data(Pt_par,
     &               pardata_send_real(:,Npar_out(procpar_mapping),
     &                                            procpar_mapping))
                call copyin_integer_data(Pt_par,
     &               pardata_send_integer(:,Npar_out(procpar_mapping),
     &                                               procpar_mapping))
              else

                do p = 0, B_nprocs-1
                  if (p .ne. procpar_mapping) then

                    call particle_in_subblock(Pt_par%xpar,Pt_par%ypar,
     &                                        p,is_in,transpar) 

                    if (is_in) then
cty                   Update transpar, ipar, jpar and kpar
                      Pt_par%idpar = -Pt_par%idpar
                      Pt_par%transpar = transpar
cty                   Upadate Npar_out, pardata_send_real and pardata_send_integer
                      Npar_out(p) = Npar_out(p) + 1
                      call copyin_real_data(Pt_par,
     &                     pardata_send_real(:,Npar_out(p),p))
                      call copyin_integer_data(Pt_par,
     &                     pardata_send_integer(:,Npar_out(p),p))
                      exit
                    end if ! (is_in)

                  end if ! (p .ne. procpar_mapping)
                end do ! p

              end if ! (is_in_mapping)

            end if ! (mark_block_id(idpar_proc_del(k)) .eq. 999)

          end if ! (ov_bgrid_block(ibpar))

        end do ! k = 1, Npar_del

cty     Exchange the number of particles that flow into or out of the present processor 
        Npar_out_total = 0
        Npar_in_total = 0
        do p = 0, nprocs-1
          procid_send = p
          procid_recv = p
          call MPI_sendrecv(Npar_out(p),1,MPI_INTEGER,procid_send,tag,
     &                      Npar_in(p),1,MPI_INTEGER,procid_recv,tag,
     &                      MPI_comm_fluid,MPI_STATUS_IGNORE,ierr)
          Npar_out_total = Npar_out_total + Npar_out(p)
          Npar_in_total = Npar_in_total + Npar_in(p)
        end do

cty     Exchange the data of particles that flow into or out of the present processor
        do p = 0, nprocs-1

          if (Npar_out(p) .eq. 0) then
            procid_send = nullproc
          else
            procid_send = p
          end if
          if (Npar_in(p) .eq. 0) then
            procid_recv = nullproc
          else
            procid_recv = p
          end if

          call MPI_sendrecv(
     &         pardata_send_real(:,1:Npar_out(p),p),
     &         Nvar_par_real*Npar_out(p),real_mp_type,
     &         procid_send,tag,
     &         pardata_recv_real(:,1:Npar_in(p),p),
     &         Nvar_par_real*Npar_in(p),real_mp_type,
     &         procid_recv,tag,
     &         MPI_comm_fluid,MPI_STATUS_IGNORE,ierr)
          call MPI_sendrecv(
     &         pardata_send_integer(:,1:Npar_out(p),p),
     &         Nvar_par_integer*Npar_out(p),MPI_INTEGER,
     &         procid_send,tag,
     &         pardata_recv_integer(:,1:Npar_in(p),p),
     &         Nvar_par_integer*Npar_in(p),MPI_INTEGER,
     &         procid_recv,tag,
     &         MPI_comm_fluid,MPI_STATUS_IGNORE,ierr)

        end do

cty     Compress the main data of Pointer_par, delete the information of
cty     the particles that flow out of the present processor, and fill it
cty     with the data of the particle at the end of the line.
        Npar_local_settled = Npar_local
        Npar_del_settled = 0
        if (Npar_local_settled-Npar_del_settled .ge. 1) then

          do k = 1, Npar_del

            do while (Npar_local_settled .ge. 1)
              if (Pointer_par(Npar_local_settled)%idpar .gt. 0) then
                exit
              else
                Npar_local_settled = Npar_local_settled - 1
                Npar_del_settled = Npar_del_settled + 1
              end if ! (Pointer_par(Npar_local_settled)%idpar .gt. 0)
            end do ! (Npar_local_settled .ge. 1)

            if (Npar_del_settled .ge. Npar_del) exit
            Pointer_par(idpar_proc_del(k)) = 
     &      Pointer_par(Npar_local_settled)
            Npar_local_settled = Npar_local_settled - 1
            Npar_del_settled = Npar_del_settled + 1
            if (Npar_del_settled .ge. Npar_local) exit

          end do ! k

        else
          Npar_local_settled = 0
        end if ! (Npar_local_settled-Npar_del_settled .ge. 1)

cty     Loads the received particle data into the master data Pointer_par 
        do p = 0, nprocs-1
          do k = 1, Npar_in(p)
            Npar_local_settled = Npar_local_settled + 1
            call copyout_real_data(
     &           pardata_recv_real(:,k,p),
     &           Pointer_par(Npar_local_settled))
            call copyout_integer_data(
     &           pardata_recv_integer(:,k,p),
     &           Pointer_par(Npar_local_settled))
          enddo 
        enddo 

cty     Update the number of particles in every processor at each time step
        Npar_local = Npar_local - Npar_del + Npar_in_total

        do k = 1, Npar_local
          Pt_par => Pointer_par(k)

cty       Update transpar, ipar, jpar and kpar
          if (Pt_par%idpar .lt. 0) then
            Pt_par%idpar = -Pt_par%idpar
            ipar_ori = Pt_par%ipar
            jpar_ori = Pt_par%jpar

            if ( (block_proc(Pt_par%procpar) .eq. 1) .and.
     &           (block_proc(procid_fluid) .eq. 1) ) then

              call update_B_grid_ij_particle
     &            (Pt_par%transpar,Pt_par%xpar,Pt_par%ypar,
     &             Pt_par%upar,Pt_par%vpar,
     &             Pt_par%procpar,procid_fluid,
     &             ipar_ori,jpar_ori,
     &             Pt_par%ipar,Pt_par%jpar)

            else if ( (block_proc(Pt_par%procpar) .eq. 2) .and.
     &              (block_proc(procid_fluid) .eq. 2) ) then

              ipar_predict = ibeg_proc(Pt_par%procpar)+ipar_ori-1
              vel_tangen = Pt_par%upar*tangen_x_global(ipar_predict)
     &                   + Pt_par%vpar*tangen_y_global(ipar_predict)
              vel_normal = Pt_par%upar*tangen_x_global(ipar_predict)
     &                   + Pt_par%vpar*tangen_y_global(ipar_predict)
              call update_O_grid_ij_particle_prediction
     &            (Pt_par%transpar,
     &             Pt_par%xpar,Pt_par%ypar,
     &             vel_tangen,vel_normal,
     &             Pt_par%procpar,procid_fluid,
     &             ipar_ori,jpar_ori,
     &             Pt_par%ipar,Pt_par%jpar)

            else if ( (block_proc(Pt_par%procpar) .eq. 1) .and.
     &                (block_proc(procid_fluid) .eq. 2) ) then

              ipar_mapping =
     &        mapping(Pt_par%procpar)%nearest_i(Pt_par%ipar,
     &                                          Pt_par%jpar)
              jpar_mapping =
     &        mapping(Pt_par%procpar)%nearest_j(Pt_par%ipar,
     &                                          Pt_par%jpar)
              nx = subblock_proc(procid_fluid)%xsize
cty           It should be notice that Pt_par%transpar has been updated before the
cty           particle is distributed to the new processor
              y = transformed_y(Pt_par%transpar,Pt_par%ypar)

              start_i = max(ipar_mapping-2,1)
              end_i = min(ipar_mapping+2,nx-1)
              step_i = 1
              start_j = jpar_mapping
              end_j = jpar_mapping-2
              step_j = -1

              outer_HO: do j = start_j, end_j, step_j
                inner_HO: do i = start_i, end_i, step_i
                  
                  call identify_O_grid_ij_particle
     &                (Pt_par%xpar,y,i,j,is_in,Pt_par%ipar,Pt_par%jpar)
                  if (is_in) then
                    Pt_par%procpar = procid_fluid 
                    exit outer_HO
                  end if

                end do inner_HO ! i
              end do outer_HO ! j 

              if (.not. is_in) then
                call update_B_O_grid_ij_particle
     &              (Pt_par%transpar,Pt_par%xpar,Pt_par%ypar,
     &               Pt_par%procpar,procid_fluid,
     &               ipar_ori,jpar_ori,
     &               Pt_par%ipar,Pt_par%jpar)
              end if ! (.not. is_in)

            else if ( (block_proc(Pt_par%procpar) .eq. 2) .and.
     &                (block_proc(procid_fluid) .eq. 1) ) then

              ipar_mapping =
     &        mapping(Pt_par%procpar)%nearest_i(Pt_par%ipar,
     &                                          Pt_par%jpar)
              jpar_mapping =
     &        mapping(Pt_par%procpar)%nearest_j(Pt_par%ipar,
     &                                          Pt_par%jpar)
              nx = subblock_proc(procid_fluid)%xsize
              ny = subblock_proc(procid_fluid)%ysize
cty           It should be notice that Pt_par%transpar has been updated before the
cty           particle is distributed to the new processor
              y = transformed_y(Pt_par%transpar,Pt_par%ypar)

              start_i = max(ipar_mapping-2,1)
              end_i = min(ipar_mapping+2,nx-1)
              step_i = 1
              start_j = max(jpar_mapping-2,1)
              end_j = min(jpar_mapping+2,ny-1)
              step_j = 1

              outer_OH: do j = start_j, end_j, step_j 
                inner_OH: do i = start_i, end_i, step_i
                  
                  call identify_B_grid_ij_particle
     &                (Pt_par%xpar,y,i,j,is_in,Pt_par%ipar,Pt_par%jpar)
                  if (is_in) then
                    Pt_par%procpar = procid_fluid 
                    exit outer_OH
                  end if

                end do inner_OH ! i
              end do outer_OH ! j 

              if (.not. is_in) then
                Pt_par%procpar = procid_fluid 
                call compute_ij_particle
     &              (Pt_par%transpar,Pt_par%xpar,Pt_par%ypar,
     &               Pt_par%ipar,Pt_par%jpar)
              end if ! (.not. is_in)

            end if

          end if ! (Pt_par%idpar .lt. 0)

cty       For pathlines mode, interpolate the updated velocities

          call compute_particles_force(Pointer_par(k)
     &                                ,c_ijk
     &                                ,Fpar
     &                                 )


          Pt_par%axpar = 0.5_wp * Pt_par%axpar + 
     &                   0.5_wp * Fpar(1)
          Pt_par%aypar = 0.5_wp * Pt_par%aypar + 
     &                   0.5_wp * Fpar(2)
          Pt_par%azpar = 0.5_wp * Pt_par%azpar + 
     &                   0.5_wp * Fpar(3)

          if (Pt_par%dpar .lt. 0) then
            Pt_par%dpar = -Pt_par%dpar
          else
            Pt_par%upar = Pt_par%upar + Pt_par%axpar*dt
            Pt_par%vpar = Pt_par%vpar + Pt_par%aypar*dt 
            Pt_par%wpar = Pt_par%wpar + Pt_par%azpar*dt
          end if

        end do ! k

        deallocate(Npar_out)
        deallocate(Npar_in)

        deallocate(pardata_send_real)
        deallocate(pardata_recv_real)
        deallocate(pardata_send_integer)
        deallocate(pardata_recv_integer)

cty     Sum Npar_local in ioid, which results in Npar_global
        call MPI_reduce(Npar_local,Npar_global,1,MPI_INTEGER,
     &                  MPI_SUM,ioid,MPI_comm_fluid,ierr)
        call capture_particles_number(time,Npar_global)

      end subroutine ! timeadvance_redistribute_particles



cty========================================================================
cty     subroutine compute_particles_force
cty========================================================================
cty
cty     Description:
cty     Calculate the force exerted on the particle by the fluid.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       Pt_par: Pointer to the particle
cty
cty     Output arguments:
cty       c_ijk: The weight cofficients 
cty       Fpar: The force on a particle
cty       Torpar: The torque on a particle
cty       Qpar: The heat on a particle
cty
cty     Notice: 
cty
cty========================================================================

      subroutine compute_particles_force(Pt_par
     &                                  ,c_ijk
     &                                  ,Fpar
     &                                  )

        use corevars, only : nstart, nit, p_deriv
        use flow_properties, only : rey, rrey
        use, intrinsic :: ieee_arithmetic

        implicit none

        type(particle), intent(in) :: Pt_par         
        real(kind=wp), intent(out) :: c_ijk(-2:2,-2:2,-2:2)
        real(kind=wp), intent(out) :: Fpar(3)

        real(kind=wp) :: Fdrag(3)

        real(kind=wp) :: uflow, vflow, wflow
        real(kind=wp) :: rhoflow, muflow

        real(kind=wp) :: relative_vel
        real(kind=wp) :: relative_velx, relative_vely, relative_velz
        real(kind=wp) :: Repar
        real(kind=wp) :: kappar 
        real(kind=wp) :: Cd 

        call interpolate_flow_information(
     &       Pt_par%transpar,
     &       Pt_par%xpar,Pt_par%ypar,
     &       Pt_par%zpar,
     &       Pt_par%procpar,Pt_par%ipar,Pt_par%jpar,Pt_par%kpar,
     &       c_ijk,
     &       uflow,vflow,wflow,
     &       rhoflow,muflow
     &       )


cty     Calculate the relative velocity
        relative_velx = uflow - Pt_par%upar
        relative_vely = vflow - Pt_par%vpar
        relative_velz = wflow - Pt_par%wpar
        relative_vel = dsqrt( relative_velx**2._wp +
     &                        relative_vely**2._wp +
     &                        relative_velz**2._wp )
cty     In case that relative_vel equals to zero
        relative_vel = relative_vel + 1.0e-16

cty     Calculate the particle Reynolds number and shear Reynolds number
        Repar = rey * Para_non_dim%diameter_par
     *              * relative_vel
     &              * rhoflow/muflow 


cty     Calculate the drag force
        kappar = 3._wp*PI*muflow
        if (Repar .le. 1000._wp) then
ct        Cd = 24._wp/Repar*(1._wp+0.15_wp*Repar**0.687_wp)
          Cd = 1._wp + 0.15_wp*Repar**0.687_wp
        else
cty       Cd = 0.44_wp
          Cd = 0.44_wp/24._wp*Repar
        end if 
        Cd = Cd / Para_non_dim%diameter_par
     &          / (rey*Para_non_dim%diameter_par) 
     &          / (1._wp/6._wp*PI*Para_non_dim%density_par)
        Fdrag(1) = kappar * Cd * relative_velx 
        Fdrag(2) = kappar * Cd * relative_vely 
        Fdrag(3) = kappar * Cd * relative_velz 

cty     Sum all the force
        Fpar(1) = Fdrag(1) 

        Fpar(2) = Fdrag(2) 

        Fpar(3) = Fdrag(3) 

      end subroutine ! compute_particles_force



cty========================================================================
cty     subroutine bounce_particles
cty========================================================================
cty
cty     Description:
cty     Update the particle after the bounce with the blade surface.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       transpar: trans_mode if (xpar, ypar) lies inside the boundaries of the y 
cty                 transformed subblock of O-type grid.
cty       ipar: Local indices of the cell where the particle locates in
cty
cty     Output arguments:
cty       xpar, ypar: Mirror coordinates of the particle
cty       upar, vpar: Mirror velocities of the particle
cty
cty     Notice: 
cty
cty========================================================================

      subroutine bounce_particles(transpar,ipar
     &                           ,xpar,ypar
     &                           ,upar,vpar
     &                            )
 
        use corevars, only : z_x, r_y
        use, intrinsic :: ieee_arithmetic

        implicit none

        real(kind=wp), intent(in) :: transpar
        integer, intent(in) :: ipar
        real(kind=wp), intent(inout) :: xpar, ypar
        real(kind=wp), intent(inout) :: upar, vpar


        call compute_mirror_particle(xpar
     &                              ,transformed_y(transpar,ypar)
     &                              ,upar,vpar
     &                              ,z_x(ipar,1),r_y(ipar,1)
     &                              ,z_x(ipar+1,1),r_y(ipar+1,1)
     &                              ,xpar,ypar
     &                              ,upar,vpar
     &                               ) 
        call transform_particles_xy_coordinates(xpar,ypar)  

      end subroutine ! bounce_particles



cty========================================================================
cty     subroutine compute_mirror_particle
cty========================================================================
cty
cty     Description:
cty     Calculate the mirror location and velocities of the particle after
cty     the particle bounces with the blade surface.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       xpar, ypar: Coordinates of the particle
cty       upar, vpar: Velocities of the particle
cty       x1, y1: First extremity of the segment
cty       x2, y2: Second extremity of the segment
cty
cty     Output arguments:
cty       x, y: Mirror coordinates of the particle
cty       u, v: Mirror velocities of the particle
cty
cty     Notice: 
cty
cty========================================================================

      subroutine compute_mirror_particle(xpar,ypar
     &                                  ,upar,vpar
     &                                  ,x1,y1,x2,y2
     &                                  ,x,y,u,v
     &                                   )

        use corevars, only : dt

        real(kind=wp), intent(in) :: xpar, ypar
        real(kind=wp), intent(in) :: upar, vpar
        real(kind=wp), intent(in) :: x1, y1, x2, y2
        real(kind=wp), intent(out) :: x, y
        real(kind=wp), intent(out) :: u, v

        real(kind=wp) :: dpx, dpy
        real(kind=wp) :: k, rec_k

        dpx = x2-x1
        dpy = y2-y1

        if (abs(dpx) .gt. abs(dpy)) then
          k = dpy/dpx
cty       Firstly, we calculate the cross point
          y = (k*(k*ypar+xpar) + (y1-k*x1)) / (1._wp+k**2)
          x = -k*(y-ypar)+xpar
cty       Subsequently, we calculate the mirror point
          y = 2._wp*y - ypar
          x = 2._wp*x - xpar
          v = (2._wp*k*upar + (-1._wp+k**2)*vpar) / (1._wp+k**2)
          u = -k*(v-vpar)+upar
        else
cty       Firstly, we calculate the cross point
          rec_k = dpx/dpy
          y = (-rec_k*(-rec_k*y1+x1) + (ypar+rec_k*xpar)) / 
     &        (1._wp+rec_k**2)
          x = rec_k*(y-y1)+x1
cty       Subsequently, we calculate the mirror point
          y = 2._wp*y - ypar
          x = 2._wp*x - xpar
          u = (2._wp*rec_k*vpar + (-1._wp+rec_k**2)*upar) / 
     &        (1._wp+rec_k**2)
          v = -rec_k*(u-upar)+vpar
        end if


      end subroutine ! compute_mirror_particle



cty========================================================================
cty     subroutine transform_particles_xyz_coordinates
cty========================================================================
cty
cty     Description:
cty     Transform the xyz coordinates of the particle based on the periodic
cty     boundary condition in the pitchwise and spanwise direction.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       xpar, ypar, zpar: Coordinates of the particle
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine transform_particles_xyz_coordinates(xpar,ypar,zpar)

        use overset, only : ov_pitch
        use corevars, only : thl

        implicit none

        real(kind=wp) :: xpar, ypar, zpar


        if (ypar .le. B_ymin) then
          ypar = ypar + ov_pitch
        else if (ypar .gt. B_ymax) then
          ypar = ypar - ov_pitch
        end if

        if (zpar .le. B_zmin) then
          zpar = zpar + thl
        else if (zpar .gt. B_zmax) then
          zpar = zpar - thl
        end if

      end subroutine ! transform_particles_xyz_coordinates



cty========================================================================
cty     subroutine transform_particles_xy_coordinates
cty========================================================================
cty
cty     Description:
cty     Transform the x and y coordinates of the particle based on the periodic
cty     boundary condition in the streamwise and pitchwise direction.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       xpar, ypar: Coordinates of the particle
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine transform_particles_xy_coordinates(xpar,ypar)

        use overset, only : ov_pitch

        implicit none

        real(kind=wp) :: xpar, ypar


        if (ypar .le. B_ymin) then
          ypar = ypar + ov_pitch
        else if (ypar .gt. B_ymax) then
          ypar = ypar - ov_pitch
        end if

      end subroutine ! transform_particle_xy_coordinates



cty========================================================================
cty     subroutine reduce_proc_list
cty========================================================================
cty
cty     Description:
cty     Delete the negative number in the array proc_list, then the same number,
cty     return the length of the updated array proc_list. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty       proc_list: Length of the updated array proc_list
cty
cty     Notice: 
cty
cty========================================================================

      subroutine reduce_proc_list(proc_list,proc_list_length)

        use mppvars, only : nullproc, nprocs

        implicit none

        integer, intent(inout) :: proc_list(3)
        integer, intent(out) :: proc_list_length

        integer :: proc_list_temp(3)     

        integer :: i, j

        proc_list_length = 0
        do i = 1, 3
          if (proc_list(i) .ne. nullproc) then
            proc_list_length = proc_list_length + 1
            proc_list_temp(proc_list_length) = proc_list(i)
          end if ! (proc_list(i) .gt. 0)
        end do ! i
        proc_list = proc_list_temp

        j = 0
        do i = 1, proc_list_length
          if (.not. any(proc_list(1:j).eq.proc_list_temp(i))) then
            j = j + 1
            proc_list(j) = proc_list_temp(i)
          end if ! (.not. any(proc_list(1:j) .eq. proc_list_temp(i)))
        end do ! i
        proc_list_length = j


      end subroutine ! reduce_proc_list



cty========================================================================
cty     subroutine obtain_proc_list
cty========================================================================
cty
cty     Description:
cty     Obtain the array proc_list and the length. 
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice: 
cty
cty========================================================================

      subroutine obtain_proc_list(procpar,ipar,jpar,
     &                            proc_list,proc_list_length)
 
        use corevars, only : nxp, nyp

        implicit none

        integer, intent(in) :: procpar
        integer, intent(in) :: ipar, jpar
        integer, intent(out) :: proc_list(3)
        integer, intent(out) :: proc_list_length

        if (ipar .le. 0.5_wp*nxp) then
          if (jpar .le. 0.5_wp*nyp) then
            proc_list = xmym_proc_list(procpar)%id
            proc_list_length = xmym_proc_list(procpar)%length
          else
            proc_list = xmyp_proc_list(procpar)%id
            proc_list_length = xmyp_proc_list(procpar)%length
          end if 
        else
          if (jpar .le. 0.5_wp*nyp) then
            proc_list = xpym_proc_list(procpar)%id
            proc_list_length = xpym_proc_list(procpar)%length
          else
            proc_list = xpyp_proc_list(procpar)%id
            proc_list_length = xpyp_proc_list(procpar)%length
          end if 
        end if 

      end subroutine ! obtain_proc_list



cty========================================================================
cty     subroutine copyin_real_data
cty========================================================================
cty
cty     Description:
cty     Copy the real(kind=wp) data of the particles into the arrays pardata.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       Pt_par: Pointer to the particle
cty
cty     Output arguments:
cty       pardata_real: Arrays contains the real(kind=wp) data of the particles
cty
cty     Notice: 
cty
cty========================================================================

      pure subroutine copyin_real_data(Pt_par,pardata_real)
       
        implicit none 

        type(particle), intent(in) :: Pt_par
        real(kind=wp), intent(out) :: pardata_real(Nvar_par_real)

        pardata_real(1) = Pt_par%idpar 
        pardata_real(2) = Pt_par%dpar
        pardata_real(3) = Pt_par%rhopar 
        pardata_real(4) = Pt_par%transpar
        pardata_real(5) = Pt_par%xpar 
        pardata_real(6) = Pt_par%ypar
        pardata_real(7) = Pt_par%zpar
        pardata_real(8) = Pt_par%upar  
        pardata_real(9) = Pt_par%vpar
        pardata_real(10) = Pt_par%wpar
        pardata_real(11) = Pt_par%axpar 
        pardata_real(12) = Pt_par%aypar
        pardata_real(13) = Pt_par%azpar

      end subroutine ! copyin_real_data



cty========================================================================
cty     subroutine copyout_real_data
cty========================================================================
cty
cty     Description:
cty     Copy the real(kind=wp) data of the arrays pardata into the particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       pardata_real: Arrays contains the real(kind=wp) data of the particles
cty
cty     Output arguments:
cty       Pt_par: Pointer to the particle
cty
cty     Notice: 
cty
cty========================================================================

      pure subroutine copyout_real_data(pardata_real,Pt_par)
       
        implicit none 

        real(kind=wp), intent(in) :: pardata_real(Nvar_par_real)
        type(particle), intent(out) :: Pt_par

        Pt_par%idpar = pardata_real(1)  
        Pt_par%dpar = pardata_real(2) 
        Pt_par%rhopar = pardata_real(3) 
        Pt_par%transpar = pardata_real(4)
        Pt_par%xpar = pardata_real(5) 
        Pt_par%ypar = pardata_real(6) 
        Pt_par%zpar = pardata_real(7) 
        Pt_par%upar = pardata_real(8) 
        Pt_par%vpar = pardata_real(9) 
        Pt_par%wpar = pardata_real(10) 
        Pt_par%axpar = pardata_real(11) 
        Pt_par%aypar = pardata_real(12) 
        Pt_par%azpar = pardata_real(13) 

      end subroutine ! copyout_real_data



cty========================================================================
cty     subroutine copyin_integer_data
cty========================================================================
cty
cty     Description:
cty     Copy the integer data of the particles into the arrays pardata.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       Pt_par: Pointer to the particle
cty
cty     Output arguments:
cty       pardata_integer: Arrays contains the integer data of the particles
cty
cty     Notice: 
cty
cty========================================================================

      pure subroutine copyin_integer_data(Pt_par,pardata_integer)
       
        implicit none 

        type(particle), intent(in) :: Pt_par
        integer, intent(out) :: pardata_integer(Nvar_par_integer)

        pardata_integer(1) = Pt_par%procpar
        pardata_integer(2) = Pt_par%ipar
        pardata_integer(3) = Pt_par%jpar
        pardata_integer(4) = Pt_par%kpar
        pardata_integer(5) = Pt_par%is_monitor

      end subroutine ! copyin_integer_data



cty========================================================================
cty     subroutine copyout_integer_data
cty========================================================================
cty
cty     Description:
cty     Copy the integer data of the arrays pardata into the particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty       pardata_integer: Arrays contains the integer data of the particles
cty
cty     Output arguments:
cty       Pt_par: Pointer to the particle
cty
cty     Notice: 
cty
cty========================================================================

      pure subroutine copyout_integer_data(pardata_integer,Pt_par)
       
        implicit none 

        integer, intent(in) :: pardata_integer(Nvar_par_integer)
        type(particle), intent(out) :: Pt_par

        Pt_par%procpar = pardata_integer(1) 
        Pt_par%ipar = pardata_integer(2) 
        Pt_par%jpar = pardata_integer(3) 
        Pt_par%kpar = pardata_integer(4) 
        Pt_par%is_monitor = pardata_integer(5) 

      end subroutine ! copyout_integer_data



cty========================================================================
cty     subroutine move_monitor_particles_file
cty========================================================================
cty
cty     Description:
cty     This routine move the file that contains particle data in the folder 
cty     Particles_continue to the new folder Particles_end.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine move_monitor_particles_file(idpar)

        implicit none

        real(kind=wp), intent(in) :: idpar

        integer :: i
        character(21) :: idpar_name
        character(100) :: particle_file, particle_file_out
        logical :: lex
        character(300) :: exec_arg

        if (Para_par%Monitor_particle_steps .gt. 0) then
        
cty       Set the file name of the particle_file and particle_file_out
          write(idpar_name,'(F9.0)') idpar
cty       Here we use 8 places to mark the total global id of the 
cty       particles, and fill the blanks with zeros.
          do i = 1, len_trim(idpar_name)
            if (idpar_name(i:i) .eq. '') then
              idpar_name(i:i) = '0' 
            else 
              exit
            end if
          end do 
          particle_file = 'Monitor_particles/Monitor_particle_'//
     &                     trim(idpar_name)//'bin'
          particle_file = trim(adjustl(particle_file))
          particle_file_out = 'Monitor_particles/out/'//
     &                        'Monitor_particle_'//
     &                         trim(idpar_name)//'bin'
          particle_file_out = trim(adjustl(particle_file_out))

          write(exec_arg,'(A,A,A,A)') 'mv -f ', particle_file, 
     &                                ' ', particle_file_out

          inquire(file=particle_file,exist=lex)
cty       Move the file to the new folder
          if (lex) then
            call system(exec_arg)
          else
            call stopnow('File Monitor_particle_'//trim(idpar_name)//
     &                   'bin does not exist!')
          end if ! (lex)

        end if ! (Para_par%Monitor_particle_steps .gt. 0)

      end subroutine ! move_monitor_particles_file



cty========================================================================
cty     subroutine add_particles
cty========================================================================
cty
cty     Description:
cty     This routine continuously add particles into the flow field.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine add_particles

        use corevars, only : nstart, nsteps, nit, time
        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_SUM, MPI_MAX,
     &                      MPI_comm_fluid, ioid, procid_fluid, 
     &                      nprocs, ioproc
        use ompvars, only : start_master, end_master

        implicit none

        integer :: k
        type(particle), pointer :: Pt_par
        real(kind=wp) :: idpar_Max_local, idpar_Max_global
        integer :: ierr

        real(kind=wp), allocatable, dimension(:,:) :: pardata
        real(kind=wp) :: transpar, xpar, ypar, zpar
        logical :: is_in

        if ( (mod(nit-nstart,Para_par%Add_particle_steps).eq.0)
     &       .and. (nit.gt.nstart) ) then

          if (start_master()) then

            idpar_Max_local = 0._wp
            do k = 1, Npar_local
              Pt_par => Pointer_par(k)
              idpar_Max_local = max(Pt_par%idpar,idpar_Max_local)
            end do ! k
cty         Get the max idpar_Max_local in ioid, which results in idpar_Max_global
            call MPI_reduce(idpar_Max_local,idpar_Max_global,
     &                      1,real_mp_type,MPI_MAX,ioid,
     &                      MPI_comm_fluid,ierr)
            call MPI_bcast(idpar_Max_global,1,real_mp_type,ioid,
     &                     MPI_comm_fluid,ierr) 

cty         Add particles from the inlet
            do k = 1, Npar_added_local 

              Npar_local = Npar_local + 1
              Pointer_par(Npar_local)%idpar = 
     &                                added_pardata_real(1,k) +
     &                                idpar_Max_global
              Pointer_par(Npar_local)%dpar = 
     &                                added_pardata_real(2,k)
              Pointer_par(Npar_local)%rhopar = 
     &                                added_pardata_real(3,k)
              Pointer_par(Npar_local)%transpar = 
     &                                added_pardata_real(4,k)
              Pointer_par(Npar_local)%xpar = 
     &                                added_pardata_real(5,k)
              Pointer_par(Npar_local)%ypar = 
     &                                added_pardata_real(6,k)
              Pointer_par(Npar_local)%zpar = 
     &                                added_pardata_real(7,k)

              Pointer_par(Npar_local)%axpar = 0._wp 
              Pointer_par(Npar_local)%aypar = 0._wp 
              Pointer_par(Npar_local)%azpar = 0._wp 


              call interpolate_flow_part_information(
     &             added_pardata_real(4,k),
     &             added_pardata_real(5,k),
     &             added_pardata_real(6,k),
     &             added_pardata_real(7,k),
     &             added_pardata_integer(1,k),
     &             added_pardata_integer(2,k),
     &             added_pardata_integer(3,k),
     &             added_pardata_integer(4,k),
     &             Pointer_par(Npar_local)%upar,
     &             Pointer_par(Npar_local)%vpar,
     &             Pointer_par(Npar_local)%wpar
     &             )


              Pointer_par(Npar_local)%procpar = 
     &                                added_pardata_integer(1,k) 
              Pointer_par(Npar_local)%ipar = 
     &                                added_pardata_integer(2,k) 
              Pointer_par(Npar_local)%jpar = 
     &                                added_pardata_integer(3,k) 
              Pointer_par(Npar_local)%kpar = 
     &                                added_pardata_integer(4,k)
              Pointer_par(Npar_local)%is_monitor = 
     &                                added_pardata_integer(5,k) 

            end do ! k


cty         Sum Npar_local in ioid, which results in Npar_global
            call MPI_reduce(Npar_local,Npar_global,1,MPI_INTEGER,
     &                      MPI_SUM,ioid,MPI_comm_fluid,ierr)

          end if ! start_master
          call end_master()

        end if ! (mod(nit-nstart,Para_par%Add_particle_steps).eq.0) 

      end subroutine ! add_particles



cty========================================================================
cty     subroutine read_added_particles_init
cty========================================================================
cty
cty     Description:
cty     This routine initializes added particles that enter into the flow 
cty     field based on the input file.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine read_added_particles_init

        use mppvars, only : MPI_INTEGER, real_mp_type, MPI_SUM, MPI_MAX,
     &                      MPI_comm_fluid, ioid, procid_fluid, 
     &                      nprocs, ioproc
        use mbvars, only : num_blocks
     &                   , ov_ogrid_block, ov_bgrid_block

        implicit none

        integer :: k, j
        integer :: lex, ierr
       
        integer :: Npar_added_init
        integer, allocatable, dimension(:) :: Npar_added_local_proc

        real(kind=wp), allocatable, dimension(:,:) :: pardata
        real(kind=wp) :: transpar, xpar, ypar, zpar
        integer :: block_id
     &           , which_block
        logical :: is_in


        inquire(file='ADDED_PARTICLES_INPUT.dat',exist=lex)
        if (lex) then

          if (ioproc) then
            open(unit=1000,file='ADDED_PARTICLES_INPUT.dat',
     &                     status='old',action='read')

cty         Notice that Npar_added_init only be read one time
            read(1000,*) Npar_added_init
            write(*,*) 'Begin to read ADDED_PARTICLES_INPUT.dat'
            write(*,'(A,I11)') 'Initial added particles number = ', 
     &                          Npar_added_init
          end if ! (ioproc) 
          call MPI_bcast(Npar_added_init,1,MPI_INTEGER,
     &                   ioid,MPI_comm_fluid,ierr)

cty       Allocate pardata arrays used for reading data and assembling particles
          allocate(added_pardata_real(7,Npar_added_init))
          allocate(added_pardata_integer(Nvar_par_integer,
     &                                   Npar_added_init))
          allocate(pardata(4,Npar_added_init))

cty       Initialize Npar_added_local for each processor
          Npar_added_local = 0

cty       Read added particles
cty       Notice the data is only be read in io processor, 
cty       then be bcasted to other processors.
          if (ioproc) then
            do k = 1, Npar_added_init

              read(1000,*) (pardata(j,k), j = 1, 4)  

            end do ! k
            close(1000) 
          end if ! (ioproc)

          call MPI_bcast(pardata,4*Npar_added_init,
     &                   real_mp_type,ioid,MPI_comm_fluid,ierr) 

          do k = 1, Npar_added_init
            xpar = pardata(2,k)
            ypar = pardata(3,k)
            zpar = pardata(4,k)

cty         Delete input particles outside the background grid previously
            if (    
     &              (xpar.lt.B_xmin .or. xpar.ge.B_xmax) 
     &         .or. (ypar.le.B_ymin .or. ypar.gt.B_ymax) 
     &         .or. (zpar.le.B_zmin .or. zpar.gt.B_zmax)
     &          ) then

              Npar_added_init = Npar_added_init - 1

              cycle

            else  

              block_id = block_proc(procid_fluid)

              if (.not. particle_in_inner_surface(xpar,ypar)) then

                do j = 1, num_blocks
                  if (ov_bgrid_block(j)) then
                    which_block = j
cty                 Here we can coarsely identify does the particle only locate in backgroun grid
                    if ( (xpar .lt. O_outer_xmin) .or.
     &                   (xpar .gt. O_outer_xmax) ) then
                      exit
                    end if
                  else if (ov_ogrid_block(j) .and. 
     &                     particle_in_outer_surface(xpar,ypar))
     &              then
                    which_block = j
                  end if
                end do

cty             Subsequently, save particles in the corresponding processor
                if (which_block .eq. block_id) then
                  call particle_in_subblock(xpar,ypar,
     &                                      procid_fluid,is_in,transpar)

                  if (is_in) then

                    Npar_added_local = Npar_added_local + 1
                    added_pardata_real(1,Npar_added_local) = 
     &                                   pardata(1,k) 
                    added_pardata_real(2,Npar_added_local) = 
     &                                   Para_non_dim%diameter_par
                    added_pardata_real(3,Npar_added_local) = 
     &                                   Para_non_dim%density_par
cty                 transpar is used to simplify the process of compute_ij_particle
                    added_pardata_real(4,Npar_added_local) = transpar 
                    added_pardata_real(5,Npar_added_local) = xpar 
                    added_pardata_real(6,Npar_added_local) = ypar
                    added_pardata_real(7,Npar_added_local) = zpar

                    added_pardata_integer(1,Npar_added_local) = 
     &                                      procid_fluid 
                    added_pardata_integer(5,Npar_added_local) = 999

                    call compute_ij_particle(
     &                   transpar,xpar,ypar, 
     &                   added_pardata_integer(2,Npar_added_local), 
     &                   added_pardata_integer(3,Npar_added_local))
                    call compute_k_particle(zpar,
     &                   added_pardata_integer(4,Npar_added_local))


                  end if ! (is_in)
                end if ! (which_block .eq. block_id)     

              else

              end if ! (.not. particle_in_inner_surface(xpar,ypar))

            end if ! outside the background grid

          end do ! k

          deallocate(pardata)

        else
          call stopnow('file ADDED_PARTICLES_INPUT.dat '//
     &                 'does not exist')
        end if   

        allocate(Npar_added_local_proc(0:nprocs-1))
cty     Collect Npar_added_local in ioid, which is save in Npar_added_local_proc(procid_fluid)
        call MPI_gather(Npar_added_local,1,MPI_INTEGER,
     &                  Npar_added_local_proc(procid_fluid),1,
     &                  MPI_INTEGER,ioid,MPI_comm_fluid,ierr)
cty     Check the particle number Npar_added_local in each processor 
        if (ioproc) then
          do k = 0, nprocs-1
            write(*,'(A,I4,A9,I4,A4,I11)') 
     &              'The added particle number in processor', k, 
     &              'of block', block_proc(k), 
     &              'is', Npar_added_local_proc(k)
          end do
        end if
        deallocate(Npar_added_local_proc)

cty     Sum Npar_added_local in ioid, which results in Npar_added_global
        call MPI_reduce(Npar_added_local,Npar_added_global,
     &                  1,MPI_INTEGER,MPI_SUM,ioid,
     &                  MPI_comm_fluid,ierr)
        call MPI_bcast(Npar_added_global,1,MPI_INTEGER,
     &                 ioid,MPI_comm_fluid,ierr) 

        Para_par%Add_particle_number = Npar_added_global

        if (ioproc) then
          if (Para_par%Add_particle_number .le. 0) then
            call stopnow
     &      ('The Add_particle_number must be positive')
          else 
            write(*,'(A,I7,A20)')
     &      ' Each time,', Para_par%Add_particle_number,
     &      'particles are added'
          end if ! (Para_par%Add_particle_number .le. 0)
          write(*,*) 'ADDED_PARTICLES_INPUT.dat has been read'
        end if

      end subroutine ! read_added_particles_init



cty========================================================================
cty     subroutine init_particles_statistics
cty========================================================================
cty
cty     Description:
cty     This routine initializes the new statistic routine for particles. 
cty     It defines the indices as well and allocates the arrays.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine init_particles_statistics

        use corevars, only : nxp, nyp, nzp
        use mbvars, only : num_blocks
     &                   , ov_bgrid_block, ov_ogrid_block
        use mppvars, only : ioproc
        use global_grid, only : grids

        implicit none

        integer :: block_id

        n_qstat_par = n_qstat_par_vel
        

        if (Para_par%Average_particle_steps .gt. 0) then

          if (ioproc) then
            write(*,*) 'For particles, '//
     &                 'the following statistcs are collected:'
            write(*,*) ' n: the number density of particles'
            write(*,*) ' u: x-velocity of particles'
            write(*,*) ' v: y-velocity of particles'
            write(*,*) ' w: z-velocity of particles'
            write(*,*) " u'u': x-Reynolds normal stress of particles" 
            write(*,*) " v'v': y-Reynolds normal stress of particles"
            write(*,*) " w'w': z-Reynolds normal stress of particles"
            write(*,*) " u'v': xy-Reynolds shear stress of particles"
            write(*,*) " u'w': xz-Reynolds shear stress of particles"
            write(*,*) " v'w': yz-Reynolds shear stress of particles"
            write(*,*) 'n_qstat_par =', n_qstat_par
          end if

cty       Allocate the array qstat_par for the calculation of particles statistics 
          allocate(qstat_par(nxp,nyp,1,n_qstat_par))
          qstat_par(:,:,:,:) = 0._wp

          if (ioproc) then
cty         Allocate the array B_qstat_par_gather and O_qstat_par_gather to
cty         gather particles statistics 
            do block_id = 1, num_blocks
              if (ov_bgrid_block(block_id)) then 
                allocate(B_qstat_par_gather(grids(block_id)%nxp_block,
     &                                      grids(block_id)%nyp_block,
     &                                      1,
     &                                      n_qstat_par))
                B_qstat_par_gather(:,:,:,:) = 0._wp
              end if

              if (ov_ogrid_block(block_id)) then
                allocate(O_qstat_par_gather(grids(block_id)%nxp_block,
     &                                      grids(block_id)%nyp_block,
     &                                      1,
     &                                      n_qstat_par))
                O_qstat_par_gather(:,:,:,:) = 0._wp
              end if
            end do ! block_id
          end if

        end if ! (Para_par%Average_particle_steps .gt. 0)

      end subroutine ! init_particles_statistics



cty========================================================================
cty     subroutine process_particles_statistics
cty========================================================================
cty
cty     Description:
cty     This routine processes the statistic of particles, which includes
cty     the subroutine collect_particles_statistics and the subroutine
cty     write_particles_statistics.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine process_particles_statistics

        use corevars, only : nstart, nit, nnstat
        use ompvars, only : start_master, end_master
 
        implicit none

        if (Para_par%Average_particle_steps .gt. 0) then

          if ( (mod(nit-nstart-1,nsprint_par) .eq. 0) .and.
     &         (nit .gt. nstart+1) ) then

cty         Collect particles data for computing statistics
            call collect_particles_statistics

            if (start_master()) then

              if (mod(nit-nstart-1,nnstat) .eq. 0) then
cty             Write statistics of particles
                call write_particles_statistics
              end if

            end if !start_master
            call end_master()

          end if ! ((mod(nit-nstart-1,nsprint_par) .eq. 0).and.(nit.gt.nstart+1))

        end if ! (Para_par%Average_particle_steps .gt. 0)

      end subroutine ! process_particles_statistics



cty========================================================================
cty     subroutine collect_particles_statistics
cty========================================================================
cty
cty     Description:
cty     This routine collect particles data for computing statistics.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine collect_particles_statistics

        implicit none
 
        integer :: num
        type(particle), pointer :: Pt_par
        integer :: ipar, jpar, kpar
        real(kind=wp) :: upar, vpar, wpar

        do num = 1, Npar_local
          Pt_par => Pointer_par(num)

          ipar = Pt_par%ipar 
          jpar = Pt_par%jpar 
          kpar = 1

          upar = Pt_par%upar 
          vpar = Pt_par%vpar 
          wpar = Pt_par%wpar

cty       ( 1) Sum the particles number
          qstat_par(ipar,jpar,kpar, 1) = 
     &    qstat_par(ipar,jpar,kpar, 1) + 1._wp
cty       ( 2) Sum the particles x-velocity upar
          qstat_par(ipar,jpar,kpar, 2) = 
     &    qstat_par(ipar,jpar,kpar, 2) + upar
cty       ( 3) Sum the particles y-velocity vpar
          qstat_par(ipar,jpar,kpar, 3) = 
     &    qstat_par(ipar,jpar,kpar, 3) + vpar 
cty       ( 4) Sum the particles z-velocity wpar
          qstat_par(ipar,jpar,kpar, 4) = 
     &    qstat_par(ipar,jpar,kpar, 4) + wpar
cty       ( 5) Sum the particles momentum flux upar*upar
          qstat_par(ipar,jpar,kpar, 5) = 
     &    qstat_par(ipar,jpar,kpar, 5) + upar*upar
cty       ( 6) Sum the particles momentum flux vpar*vpar
          qstat_par(ipar,jpar,kpar, 6) = 
     &    qstat_par(ipar,jpar,kpar, 6) + vpar*vpar
cty       ( 7) Sum the particles momentum flux wpar*wpar
          qstat_par(ipar,jpar,kpar, 7) = 
     &    qstat_par(ipar,jpar,kpar, 7) + wpar*wpar
cty       ( 8) Sum the particles momentum flux upar*vpar
          qstat_par(ipar,jpar,kpar, 8) = 
     &    qstat_par(ipar,jpar,kpar, 8) + upar*vpar
cty       ( 9) Sum the particles momentum flux upar*wpar
          qstat_par(ipar,jpar,kpar, 9) = 
     &    qstat_par(ipar,jpar,kpar, 9) + upar*wpar
cty       (10) Sum the particles momentum flux vpar*wpar
          qstat_par(ipar,jpar,kpar,10) = 
     &    qstat_par(ipar,jpar,kpar,10) + vpar*wpar

        end do ! num

      end subroutine ! collect_particles_statistics



cty========================================================================
cty     subroutine write_particles_statistics
cty========================================================================
cty
cty     Description:
cty     This routine write the statistic of particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine write_particles_statistics

        use corevars, only : nxp, nyp, nzp, nit
        use mbvars, only : num_blocks
     &                   , ov_bgrid_block, ov_ogrid_block
        use mppvars, only : ioproc
        use global_grid, only : grids

        implicit none

        character(15) :: id_name
        character(300) :: particle_file 
    
        integer :: var, i, j, k
        integer :: block_id

        call gather_particles_statistics

        if (ioproc) then
        
          do block_id = 1, num_blocks

cty         Set the file name of the particle_file
            write(id_name,'(I15)') nit-1
            write(particle_file,'(A,I0)') 
     &                          'Average_particles/Average_particle_b',
     &                           block_id
            particle_file = trim(particle_file)//'_'//
     &                      trim(adjustl(id_name))//'.bin'
            particle_file = trim(adjustl(particle_file))

cty         Open the particle_file
            open(unit=999,file=particle_file,
     &           form='unformatted',status='replace',
     &           action='write')  

            write(999) grids(block_id)%nxp_block,
     &                 grids(block_id)%nyp_block,
     &                 1,
     &                 n_qstat_par

            if (ov_bgrid_block(block_id)) then 

              do var = 1, n_qstat_par
                do k = 1, 1
                  do j = 1, grids(block_id)%nyp_block
                    do i = 1, grids(block_id)%nxp_block
                      write(999) B_qstat_par_gather(i,j,k,var)
                    end do ! i
                  end do ! j
                end do ! k
              end do ! var

            end if ! (ov_bgrid_block(block_id))

            if (ov_ogrid_block(block_id)) then

              do var = 1, n_qstat_par
                do k = 1, 1
                  do j = 1, grids(block_id)%nyp_block
                    do i = 1, grids(block_id)%nxp_block
                      write(999) O_qstat_par_gather(i,j,k,var)
                    end do ! i
                  end do ! j
                end do ! k
              end do ! var

            end if ! (ov_ogrid_block(block_id))

            close(999)

          end do ! block_id

        end if

        qstat_par(:,:,:,:) = 0._wp
        B_qstat_par_gather(:,:,:,:) = 0._wp
        O_qstat_par_gather(:,:,:,:) = 0._wp

      end subroutine ! write_particles_statistics



cty========================================================================
cty     subroutine gather_particles_statistics
cty========================================================================
cty
cty     Description:
cty     This routine gather the data of particles.
cty     -------------------------------------------------------------------
cty
cty     Input arguments:
cty
cty     Output arguments:
cty
cty     Notice:
cty
cty========================================================================

      subroutine gather_particles_statistics

        use corevars, only : nxp, nyp, nzp
        use mbvars, only : ov_bgrid_block, ov_ogrid_block
        use mppvars, only : real_mp_type, MPI_status_size, 
     &                      MPI_comm_fluid, ioid, procid_fluid, 
     &                      nprocs, ioproc
        use global_grid, only : grids

        implicit none

        integer :: block_id

        integer :: l, data_number, data_number_send, data_number_recv
        real(kind=wp), allocatable, dimension(:) :: data_send, data_recv

        integer :: var, i, j, k
        integer :: ierr
        integer :: status(MPI_status_size)

        integer :: p

cty     For the io id, we can assemble the qstat_par into B_qstat_par_gather directly
        if (ioproc) then
          do var = 1, n_qstat_par
            do k = 1, 1
              do j = 1, nyp
                do i = 1, nxp

                  B_qstat_par_gather(i,j,k,var) = 
     &              qstat_par(i,j,k,var)

                end do ! i
              end do ! j
            end do ! k
          end do ! var
        end if

        block_id = block_proc(procid_fluid)
cty     Set data_send and data_recv, which are large enough
        data_number = grids(block_id)%nxp_block*
     &                grids(block_id)%nyp_block*
     &                1*
     &                n_qstat_par
        allocate(data_send(data_number))
        allocate(data_recv(data_number))

        if (.not. ioproc) then

cty       For other ids, we need to assemble the qstat_par into data_send
          l = 0 
          do var = 1, n_qstat_par
            do k = 1, 1
              do j = 1, nyp
                do i = 1, nxp
                  l = l + 1
                  data_send(l) = qstat_par(i,j,k,var)
                end do ! i
              end do ! j
            end do ! k
          end do ! var

          data_number_send = nxp*nyp*1*n_qstat_par
cty       Send the data_send to the io id
          call MPI_send(data_send,data_number_send,real_mp_type,
     &                  ioid,procid_fluid,MPI_comm_fluid,ierr)

        else

          do p = 1, nprocs-1

            block_id = block_proc(p)

            data_number_recv = (iend_proc(p)-ibeg_proc(p)+1)*
     &                         (jend_proc(p)-jbeg_proc(p)+1)*
     &                          1*
     &                          n_qstat_par
cty         Receive the data_send in data_recv 
            call MPI_recv(data_recv,data_number_recv,real_mp_type,
     &                    p,p,MPI_comm_fluid,status,ierr)

            if (ov_bgrid_block(block_id)) then 

cty           Assemble the data_recv in 
              l = 0 
              do var = 1, n_qstat_par
                do k = 1, 1
                  do j = jbeg_proc(p), jend_proc(p)
                    do i = ibeg_proc(p), iend_proc(p)
                      l = l + 1

                        B_qstat_par_gather(i,j,k,var) = 
     &                    data_recv(l)

                    end do ! i
                  end do ! j
                end do ! k
              end do ! var

            end if

            if (ov_ogrid_block(block_id)) then

cty           Assemble the data_recv in 
              l = 0 
              do var = 1, n_qstat_par
                do k = 1, 1
                  do j = jbeg_proc(p), jend_proc(p)
                    do i = ibeg_proc(p), iend_proc(p)
                      l = l + 1

                      O_qstat_par_gather(i,j,k,var) = 
     &                  data_recv(l)

                    end do ! i
                  end do ! j
                end do ! k
              end do ! var

            end if

          end do ! p

        end if

        deallocate(data_send)
        deallocate(data_recv)

      end subroutine ! gather_particles_statistics



      end module

