program test
use differential_evolution
implicit none

real(rp), parameter :: bohr2angstrom      = 0.52917721067_rp
real(rp), parameter :: angstrom2bohr      = 1._rp/bohr2angstrom
real(rp), parameter :: hartree2kcal       = 627.509_rp
!https://de.wikipedia.org/wiki/Van-der-Waals-Radius
real(rp), parameter :: vdW_scaling = 1.0_rp/3.0_rp ! r must be smaller than scaling*vdW_radius to be feasible
real(rp), dimension(55), parameter :: vdW_radius = (/ &
                                     1.20_rp*angstrom2bohr, & ! H
                                     1.40_rp*angstrom2bohr, & ! He
                                     1.82_rp*angstrom2bohr, & ! Li
                                     1.53_rp*angstrom2bohr, & ! Be
                                     1.92_rp*angstrom2bohr, & ! B
                                     1.70_rp*angstrom2bohr, & ! C
                                     1.55_rp*angstrom2bohr, & ! N
                                     1.52_rp*angstrom2bohr, & ! O
                                     1.47_rp*angstrom2bohr, & ! F
                                     1.54_rp*angstrom2bohr, & ! Ne
                                     2.27_rp*angstrom2bohr, & ! Na
                                     1.73_rp*angstrom2bohr, & ! Mg
                                     1.84_rp*angstrom2bohr, & ! Al
                                     2.10_rp*angstrom2bohr, & ! Si
                                     1.80_rp*angstrom2bohr, & ! P
                                     1.80_rp*angstrom2bohr, & ! S
                                     1.75_rp*angstrom2bohr, & ! Cl
                                     1.88_rp*angstrom2bohr, & ! Ar
                                     2.75_rp*angstrom2bohr, & ! K
                                     2.31_rp*angstrom2bohr, & ! Ca
                                     0.00_rp*angstrom2bohr, & ! Sc
                                     0.00_rp*angstrom2bohr, & ! Ti
                                     0.00_rp*angstrom2bohr, & ! V
                                     0.00_rp*angstrom2bohr, & ! Cr
                                     0.00_rp*angstrom2bohr, & ! Mn
                                     0.00_rp*angstrom2bohr, & ! Fe
                                     0.00_rp*angstrom2bohr, & ! Co
                                     1.63_rp*angstrom2bohr, & ! Ni
                                     1.40_rp*angstrom2bohr, & ! Cu
                                     1.39_rp*angstrom2bohr, & ! Zn
                                     1.87_rp*angstrom2bohr, & ! Ga
                                     2.11_rp*angstrom2bohr, & ! Ge
                                     1.85_rp*angstrom2bohr, & ! As
                                     1.90_rp*angstrom2bohr, & ! Se
                                     1.85_rp*angstrom2bohr, & ! Br
                                     2.02_rp*angstrom2bohr, & ! Kr
                                     3.03_rp*angstrom2bohr, & ! Rb
                                     2.49_rp*angstrom2bohr, & ! Sr
                                     0.00_rp*angstrom2bohr, & ! Y
                                     0.00_rp*angstrom2bohr, & ! Zr
                                     0.00_rp*angstrom2bohr, & ! Nb
                                     0.00_rp*angstrom2bohr, & ! Mo
                                     0.00_rp*angstrom2bohr, & ! Tc
                                     0.00_rp*angstrom2bohr, & ! Ru
                                     0.00_rp*angstrom2bohr, & ! Rh
                                     1.63_rp*angstrom2bohr, & ! Pd
                                     1.72_rp*angstrom2bohr, & ! Ag
                                     1.58_rp*angstrom2bohr, & ! Cd
                                     1.93_rp*angstrom2bohr, & ! In
                                     2.17_rp*angstrom2bohr, & ! Sn
                                     2.06_rp*angstrom2bohr, & ! Sb
                                     2.06_rp*angstrom2bohr, & ! Te
                                     1.98_rp*angstrom2bohr, & ! I
                                     2.16_rp*angstrom2bohr, & ! Xe
                                     0.00_rp*angstrom2bohr  &
                                     /)

! program parameters
logical  :: verbose           = .false. ! toggles printing mode
logical  :: use_greedy_fit    = .false. ! greedily fit individual multipoles first
logical  :: greedy_only_multi = .false. ! stops greedy fit after fitting atomic multipoles
!logical  :: use_symmetry      = .false. ! toggles symmetry constraint mode
logical  :: fit_multipoles    = .false. ! fit multipoles instead of charges
logical  :: generate_atomic   = .false. ! for visualizing charge fits
logical  :: refine_solution   = .false. ! when used to refine a solution
integer, parameter :: lmax    = 5       ! where is the multipole expansion truncated? maximum = 5!
integer  :: lstart            = 0       ! where do we start the fitting of the multipole expansion?
integer  :: lstop             = lmax    ! where do we start the fitting of the multipole expansion?
integer  :: lcur              = 0       ! where is the multipole expansion truncated? maximum = 5!
integer  :: num_charges       = 1       ! how many charges to fit (current)
integer  :: num_charges_min   = 2       ! maximum number of charges to fit
integer  :: num_charges_max   = 2       ! maximum number of charges to fit
integer  :: num_charges_min_multipole = 1 
integer  :: num_charges_max_multipole = 5
integer  :: num_trials        = 1       ! maximum number of trials per number of charges
real(rp) :: total_charge      = 0._rp   ! total charge
real(rp) :: total_charge2     = 0._rp   ! total charge
real(rp) :: max_charge        = 1._rp  ! maximum allowed absolute value of a charge (in search range)
real(rp) :: max_extend        = 5._rp   ! gets calculated automatically from vdW radii

character(len=1024) :: input_esp_cubefile = '', &
                       compare_esp_cubefile = '', &
                       input_multipolefile = '', &
                       input_xyzfile = '', &
                       input_density_cubefile = '', &
                       prefix = '', &
                       dummystring = '', dummystring2 = '', dummystring3 = '' ! input files
integer :: ppos !for use with scan()

! other program "modes", only used for analysis or generating cube files at better resolution
logical :: analysis_mode = .false. ! analysis mode: compare two cube files and compute RMSDs
logical :: generate_mode = .false. ! for generating esp cube file from charges (at higher resolution)


!how to "prune" grid points
logical  :: use_vdW_grid_cutoff = .true.
real(rp) :: vdw_grid_min_cutoff = 1.20_rp      !radius smaller than this*vdw_radius gets ignored
real(rp) :: vdw_grid_max_cutoff = 2.20_rp      !radius larger than this*vdw_radius gets ignored
logical  :: use_density_grid_cutoff = .false.
real(rp) :: density_grid_min_cutoff = 3.162277660e-4_rp !density smaller than this is ignored
real(rp) :: density_grid_max_cutoff = 2.0e-3 ! 1.0e-3_rp !density larger than this is ignored


! information about the ESP
integer :: Natom ! number of atoms
real(rp), dimension(3)                :: origin, axisX, axisY, axisZ   ! coordinate system of grid
integer                               :: NgridX, NgridY, NgridZ, Ngrid ! number of grid points in each direction
real(rp)                              :: Ngridr                        ! total number of grid points as real (for mean etc.)
real(rp), dimension(:),   allocatable :: esp_grid, esp_grid2           ! values of the electrostatic potential
real(rp), dimension(:,:), allocatable :: gridval                       ! stores at which gridpoints the ESP is interesting (outside atoms)
integer,  dimension(:),   allocatable :: atom_num                      ! stores the atomic numbers of the atoms
real(rp), dimension(:,:), allocatable :: atom_pos                      ! stores the atomic positions 
real(rp), dimension(3)                :: atom_com                      ! stores the atomic center of mass (for checking charge symmetry)
real(rp), dimension(:), allocatable :: multipole                       ! stores multipole parameters

real(rp), dimension(:,:,:), allocatable :: multipole_solutions         ! stores charge positions for the multipole fits (num_charges_max_multipole*4,num_charges_max_multipole,Natom)
real(rp), dimension(:,:), allocatable   :: multipole_solutions_rmse    ! stores RMSE of the multipole solutions (needed for greedy mode)
real(rp), dimension(:), allocatable :: multipole_best, multipole_save  ! current best charge fit in the multipole fit
real(rp), dimension(:,:,:), allocatable, save :: symmetry_ops          ! stores the symmetry operations in matrix form

! these are only used for visualizing the grid with R in the end
real(rp), dimension(:,:), allocatable :: sliceXY,sliceXZ,sliceYZ       ! cuts along planes
real(rp), dimension(:,:), allocatable :: sliceXY2,sliceXZ2,sliceYZ2    ! cuts along planes
logical,  dimension(:,:), allocatable :: usedXY, usedXZ, usedYZ        ! is the grid point at that cut used in the fit?

! fitting variables
integer ::  natmfit = 0                               ! defines number of atoms included in fit
integer,  dimension(:),   allocatable :: fitatoms     ! contains list of atom indices to be fitted
real(rp), dimension(:),   allocatable :: charges      ! contains position (x,y,z) and magnitude (q) of charges
                                                      ! x1: charge(1), y1: charge(2), z1: charge(3), q1: charge(4)
                                                      ! x2: charge(5), ...
real(rp), dimension(:),   allocatable :: bestcharges  ! best solution found so far
real(rp), dimension(:,:), allocatable :: search_range ! has a minimum and a maximum value for every entry of charges


! for error handling
integer :: ios

integer :: i,j,k,l,a,b,g,try,qdim,lcheck,n_steps !current dimensionality of charge

real(rp) :: deriv, tmp, tmp2, step_size, stepsize, learning_rate, RMSE_best, RMSE_tmp, MAE_tmp, maxAE_tmp, RMSE_a1, RMSE_a2

integer :: cmd_count
character(len=1024) :: arg, bla

call system("mkdir -p slices")

! read command line arguments
! call read_command_line_arguments()

!make a folder to store output for cleanliness
if(trim(prefix) /= '') then
    call execute_command_line("mkdir -p "//trim(prefix))
end if 

! ! initialize RNG (for drawing the population)
! call init_random_seed(0)

! ! read in the cube file
! call read_cube_file(trim(input_esp_cubefile),trim(input_density_cubefile))

! initialize list of atoms to be fitted if nothing was defined with -atom flag
if(natmfit == 0) then
  natmfit = Natom
  allocate(fitatoms(natmfit))
  do i=1,Natom
    fitatoms(i)=i
  enddo
endif

stepsize = 0.1
learning_rate = 1 ! https://rosettacode.org/wiki/Gradient_descent#Fortran
n_steps = 1000

!read argument count
cmd_count = command_argument_count()

! loop through command line arguments
do i = 1,cmd_count
    call get_command_argument(i, arg, l)
    write(*,'(A)'), arg(1:l)

        ! input xyz cube file
    if(arg(1:l) == '-xyz') then
        call get_command_argument(i+1, arg, l)
        read(arg,'(A)',iostat = ios) input_xyzfile
        if(ios /= 0) call throw_error('Could not read command line argument "-xyz"')
        call get_command_argument(i+2, arg, l)
        read(arg,'(A)',iostat = ios) input_multipolefile
    end if

    if(arg(1:l) == '-stepsize') then
        call get_command_argument(i+1, arg, l)
        read(arg,*,iostat = ios) stepsize
    end if
    
    if(arg(1:l) == '-n_steps') then
        call get_command_argument(i+1, arg, l)
        read(arg,*,iostat = ios) n_steps
    end if
    
    if(arg(1:l) == '-learningrate') then
        call get_command_argument(i+1, arg, l)
        read(arg,*,iostat = ios) learning_rate
    end if

    ! input dens cube file
    if(arg(1:l) == '-dens') then
        call get_command_argument(i+1, arg, l)
        read(arg,'(A)',iostat = ios) input_density_cubefile
        if(ios /= 0) call throw_error('Could not read command line argument "-dens"')
    end if

    ! input esp cube file
    if(arg(1:l) == '-esp') then
        call get_command_argument(i+1, arg, l)
        read(arg,'(A)',iostat = ios) input_esp_cubefile
        if(ios /= 0) call throw_error('Could not read command line argument "-esp"')
    end if

end do

vdw_grid_min_cutoff = 1.2_rp
vdw_grid_max_cutoff = 100._rp
!
! Read cube files
call read_cube_file(trim(input_esp_cubefile),trim(input_density_cubefile))

open(30, file=trim(input_xyzfile), status="old", action="read", iostat = ios)
if(ios /= 0) call throw_error('Could not open "'//trim(input_xyzfile)//'" for reading')
read(30,*,iostat=ios) num_charges
allocate(charges(4*num_charges-1), bestcharges(4*num_charges-1), search_range(2,4*num_charges-1), stat=ios)
qdim = 4*num_charges-1
if((ios /= 0).or.(num_charges < 1)) call throw_error('"'//trim(input_xyzfile)//'" has the wrong format.')
read(30,*,iostat=ios) !skip comment line
total_charge  = 0._rp
!read charge coordinates and magnitude (also automatically determines total charge)
do i = 1,qdim,4
    if(i+3 < qdim) then
        read(30,*) dummystring, charges(i:i+3)
        total_charge  = total_charge  + charges(i+3)
    else
        read(30,*) dummystring, charges(i:i+2), tmp
        total_charge = total_charge + tmp
    end if
    charges(i:i+2) = charges(i:i+2)*angstrom2bohr
    ! Print the charges
    print*, charges(i:i+2)
end do

read(30,*,iostat=ios) !skip blank line
!    read(30,*) dummystring, RMSE_best !read RMSE (commented as we may be using
!    fragments, in which case RMSE in reconstructed file is incorrect)
close(30)


write(*, '(A)') "First" 
RMSE_tmp = rmse_qtot(charges(1:qdim))
write(*,'(A30,2ES23.9,I10)') "Total", RMSE_tmp*hartree2kcal !sqrt(rmse_tot/npts)*hartree2kcalmol
RMSE_best = RMSE_tmp

write(*, '(A,F10.2)') "Step size: ", stepsize
step_size = stepsize * angstrom2bohr
write(*, '(A,F10.2)') "Step size in Bohr: ", step_size

write(*, '(A,F10.2)') "Learning Rate: ", learning_rate
write(*, '(A,I10)') "# Steps: ", n_steps

! Loop gradient descent for n steps
do g=1,n_steps,1
    do i = 1,qdim,1   
        if (mod(i, 4) > 0) then
        ! print*, i, mod(i, 4)

        ! Numerical gradient
        charges(i) = charges(i) + step_size ! take a step forward
        RMSE_a1 = rmse_qtot(charges(1:qdim)) ! calculate RMSD
        charges(i) = charges(i) - 2 * step_size ! take two steps backwards
        RMSE_a2 = rmse_qtot(charges(1:qdim)) ! calculate RMSD
        charges(i) = charges(i) + step_size ! reset position
        deriv = (RMSE_a1 - RMSE_a2)/2*step_size ! change in Loss function versus step size (two steps)
        
        !write(*,'(A30,2ES23.9,I10)') "E dif/step", deriv !hartree2kcal !sqrt(rmse_tot/npts)*hartree2kcalmol
        !write(*,'(A30,2ES23.9,I10)') "x_start", charges(i)*bohr2angstrom
        charges(i) = charges(i) - learning_rate * deriv
        !write(*,'(A30,2ES23.9,I10)') "x_final", charges(i)*bohr2angstrom
        !RMSE_a1 = rmse_qtot(charges(1:qdim))
        !write(*,'(A30,2ES23.9,I10)') "Error", RMSE_a1*hartree2kcal
        end if
    end do
    RMSE_a1 = rmse_qtot(charges(1:qdim))
    write(*,'(A30,I3,2ES23.9,I10)') "Error", g, RMSE_a1*hartree2kcal
end do

!  "Total energy"
RMSE_a1 = rmse_qtot(charges(1:qdim))
write(*,'(A30,2ES23.9,I10)') "Error", RMSE_a1*hartree2kcal

call write_xyz_file(charges(1:qdim),filename="refined.xyz")


contains

!
!    Functions
!    F  
!    F
!    FFFFFFFFF       U  N C T I O N S 
!    F
!    F
!    F


!-------------------------------------------------------------------------------
! computes root mean squared error of the current fit to the true esp
real(rp) function rmse(q)
    implicit none
    real(rp), dimension(:) :: q ! input charges
    integer :: idx
    rmse = 0._rp
    do idx = 1,Ngrid
        rmse = rmse + (coulomb_potential(gridval(:,idx),q) - esp_grid(idx))**2
    end do
    rmse = sqrt(rmse/Ngridr)
end function rmse
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! computes Coulomb potential for the given charges q at position x
real(rp) function coulomb_potential(x,q)
    implicit none
    real(rp), dimension(3), intent(in) :: x ! position
    real(rp), dimension(:) :: q ! charges
    real(rp) :: r
    integer :: i
    coulomb_potential = 0._rp
    do i=1,size(q,dim=1),4
        r = sqrt(sum((q(i:i+2)-x)**2))  ! compute distance
        if(r < 1.e-9_rp) r = 1.e-9_rp ! prevent division by 0
        coulomb_potential = coulomb_potential + q(i+3)/r
    end do
end function coulomb_potential
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! computes root mean squared error of the current fit to the true esp, using constraint charges
! (this means that the total charge must add up to a specific value)
real(rp) function rmse_qtot(qin)
    implicit none
    real(rp), dimension(:) :: qin ! input charges
    real(rp), dimension(size(qin,dim=1)+1) :: q   ! complete charges
    real(rp) :: qsum
    integer :: i
    qsum = 0._rp
    do i = 1,size(qin,dim=1)-3,4
        qsum = qsum + qin(i+3)
    end do
    q = qin
    q(size(q,dim=1)) = total_charge - qsum
    rmse_qtot = rmse(q)
end function rmse_qtot
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read Gaussian cube files
subroutine read_cube_file(filepath,density_filepath)
    implicit none
    character(len=*), intent(in)           :: filepath
    character(len=*), intent(in), optional :: density_filepath !for density cutoff
    character(len=128) :: ctmp ! dummy character variable
    real(rp) :: x(3) ! dummy vector
    real(rp) :: rtmp,rtmp2 ! dummy real variable
    integer :: ios ! keeps track of io status
    integer :: i,j,k,idx,lcount
    
    if(use_vdW_grid_cutoff .and. use_density_grid_cutoff) &
        call throw_error('Only one grid cutoff scheme can be selected at once.')
        
    if(use_density_grid_cutoff.and..not.present(density_filepath)) &
        call throw_error('Density grid cutoff scheme was selected, but no density cube file was given.')
    
    
    
    if(verbose) write(*,'(A)') 'Reading "'//trim(filepath)//'"...'
    if(verbose) write(*,*)
    
    !open cube files
    open(30, file=trim(filepath), status= "old", action= "read", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "'//trim(filepath)//'".')
    if(use_density_grid_cutoff) then
        open(31, file=trim(density_filepath), status= "old", action= "read", iostat = ios)
        if(ios /= 0) call throw_error('Could not open "'//trim(density_filepath)//'".')
    end if
    
    ! skip the title in the header
    read(30,'(A128)',iostat = ios) ctmp
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Bad Format.')
    if(verbose) write(*,*) trim(ctmp)
    read(30,'(A128)',iostat = ios) ctmp
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Bad Format.')
    if(verbose) write(*,*) trim(ctmp)
    if(use_density_grid_cutoff) then
        read(31,'(A128)',iostat = ios) ctmp
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Bad Format.')
        read(31,'(A128)',iostat = ios) ctmp
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Bad Format.')
    end if
    
    ! read information about coordinate system and the number of atoms
    read(30,*,iostat = ios) Natom, origin(:)
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Missing Header?')
    read(30,*,iostat = ios) NgridX, axisX(:)
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Missing Header?')
    read(30,*,iostat = ios) NgridY, axisY(:)
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Missing Header?')
    read(30,*,iostat = ios) NgridZ, axisZ(:)
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Missing Header?')
    if(use_density_grid_cutoff) then ! skip this information in density cube file
        read(31,*,iostat = ios) 
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Missing Header?')
        read(31,*,iostat = ios) 
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Missing Header?')
        read(31,*,iostat = ios) 
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Missing Header?')
        read(31,*,iostat = ios) 
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Missing Header?')
    end if
    

    if(verbose) write(*,'(I5,3F12.6)') Natom, origin(:)
    if(verbose) write(*,'(I5,3F12.6)') NgridX, axisX(:)
    if(verbose) write(*,'(I5,3F12.6)') NgridY, axisY(:)
    if(verbose) write(*,'(I5,3F12.6)') NgridZ, axisZ(:) 
    
    ! allocate memory to store slices through the potential (for visualization with R)
    if(.not.allocated(sliceXY)) allocate(sliceXY(NgridX,NgridY), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for sliceXY.')
    if(.not.allocated(sliceXZ)) allocate(sliceXZ(NgridX,NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for sliceXZ.')
    if(.not.allocated(sliceYZ)) allocate(sliceYZ(NgridY,NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for sliceYZ.')
    if(.not.allocated(usedXY))  allocate(usedXY (NgridX,NgridY), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for usedXY.')
    if(.not.allocated(usedXZ))  allocate(usedXZ (NgridX,NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for usedXZ.')
    if(.not.allocated(usedYZ))  allocate(usedYZ (NgridY,NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for usedYZ.')
    usedXY = .true.; usedXZ = .true.; usedYZ = .true.; !initialize to true
    
    
    ! allocate memory to store atom information
    if(.not.allocated(atom_num)) allocate(atom_num(Natom),   stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    if(.not.allocated(atom_pos)) allocate(atom_pos(3,Natom), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    
    ! read atom information
    do i = 1,Natom
        read(30,*,iostat=ios) atom_num(i), rtmp, atom_pos(:,i)
        if(ios /= 0) call throw_error('Could not read atom information.')
        if(verbose) write(*,'(I5,4F12.6)') atom_num(i), real(atom_num(i),rp), atom_pos(:,i)
    end do
    if(use_density_grid_cutoff) then ! skip atom information in density cube file
        do i = 1,Natom
            read(31,*,iostat=ios) 
            if(ios /= 0) call throw_error('Density cube file has unmatching format.')
        end do
    end if
        
    ! allocate memory to store grid information
    if(.not.allocated(esp_grid)) allocate(esp_grid(NgridX*NgridY*NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    if(.not.allocated(gridval))  allocate(gridval(3,NgridX*NgridY*NgridZ),stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    
    
    ! find the "interesting" grid points
    if(use_vdW_grid_cutoff) then !based on vdW radii                
        ! read grid information (the grid is stored in a flat array)
        ! for cache friendliness (looping over these values is
        ! performance critical for cost function evaluation)
        idx = 0
        lcount = 0
        do i = 1,NgridX
            do j = 1,NgridY
                do k = 1,NgridZ
                    lcount = lcount + 1
                    x = origin + (i-1)*axisX + (j-1)*axisY + (k-1)*axisZ
                    if(in_interaction_belt(x,vdw_grid_min_cutoff, &
                       vdw_grid_max_cutoff)) then !value gets added
                        idx = idx + 1
                        read(30,'(ES13.5)',advance='no',iostat = ios) esp_grid(idx)
                        if(ios /= 0) call throw_error('Could not read ESP grid.')
                        gridval(:,idx) = x
                    else
                        read(30,'(ES13.5)',advance='no',iostat = ios) rtmp
                        if(ios /= 0) call throw_error('Could not read ESP grid.')
                    end if
                    
                    !to accomodate weird format, we sometimes have to skip lines
                    if(mod(lcount,6) == 0 .or. mod(lcount,NgridZ) == 0) then
                        read(30,*) ! line break
                        if(mod(lcount,NgridZ) == 0) lcount = 0
                    end if
                    
                    ! everything below is for later visualization with R (very useful to assess quality)
                    if( i == NgridX/2 ) then  
                        if(.not.in_interaction_belt(x,vdw_grid_min_cutoff, &
                           vdw_grid_max_cutoff)) then !value is NOT used for fitting
                            usedYZ(j,k)  = .false.
                            sliceYZ(j,k) = rtmp
                        else
                            sliceYZ(j,k) = esp_grid(idx)
                        end if
                    end if
                    if( j == NgridY/2 ) then
                        if(.not.in_interaction_belt(x,vdw_grid_min_cutoff, &
                           vdw_grid_max_cutoff)) then !value is NOT used for fitting
                            usedXZ(i,k)  = .false.
                            sliceXZ(i,k) = rtmp
                        else
                            sliceXZ(i,k) = esp_grid(idx)
                        end if
                    end if
                    if( k == NgridZ/2 ) then
                        if(.not.in_interaction_belt(x,vdw_grid_min_cutoff, &
                           vdw_grid_max_cutoff)) then !value is NOT used for fitting
                            usedXY(i,j)  = .false.
                            sliceXY(i,j) = rtmp
                        else
                            sliceXY(i,j) = esp_grid(idx)
                        end if
                    end if
                end do
            end do
        end do  
    else if(use_density_grid_cutoff) then ! based on density cube file        
        idx = 0
        lcount = 0
        do i = 1,NgridX
            do j = 1,NgridY
                do k = 1,NgridZ
                    lcount = lcount + 1
                    read(31,'(ES13.5)',advance='no',iostat = ios) rtmp
                    if(ios /= 0) call throw_error('Density cube file has wrong format.')
                    if((rtmp < density_grid_max_cutoff) .and. (rtmp > density_grid_min_cutoff)) then !add point
                        x = origin + (i-1)*axisX + (j-1)*axisY + (k-1)*axisZ
                        idx = idx + 1
                        read(30,'(ES13.5)',advance='no',iostat = ios) esp_grid(idx)
                        if(ios /= 0) call throw_error('Could not read ESP grid.')
                        gridval(:,idx) = x
                    else !skip point
                        read(30,'(ES13.5)',advance='no',iostat = ios) rtmp2
                        if(ios /= 0) call throw_error('Could not read ESP grid.')
                    end if
                    !to accomodate weird format, we sometimes have to skip lines
                    if(mod(lcount,6) == 0 .or. mod(lcount,NgridZ) == 0) then
                        read(30,*) ! line break
                        read(31,*)
                        if(mod(lcount,NgridZ) == 0) lcount = 0
                    end if
                    
                    if( i == NgridX/2 ) then  
                        if(.not.((rtmp < density_grid_max_cutoff) .and. (rtmp > density_grid_min_cutoff))) then !value is NOT used for fitting
                            usedYZ(j,k)  = .false.
                            sliceYZ(j,k) = rtmp2
                        else
                            sliceYZ(j,k) = esp_grid(idx)
                        end if
                    end if
                    if( j == NgridY/2 ) then
                        if(.not.((rtmp < density_grid_max_cutoff) .and. (rtmp > density_grid_min_cutoff))) then !value is NOT used for fitting
                            usedXZ(i,k)  = .false.
                            sliceXZ(i,k) = rtmp2
                        else
                            sliceXZ(i,k) = esp_grid(idx)
                        end if
                    end if
                    if( k == NgridZ/2 ) then
                        if(.not.((rtmp < density_grid_max_cutoff) .and. (rtmp > density_grid_min_cutoff))) then !value is NOT used for fitting
                            usedXY(i,j)  = .false.
                            sliceXY(i,j) = rtmp2
                        else
                            sliceXY(i,j) = esp_grid(idx)
                        end if
                    end if
                    
                end do
            end do
        end do  
        close(31)
    end if
    !store number of interesting grid points
    Ngrid = idx  
    Ngridr = real(Ngrid,rp) 
    
    if(verbose) write(*,'(6ES13.5)') esp_grid(1:6)
    if(verbose) write(*,'(6(A8,5X))') ".",".",".",".",".","."
    if(verbose) write(*,'(6(A8,5X))') ".",".",".",".",".","."
    if(verbose) write(*,'(6(A8,5X))') ".",".",".",".",".","."
    if(verbose) write(*,'(6ES13.5)') esp_grid(Ngrid-5:Ngrid)
    if(verbose) write(*,*)
    if(verbose) write(*,'(I0,A,I0,A,F4.1,A)') Ngrid, ' out of ', NgridX*NgridY*NgridZ, &
            ' gridpoints are considered (',100*Ngrid/real(NgridX*NgridY*NgridZ,rp),'%).'
    
    close(30)
    
    if(verbose) write(*,*)
    if(verbose) write(*,'(A)') '...done!'
    if(verbose) write(*,*)

    return
end subroutine read_cube_file
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! deallocate memory
subroutine dealloc()
    implicit none
    if(allocated(atom_num))     deallocate(atom_num)
    if(allocated(atom_pos))     deallocate(atom_pos)
    if(allocated(esp_grid))     deallocate(esp_grid)
    if(allocated(charges))      deallocate(charges)
    if(allocated(bestcharges))  deallocate(bestcharges)
    if(allocated(search_range)) deallocate(search_range)
    if(allocated(sliceXY))      deallocate(sliceXY)
    if(allocated(usedXY))       deallocate(usedXY)
    if(allocated(sliceXZ))      deallocate(sliceXZ)
    if(allocated(usedXZ))       deallocate(usedXZ)
    if(allocated(sliceYZ))      deallocate(sliceYZ)
    if(allocated(usedYZ))       deallocate(usedYZ)   
end subroutine dealloc
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! throws an error message and terminates the code
subroutine throw_error(message)
    implicit none
    character(len=*), intent(in) :: message
    write(*,'(A)') "ERROR: "//message
    call dealloc()
    stop   
end subroutine throw_error
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! checks whether a point is inside the interaction belt as defined by the vdW radii
! Atomic Multipoles: Electrostatic Potential Fit, Local Reference Axis Systems,
! and Conformational Dependence
logical function in_interaction_belt(q,mincut,maxcut)
    implicit none
    real(rp), dimension(:) :: q
    real(rp) :: r, rmin ! shortest distance to atom
    real(rp) :: mincut, maxcut ! minimum and maximum grid cut-offs passed
    integer  :: i, a
    ! loop over the charge coordinates 
    do i = 1,size(q,dim=1),4
        !find atom with minimal relative distance
        rmin = huge(0._rp) 
        do a = 1,Natom
            r = sqrt(sum((atom_pos(:,a)-q(i:i+2))**2))/(vdW_radius(atom_num(a)))
            !print*, atom_num(a), vdW_radius(atom_num(a))*bohr2angstrom
            if(r < rmin) rmin = r
        end do
        ! this means the radius is not in the defined interaction cutoff
        if(.not.((rmin >= mincut) .and. (rmin <= maxcut))) then
                in_interaction_belt = .false.
                !print*, r
                return
        end if
    end do
    in_interaction_belt = .true.
    return
end function in_interaction_belt
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! writes a xyz file containing the results
subroutine write_xyz_file(charges,a,filename)
    implicit none 
    integer :: ios
    integer, optional :: a
    real(rp), dimension(qdim), intent(in) :: charges
    character(len=*), intent(in), optional :: filename
    character(len=1024) :: outfile, dummy
    
    write(outfile,'(I0)') num_charges
    if(present(a)) then
        write(dummy,'(I0)') a
        outfile = "multipole"//trim(dummy)//"_"//trim(outfile)//"charges.xyz"
    else
        outfile = trim(outfile)//"charges.xyz"
    end if
    if(trim(prefix) /= '') outfile = trim(prefix)//"/"//trim(outfile)
    
    if(present(filename)) outfile = filename
    
    
    open(30, file=trim(outfile), status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "'//trim(outfile)//'" for writing')
    write(30,'(I0)') num_charges
    write(30,'(A,4A26)')  "s","x[A]","y[A]","z[A]","q[e]"
    tmp = 0._rp
    do i = 1,qdim,4
        if(i+3 <= qdim) then
            if(charges(i+3) > 0._rp) then
                write(30,'(A,1X,4(F25.16,1X))') "N",charges(i:i+2)*bohr2angstrom,charges(i+3)
            else
                write(30,'(A,1X,4(F25.16,1X))') "O",charges(i:i+2)*bohr2angstrom,charges(i+3)
            end if
            tmp = tmp + charges(i+3)
        else
            tmp = total_charge-tmp
            if(tmp > 0._rp) then
                write(30,'(A,1X,4(F25.16,1X))') "N",charges(i:i+2)*bohr2angstrom,tmp
            else
                write(30,'(A,1X,4(F25.16,1X))') "O",charges(i:i+2)*bohr2angstrom,tmp 
            end if
        end if
    end do
    write(30,*)
    write(30,'(A,ES23.9,A)') "        RMSE ", &
                      RMSE_tmp*hartree2kcal," kcal/mol"
    ! write(30,'(A,ES23.9,A)') "         MAE ", &
    !                   MAE_tmp*hartree2kcal," kcal/mol"
    ! write(30,'(A,ES23.9,A)') "     max. AE ", &
    !                   maxAE_tmp*hartree2kcal," kcal/mol"
    write(30,*)
    write(30,'(A)') "Coordinates in bohr"
    write(30,'(A,4A26)')  "s","x[bohr]","y[bohr]","z[bohr]","q[e]"
    tmp = 0._rp
    do i = 1,qdim,4
        if(i+3 <= qdim) then
            if(charges(i+3) > 0._rp) then
                write(30,'(A,1X,4(F25.16,1X))') "+",charges(i:i+2),charges(i+3)
            else
                write(30,'(A,1X,4(F25.16,1X))') "-",charges(i:i+2),charges(i+3)
            end if
            tmp = tmp + charges(i+3)
        else
            tmp = total_charge-tmp
            if(tmp > 0._rp) then
                write(30,'(A,1X,4(F25.16,1X))') "+",charges(i:i+2),tmp
            else
                write(30,'(A,1X,4(F25.16,1X))') "-",charges(i:i+2),tmp 
            end if
        end if
    end do
    close(30)
    end subroutine write_xyz_file
    !-------------------------------------------------------------------------------
   

end program test
