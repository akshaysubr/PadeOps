module RoutinesRotation
    use kind_parameters, only: rkind, mpirkind
    use exits, only: GracefulExit
    use constants, only : one, two
    use decomp_2d
    use decomp_2d_io
    use mpi 
    
    implicit none
    contains

    subroutine get_xy_meanC_from_fC(f, fmean, nx, ny, nz, ftemp, ftemp_in_Y, ftemp_in_Z, gpC)
        real(rkind), dimension(:,:,:), intent(in)    :: f
        real(rkind), dimension(:),     intent(out)   :: fmean
        integer, intent(in) :: nx, ny, nz
        real(rkind), dimension(:,:,:), intent(inout) :: ftemp, ftemp_in_Y, ftemp_in_Z
        type(decomp_info) :: gpC

        integer :: ierr, j, k
        real(rkind) :: avgFact

        avgFact = 1.d0/real(nx*ny,rkind)

        ! take average in x
        do k = 1, size(f,3)
          do j = 1, size(f,2)
              ftemp(1,j,k) = sum(f(1:nx,j,k))
          enddo
        enddo

        ! transpose x to y
        call transpose_x_to_y(ftemp, ftemp_in_Y, gpC)

        ! take average in y
        do k = 1, size(ftemp_in_Y,3)
          do j = 2, size(ftemp_in_Y,2)
              ftemp_in_Y(1,1,k) = ftemp_in_Y(1,1,k) + ftemp_in_Y(1,j,k)
          enddo
        enddo

        ! transpose y to z
        call transpose_y_to_z(ftemp_in_Y, ftemp_in_Z, gpC)

        if (nrank == 0) then
            fmean = ftemp_in_Z(1,1,:)*avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 
        
        call mpi_bcast(fmean, nz, mpirkind, 0, mpi_comm_world, ierr)
        !end if

    end subroutine 
    
end module



program rotateFields
    use kind_parameters, only: rkind, clen
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use constants, only: deg_to_radians
    use mpi 
    use timer, only: tic, toc
    use decomp_2d_io
    use RoutinesRotation
   
    implicit none

    character(len=clen) :: inputfile 
    character(len=clen) :: outputdir, inputdir
    integer :: ioUnit, nx, ny, nz, ierr, inputFile_TID, inputFile_RID, outputFile_TID, outputFile_RID, i1, i2, i, j, k
    type(decomp_info) :: gpC, gpE
    real(rkind), dimension(:,:,:), allocatable :: u, v, urot, vrot, ftemp, ftemp_in_Y, ftemp_in_Z
    real(rkind), dimension(:),     allocatable :: zline, umean, vmean, urotmean, vrotmean
    character(len=clen) :: tempname, fname
    real(rkind) :: tsim, zH, umean_at_hub, vmean_at_hub, angle_at_hub, angle_rotate
    real(rkind) :: angle_with_x_deg, cos_rotangle, sin_rotangle, Lz, dz 
    namelist /INPUT/ nx, ny, nz, inputdir, outputdir, inputFile_TID, inputFile_RID, &
                     outputFile_TID, outputFile_RID, zH, Lz, angle_with_x_deg

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)

    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)
  

    !!!!!!! READ HEADER !!!!!!!!!!!
    write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",inputFile_RID, "_info.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=10,file=fname,access='sequential',form='formatted')
    read (10, *)  tsim
    close(10)
    call message(0, "Upsampling File data dumped at tSIM=", tsim)

    !!!!!!!!!!!!! CELL FIELDS !!!!!!!!!!!!!!!
    allocate(u(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(v(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(urot(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(vrot(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))

    allocate(ftemp(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(ftemp_in_Y(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
    allocate(ftemp_in_Z(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))
    allocate(zline(gpC%zsz(3)))
    allocate(umean(gpC%zsz(3)))
    allocate(vmean(gpC%zsz(3)))
    allocate(urotmean(gpC%zsz(3)))
    allocate(vrotmean(gpC%zsz(3)))
   
    ! Set up z grid
    dz = Lz / real(nz, rkind)
    do i = 1, nz
      zline(i) = (i-0.5d0) * dz
    enddo
 
    ! Read u and v fields
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_u.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,u,fname, gpC)

    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_v.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,v,fname, gpC)

    !! --- STEP 1 :: Calculate horizontal averages of u and v fields --- 
    call get_xy_meanC_from_fC(u, umean, nx, ny, nz, ftemp, ftemp_in_Y, ftemp_in_Z, gpC)
    call get_xy_meanC_from_fC(v, vmean, nx, ny, nz, ftemp, ftemp_in_Y, ftemp_in_Z, gpC)

    !! --- STEP 2 :: Calculate the angle at desired height ---
    !zH = 100.d0
    i1 = minloc(abs(zline-zH),1)
    if(zline(i1) > zH) then
      i1 = i1-1
    endif
    if(i1 < 1) i1 = 1
    
    i2 = i1+1
    if(i2 > nz) i2 = i2-1

    print *, i1, i2
    print *, zline(i1), zline(i2), zH

    umean_at_hub = umean(i1) + (zH-zline(i1))/dz*(umean(i2)-umean(i1))
    vmean_at_hub = vmean(i1) + (zH-zline(i1))/dz*(vmean(i2)-vmean(i1))

    angle_at_hub = atan2(vmean_at_hub, umean_at_hub)

    !! --- STEP 3 :: rotate the fields by -angle_at_hub ---
    angle_rotate = angle_with_x_deg * deg_to_radians - angle_at_hub
    sin_rotangle = sin(angle_rotate)
    cos_rotangle = cos(angle_rotate)

    do k = 1, size(u, 3)
     do j = 1, size(u, 2)
      do i = 1, size(u, 1)
          urot(i,j,k) = u(i,j,k) * cos_rotangle - v(i,j,k) * sin_rotangle 
          vrot(i,j,k) = u(i,j,k) * sin_rotangle + v(i,j,k) * cos_rotangle 
      enddo
     enddo
    enddo

    !! diagnose the rotated fields
    call get_xy_meanC_from_fC(urot, urotmean, nx, ny, nz, ftemp, ftemp_in_Y, ftemp_in_Z, gpC)
    call get_xy_meanC_from_fC(vrot, vrotmean, nx, ny, nz, ftemp, ftemp_in_Y, ftemp_in_Z, gpC)
    if(nrank==0) then
        print '(a,e19.12)', 'umean_at_hub = ', umean_at_hub
        print '(a,e19.12)', 'vmean_at_hub = ', vmean_at_hub
        print '(a,e19.12)', 'angle_at_hub = ', angle_at_hub

        open(10, file='original_rotated_fields.dat')
        do i = 1, nz
          write(10,'(5(e19.12,1x))') zline(i), umean(i), vmean(i), urotmean(i), vrotmean(i)
        enddo
        close(10)
    endif

    !! --- STEP 4 :: write out urot and vrot
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_u.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    call decomp_2d_write_one(1, urot, fname, gpC)

    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_v.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    call decomp_2d_write_one(1, vrot, fname, gpC)

    !!!!!! WRITE HEADER !!!!!!!!!!!
    if (nrank == 0) then
        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",outputFile_RID, "_info.",outputFile_TID
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        OPEN(UNIT=10, FILE=trim(fname))
        write(10,"(100g15.5)") tsim
        close(10)
    end if 

    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program 
