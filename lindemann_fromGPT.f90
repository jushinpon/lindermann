program computelindemann
implicit real*8(a-h,o-z)
integer natom, totalfile, counter, Denominator
real*8 xtotal, xone, xtwo, ytotal, yone, ytwo, ztotal, zone, ztwo, distance, lindermann
real*8, allocatable :: x(:), y(:), z(:)
real*8, allocatable :: rij2_avg(:,:), rij_avg(:,:)
integer, allocatable :: atype(:)
integer startStep, endStep, incStep, startT, endT, incT 
character*128 filename
real*8 part1

open(22, file='./lindemann_input.dat', status='old')
read(22, *) 
read(22, *) natom
read(22, *) 
read(22, *) totalfile
read(22, *) 
read(22, *) startT
read(22, *) 
read(22, *) endT
read(22, *) 
read(22, *) incT
read(22, *) 
read(22, *) startStep
read(22, *) 
read(22, *) endStep
read(22, *) 
read(22, *) incStep
close(22)

Denominator = (natom - 1)
part1 = 1.0 / dble(Denominator)

open(101, file='Temp_lindemann.dat', status='unknown')
write(101, *) "Total atom number: ", natom

allocate(x(natom), y(natom), z(natom))
allocate(atype(natom))
allocate(rij2_avg(natom, natom))
allocate(rij_avg(natom, natom))

do i = startT, endT, incT
    write(101, *) "T= ", i

    rij2_avg(:,:) = 0.0
    rij_avg(:,:) = 0.0
    totalfile = 0

    do j = startStep, endStep, incStep
        write(filename, '(a,a,i5.5,a,i5.5,a)') './temp/', 'Lind_', i, '_', j, '.cfg'
        open(1111, file=filename, status='old')
        read(1111,*)  ! ITEM: TIMESTEP
        read(1111,*)  ! timestep
        read(1111,*)  ! ITEM: NUMBER OF ATOMS
        read(1111,*)  ! natom
        read(1111,*)  ! ITEM: BOX BOUNDS
        read(1111, *) xlo, xhi
        read(1111, *) ylo, yhi
        read(1111, *) zlo, zhi
        read(1111, *) ! ITEM: ATOMS id type x y z

        do k = 1, natom
            read(1111, *) id, atype(id), x(id), y(id), z(id)
        end do
        close(1111)

        do l = 1, natom - 1
            do m = l + 1, natom
                dx = x(l) - x(m)
                dy = y(l) - y(m)
                dz = z(l) - z(m)
                dist2 = dx*dx + dy*dy + dz*dz
                dist = sqrt(dist2)

                rij2_avg(l,m) = rij2_avg(l,m) + dist2
                rij_avg(l,m)  = rij_avg(l,m) + dist

                rij2_avg(m,l) = rij2_avg(l,m)
                rij_avg(m,l)  = rij_avg(l,m)
            end do
        end do

        totalfile = totalfile + 1
    end do

    ! Normalize
    do l = 1, natom
        lindermann = 0.0
        do m = 1, natom
            if (l == m) cycle
            r_avg = rij_avg(l,m) / dble(totalfile)
            r2_avg = rij2_avg(l,m) / dble(totalfile)
            lindermann = lindermann + sqrt(r2_avg - r_avg*r_avg) / r_avg
        end do
        lindermann = lindermann * part1
        write(101, *) l, atype(l), lindermann
    end do

    write(*,*) "Atomic Lindemann for T = ", i, " Done!"
end do

close(101)
end program computelindemann
