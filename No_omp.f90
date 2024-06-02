program computelindemann
!use omp_lib
implicit real*8(a-h,o-z)
integer natom, totalfile, counter, Denominator, nthreads, tid
real*8 xtotal, xone, xtwo, ytotal, yone, ytwo, ztotal, zone, ztwo, distance,lindermann
real*8 part1
real*8, allocatable::x(:,:), y(:,:), z(:,:), omptotalfile(:)
real*8, allocatable::rij2(:,:,:), rij(:,:,:)
integer, allocatable::atype(:)
integer startStep, endStep, incStep, startT, endT, incT 
character*128 filename

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

nthreads = totalfile

!call omp_set_num_threads(nthreads)

allocate(x(0:totalfile-1, natom))
allocate(y(0:totalfile-1, natom))
allocate(z(0:totalfile-1, natom))
allocate(atype(natom))

allocate(omptotalfile(0:totalfile-1))
!omptotalfile(:) = dble(totalfile)
allocate(rij2(0:totalfile-1,natom, natom))
allocate(rij(0:totalfile-1,natom, natom))

do 1 i = startT, endT, incT
    
    write(101, *) "T= ", i
    counter = 0
    do 2 j = startStep, endStep, incStep
        write(filename, '(a,a,i5.5,a,i5.5,a)') './temp/', 'Lind_', i, '_', j, '.cfg'
        open(1111, file=filename, status='old')
        read(1111,*) !ITEM: TIMESTEP	
	    read(1111,*) ! timestep
	    read(1111,*) !ITEM: NUMBER OF ATOMS
	    read(1111,*) !natom
	    read(1111,*) !ITEM: BOX BOUNDS pp pp pp
        read(1111, *) xlo, xhi
        read(1111, *) ylo, yhi
        read(1111, *) zlo, zhi
        read(1111, *) !ITEM: ATOMS id type x y z
        
        do 1001 k = 1, natom
            !write(*, *) k
            read(1111, *) id, atype(id), xcoor, ycoor, zcoor
            x(counter, id) = xcoor
            y(counter, id) = ycoor
            z(counter, id) = zcoor
        1001 continue
        close(1111)
        counter = counter + 1
    2 continue
    !write(*,*)"Before omp"
    !pause
    omptotalfile(:) = dble(totalfile)
    rij2 = 0.0
    rij = 0.0   
    do tid = 0, nthreads -1
    do l = 1, natom - 1
        do m = l + 1, natom
          !  if (l == m) cycle
            !write(*,*)"tid = ",tid
            xone = x(tid, l)
            yone = y(tid, l)
            zone = z(tid, l)
            xtwo = x(tid, m)
            ytwo = y(tid, m)
            ztwo = z(tid, m)

            xtotal = (xone - xtwo) * (xone - xtwo)
            ytotal = (yone - ytwo) * (yone - ytwo)
            ztotal = (zone - ztwo) * (zone - ztwo)
            temp = xtotal + ytotal + ztotal
            distance = SQRT(temp)

            rij2(tid,l, m) = temp
            rij(tid,l, m) = distance

            rij2(tid,m,l) = temp
            rij(tid,m,l) = distance
        end do
    end do
    end do
   
    !r2_sum = r2_sum/dble(nthreads)
    !r_sum = r_sum/dble(nthreads)

    !avetotaldis2_sum = 0.0
    do l = 1, natom
        lindermann = 0.0
        !sum_r = 0.0
        !sum_r2 = 0.0
        !get local lindermann below
       

        do m = 1, natom
            if (l == m) cycle
            averij = 0.0
            averij2 = 0.0
            do k = 0,nthreads -1 
                averij = averij + rij(k,l,m)
                !averij2 = averij2 + rij(k,l,m) * rij(k,l,m)
                averij2 = averij2 + rij2(k,l,m)
            enddo

            averij = averij/dble(nthreads)
            averij2 = averij2/dble(nthreads)
             !/ nthreads
           !  if(m == 1034)then
                 !write(*,*)m,": ",averij2,averij * averij,averij
                 !write(*,*)averij2 - averij * averij
            !endif
            lindermann = lindermann  + SQRT(averij2 - averij * averij)/averij
            !sum_r = sum_r + r_total
            !sum_r2 = sum_r2 + r2_total
        end do
        !write(*, *)l,atype(l),lindermann
        !write(*, *)l,atype(l),lindermann*part1
        lindermann = lindermann*part1
        write(101, *)l,atype(l),lindermann
        !sum_r = sum_r / Denominator
        !sum_r2 = sum_r2 / Denominator
        !delta_i = sqrt(sum_r2 - sum_r ** 2) / sum_r
        !avetotaldis_sum = avetotaldis_sum + sum_r
        !avetotaldis2_sum = avetotaldis2_sum + sum_r2
    end do

    !squareave = avetotaldis_sum * avetotaldis_sum
    !part2 = SQRT(avetotaldis2_sum - squareave) / avetotaldis_sum
!
    !lindemann = part1 * part2
    write(*, *) "Atomic Lindermann for T = ", i, ' Done!'
1 continue

close(101)
end
