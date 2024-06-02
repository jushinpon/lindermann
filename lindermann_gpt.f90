program computelindemann
use omp_lib
implicit real*8(a-h,o-z)
integer natom, totalfile, counter, Denominator, nthreads, tid
real*8 xtotal, xone, xtwo, ytotal, yone, ytwo, ztotal, zone, ztwo, distance
real*8 part1, lindemann, avetotaldis_sum, avetotaldis2_sum, squareave, part2
real*8, allocatable::x(:,:), y(:,:), z(:,:), omptotalfile(:)
real*8, allocatable::totaldis2(:), dis_square(:), avetotaldis2(:), totaldis(:), avetotaldis(:)
integer, allocatable::atype(:,:)
integer startStep, endStep, incStep, startT, endT, incT 
character*128 filename

open(22, file='./linder_input.dat', status='old')  !trim(filename)
read(22, *) !total atoms in a cfg file:
read(22, *) natom
read(22, *) !total cfg files at a temperature:
read(22, *) totalfile
read(22, *) !starting temperature:
read(22, *) startT
read(22, *) !ending temperature:
read(22, *) endT
read(22, *) !temperature increment:
read(22, *) incT
read(22, *) !starting step:
read(22, *) startStep
read(22, *) !ending step:
read(22, *) endStep
read(22, *) !step increment:
read(22, *) incStep
close(22)

Denominator = (natom - 1)
part1 = 1.0 / dble(Denominator)
part1_sq = part1 * part1

open(101, file='Temp_linder.dat', status='unknown')
write(101, *) "Total atom number: ", natom

nthreads = totalfile

call omp_set_num_threads(nthreads)

allocate(x(0:totalfile-1, natom))
allocate(y(0:totalfile-1, natom))
allocate(z(0:totalfile-1, natom))
allocate(atype(0:totalfile-1, natom))

allocate(totaldis2(0:totalfile-1))
allocate(dis_square(0:totalfile-1))
allocate(avetotaldis2(0:totalfile-1))
allocate(totaldis(0:totalfile-1))
allocate(avetotaldis(0:totalfile-1))
allocate(omptotalfile(0:totalfile-1))

do 1 i = startT, endT, incT
    counter = 0 !for different thread id
    write(101, *) "T= ", i

    do 2 j = startStep, endStep, incStep
        write(filename, '(a,a,i5.5,a,i5.5,a)') './temp/', 'Lind_', i, '_', j, '.cfg'
        open(1111, file=filename, status='old')  !trim(filename)
        read(1111, *) !ITEM: TIMESTEP	
        read(1111, *) ! timestep
        read(1111, *) !ITEM: NUMBER OF ATOMS
        read(1111, *) !natom
        read(1111, *) !ITEM: BOX BOUNDS pp pp pp
        read(1111, *) xlo, xhi
        read(1111, *) ylo, yhi
        read(1111, *) zlo, zhi
        read(1111, *) !ITEM: ATOMS id type x y z

        do 1001 k = 1, natom
            read(1111, *) id, atype(counter, id), xcoor, ycoor, zcoor
            x(counter, id) = xcoor
            y(counter, id) = ycoor
            z(counter, id) = zcoor
        1001 continue
        close(1111)
        counter = counter + 1
    2 continue

    omptotalfile(:) = dble(totalfile)

    do 1000 l = 1, natom
        !$OMP PARALLEL PRIVATE(tid, n, xone, xtwo, xtotal, yone, ytwo, ytotal, zone, ztwo, ztotal, distance)
        tid = omp_get_thread_num()
        !$OMP DO
        do 1003 n = 1, nthreads
            distance = 0.0
            totaldis(tid) = 0.0
            totaldis2(tid) = 0.0

            do 1002 m = 1, natom
                if (l == m) cycle

                xone = x(tid, l)
                yone = y(tid, l)
                zone = z(tid, l)
                xtwo = x(tid, m)
                ytwo = y(tid, m)
                ztwo = z(tid, m)

                xtotal = (xone - xtwo) * (xone - xtwo)
                ytotal = (yone - ytwo) * (yone - ytwo)
                ztotal = (zone - ztwo) * (zone - ztwo)
                distance = SQRT(xtotal + ytotal + ztotal)

                totaldis(tid) = totaldis(tid) + distance
                totaldis2(tid) = totaldis2(tid) + (xtotal + ytotal + ztotal)
            1002 continue
        1003 continue
        !$OMP END DO
        !$OMP END PARALLEL

        dis_square(:) = totaldis2(:)
        avetotaldis2(:) = dis_square(:) / omptotalfile(:)
        avetotaldis(:) = totaldis(:) / omptotalfile(:)
        avetotaldis_sum = 0.0
        avetotaldis2_sum = 0.0
        do 104 n = 0, nthreads-1
            avetotaldis_sum = avetotaldis_sum + avetotaldis(n)
            avetotaldis2_sum = avetotaldis2_sum + avetotaldis2(n)
        104 continue

        squareave = avetotaldis_sum * avetotaldis_sum
        part2 = SQRT(avetotaldis2_sum - squareave) / avetotaldis_sum
        write(*,*)avetotaldis2_sum,squareave
         write(*,*)avetotaldis2_sum,avetotaldis_sum
         write(*,*)avetotaldis2_sum,avetotaldis_sum
         write(*,*)part1,part2
         pause
        lindemann = part1 * part2
        write(101, *) lindemann
    1000 continue

    write(*, *) "Atomic Lindermann for T = ", i, ' Done!'
1 continue

close(101)
end
