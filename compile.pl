#foreach my $i (@files){
#unlink "$i";
#print "Kill $i file\n";
#}
#`rm -rf output`;
unlink "./lindermann.x";
#system("gfortran -fopenmp -Wall -Wno-tabs -o lindermann.x ./lindemann_fortran_new.f90");
#20240602:openmp version has some bugs, so only single thread version works good now.
system("gfortran -o lindemann.x ./No_omp.f90");
`chmod 755 ./lindemann.x`;
#system("gfortran -o maxent.exe 00MAXENT_main.f90");
#$temp = system("gfortran -O3 -o maxent.exe 00MAXENT_main.f90");
#$temp = system("gfortran -O3 -o maxent.x 00MAXENT_main.f95");
#$temp = system("gfortran -fopenmp -O3 -o maxent.x 00MAXENT_main.f95");
#$temp = system("gfortran -g -fcheck=all -Wall 00MAXENT_main.f90");

#sleep(1);
#$temp = './maxent.exe';#.' > 00printout.txt';
#print $temp;
#system($temp);