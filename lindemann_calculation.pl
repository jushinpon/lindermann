=b
Put all folders you want to get Temperature-Lindermann profiles under this folder
=cut
use warnings;
use strict;
use Cwd;

my $currentPath = getcwd();# dir for all scripts
chdir("..");
my $mainPath = getcwd();# main path of Perl4dpgen dir
chdir("$currentPath");

#my $lindermannX = "./lindemann.x"; #lindermann executable path
my $lindermannX = "/opt/lindermann/lindemann.x"; #lindermann executable path

`rm -rf $currentPath/temp`;#remove old temp folder first
my @all_folders = `find $currentPath -mindepth 1 -maxdepth 1 -type d|grep -v .git `;
map { s/^\s+|\s+$//g; } @all_folders;
my $job_count = 0;
my $total_jobs = @all_folders;
for my $fold (@all_folders){
    print "$fold\n";
    $job_count++;
    print "*** Calculating Temperature-Lindermann profile for $fold ($job_count/$total_jobs)\n";
    my @temp = `find $fold -mindepth 1 -maxdepth 1 -type f -name "*_*_*.cfg"|sort `;
    map { s/^\s+|\s+$//g; } @temp;
    `rm -rf $currentPath/temp`;
    `mkdir $currentPath/temp`;
    my %T;
    for my $cfg (@temp){
        $cfg =~ /.+?\/(\w+)_(\d+)_(\d+)\.cfg/;
         # Check if any of the variables are empty and die if true
        die "Error: Prefix is empty" unless defined $1 && $1 ne '';
        die "Error: First Number is empty" unless defined $2 && $2 ne '';
        die "Error: Second Number is empty" unless defined $3 && $3 ne '';
        push @{$T{$2}},$3;
        my $link = "$currentPath/temp/Lind_".sprintf("%05d", $2)."_".sprintf("%05d", $3).".cfg";
        system("ln -s $cfg $link");
        #print "$cfg\n $link\n";
        #die;
    }
    my $ref_cfg = $temp[0];
    my $natom = `grep -v '^[[:space:]]*\$' $ref_cfg|grep -A 1 "ITEM: NUMBER OF ATOMS"|grep -v "ITEM: NUMBER OF ATOMS"|grep -v -- '--'`;
    $natom =~ s/^\s+|\s+$//g;# for open a dynamic array in fortran
        
    my @all_T = sort { $a <=> $b } keys %T;
    my $startT = $all_T[0];
    my $endT = $all_T[-1];
    my $incT = $all_T[1] - $all_T[0];

    my @all_step = sort { $a <=> $b } @{$T{$all_T[0]}};
    my $tot_cfg4T = @all_step; 
    my $startStep = $all_step[0];
    my $endStep = $all_step[-1];
    my $incStep = $all_step[1] - $all_step[0];

    `rm -f lindemann_input.dat`;
#print "myfile: $file\n";
my $lind_input = <<"END_MESSAGE";
total atoms in a cfg file:
$natom
total cfg files at a temperature:
$tot_cfg4T
starting temperature:
$startT
ending temperature:
$endT
temperature increment:
$incT
starting step:
$startStep 
ending step:
$endStep
step increment:
$incStep

END_MESSAGE

    open(FH, '>', "lindemann_input.dat") or die $!;
    print FH $lind_input;
    close(FH);  
    sleep(1);
    system("$lindermannX");
    sleep(1);
    system("perl lindemann_analysis.pl");
    sleep(1);
   # my $data_path = `dirname $id`;
   # $data_path =~ s/^\s+|\s+$//g;
    my $foldername = `basename $fold`;
    $foldername =~ s/^\s+|\s+$//g;
    `mv Temp_lindemann.dat $foldername-local_lindemann.dat`;
    `mv lindemann_results.csv $foldername-lindemann_results.csv`;
    `rm -rf $currentPath/temp`;
    `rm -f $currentPath/lindemann_input.dat`;
}
