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
my $sour = "Temp_lindemann.dat";
my $natom = `grep "Total atom number:" ./$sour|awk '{print \$4}'`;
$natom =~ s/^\s+|\s+$//g;

my @all_Temperatures = `grep "T=" ./$sour|awk '{print \$2}'`;
map { s/^\s+|\s+$//g; } @all_Temperatures;

my @all_lindemann = `grep -A $natom "T=" ./$sour|grep -v "T="|grep -v -- '--'`;
map { s/^\s+|\s+$//g; } @all_lindemann;

my $total_lin = $natom * @all_Temperatures;
my $by_grep = @all_lindemann;

die "totoal lindemann number is not correct! (grep: $by_grep, natom* all T number: $total_lin)\n" unless($total_lin == $by_grep);

my %T2Lindemann;
for my $t (0 .. $#all_Temperatures){
    my $temperature = $all_Temperatures[$t];
    $T2Lindemann{$temperature} = [@all_lindemann[$t * $natom .. ($t + 1) * $natom - 1]];
}

#get all atomtypes
my %types;
my $lowest_T = $all_Temperatures[0];
my @temp = @{$T2Lindemann{$lowest_T}};#use the lowest temperature to find all atom types
for (@temp){
    my @temp = split(/\s+/,$_);
    if(exists $types{$temp[1]}){
        $types{$temp[1]}++;
    }
    else{
        $types{$temp[1]} = 1;
    }
    #print "$_\n";
};

my $types = join(",",sort {$a <=> $b} keys %types);
my $columns = "Temperature,All,"."$types";

open(FH, '>', "lindemann_results.csv") or die $!;

print FH "$columns\n";

for my $t (@all_Temperatures){

    my $lind4all = 0.0;
    my %lind4types;
    @lind4types{keys %types} = 0.0;
    
    for my $id (@{$T2Lindemann{$t}}){
        #1 2   4.2110094767136407E-003
        my @temp = split(/\s+/,$id);
        $lind4all += $temp[2];
        $lind4types{$temp[1]} += $temp[2];
    }
    $lind4all = $lind4all/$natom;
    my @aveLin;
    for my $ty (sort {$a <=> $b} keys %types){
        $lind4types{$ty} = $lind4types{$ty}/$types{$ty}; 
        push @aveLin, $lind4types{$ty};  
    }

    my $aveLin = join(",",@aveLin);
    print FH "$t,$lind4all,$aveLin\n";
}
close(FH);


#
#for (sort {$a <=> $b} keys %types){
#    print "$_,$types{$_}\n";
#}
#    open(FH, '>', "linder_input.dat") or die $!;
#    print FH $lind_input;
#    close(FH);
#   die;
#    sleep(1);
#    system("$lindermannX/lindermann.x");
#   # my $data_path = `dirname $id`;
#   # $data_path =~ s/^\s+|\s+$//g;
#    my $foldername = `basename $fold`;
#    $foldername =~ s/^\s+|\s+$//g;
#    `mv Temp_linder.csv Temp_linder_$foldername.csv`
#}
#`rm -rf $currentPath/temp`;
#`rm -f $currentPath/linder_input.dat`;
