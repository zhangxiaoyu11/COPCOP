#!/usr/bin/perl -w
use strict;
use warnings;
use File_IO;        # File_IO.pm file
use Core_algorithm;     # Core_algorithm.pm file
use Extra_information;      # Extra_information.pm file

srand(time|$$);

##----------INITIALIZATION----------##

my($header) =
    "*************************************************************

  COPCOP: detecting signals of CO-Positive and CO-Purifying

  Xiaoyu Zhang
  School of Computer Science
  National University of Defense Technology
  Changsha, China
  zhangxiaoyu11\@nudt.edu.cn

*************************************************************\n";
print "\n".$header;

print "\n[INITIALIZATION]\n";

my($control, @parameter);
$control = 'COPCOP.ctl';
File_IO::Control_reader($control, \@parameter);

print "\nReading control file options...\n\n";

my($print_parameter) = "Input file1: ".$parameter[0]."\n".
    "Input file2: ".$parameter[1]."\n".
    "Output file: ".$parameter[2]."\n".
    "Type of data 1: ".$parameter[3]."\n".
    "Type of data 2: ".$parameter[4]."\n".
    "Reference sequence 1: ".$parameter[5]."\n".
    "Reference sequence 2: ".$parameter[6]."\n".
    "Unit width: ".$parameter[7]."\n".
    "Significance test: ".$parameter[8]."\n".
    "Threshold R-value:  ".$parameter[9]."\n".
    "Threshold alpha-value: ".$parameter[10]."\n".
    "Random sampling number: ".$parameter[11]."\n".
    "Minimum R: ".$parameter[12]."\n".
    "Distance function: ".$parameter[13]."\n";
print $print_parameter;

my($out_file);
$out_file = File_IO::Output_file($parameter[2]);

print $out_file $header;
print $out_file "\n[INPUT PARAMETERS]\n";
print $out_file "-------------------------------------------------------------\n";
print $out_file $print_parameter;
print $out_file "-------------------------------------------------------------\n";

print "\nReading sequence files...\n\n";

my(@sequences1,@seq_names1,@sequences2,@seq_names2);
File_IO::Sequences_reader($parameter[0],\@sequences1,\@seq_names1);
File_IO::Sequences_reader($parameter[1],\@sequences2,\@seq_names2);
if($parameter[3] == 1){
    File_IO::Codon2aa(\@sequences1);
}
if($parameter[4] == 1){
    File_IO::Codon2aa(\@sequences2);
}
my($ref_seq_no1, $ref_seq_no2) = (0,0);
$ref_seq_no1 = $parameter[5]-1;
$ref_seq_no2 = $parameter[6]-1;
my($seq_number) = scalar(@seq_names1);

if($ref_seq_no1 >= $seq_number)
{
    print "Error: Please re-check Reference sequence 1 in the control file!\n";
    exit;
}
if($ref_seq_no2 >= $seq_number)
{
    print "Error: Please re-check Reference sequence 2 in the control file!\n";
    exit;
}
File_IO::Check_seq_len(@sequences1);
File_IO::Check_seq_len(@sequences2);


##----------CORE_ALGORITHM----------##
print "[CO-EVOLUTION ANALYSIS]\n\n";

print $out_file "\n[CO-EVOLUTION ANALYSIS]\n";

print "Computing gene distances...\n\n";
my(@gene_distances1,@gene_distances2);
Core_algorithm::Get_gene_dists(\@sequences1,$parameter[13],\@gene_distances1);
Core_algorithm::Get_gene_dists(\@sequences2,$parameter[13],\@gene_distances2);

print "Computing variances...\n\n";
my($unit_width) = $parameter[7];
my(@D1,@D2);
Core_algorithm::Get_D(\@sequences1,\@gene_distances1,$unit_width,$parameter[13],\@D1);
Core_algorithm::Get_D(\@sequences2,\@gene_distances2,$unit_width,$parameter[13],\@D2);

print "Computing Pearson correlation coefficient...\n\n";
my($seq_pair_number) = $seq_number*($seq_number-1)/2;
my(@Correlations,@Correlations_matrix);
Core_algorithm::Get_Correlation(\@D1,\@D2,$seq_pair_number,\@Correlations,\@Correlations_matrix);

print "Computing the threshold of correlation coefficient...\n\n";
my($threshold) = 0;
if($parameter[8] == 0){
    $threshold = $parameter[9];
}elsif($parameter[8] == 1){
    Core_algorithm::Get_threshold($parameter[10],$parameter[11],\$threshold,\@Correlations);
}

print "Looking for target unit pairs...\n\n";
my(@Selected_unit_pairs,@Selected_pairs_coefficient);
Core_algorithm::Find_unit_pairs($threshold,\@Correlations_matrix,\@Selected_unit_pairs,\@Selected_pairs_coefficient);


##----------OUTPUT_CO-EVOLUTION_INFORMATION----------##
print $out_file "Threshold: ".$threshold."\n";

print $out_file "-------------------------------------------------------------\n";
print $out_file "Position1\t\tPosition2\t\t\tCorrelation\n";
print $out_file "-------------------------------------------------------------\n";
File_IO::Print_Unit_Pairs($out_file,\@Selected_unit_pairs,\@Selected_pairs_coefficient);


##----------OUTPUT_EXTRA_INFORMATION----------##
##File_IO::Print2file_array(\@Correlations,"Cor.out");
##File_IO::Print2file_matrix(\@Correlations_matrix,"Cor_matrix.out");

##my(@con_degree1, @con_degree2);
##Extra_information::Conserved_degree(\@sequences1,\@con_degree1);
##Extra_information::Conserved_degree(\@sequences2,\@con_degree2);
##File_IO::Print2file_array(\@con_degree1,"Con_degree1.out");
##File_IO::Print2file_array(\@con_degree2,"Con_degree2.out");

##File_IO::Print2file_matrix(\@Selected_unit_pairs,"Unit_pairs.out");

print "[FINISHED]\n\n";






