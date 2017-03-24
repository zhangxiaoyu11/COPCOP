package Core_algorithm;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
use Distances;
use Statistics;
use File_IO;

sub Get_gene_dists{
    my($sequences, $dist_func, $gene_distances) = @_;

    my($BLOSUM_file_name) = 'Blosum';
    my(@BLOSUM_matrix, @Blosum_aa_order);
    File_IO::Read_blosum(\$BLOSUM_file_name, \@BLOSUM_matrix, \@Blosum_aa_order);

    my($seq1_m, $seq2_n, $gene_dists_i) = (0,1,0);
    while($seq1_m <= scalar(@$sequences) - 2){
        while($seq2_n <= scalar(@$sequences) - 1){

            if($dist_func == 0){
                $$gene_distances[$gene_dists_i] = Distances::Poisson_dist($$sequences[$seq1_m],$$sequences[$seq2_n]);
            }

            if($dist_func == 1){
                $$gene_distances[$gene_dists_i] = Distances::BLOSUM_dist($$sequences[$seq1_m],$$sequences[$seq2_n],\@BLOSUM_matrix,\@Blosum_aa_order);
            }

            $gene_dists_i += 1;
            $seq2_n += 1;
        }
        $seq1_m += 1;
        $seq2_n = $seq1_m + 1;
    }
}

sub Get_D{
    my($sequences, $gene_distances, $unit_width, $dist_func, $D) = @_;

    my($BLOSUM_file_name) = 'Blosum';
    my(@BLOSUM_matrix, @Blosum_aa_order);
    File_IO::Read_blosum(\$BLOSUM_file_name, \@BLOSUM_matrix, \@Blosum_aa_order);

    my($site_p, $seq1_m, $seq2_n, $seq_pair_i) = (0,0,1,0);
    my($unit_distance) = 0;
    my($unit1, $unit2);
    ## need to consider the unit width
    while($site_p <= length($$sequences[0]) - $unit_width){
        while($seq1_m <= scalar(@$sequences) - 2){
            while($seq2_n <= scalar(@$sequences) - 1){
                $unit1 = substr($$sequences[$seq1_m],$site_p,$unit_width);
                $unit2 = substr($$sequences[$seq2_n],$site_p,$unit_width);

                if($dist_func == 0){
                    $unit_distance = Distances::Poisson_dist($unit1,$unit2);
                }

                if($dist_func == 1){
                    $unit_distance = Distances::BLOSUM_dist($unit1,$unit2,\@BLOSUM_matrix,\@Blosum_aa_order);
                }

                $$D[$site_p][$seq_pair_i] = ($unit_distance-$$gene_distances[$seq_pair_i])**2;
                $seq_pair_i += 1;
                $seq2_n += 1;
            }
            $seq1_m += 1;
            $seq2_n = $seq1_m + 1;
        }
        ($seq1_m, $seq2_n, $seq_pair_i) = (0,1,0);
        $site_p += 1;
    }
}

sub Get_Correlation{
    my($D1, $D2, $seq_pair_number, $Correlations, $Correlations_matrix) = @_;

    my($site1_i, $site2_j, $seq_pair1_m, $seq_pair2_n, $cor_number_p) = (0, 0, 0, 0, 0);

    while($site1_i <= scalar(@$D1) - 1)
    {
        my(@D_column1);
        $seq_pair1_m = 0;
        while($seq_pair1_m < $seq_pair_number){
            $D_column1[$seq_pair1_m] = $$D1[$site1_i][$seq_pair1_m];
            $seq_pair1_m += 1;
        }
        while($site2_j <= scalar(@$D2) - 1)
        {
            my(@D_column2);
            $seq_pair2_n = 0;
            while($seq_pair2_n < $seq_pair_number){
                $D_column2[$seq_pair2_n] = $$D2[$site2_j][$seq_pair2_n];
                $seq_pair2_n += 1;
            }
            $$Correlations[$cor_number_p] = Statistics::Correlation(\@D_column1, \@D_column2);
            $$Correlations_matrix[$site1_i][$site2_j] = $$Correlations[$cor_number_p];
            $cor_number_p += 1;
            $site2_j += 1;
        }
        $site1_i += 1;
        $site2_j = 0;
    }
}

sub Get_threshold{
    my($alpha, $sample_number, $threshold, $correlations) = @_;

    my(@numbers);
    my($i, $Mean, $SE) = (0,0,0);
    while($i < $sample_number)
    {
        $numbers[$i] = $$correlations[int(rand(scalar(@$correlations)))];
        $i += 1;
    }
    $Mean = Statistics::Average(@numbers);
    $SE = Statistics::SE(\@numbers, $Mean);
    if($alpha == 0.05)
    {
        $$threshold = ((1.69 * $SE) + $Mean);
    }
    if(($alpha < 0.05) && ($alpha > 0.025))
    {
        $$threshold = ((1.95 * $SE) + $Mean);
    }
    if(($alpha <= 0.025) && ($alpha > 0.01))
    {
        $$threshold = ((1.97 * $SE) + $Mean);
    }
    if(($alpha <= 0.01) && ($alpha > 0.001))
    {
        $$threshold = ((2.33 * $SE) + $Mean);
    }
    if($alpha <= 0.001)
    {
        $$threshold = ((3.1 * $SE) + $Mean);
    }
}

sub Find_unit_pairs{
    my($threshold, $Correlations_matrix, $Selected_unit_pairs, $Selected_pairs_coefficient) = @_;

    my($pair_number_m) = 0;

    for(my $i = 0; $i < scalar(@$Correlations_matrix); $i++){
        for(my $j = 0; $j < scalar(@{$$Correlations_matrix[0]}); $j++){
            if($$Correlations_matrix[$i][$j] >= $threshold || $$Correlations_matrix[$i][$j] <= -$threshold){
                $$Selected_unit_pairs[$pair_number_m][0] = $i;
                $$Selected_unit_pairs[$pair_number_m][1] = $j;
                $$Selected_pairs_coefficient[$pair_number_m] = $$Correlations_matrix[$i][$j];
                $pair_number_m++;
            }
        }
    }
}

sub Find_real_positions{
    my($sequences1,$sequences2,$reference_No1,$reference_No2,$real_positions1,$real_positions2) = @_;

    my(@gaps1,@gaps2)

}

1;