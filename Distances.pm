package Distances;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
use Statistics;

sub Poisson_dist
{
    my($seq1, $seq2) = @_;
    my($i) = 0;
    my($dist) = 0;
    while($i <= length($seq1) - 1)
    {
        my($aa1) = substr($seq1, $i, 1);
        my($aa2) = substr($seq2, $i, 1);
        if(($aa1 ne $aa2) && ($aa1 ne '-') && ($aa2 ne '-'))
        {
            $dist += 1;
        }
        $i += 1;
    }
    return Poisson($dist, length($seq1));
}

sub Poisson_dist1
{
    my($seq1, $seq2) = @_;
    my($i) = 0;
    my($dist) = 0;
    while($i <= length($seq1) - 1)
    {
        my($aa1) = substr($seq1, $i, 1);
        my($aa2) = substr($seq2, $i, 1);
        if(($aa1 ne $aa2) && ($aa1 ne '-') && ($aa2 ne '-'))
        {
            $dist += 1;
        }
        $i += 1;
    }
    return Poisson($dist, length($seq1));
}

sub Poisson_dist2
{
    my($seq1, $seq2) = @_;
    my($i) = 0;
    my($dist) = 0;
    while($i <= length($seq1) - 1)
    {
        my($aa1) = substr($seq1, $i, 1);
        my($aa2) = substr($seq2, $i, 1);
        if($aa1 ne $aa2)
        {
            $dist += 1;
        }
        $i += 1;
    }
    return Poisson($dist, length($seq1));
}

sub Poisson_dist3
{
    my($seq1, $seq2) = @_;
    my($i) = 0;
    my($dist) = 0;
    while($i <= length($seq1) - 1)
    {
        my($aa1) = substr($seq1, $i, 1);
        my($aa2) = substr($seq2, $i, 1);
        if(($aa1 ne $aa2) || ($aa1 eq '-') || ($aa2 eq '-'))
        {
            $dist += 1;
        }
        $i += 1;
    }
    return Poisson($dist, length($seq1));
}

sub Poisson
{
    my($distance, $long) = @_;
    my($poisson_d);
    if((1 - ($distance/$long)) != 0)
    {
        $poisson_d = -log(1 - ($distance/$long));
    }else
    {
        $poisson_d = 0;
    }
    return $poisson_d;
}

sub BLOSUM_dist
{
    my($seq1, $seq2, $BLOSUM_matrix, $Blosum_aa_order) = @_;

    my($i) = 0;
    my(@BLOSUM_seq);
    while($i <= length($seq1) - 1){
        my($m, $n) = (0,0);
        my($aa1) = substr($seq1, $i, 1);
        my($aa2) = substr($seq2, $i, 1);
        while(($aa1 ne $$Blosum_aa_order[$m]) && ($m < scalar(@$Blosum_aa_order) - 1))
        {
            $m += 1;
        }
        while(($aa2 ne $$Blosum_aa_order[$n]) && ($n < scalar(@$Blosum_aa_order) - 1))
        {
            $n += 1;
        }
        $BLOSUM_seq[$i] = $$BLOSUM_matrix[$m][$n];
        $i += 1;
    }
    my($average_BLOSUM_dist);
    $average_BLOSUM_dist = Statistics::Average(@BLOSUM_seq);

    return $average_BLOSUM_dist;
}

1;