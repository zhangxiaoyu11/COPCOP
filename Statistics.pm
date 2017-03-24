package Statistics;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;

sub Average
{
    my(@array) = @_;
    my($sum) = 0;
    for(my $i = 0; $i <= scalar(@array) - 1; $i++)
    {
        $sum += $array[$i];
    }
    return ($sum/scalar(@array));
}

sub Correlation
{
    my($sampleA, $sampleB) = @_;
    my(@A, @B, $meanA, $meanB) = (0,0,0,0);
    for(my $k = 0; $k <= scalar(@$sampleA) - 1; $k++)
    {
        $A[$k] = $$sampleA[$k];
    }
    for(my $k = 0; $k <= scalar(@$sampleB) - 1; $k++)
    {
        $B[$k] = $$sampleB[$k];
    }
    $meanA = Average(@A);
    $meanB = Average(@B);
    my($i) = 0;
    my($nominator) = 0;
    my($X) = 0;
    my($Y) = 0;
    my($correl);
    while($i <= scalar(@$sampleA) - 1)
    {
        if($$sampleB[$i])
        {
            $nominator += ($$sampleA[$i] - $meanA) * ($$sampleB[$i] - $meanB);
            $Y += (($$sampleB[$i] - $meanB)**2);
        }else
        {
            $nominator += 0;
            $Y += 0;
        }
        $X += (($$sampleA[$i] - $meanA)**2);

        $i += 1;
    }
    if(($X > 0) && ($Y > 0))
    {
        $correl = ($nominator/sqrt($X*$Y));
    }else
    {
        $correl = 0;
    }
    return $correl;
}

sub SE
{
    my($numbers, $media) = @_;
    my($Sum_squares) = 0;
    my($i);
    for($i = 0; $i <= scalar(@$numbers) - 1; $i++)
    {
        $Sum_squares += (($$numbers[$i] - $media)**2);
    }
    my($var) = sqrt($Sum_squares/($i + 1));
    return $var;
}

1;