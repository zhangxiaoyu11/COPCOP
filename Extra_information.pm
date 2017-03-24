package Extra_information;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;

sub Conserved_degree
{
    my($sequences, $con_degree) = @_;
    my($i, $j) = (0,0);
    while($i <= length($$sequences[0]) - 1)
    {
        $j = 0;
        my(@aa_array);
        $$con_degree[$i] = 0;
        while($j <= scalar(@$sequences) - 1)
        {
            $aa_array[$j] = substr($$sequences[$j], $i, 1);
            if($j > 0)
            {
                my($find, $search) = (0,0);
                while($search <= scalar(@aa_array) - 2)
                {
                    if($aa_array[$j] eq $aa_array[$search])
                    {
                        $find += 1;
                        last;
                    }
                    $search += 1;
                }
                if($find == 0)
                {
                    $$con_degree[$i] += 1;
                }
            }
            $j += 1;
        }
        $i += 1;
    }
}

1;