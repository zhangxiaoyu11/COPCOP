package File_IO;
require Exporter;
@ISA = qw(Exporter);

use strict;
use warnings;
use Seq_manag;

sub Control_reader
{
    my($control, $information) = @_;
    my($Ctl) = Input_file($control);
    my($i) = 0;
    while(<$Ctl>)
    {
        $_ =~ s/ //g;
        my(@line) = split(/\:/, $_);
        my(@line2) = split(/\*/, $line[1]);
        $$information[$i] = $line2[0];
        $$information[$i] =~ s/\s+//g;
        $i += 1;
    }
}

sub Input_file
{
    my($infile) = @_;
    my($input_file);
    unless(open($input_file, $infile))
    {
        print STDERR "Cannot find $infile!\n";
        exit;
    }
    return $input_file;
}

sub Output_file
{
    my($outfile) = @_;
    my($output_file);
    unless(open($output_file, ">$outfile"))
    {
        print STDERR "Cannot open $outfile!\n";
        exit;
    }
    return $output_file;
}

sub Sequences_reader
{
    my($in, $Seqs, $names) = @_;
    my($infile) = Input_file($in);
    my($j, $longitud, $k, $controlador) = (0,0,0,0);
    my(@file_information);
    while(<$infile>)
    {
        chomp $_;
        $file_information[$j] = $_;
        $j += 1;
    }
    $j = 0;

    ##Reading files in Phylip format
    if($file_information[0] =~ /^\d/)				##This will loop would read phylip format
    {
        while($j <= scalar(@file_information) - 1)
        {
            if(length($file_information[$j]) == 0)
            {
                $j += 1;
                next;
            }
            my(@array) = split(/\s+/, $file_information[$j]);
            if($j == 0)
            {
                $longitud = $array[1];
            }else
            {
                if($controlador == 0)
                {

                    $$names[$k] = $array[0];
                    $controlador += 1;
                    if(scalar(@array) > 1)
                    {
                        $$Seqs[$k] = $array[1];
                    }else
                    {
                        $$Seqs[$k] = '';
                    }
                }else
                {							#line 100
                    if(length($$Seqs[$k]) == $longitud)
                    {
                        $k += 1;
                        $$names[$k] = $array[0];
                        if(scalar(@array) > 1)
                        {
                            $$Seqs[$k] = $array[1];
                        }else
                        {
                            $$Seqs[$k] = '';
                        }
                    }else
                    {
                        $$Seqs[$k] .= $array[0];
                    }
                }
            }
            $j += 1;
        }
    }else ##Reading files in Mega or fasta format
    {
        while($j <= scalar(@file_information) - 1)
        {
            if(length($file_information[$j]) == 0)
            {
                $j += 1;
                next;
            }
            my(@array) = split(/\s+/, $file_information[$j]);
            if(($array[0] =~ /^(#mega)/i) || ($array[0] =~ /^(title:)/i))
            {
                $j += 1;
                next;
            }
            if($controlador == 0)
            {
                $$names[$k] = substr($array[0], 1, (length($array[0]) - 1));
                $controlador += 1;
                if(scalar(@array) > 1)
                {
                    $$Seqs[$k] = $array[1];
                }else
                {
                    $$Seqs[$k] = '';
                }
            }else
            {
                if(($array[0] =~ /^>/) || ($array[0] =~ /^#/))
                {
                    $k += 1;
                    $$names[$k] = substr($array[0], 1, (length($array[0]) - 1));
                    if(scalar(@array) > 1)
                    {
                        $$Seqs[$k] = $array[1];
                    }else
                    {
                        $$Seqs[$k] = '';
                    }
                }else
                {
                    $$Seqs[$k] .= $array[0];
                }
            }
            $j += 1;
        }
    }
    my($N_seqs) = 0;
    while($N_seqs <= scalar(@$Seqs) - 1)
    {
        chomp $$Seqs[$N_seqs];
        $$Seqs[$N_seqs] =~ s/\r//g;
        $N_seqs += 1;
    }
}

sub Codon2aa
{
    my($sequences) = @_;
    my($i) = 0;
    while($i <= scalar(@$sequences) - 1)
    {
        my($seq) = $$sequences[$i];
        $$sequences[$i] = Seq_manag::Nuc_amino($seq);
        $i += 1;
    }
}

sub Check_seq_len
{
    my(@sequences) = @_;
    my($i) = 0;
    while($i <= scalar(@sequences) - 1)
    {
        if(length($sequences[$i]) != length($sequences[0]))
        {
            print "Error: Sequence " . ($i + 1) . " has different length!\n";
            exit;
        }
        $i += 1;
    }
    return 0;
}

sub Print2file_array
{
    my($array,$out_filename) = @_;
    my($out_stream);
    $out_stream = Output_file($out_filename);
    for(my $i = 0; $i < scalar(@$array); $i++){
        print $out_stream $$array[$i]."\n";
    }
}

sub Print2file_matrix
{
    my($matrix,$out_filename) = @_;
    my($out_stream);
    $out_stream = Output_file($out_filename);
    for(my $i = 0; $i < scalar(@$matrix); $i++){
        for(my $j = 0; $j < scalar(@{$$matrix[0]}); $j++){
            print $out_stream $$matrix[$i][$j]."\t";
        }
        print $out_stream "\n";
    }
}

sub Read_blosum
{
    my($file, $B_array, $B_aa) = @_;

    my($j, $i) = (0, 1);
    my(@array, $in);
    $in = Input_file($$file);
    while(<$in>)
    {
        if($_ =~ /^(BLOSUM62)/)
        {
            next;
        }
        $_ =~ s/ {2}/ /g;
        @array = split(/\s/, $_);
        $$B_aa[$j] = $array[0];
        $$B_aa[$j] =~ s/\s+//g;
        for($i = 1; $i <= 21; $i++)
        {
            $$B_array[$j][$i -1] = $array[$i];
        }
        $j += 1;
    }
}

sub Print_Unit_Pairs
{
    my($out_file, $Selected_unit_pairs, $Selected_pairs_coefficient) = @_;

    for(my $i = 0; $i < scalar(@$Selected_pairs_coefficient); $i++){
        print $out_file $$Selected_unit_pairs[$i][0]."\t\t\t\t".$$Selected_unit_pairs[$i][1]."\t\t\t\t\t".$$Selected_pairs_coefficient[$i]."\n";
    }
}

1;