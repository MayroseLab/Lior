use strict;
use Bio::EnsEMBL::Registry;
use List::Util qw(min);

my ($species, $group, $out_bed) = @ARGV;
open(FH, '>', $out_bed) or die $!;

my $registry = 'Bio::EnsEMBL::Registry';

if ($group eq "vertebrates"){
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous');
}
else{
$registry->load_registry_from_db(
    -host => 'mysql-eg-publicsql.ebi.ac.uk', # alternatively 'useastdb.ensembl.org'
    -port => 4157,
    -user => 'anonymous');
}

my $W = 10000000;

my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice',  );
my @slices = @{ $slice_adaptor->fetch_all('toplevel') };
foreach my $slice (@slices) {
    my $chrom = $slice->seq_region_name();
    my $chrom_length = $slice->length();
    my $s = 1;
    my $e = min($W, $chrom_length);
    while ($e <= $chrom_length){
        my $subslice = $slice->sub_Slice($s,$e);
        my @repeats = @{ $subslice->get_all_RepeatFeatures() };
        foreach my $repeat (@repeats) {
            if ($repeat->start() < 1 or $repeat->end() < 1){
                next;
            }
            my $start = $repeat->start() + $s - 2;
            my $end = $repeat->end() + $s - 1;
            my $rc = $repeat->repeat_consensus();
            my $rep_class = $rc->repeat_class;
            my $rep_type = $rc->repeat_type;
            #if ($rep_class =~ /Helitron/){
            #    $rep_type = "Type II Transposons";
            #}
            #elsif ($repeat->display_id() =~ /LTR/){
            #    $rep_type = "LTRs";
            #}
            if ($rep_type eq "Dust"){
                $rep_type = "Low complexity regions";
            }
            elsif ($rep_type eq "repeatdetector"){
                $rep_type = "Unknown";
            }            
            $rep_type =~ s/ /_/g;
            $rep_type =~ s/\//_/g;
            printf FH ( "%s\t%d\t%d\t%s\t%s\t%s\n",
                    $chrom, $start, $end, $repeat->display_id(), $rep_class, $rep_type);
        }
        last if ($e == $chrom_length);
        $s = $e + 1;
        $e = min($s + $W, $chrom_length);
    }
}
close(FH);
