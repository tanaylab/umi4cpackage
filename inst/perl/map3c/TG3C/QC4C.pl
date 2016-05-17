use strict;
use FindBin;
use lib ("$FindBin::Bin/..", "$FindBin::Bin/../..");
use Cwd 'realpath';
use Getopt::Long;
use map3c::ConfOpt;
require map3c::TG3C::BaitsTab4C;
require map3c::TG3C::SamplesTab4C;

my ($workdir, $conf_file, $sample_id);
my $out_fn = undef;
my $pipeline_home = realpath("$FindBin::Bin/../../..");

die "usage perl $0 -workdir workdir -conf_file conf_file -sample_id [ID] [-output ..] \n" if @ARGV < 2;

GetOptions("workdir=s"   => \$workdir,
           "conf_file=s" => \$conf_file,
           "output=s"    => \$out_fn,
           "sample_id=i" => \$sample_id)
	or die "not all args supplied\n";


#init conf object and read conf file
my($opt) = ConfOpt::new("ConfOpt");
$opt->read_conf($conf_file);

my($bait_tab_fn)  = "$pipeline_home/" . $opt->get_opt("TG3C.baits_tab");
my($baits) = TG3C::BaitsTab4C::new("TG3C::BaitsTab4C");
$baits->read_tab($bait_tab_fn);

my($samp_tab_fn)  = "$pipeline_home/" . $opt->get_opt("TG3C.samples_tab");
my($samps) = TG3C::SamplesTab4C::new("TG3C::SamplesTab4C");
$samps->read_tab($samp_tab_fn);

my($sample) = $samps->get_sample($sample_id);
my $strip_fn = "$workdir/strip_R1.fastq";
my $bowtie_log = "$workdir/logs/bt.log";
my $umis_log   = "$workdir/logs/bait_umis.log";
my $fc_fn = "$workdir/fendchain";
my $adj_full_fn = "$workdir/adj.full.coord";

my $strip_stats = strip_stats($sample, $baits);
my $umi_counts  = parse_umi_counts($umis_log);
my ($overall_map, $unique_map) = parse_bowtie_log($bowtie_log);
my ($fc_lengths, $non_dig, $non_dig_frac, $fc_reads) = calc_nodig($fc_fn);
my (%cis_decay) = cis_decay_stats($adj_full_fn);

#Write stats output
open(OUTPUT, "> " .($out_fn || '-')) || die "cannot write OUTPUT\n";

#strip stats
print OUTPUT "-" x 20 . "extracted bait reads" . "-" x 20;
print OUTPUT "\nbait_name\tpref_reads\tpad_reads\tspecificity(% of total)\tUMIs\n";
foreach my $name (sort keys %{$strip_stats}) {
    next if $name eq "tot_reads";
    my $count     = $strip_stats->{$name}->{count};
    my $pad_count = $strip_stats->{$name}->{pad_count};
    my $cur_umis  = $umi_counts->{$name};
    printf OUTPUT "%s\t%i\t%i\t(%.02f%%)\t%i\n", $name, $count, $pad_count, eval{100 * $pad_count / $count}, $cur_umis;
}
print OUTPUT "total_reads\t$strip_stats->{'tot_reads'}(100%)\n\n";

#bowtie log
print OUTPUT "-" x 20 . "mapping stats" . "-" x 20 . "\n";
print OUTPUT "overall fragments mapped: $overall_map\n";
print OUTPUT "unique mapped fragments: $unique_map\n\n";

#fendchain stats
print  OUTPUT "-" x 20 . "library stats" . "-" x 20 . "\n";
printf OUTPUT "reads from undigested (or self-ligated) frags: %.i (%.02f%%)\n", $non_dig, $non_dig_frac*100;
print  OUTPUT "\n---fragment coverage:\n";
print  OUTPUT "frags_covered\tnum_reads(%)\n";
foreach (sort {$a <=> $b} keys %{$fc_lengths}) {
    printf OUTPUT "%i\t%i (%.02f%%)\n", $_, ${$fc_lengths}{$_}, 100*${$fc_lengths}{$_} / $fc_reads;
}

#cis decay stats
print  OUTPUT "\n---cis decay\n";
print  OUTPUT "bin\tmolecules_in_bin\n";
printf OUTPUT "0-1kb\t%i (%.02f%%)\n", $cis_decay{'1k'}, 100*$cis_decay{'1k'}/$cis_decay{'tot'};
printf OUTPUT "1-10kb\t%i (%.02f%%)\n", $cis_decay{'10k'}, 100*$cis_decay{'10k'}/$cis_decay{'tot'};
printf OUTPUT "10-100kb\t%i (%.02f%%)\n", $cis_decay{'100k'}, 100*$cis_decay{'100k'}/$cis_decay{'tot'};
printf OUTPUT "100-1000kb\t%i (%.02f%%)\n", $cis_decay{'1000k'}, 100*$cis_decay{'1000k'}/$cis_decay{'tot'};
printf OUTPUT ">1000kb\t%i (%.02f%%)\n", $cis_decay{'1000k+'}, 100*$cis_decay{'1000k+'}/$cis_decay{'tot'};
printf OUTPUT "trans\t%i (%.02f%%)\n", $cis_decay{'trans'}, 100*$cis_decay{'trans'}/$cis_decay{'tot'};
printf OUTPUT "\ntotal complexity\t%i (%.02f%%)\n", $cis_decay{'tot'}, 100*$cis_decay{'tot'}/$cis_decay{'tot'};



#Calc number of bait+pad reads
sub strip_stats {
    my ($sample, $baits) = (@_);
    my $bait_counts = {};

    my $fastq_dir = $sample->{fastqs_dir};
    my $fastq_re = $sample->{fastqs_regex};
    my @fastqs   = <$fastq_dir/*$fastq_re*R1*fastq>;
    if ($#fastqs < 0) {
        die "Can't find fastqs in $fastq_dir make sure they contain R1 in their filename\n";
    }

    my @bait_ids = @{$sample->{Bait_IDs}};
    #read bait table
    foreach my $id (@bait_ids) {
        my $cur_bait  = $baits->get_bait($id);
        my $bait_name = $cur_bait->{Bait_name};
        my $bait_seq  = $cur_bait->{Bait_seq};
        my $bait_pad  = $cur_bait->{Bait_pad};

        $bait_counts->{$bait_name}->{seq}     = $bait_seq;
        $bait_counts->{$bait_name}->{seq_pad} = $bait_seq.$bait_pad;

        $bait_counts->{$bait_name}->{count}       = 0;
        $bait_counts->{$bait_name}->{pad_count}   = 0;
    }

    $bait_counts->{"tot_reads"} = 0;

    foreach my $fn (@fastqs) {
        print STDERR "will work on $fn\n";
        my @baits = split(/,/, $baits);
        open(my $fh, "<", $fn);

        while(<$fh>) {
            my $l1 = $_;
            my $l2 = <$fh>;
            my $l3 = <$fh>;
            my $l4 = <$fh>;

            $bait_counts->{"tot_reads"}++;

            #count for each bait number of appearances
            foreach my $bait_nm (keys %{$bait_counts}) {
                if ($bait_nm eq "tot_reads") {
                    next;
                }
                my $bait_seq = $bait_counts->{$bait_nm}->{seq};
                my $seq_pad  = $bait_counts->{$bait_nm}->{seq_pad};

                if ($l2 =~ /^$bait_seq/i) {
                    $bait_counts->{$bait_nm}->{count}++;
                    if ($l2 =~ /^$seq_pad/i) {
                        $bait_counts->{$bait_nm}->{pad_count}++;
                    }
                last;
                }
            }
        }
    }

    return($bait_counts);
}

#parse per-bait umi counts
sub parse_umi_counts {
    my $fn = $_[0];
    my %umis;

    open(my $fh, "<", $fn) || warn "couldn't open $fn";
    while (<$fh>) {
        chomp;
        my @f = split(/\t/, $_);
        $umis{$f[0]} = $f[1];
    }
    return(\%umis);
}

#parse bowtie log file
sub parse_bowtie_log {
    my $fn = $_[0];
    open(my $fh, "<", $fn) || warn("couldn't find $fn\n");

    my @lines = <$fh>;

    my @overall = grep(/overall/, @lines);
    my ($overall) = ($overall[0]) =~ /(\d+\.\d+%).+/;

    my @unique = grep(/exactly\s1/, @lines);
    my ($unique) = ($unique[0]) =~ /\((\d+\.\d+%)\)/;

    return($overall, $unique);

}

#calc none digested reads
sub calc_nodig {
    my $fc_fn = $_[0];
    open(my $fh, "<", $fc_fn) || warn("couldn't find $fc_fn\n");

    #calc distrubtion of fc lengths
    my %fc_lengths;
    my $non_dig = 0;
    my $tot_reads = 0;

    while (<$fh>) {
        chomp;
        $tot_reads++;
        my @f = split(/\t/, $_);
        $fc_lengths{$f[1]}++;

        if ($f[1] < 2) {
            next;
        }

        my ($fend1) = ($f[2]) =~ /^(\d+):.+/;
        my ($fend2) = ($f[3]) =~ /^(-?\d+):.+/;
        if (abs ($fend1 - $fend2) == 1) {
            $non_dig++;
        }
    }

    return(\%fc_lengths, $non_dig, $non_dig / $tot_reads, $tot_reads);
}

#calc cis decay stats
sub cis_decay_stats{
    my $adj_fn = $_[0];
    open(my $fh, "<", $adj_fn) || warn("couldn't find $adj_fn\n");
    my %binned = ("1k"      => 0,
                  "10k"     => 0,
                  "100k"    => 0,
                  "1000k"   => 0,
                  "1000k+"  => 0,
                  "trans"   => 0,
                  "tot"     => 0);

    while (<$fh>) {
        chomp;
        if ($_ =~ /^fend1/) {
            next;
        }

        my @f = split(/\t/, $_);
        my ($chrom1, $start1, $chrom2, $start2, $mols) = (@f[1..2], @f[4..6]);
        my $offset = abs $start1 - $start2;

        $binned{'tot'} += $mols;

        if ($chrom1 ne $chrom2) {
            $binned{'trans'} += $mols;
            next;
        } elsif ($offset > 0 && $offset <= 1000) {
            $binned{'1k'} += $mols;
            next;
        } elsif ($offset > 1000 && $offset <= 10000) {
            $binned{'10k'} += $mols;
            next;
        } elsif ($offset > 10000 && $offset <= 100000) {
            $binned{'100k'} += $mols;
            next;
        } elsif ($offset > 100000 && $offset <= 1000000) {
            $binned{'1000k'} += $mols;
            next;
        } elsif ($offset > 1000000) {
            $binned{'1000k+'} += $mols;
            next;
        }
    }
    return(%binned);
}
