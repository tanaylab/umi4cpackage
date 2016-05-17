# Split adj file by baits when multiplexed

use strict;
use FindBin;
use Getopt::Long;
use Cwd 'realpath';
use lib ("$FindBin::Bin/..", "$FindBin::Bin/../..");
use File::Path("make_path");
require map3c::ConfOpt;
require map3c::ReDB;
require map3c::TG3C::BaitsTab4C;
require map3c::TG3C::SamplesTab4C;

if (scalar @ARGV < 4) {
    die "usage $0 --adj <adj> --samples_tab <tab> --baits_tab <tab> --sample_id <samples_id> --redb_dir <> --re_seq <> --log_file <>\n";
}

my ($conf_file, $adj_fn, $sample_id, $redb_dir, $re_seq, $bait_tab_fn, $samp_tab_fn);
my $log_file = undef;
GetOptions(
           "baits_tab=s" => \$bait_tab_fn,
           "samples_tab=s" => \$samp_tab_fn,
           "adj=s"       => \$adj_fn,
           "sample_id=i" => \$sample_id,
           "redb_dir=s"  => \$redb_dir,
           "log_file=s"  => \$log_file,
           "re_seq=s"    => \$re_seq
           )
	or die "not all args supplied\n";

print STDERR "adj_fn = $adj_fn\n";
my $redb = ReDB::new("ReDB");
$redb->init_from_tabs($redb_dir."/$re_seq.frags");

my($baits) = TG3C::BaitsTab4C::new("TG3C::BaitsTab4C");
$baits->read_tab($bait_tab_fn);

my($samps) = TG3C::SamplesTab4C::new("TG3C::SamplesTab4C");
$samps->read_tab($samp_tab_fn);

my($sample) = $samps->get_sample($sample_id);

# Retrieve bait fends
my ($bait_fends, $bait_names) = get_bait_fends($sample, $baits, $redb);
my %bait_fends = %{$bait_fends};
my @bait_names = @{$bait_names};

#create fhs
my %baits_out_fhs;
my %baits_umi_count; # initialize umi counter
foreach my $bait_name (values %bait_fends) {
    print STDERR "\ncreating fh for $bait_name\n";
    if ( -f "$adj_fn.$bait_name") {
        system ("rm -f $adj_fn.$bait_name");
    }
    open (my $fh, ">>", "$adj_fn.$bait_name") || die "Cannot creat $adj_fn.$bait_name";
    print $fh "fend1\tfend2\tcount\n";
    $baits_out_fhs{$bait_name} = $fh;
    $baits_umi_count{$bait_name} = 0;
}
if ( -f "$adj_fn.err") {
    system("rm -f $adj_fn.err");
}
open (my $err_fh, ">>", "$adj_fn.err") || die "Cannot create $adj_fn.err\n";


open (my $adj_fh, "<", $adj_fn) || die "Cannot open adj $adj_fn\n";
while (<$adj_fh>) {
    chomp;
    next if ($_ =~ /fend/);

    my $line = $_;

    my @f = split(/\t/, $_);
    my $fend1 = $f[0];
    print_to_bait_fh(\%bait_fends, \%baits_out_fhs, $fend1, $line);
}

# Print out umi counts
open(LOG, "> " .($log_file || '-')) || die "cannot write $log_file\n";
print LOG "Bait_name\tUMIs\n";
foreach my $name (@bait_names) {
    print LOG "$name\t$baits_umi_count{$name}\n";
}

sub print_to_bait_fh {
    my ($bait_fends, $baits_out_fhs, $fend1, $line) = (@_);
    my $cur_bait_name;
    my $matched = 0;

    #step one fend1 before and after to map
    for (my $cur_f = $fend1-3; $cur_f <= $fend1+3; $cur_f++) {

        if (exists $bait_fends->{$cur_f}) {
            $cur_bait_name = $bait_fends->{$cur_f};
            my $cur_fh = $baits_out_fhs->{$cur_bait_name};

            my @f = split("\t", $line);
            print $cur_fh "$line\n";
            $matched++;
            if ($f[2] > 0) {
                $baits_umi_count{$cur_bait_name} += $f[2];
            }
            last;
         }
    }
    if ($matched == 0) {
        print $err_fh "$line\n";
    }
}

sub get_bait_fends {
     my ($sample, $baits, $redb) = (@_);
     my @bait_ids = @{$sample->{Bait_IDs}};
     my %bait_fends;
     my @bait_names;

     foreach my $id (@bait_ids) {

         my $cur_bait    = $baits->get_bait($id);

         my $bait_name   = $cur_bait->{Bait_name};
         my $bait_seq    = $cur_bait->{Bait_seq};
         my $bait_pad    = $cur_bait->{Bait_pad};
         my $bait_chr    = $cur_bait->{Bait_chr};
         my $bait_coord  = $cur_bait->{Bait_coord};
         # Estimate fend (ignore strandness)
         my ($bait_fid, $bait_fid_off, $bait_fid_side) = $redb->map_frag($bait_chr, $bait_coord, 1500);
         my $bait_fend;
         if ($bait_fid_side == 3) {
            $bait_fend = $bait_fid*2+1;
         } else {
            $bait_fend = $bait_fid * 2;
         }
         $bait_fends{$bait_fend} = $bait_name;
         push @bait_names, $bait_name;
         print STDERR "\nWorking on bait name $bait_name fend: $bait_fend chr: $bait_chr coord $bait_coord\n";
     }
     return (\%bait_fends, \@bait_names);
}
