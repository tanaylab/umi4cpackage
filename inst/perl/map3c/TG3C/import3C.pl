use strict;
use FindBin;
use Cwd 'realpath';
use lib ( realpath("$FindBin::Bin/.."), realpath("$FindBin::Bin/../..") );
use File::Path("make_path");
require map3c::ConfOpt;
require map3c::TG3C::BaitsTab4C;
require map3c::TG3C::SamplesTab4C;

# Validate configuration option names

my ($opt) = ConfOpt::new("ConfOpt");

$opt->parse_argv_conf( \@ARGV );

my ($done) = "";
my ($ret)  = 0;

############################
# Construct workdir
############################
my $dir_lscripts = $opt->get_opt("TG3C.perl_scripts_path") . "/map3c/";
my ($wd) = $opt->get_opt("TG3C.workdir");

my ($bait_tab_fn) = $opt->get_opt("TG3C.baits_tab");
my ($baits)       = TG3C::BaitsTab4C::new("TG3C::BaitsTab4C");
$baits->read_tab($bait_tab_fn);

my ($samp_tab_fn) = $opt->get_opt("TG3C.samples_tab");
my ($samps)       = TG3C::SamplesTab4C::new("TG3C::SamplesTab4C");
$samps->read_tab($samp_tab_fn);

my ($sample_id) = $opt->get_opt( "TG3C.sample_id", "" );

my ($sample) = $samps->get_sample($sample_id);
if ( !defined($sample) ) {
    die "Trying to run on non existing sample, id $sample_id\n";
}

my ($exp_nm)    = $sample->{Experiment_name};
my ($sample_nm) = $sample->{Sample_name};

if ( -e "$wd/done_ok" ) {
    unlink("$wd/done_ok");
}

if ( !-e "$wd/$exp_nm.$sample_nm" ) {
    make_path("$wd/$exp_nm.$sample_nm");
}
$wd = "$wd/$exp_nm.$sample_nm";

if ( !-e "$wd/logs" ) {
    make_path("$wd/logs");
}
if ( !-e "$wd/conf" ) {
    make_path("$wd/conf");
}

############################
# Parse only the baits' reads and 
# filter out non-specific reads.
############################
if ( $opt->get_opt( "TG3C.do_strip", 1 ) ) {
    print STDERR "\n3CPipe: will do stripping\n";
    
    my ($fastq_dir) = $sample->{fastqs_dir};
    if ($fastq_dir !~ /\/$/) { #Check if dir ends with backslash and add it if not.
        $fastq_dir .= "/";
    }
    
    my ($fastq_fn_regexp) = $sample->{fastqs_regex};

    #Check for gz files and extract them
    my @dir_files = <$fastq_dir/*fastq.gz>;
    if ( $#dir_files >= 0 ) {
        print STDERR "found gz fastqs - will extract them\n";
        foreach my $fn (@dir_files) {
            my $cmd = "gzip -d $fn";
            $ret = system($cmd);
            die "\nERROR:\n===================\ngzip extract failed\n"
              if ($ret);
        }
    }

    my ($sample_baits) = $sample->{Bait_IDs};
    my (@seq_prefs);
    for ( my ($i) = 0 ; $i <= $#$sample_baits ; $i++ ) {
        push( @seq_prefs,
                $baits->get_bait( $sample_baits->[$i] )->{Bait_seq}
              . $baits->get_bait( $sample_baits->[$i] )->{Bait_pad} );
    }
    my ($seq_pref_regex) = "(^" . join( "|^", @seq_prefs ) . ")(.+)";

    open( CONF, ">$wd/strip.conf" );
    print CONF "inp_base=$fastq_dir\n";
    print CONF "inp_re=$fastq_fn_regexp\n";
    print CONF "out_fn1=$wd/strip_R1.fastq\n";
    print CONF "inp_r1_re=R1\n";
    print CONF "re1=\$R1=~/$seq_pref_regex/i\n";
    print CONF "O1=\$O1=\$S[0].\$S[1]\n";

    print CONF "out_fn2=$wd/strip_R2.fastq\n";
    print CONF "inp_r2_re=R2\n";
    print CONF "r2=1\n";
    print CONF "O2=\$O2=\$R2\n";
    close CONF;

    system("perl $dir_lscripts/TGPipeLib/fastq_splicer.pl \@$wd/strip.conf");

    if ( !-e "$wd/strip_R1.fastq" ) {
        die
"\n======================\nERROR: stripping failed, output does not exist\n";
    }
    if ( -z "$wd/strip_R1.fastq" ) {
        die
          "\n======================\nERROR: stripping failed - output empty\n";
    }
}

#############################
# Segmenet reads by RE seq
#############################
if ( $opt->get_opt( "TG3C.do_seg", 1 ) ) {
    print STDERR "\n3CPipe: will do segmentation\n";

    my ($first_end_code)  = $opt->get_opt( "TG3C.fastq_first_code",  "R1" );
    my ($second_end_code) = $opt->get_opt( "TG3C.fastq_second_code", "R2" );

    my $fastq_dir       = "$wd/";
    my $fastq_fn_regexp = "strip";

    print STDERR "seg $fastq_dir regexp $fastq_fn_regexp\n";
    my ($RE_seq)  = $opt->get_opt("TG3C.RE_seq");
    my ($min_len) = $opt->get_opt("TG3C.segment_min_len");
    my $min_avg_read_qual = $opt->get_opt( "TG3C.min_avg_read_qual", 20 );
    my $min_base_qual     = $opt->get_opt( "TG3C.min_base_qual",     30 );

    my $cmd =
        "perl $dir_lscripts/TG3C/split_fastq_on_re.pl "
      . "--scope_dir $fastq_dir --regexp $fastq_fn_regexp "
      . "--out_fn $wd/split.fastq --re_seq $RE_seq "
      . "--min_frag_len $min_len --min_base_qual $min_base_qual "
      . "--reads_to_sample 10000 --min_len 30 "
      . "--min_read_qual $min_avg_read_qual 2>$wd/logs/seg.log";

    print STDERR "\n$cmd\n";
    $ret = system($cmd);
    if (   $ret
        || !-e "$wd/split.fastq"
        || -z "$wd/split.fastq" )
    {
        die "\n======================\nERROR: segmentation on RE seq failed\n";
    }
}

############################
# Map with bowtie2
############################
if ( $opt->get_opt( "TG3C.do_map", 1 ) ) {
    print STDERR "3CPipe: bowtie mapping\n";

    my ($fq_lst);
    my ($fq_lst_r2);

    my (@fns) = <$wd/split.fastq>;
    $fq_lst = join( ",", @fns );

    my ($bt2_bin) = $opt->get_opt("TG3C.bowtie2_bin");
    if ( !-f $bt2_bin ) {
        die
          "Cannot find bowtie2 binary. Check path to it is set in paths.conf\n";
    }
    my ($bt2_ndx) = $opt->get_opt("TG3C.bowtie2_ndx");
    my ($nthreads) = $opt->get_opt( "TG3C.bowtie2_threads", 2 );
    my $align_mode = $opt->get_opt( "TG3C.bowtie2_align_mode", "end-to-end" );

    my ($sam_fn1) = "$wd/read1.sam";
    print STDERR "invoking bowtie2\n";
    my ($cmd) =
"$bt2_bin -p $nthreads --reorder --$align_mode -x $bt2_ndx -U $fq_lst -S $sam_fn1 2>$wd/logs/bt.log";

    $ret = system($cmd);
    if ($ret) { die "\nERROR:\n===================\nbowtie mapping failed\n"; }
}

############################
# sam->chain->fendchain 
############################
if ( $opt->get_opt( "TG3C.do_fendchain", 1 ) ) {
    print STDERR "3CPipe: sam to fendchain\n";
    my ($dir_redb) = $opt->get_opt("TG3C.redb");
    my ($re)       = $opt->get_opt("TG3C.RE_seq");
    my ($min_qual) = $opt->get_opt( "TG3C.min_map_qual", 30 );
    my ($min_len)  = $opt->get_opt("TG3C.segment_min_len");
    print STDERR "\tchaining\n";
    $ret = system(
"perl $dir_lscripts/TG3C/sam2chain.pl $wd/read1.sam $wd/chain $min_qual $min_len 2>$wd/logs/sam2chain.log"
    );
    if ($ret) { die "\nERROR:\n===================\nsam2chain.pl failed\n"; }
    print STDERR "\tchain to fends\n";
    $ret = system(
"perl $dir_lscripts/TG3C/chain2fendchain.pl $dir_redb $re $wd/chain $wd/fendchain 2>$wd/logs/fendchain.log"
    );

    if ($ret) {
        die "\nERROR:\n===================\nchain2fendchain.pl failed\n";
    }
}

############################
# fendchain to adj table
############################
if ( $opt->get_opt( "TG3C.do_adj", 1 ) ) {
    print STDERR "3CPipe: fendchain to adj\n";
    my ($switch_ratio) = $opt->get_opt("TG3C.switch_ratio");

    my ($fendchain_in) = "$wd/fendchain";

    my ($filter_near_sonic) = $opt->get_opt("TG3C.filter_near_sonic");

    my $cmd =
        "perl $dir_lscripts/TG3C/fendchain2adj.pl "
      . "--fendchain_fn $fendchain_in "
      . "--out_fn $wd/adj "
      . "--switch_ratio $switch_ratio "
      . "--sonication_distance $filter_near_sonic "
      . "2>$wd/logs/fendchain2adj.log";
    $ret = system($cmd);
    if ($ret) {
        die "\nERROR:\n===================\ffendchain2adj.pl failed\n";
    }
}

if ( $opt->get_opt( "TG3C.do_adj_coord", 0 ) ) {
    print STDERR "3CPipe: adj 2 coord\n";
    my ($dir_redb) = $opt->get_opt("TG3C.redb");
    my ($re)       = $opt->get_opt("TG3C.RE_seq");
    if ( -e "$wd/adj.full" ) {
        $ret = system(
"perl $dir_lscripts/TG3C/adj2coords.pl $dir_redb $re $wd/adj,$wd/adj.full 2>$wd/logs/adj2coord.log"
        );
    }
    else {
        $ret = system(
"perl $dir_lscripts/TG3C/adj2coords.pl $dir_redb $re $wd/adj 2>$wd/logs/adj2coord.log"
        );
    }
    if ($ret) { die "\nERROR:\n===================\nadj2coords.pl failed\n"; }
}

############################
# Split adj by baits
############################
if ( $opt->get_opt( "TG3C.split_tracks_by_baits", "TRUE" ) == "TRUE" ) {
    print STDERR "3CPipe: will split tracks to baits\n";
    my ($dir_redb) = $opt->get_opt("TG3C.redb");
    my ($re)       = $opt->get_opt("TG3C.RE_seq");
    my @conf_files = grep( /^\@/, @ARGV );
    my $conf_file = substr( $conf_files[0], 1 );
    my $cmd =
        "perl $dir_lscripts/TG3C/split_mplex_adj.pl "
      . "-adj $wd/adj "
      . "-baits_tab $bait_tab_fn "
      . "-samples_tab $samp_tab_fn "
      . "-sample_id $sample_id "
      . "-redb_dir $dir_redb "
      . "-re_seq $re "
      . "-log_file $wd/logs/bait_umis.log";
    my $ret = system($cmd);
    if ($ret) { warn "Couldn't split adj by baits\n"; }
}

############################
# Get stats
############################
if ( $opt->get_opt( "TG3C.do_QC", 1 ) ) {
    print STDERR "\n3CPipe: will write 3C_stats.txt file\n";
    my @conf_files = grep( /^\@/, @ARGV );
    my $conf_file = substr( $conf_files[0], 1 );
    my $cmd =
        "perl $dir_lscripts/TG3C/QC4C.pl -workdir $wd "
      . " -samples_tab $samp_tab_fn -baits_tab $bait_tab_fn -output "
      . "$wd/3C_stats.txt -sample_id $sample_id";
    my $ret = system($cmd);
    if ($ret) { warn "\nERROR:\n===================\nQC4C.pl failed\n"; }
}

if ( -e "$wd/adj" ) {
    open( OK, ">$wd/done_ok" );
    print OK "ok";
    close OK;
}
