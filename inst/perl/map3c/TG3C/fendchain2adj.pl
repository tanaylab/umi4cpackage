use strict;
use Getopt::Long;


sub edit_d {
        my($s1, $s2) = @_;

        my($d) = 0;
        for(my($i) = 0; $i < length($s1); $i++) {
                my($c1) = substr($s1, $i, 1);
                my($c2) = substr($s2, $i, 1);
                if($c1 ne "N" && $c2 ne "N" && $c1 ne $c2) {
                        $d++;
                }
        }
        return($d);
}

sub is_umi_overlap 
{
    my($umis_a, $umis_b, $max_filter_hamming, $max_sonic_d) = @_;

    my($umis_a1, $umis_a2) = split("#", $umis_a);
    my($umis_b1, $umis_b2) = split("#", $umis_b);

    if(edit_d($umis_a1, $umis_b1) < $max_filter_hamming) {
        $::tot_filt_edit++;
        return(1);
    }

    if($umis_a2 eq $umis_b2) {
        $::tot_filt_umi2++;
        return(1);
    }

    my(@f_a2) = split(":", $umis_a2);
    my(@f_b2) = split(":", $umis_b2);
    my($end_point_a2) = $f_a2[$#f_a2]; 
    my($end_point_b2) = $f_b2[$#f_b2]; 
    my($fend_a2) = $f_a2[$#f_a2-5]; 
    my($fend_b2) = $f_b2[$#f_b2-5]; 
    if($fend_a2 == $fend_b2
    && abs($end_point_a2 - $end_point_b2) <= $max_sonic_d) {
        $::tot_filt_sonic++;
        return(1);
    }
    return(0);
}

###############Code

if($#ARGV < 1) {
    die "usage $0 --fendchain_fn fc_fn [--out_fn adj] --switch_ratio FLOAT --max_hamming_distance INT --sonication_distance INT\n";
}

my($fc_tab);
my($switch_ratio) = 0.1;
my($max_filter_hamming) = 3;
my($max_sonic_d) = 1;
my($out_fn) = "adj";

GetOptions("fendchain_fn=s"         => \$fc_tab,
           "out_fn:s"               => \$out_fn,
           "switch_ratio:f"         => \$switch_ratio,
           "max_hamming_distance:i" => \$max_filter_hamming,
           "sonication_distance:i"  => \$max_sonic_d)
or die "not all args supplied\n";

my(%adj_umis);

open(FC, $fc_tab) || die "cannot open fendchain tab $fc_tab\n";
open(my $out_fh, ">", $out_fn) or die "cannot write output to $out_fn\n";
open(my $out_full_fh, ">", "$out_fn.full") or die "cannot write output to $out_fn.full\n";

print $out_fh "fend1\tfend2\tcount\n";
print $out_full_fh "fend1\tfend2\tmol\treads\tpair_max_r\n";

my($count) = 0;
while(<FC>) {
	chop;
	my(@f) = split("\t", $_);
	if($f[1] == 1) {
		next;
	}

        if($count % 100000 == 0) {
                print STDERR "count ".int($count/100000)."M..";
        }
        $count++;

#HISEQ:171:HALHWADXX:1:1101:2949:2211:TCCCCGTGGG	2	3392655:1+:30:744:774:1:42	3392656:g+:221:-1:220:1565:42

    #trimming quality tags from fc entries
        
        for(my($i) =2; $i <= $#f; $i++) {
                $f[$i] =~s/:\d+$//;
        }
        
        my(@head_f) = split(":", $f[0]);
        my($umi1) = $head_f[$#head_f];
        my($umi2) = join("\t", @f[2..$#f]);

        my(@fe1) = split(":", $f[2]);
        my(@fe2) = split(":", $f[3]);
	my($fe_id1) = $fe1[0]*2 + (($fe1[1]=~/\+/) ? 1 : 0);
	my($fe_id2) = $fe2[0]*2 + (($fe2[1]=~/\+/) ? 0 : 1);

        if(!exists($adj_umis{"$fe_id1\t$fe_id2"})) {
            $adj_umis{"$fe_id1\t$fe_id2"} = {};
        }
        $adj_umis{"$fe_id1\t$fe_id2"}->{"$umi1#$umi2"}++;
}

my($adj);

foreach $adj (keys %adj_umis) {
    my($umis_r) = $adj_umis{$adj};
    my(@umis) = keys %$umis_r;

    my(@sumis) = sort { $umis_r->{$b} <=> $umis_r->{$a} } @umis;

    my(@umi_mark);

    my($tot_reads) = 0;
    my($tot_mols) = 0;

    my(@mol_reads);
    my($max_mol_reads) = 0;

    for(my($i) = 0; $i <= $#sumis; $i++) {
        if(defined($umi_mark[$i])) {
            next;
        }
        my($mreads) = $umis_r->{$sumis[$i]};
        for(my($j) = $i+1; $j <= $#sumis; $j++) {
            if(defined($umi_mark[$j])) {
                next;
            }
            if(is_umi_overlap($sumis[$i], $sumis[$j], $max_filter_hamming, $max_sonic_d)) {
                $umi_mark[$j] = $i;
                $mreads += $umis_r->{$sumis[$j]};
            }
        }
        push(@mol_reads, $mreads);
        if($mreads > $max_mol_reads) {
            $max_mol_reads = $mreads;
        }
        $tot_reads += $mreads;
    }

    for(my($mol) = 0; $mol <= $#mol_reads; $mol++) {
        if($mol_reads[$mol] > $max_mol_reads * $switch_ratio) {
            $tot_mols++;
        } else {
            $::tot_filt_switch++;
        }
    }

    print $out_full_fh "$adj\t$tot_mols\t$tot_reads\t$max_mol_reads\n";
    print $out_fh "$adj\t$tot_mols\n";
}

print STDERR "\n#########Done, filt stat######\n";
print STDERR "Sonic (thres=$max_sonic_d) ".($::tot_filt_sonic)."\n";
print STDERR "Switch (thres=$switch_ratio) ".($::tot_filt_switch)."\n";
print STDERR "Edit (hamming dist=$max_filter_hamming) ".($::tot_filt_edit)."\n";
print STDERR "Umi2 ".($::tot_filt_umi2)."\n";
