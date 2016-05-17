use strict;
use FindBin;
use lib ("$FindBin::Bin/..", "$FindBin::Bin/../..");

require map3c::ReDB;

#read re db from tab


if($#ARGV != 3) {
	die "usage chain2fendchain.pl re_db re chain out_fendchain\n";
}

my($re) = $ARGV[1];
my(@chain_fns) = <$ARGV[2]>;

if($#chain_fns == -1) {
	die "no input chains in $ARGV[2]\n";
}

my($redb) = ReDB::new("ReDB");
print STDERR "will init re db\n";
$redb->init_from_tabs($ARGV[0]."/$re.frags");
print STDERR "done read tab\n";

my($search_scope) = 1500;

open(OUT, ">$ARGV[3]") || die "cannot write fend chain tab\n";
open(ERR, ">$ARGV[3].err");

my($chain_fn);

foreach $chain_fn (@chain_fns) {

	open(CHAIN, $chain_fn) || die "cannot open chain tab\n";

	print STDERR "fend transoforming $chain_fn\n";

	while(<CHAIN>) {
		chop;
		my($id, $num, @f) = split("\t", $_);
		if($num =~/ERR/) {
			next;
		}
		my($suf) = "$num";
		my(@c_fids);
		my(@c_pos);
		my(@c_maxpos);
		my(@c_end);
		my(@c_strand);
		my(@c_qual);
		my($err) = "";
		for(my($i) = 0; $i <= $#f; $i++) {
			my($rfr, $rto, $chr, $pos, $strand, $qual, $end) = $f[$i]=~/R([\-]*\d+)-([\-]*\d+):(.+):(\d+):(.):(\d+):(.)/;
			if(!defined($rfr)) {
				die "bad chain format, $f[$i], out of here\n";
			} 
			#map 
			my($maxpos) = $pos + abs($rto-$rfr);
			my($fid1) = $redb->map_frag($chr, $pos+15, $search_scope);
			my($fid2) = $redb->map_frag($chr, $maxpos-15, $search_scope);
			if($fid1 != $fid2) {
				my($off5) = $redb->frag_5_dist($fid1, $pos);
				my($off3) = $redb->frag_3_dist($fid1, $maxpos);
				my($off52) = $redb->frag_5_dist($fid2, $pos);
				my($off32) = $redb->frag_3_dist($fid2, $maxpos);
				print STDERR "incompatible frags $fid1 $fid2 at $chr $pos $maxpos offs $off5 $off3 offs2 $off52 $off32 end $end and fi $f[$i]\n";
				$qual = "missRE";
		#		next;
			}
			
			push(@c_fids, $fid1);
			push(@c_pos, $pos);
			push(@c_maxpos, $maxpos);
			push(@c_end, $end);
			push(@c_strand, $strand);
			push(@c_qual, $qual);
		}
		my($cf) = "";
		for(my($i) = 0; $i <= $#c_fids; $i++) {
			while($i != $#c_fids
			&& $c_fids[$i] == $c_fids[$i+1]) {
				if($c_strand[$i] == "+"
				&& $c_strand[$i+1] == "-"
				&& $c_maxpos[$i] <= $c_pos[$i+1]+5) {
					$c_strand[$i+1] = "+";
					$c_pos[$i+1] = $c_pos[$i];
					$c_end[$i+1] = "g";
					$num--;
					$i++;
				} elsif($c_strand[$i] == "-"
				&& $c_strand[$i+1] == "+"
				&& $c_pos[$i] >= $c_maxpos[$i+1]-5) {
					$c_strand[$i+1] = "-";
					$c_maxpos[$i+1] = $c_maxpos[$i];
					$c_end[$i+1] = "g";
					$num--;
					$i++;
				} else {
					$i++;
					$err = "FID_Overlap";
				}
			}
			
			my($off1);
			my($off12);
			my($off2);
			if($c_strand[$i] eq "+") {
				$off1 = $redb->frag_5_dist($c_fids[$i], $c_pos[$i]);
				$off12 = $redb->frag_5_dist($c_fids[$i], $c_maxpos[$i]);
				$off2 = $redb->frag_3_dist($c_fids[$i], $c_maxpos[$i]);
			} else {
				$off2 = $redb->frag_5_dist($c_fids[$i], $c_pos[$i]);
				$off1 = $redb->frag_3_dist($c_fids[$i], $c_maxpos[$i]);
				$off12 = $redb->frag_3_dist($c_fids[$i], $c_pos[$i]);
			}
			$cf .= "\t$c_fids[$i]:$c_end[$i]$c_strand[$i]:";
			$cf .= ($c_maxpos[$i] - $c_pos[$i]);
			$cf .= ":$off1:$off12:$off2:$c_qual[$i]";
		}
		if($err eq "") {
			print OUT "$id\t$num$cf\n";
		} else {
			print ERR "$id\t$num$cf\n";
		}
	}
}
