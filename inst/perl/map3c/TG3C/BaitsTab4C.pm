use strict;

package TG3C::BaitsTab4C;

#use diagnostics;

1;

#Bait_ID	Bait_name	Bait_seq	Bait_pad	Bait_chr	Bait_coord
#1	ANK1_TSS	AGTTTAGCAGACTCAAAGGAAAGC	CTCTAAGATC	8	41654693
#
sub new($) {
	my($clas) = @_;

	my($self) = {};
	
	bless($self,$clas);

        $self->{baits} = [];

        return($self);
}

sub read_tab($$) {
    my($self, $tab_fn) = @_;
    open(TAB, $tab_fn) || die "cannot open bait tab $tab_fn\n";

    my($h);
    $h = <TAB>;
    chop $h;
    my(@head) = split("\t", $h);

    while(<TAB>) {
        chop;
        my(@f) = split("\t", $_);
        
        my($id) = $f[0];
        my($b) = {};

        for(my($i) = 1; $i <= $#head; $i++) {
            if (($head[$i] eq "Bait_name") && 
                $f[$i] =~ m/[^a-zA-Z0-9_]/) {
                die("Error in $tab_fn: value was $f[$i]:
                    bait name must contain only alphanumeric characters and \"_\"\n");
            } if ($f[$i] =~ m/\s/) {
                die("Error in $tab_fn: in $f[$i]: values cannot contain whitespaces!\n");
            }
            $b->{$head[$i]} = $f[$i];
        } 
        $self->{baits}->[$id] = $b;
    }
}

sub get_bait($$) {
    my($self, $id) = @_;
    return($self->{baits}->[$id]);
}

