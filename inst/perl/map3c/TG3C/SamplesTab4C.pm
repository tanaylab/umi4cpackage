use strict;

package TG3C::SamplesTab4C;

#use diagnostics;

1;

#Sample_ID	Sample_name	Experiment_name	fastqs_regex	fastqs_dir	Bait_IDs
#1	CMK	example	CMK	example_data/	1	

sub new($) {
	my($clas) = @_;

	my($self) = {};
	
	bless($self,$clas);

        $self->{samp} = [];
        
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
            if (($head[$i] eq "Sample_name" || $head[$i] eq "Experiment_name") && 
                $f[$i] =~ m/[^a-zA-Z0-9_]/) {
                die("\nError in $tab_fn: value was $f[$i]:
                    Sample name and Experiment name must contain only alphanumeric characters and \"_\"\n\n");
            } if ($f[$i] =~ m/\s/) {
                die("Error in $tab_fn: in $f[$i]: values cannot contain whitespaces!\n");
            }
            if($head[$i] eq "Bait_IDs") {
                my($bait_ids) = [];
                @$bait_ids = split(",", $f[$i]);
                $b->{$head[$i]} = $bait_ids;
            } else {
                $b->{$head[$i]} = $f[$i];
            }
        } 
        $self->{samp}->[$id] = $b;
    }
}

sub get_sample($$) {
    my($self, $id) = @_;

    if(!defined($self->{samp}->[$id])) {
        return(undef);
    } else {
        return($self->{samp}->[$id]);
    }
}

