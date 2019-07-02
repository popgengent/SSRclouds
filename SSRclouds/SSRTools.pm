#Perl -Sx "{0}"; Exit {Status}
#!perl

package SSRTools;
use strict;
use warnings;
use Data::Dumper;

require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '0.01';
@ISA         = qw(Exporter);
@EXPORT_OK   = qw(makeMotifs makeOligos getAlts concatenate);

my $silent = 0; my $quiet = 1; my $murmur = 2; my $talk = 3; my $loud = 5;
my $read = "<"; my $write = ">"; my $append = ">>"; my $read_append = "+>>"; my $no_record = 1;

my @nucs = ("A", "C", "T", "G");  # nucleotides

# subs will go beneath here
sub makeMotifs{ my ($motifs, $IDs, $motifmax, $revcom) = @_;
    for (my $motiflen = 1; $motiflen <= $motifmax; $motiflen++) {  # length is inclusive
        if ($motiflen eq 1){ getAltMotifs("", $motifs, $IDs, $motifmax, $revcom); }
        else { foreach my $submotif (keys %{$motifs}) {
            if ( length $submotif  eq ( $motiflen - 1 )) { 
                getAltMotifs( $submotif, $motifs, $IDs, $motifmax, $revcom); 
            }}
        }
    }
}

# makes concatemers of a supplied motif of the length specified by oligolenth
sub makeOligos{  my ($motifs, $oligolength )=@_;
	$motifs->{alloligos} = hash_ref();
	foreach my $motif (keys %$motifs) {
		$motifs->{$motif}->{oligos} = hash_ref();
		$motifs->{$motif}->{oligofams} = hash_ref();
		my $motifam = $motifs->{$motif};
		foreach my $alt ( keys %{$motifam->{alts}}) { 
			my $oligo = concatenate($alt, $oligolength, length $alt);
			if( !$motifam->{oligos}->{rc($oligo)} ) {      
				$motifam->{oligos}->{$oligo} = rc ($oligo);
			} else {}
			$motifam->{oligofams}->{$oligo} = 1;
			$motifs->{alloligos}->{$oligo} = 1; $motifs->{alloligos}->{rc($oligo)} = 1;
		} 
	}
}


# creates new motifs by adding individual nucleotides to submotifs provided in make Motifs, assign len and 
# alternative motifs by passing newly created motifs to getAlts
sub getAltMotifs{ my ($submotif, $newmotifs, $IDs, $motifmax, $revcom) = @_;
    foreach my $nuc (@nucs) {
        my $motif = $submotif.$nuc;
        if( !$newmotifs->{$motif} ) {  # don't redo ones that have already been found
            $newmotifs->{$motif} = hashifempty ($newmotifs, $motif); 
            my $motifam = $newmotifs->{$motif}; 
            $motifam->{len} = length $motif; 
            $motifam->{alts} = getAlts($motif, $motifam->{len}, $revcom); # anonymous hash of alternative motifs 
			$motifam->{id} = (sort (keys(%{$motifam->{alts}})))[0];
			push (@$IDs, $motifam->{id});
            my @concats;  makeConcats($motifam->{alts}, $motifam->{len}, $motifmax, \@concats);
            foreach my $newalt (@concats) { $motifam->{alts}->{$newalt}= 1; }
            foreach my $alt (keys %{$motifam->{alts}}) { $newmotifs->{$alt} = $newmotifs->{$motif}; }
       } else { } # already exists, don't recreate it
    }
}


 # creates concatemers of whole motifs until concatemer exceeds specified length
sub makeConcats{ my ($alts, $len, $max, $concats) = @_; 
    foreach my $alt (keys %{$alts}){
        for (my $multiple = 2; $len * $multiple <= $max; $multiple++){
            push(@$concats, concatenate($alt, $len * $multiple, $len)); 
        }
    }
}

# return translations and reverse complements of input motif
sub getAlts { my ($motif, $len, $revcom) = @_; 	
	my $offset = -$len + 1;
	my $altmotifs =  hash_ref();
	for (my $i = 0; $i < $len; $i++) {
		my $newmotif = substr ($motif, $offset + $i).substr ($motif, 0, $offset + $i);
        $altmotifs->{$newmotif}=1; if ( $revcom ) { $altmotifs->{rc($newmotif)}=1; }
	}
	return $altmotifs; # need to remember this to avoid memory leak
}

# concatenate motif to return an oligo of input length
sub concatenate {	my $motif = shift; my $oligolen = shift; my $len = length $motif; 
	my $repeats = int $oligolen/$len;
	my $oligo = "";
	for (my $i = 0; $i <= $repeats; $i++) { $oligo .= $motif; }
	return substr ($oligo, 0, $oligolen);
}

# takes supplied sequence and makes a reverse complement of it
sub rc { my $motif = shift; 
	my $revcom = reverse $motif;
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	return $revcom;
}

sub hashifempty{ 
	my $hash = shift; my $subhash = shift;
	if (!$hash->{$subhash}) { $hash->{$subhash} = hash_ref(); }
	return $hash->{$subhash};
}

sub scalar_ref { my $scalar; return \$scalar; }
sub array_ref { my @array; return \@array; }
sub hash_ref { my %hash; return \%hash; }
# subs go above here
1;
