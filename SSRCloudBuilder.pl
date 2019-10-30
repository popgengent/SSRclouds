#Perl -Sx "{0}"; Exit {Status}
#!perl

# There are a bunch of subroutines in here that need will need to be cleaned out if not in use.

# v.02 has been modified from v.01 in the following ways:
#1- determination of the core count is now dependent on the number of annotations sent by SSRFinderv3.
#2- during cloud formation, clouds are expanded only from the most represented oligo and its reverse
# complement. v.01 used every oligo above the core setting for expansion.
#3- minimum locus length has also been eliminated

use strict;
use warnings;

#use merge_clouds_bed qw(merge_locs);
use File::Glob ':glob';
use POSIX;
# use Math::Random;
use Data::Dumper;
use Cwd;
use SSRTools qw(getAlts concatenate);

#globals
my %program =(
	     id => 1,
	     name => "SSRCloudBuilder.pl",
          version => ".02",
	     nickname => "CloudBuilder",
	     authors => "Jonathan Shortt, David D. Pollock, Corey Cox",
	     began => "06/27/2014",
	     modified => "07/14/2015", 	#updateable
	     uses => "Find locations of SSRs and SSR-derived sequences- used in conjunction with SSR_Finderv3",
	     runrecord => "RunRecord_SSR_finder.txt",
	     computer => "Talisker",	#updateable
	    );

#Settings for using control/default/factory system
my $factoryfile = "factory"; 	# factory file setting is meant to be an example, or fixed if *

#main memory
my $silent = 0; my $quiet = 1; my $murmur = 2; my $talk = 3; my $loud = 5; 
my $read = "<"; my $write = ">"; my $append = ">>"; my $read_append = "+>>"; my $no_record = 1;


my $flags = {
    kmerbits => 0,
    recordpos => 1,  ### currently positions are needed to build clouds!!!
    cloudbits => 0,
    cloudid => 1,    ### currently required for proper export and position annotation
    consensus => 1,
    old_himp_clouds => 0,
    rev_com_clouds => 1
};

### Your variables that will be used through the program (main memory) should go here.
my $seqs = {}; my $ptrs = {}; my $kseqs = {}; my $clouds = {}; my $kmers = {};
my $search_kmers = {}; my $locs = []; my $cloud_info = [];

my $num_stringency_levels = 5; #change R script if this is changed
 
use constant KEYNUCS => ('A', 'C', 'G', 'T');

# Read in args
my $fasta_file= $ARGV[0]; my $cloud_file = $ARGV[1]; 
my $most_stringent_core_threshold = $ARGV[2]; 
my $kmer_len = $ARGV[3]; my $shell_factor = $ARGV[4]; 
my $rev_com_clouds = $ARGV[5]; my $cloud_info_file = $ARGV[6]; my $motif = $ARGV[7];

# Main

my $fp_out = openup ($cloud_info_file, $append, "recording information to $cloud_info_file");
readFasta($fasta_file, "building P-Clouds", $seqs );
countKmersInSeqs ($seqs, $kmer_len, $kmers );
my $short_perfs = getKLenMotifs ( $motif, $rev_com_clouds, 12 ); # put perfect 12-mers of repeats in clouds here.
foreach my $short_perf (@$short_perfs) {
	# print STDERR "$short_perf\n";
	$kmers->{$short_perf} = hash_ref();
	my $kinfo = $kmers->{$short_perf};
	$kinfo->{'id'} = $short_perf;
	$kinfo->{'count'} = 999999;
	# select STDOUT; print Dumper ($kmers);
	addKmer2cloud($kinfo, 1);
}
# print STDERR "------------------------------\n";
# select STDOUT; print Dumper ($kmers);
# printx(1, \&exportPClouds, $cloud_file, "exporting PClouds", $kmers, $motif);
# exit;
my $perf_oligos = getKLenMotifs ( $motif, $rev_com_clouds, $kmer_len );
my ($cores, $core_reduction_factor) = getCloudCores ( $most_stringent_core_threshold, $num_stringency_levels, $fp_out, $motif );
putKmersInClouds ( $perf_oligos, $cores, $kmers, $kmer_len );
if (check_for_print ($kmers)) {
	foreach my $oligo_count ( @$cloud_info ) {
		if ( $oligo_count) { print $fp_out "$oligo_count\t"; }
		else { print $fp_out "0\t"; }
	}
	print $fp_out "\n";
	printx(1, \&exportPClouds, $cloud_file, "exporting PClouds", $kmers, $motif);
	# graphClouds($cloud_file);
} 

exit;


# Subs
sub check_for_print { my $kmers = shift;
	foreach my $kmer (keys %$kmers) {
		my $kinfo = $kmers->{$kmer};
		if ($kinfo->{cloud}) { return 1; }
	}
 	return 0;
}

sub getKLenMotifs { my ($motif, $rev_com_clouds, $kmer_len) = @_; my $perf_oligos = [];
	foreach my $alt_motif ( keys %{getAlts($motif, length $motif, $rev_com_clouds)} ) {
		my $klen_motif = concatenate ($alt_motif, $kmer_len);
		push ( @$perf_oligos, $klen_motif );
	}
	return $perf_oligos;
}

sub getCloudCores { my ( $most_stringent_core_threshold, $num_cores, $fp_out, $motif) = @_; my $cores = [];
	my $core_reduction_factor = $most_stringent_core_threshold**(1/$num_cores);
	my $next_core_count = $most_stringent_core_threshold;
	print $fp_out "$motif\t$core_reduction_factor\t";
	for (my $i = 0; $i < $num_cores; $i++) {
		push (@$cores, $next_core_count);
		print $fp_out "$next_core_count\t";
		$next_core_count = (($next_core_count)/($core_reduction_factor));
	}
	return ($cores, $core_reduction_factor);
}

sub graphClouds { my $data = shift;
	my $dir = cwd(); my $datafile = $dir."/".$data;
	my $output1 = $datafile."_NumKmersByCloud"; my $output2 = $datafile."_countsOfKmersByCloud";
	system "Rscript cloudStats.R $datafile $output1 $output2";

}

sub printx { my $flag = shift; my $sub_ref = shift; my $target=shift; my $blurb = shift; my @args = @_;
	print STDERR "Printing Generic $blurb ($flag, $sub_ref, $target, $blurb, @_, etc.)\n\n";
	my $file_handle; my $handle;
	if ($flag > 0) {
		if (!$target) { $handle = \*STDERR; }
        if (!$blurb) { $blurb = "unknown information"; }
		#if ($target eq \*OUTFILE || $target eq \*STDERR || $target eq \*STDOUT) {
           # print STDERR "standard output ID: file pointer is $target\n"; $handle = $target;
       # }
		else {
            $file_handle = openup($target, $write, "for printing via printx");
            print STDERR "output in form of file name:  $target\n\n";  $handle = $file_handle;
        }
		if ($handle) { $sub_ref->($handle, @args);  return; }# call print function and exit sub
	}
    if ($file_handle) { close($file_handle); }
    print STDERR "\nNothing printed for $blurb, flag was off ($flag) or can't open file ($target).\n";
}

sub hashifempty{ 
	my $hash = shift; my $subhash = shift;
	if (!$hash->{$subhash}) { $hash->{$subhash} = hash_ref(); }
	return $hash->{$subhash};
}

sub scalar_ref { my $scalar; return \$scalar; }
sub array_ref { my @array; return \@array; }
sub hash_ref { my %hash; return \%hash; }
sub skip {}

##
##		Basic subroutines for program communication and recording
##
sub openup { my $file = shift; my $mode = shift; my $blurb = shift; my $no_record = shift; if (! $blurb) { $blurb = "unknown"; }
    open (my $file_handle, $mode, $file) or warn "Couldn't open file: $file. \n";
    if (!$no_record) { print STDERR "Opened $file as $blurb.\n"; }
	return $file_handle;
}

sub readFasta { my $target = shift;my $blurb = shift; my $seqs = shift; my $name = ""; my $seq = ""; my $fp;
    if (eval {fileno $target}) { $fp = $target; }
    else { $fp = openup($target, $read, "fasta file for ".$blurb); }
    while (<$fp>) { if ($_ =~ s/>//) { chomp; $name = $_; }; last; } # get first sequence name
    while (<$fp>) { chomp;
        if ($_ =~ s/>//) { $seqs->{$name} = $seq; $name = $_; $seq = ""; }
        else { $seq .= uc($_); }
    }
    close $fp;
}

# This seems like it has some memory liabilities but it is needed for how we currently build clouds.
sub makeSeqHash { my ($name, $seq, $seqs, $ptrs) = @_;
    my $seqinfo = {}; $seqinfo->{name} = $name;
    $seqs->{$seqinfo.""} = $seq; $ptrs->{$seqinfo.""} = $seqinfo;
}

sub putKmersInClouds { my ( $perf_oligos, $cores, $kmers, $kmer_len ) = @_; 
	foreach my $kmer (@$perf_oligos ) {
		if ( $kmers->{$kmer} ) {
			my $kinfo = $kmers->{$kmer};
			if ( checkifeligible ($kinfo, $cores) ) {
				addKmer2cloud ($kinfo, 1);
				$kinfo->{count} = 99999;
				mutateKmer ($kmers, $kmer, $cores, $kmer_len );
			}
		}
	}		
}

sub mutateKmer { my ( $kmers, $kmer, $cores, $kmer_len ) = @_;
	no warnings 'recursion';
	for (my $i = 0; $i < $kmer_len; $i++) {
		foreach my $nuc (KEYNUCS) {
			my $mutant_kmer = $kmer; substr($mutant_kmer, $i, 1, $nuc);
			if ($mutant_kmer eq $kmer) { next; }
			if ($kmers->{$mutant_kmer} ) {					#this line just checks if it exists, wouldn't it be better to do some explicit test like check its count?
				my $kinfo = $kmers->{$mutant_kmer};
				if (checkifeligible ($kinfo, $cores) ) {
					if (!$kinfo->{cloud} || $kinfo->{cloud} > $kmers->{$kmer}->{cloud}) {
                        # print "\tExpanding from $kmer to $mutant_kmer\n";
						my $cloud_number = findShell($kinfo, $cores, $kmers->{$kmer}->{cloud});
						addKmer2cloud ($kmers->{$mutant_kmer}, $cloud_number);
						mutateKmer ($kmers, $mutant_kmer, $cores, $kmer_len );
					}
				}	
			}
            $mutant_kmer = revComp($mutant_kmer);
            if ($mutant_kmer eq $kmer) { next; }
			if ($kmers->{$mutant_kmer} ) {					#this line just checks if it exists, wouldn't it be better to do some explicit test like check its count?
				my $kinfo = $kmers->{$mutant_kmer};
				if (checkifeligible ($kinfo, $cores) ) {
					if (!$kinfo->{cloud} || $kinfo->{cloud} > $kmers->{$kmer}->{cloud}) {
                        # print "\tExpanding from $kmer to $mutant_kmer through rc\n";
						my $cloud_number = findShell($kinfo, $cores, $kmers->{$kmer}->{cloud});
						addKmer2cloud ($kmers->{$mutant_kmer}, $cloud_number);
						mutateKmer ($kmers, $mutant_kmer, $cores, $kmer_len );
					}
				}	
			}
		}
	}
}

sub findShell { my ($kinfo, $cores, $previous_cloud_number) = @_;
	my $count = $kinfo->{count};
	for (my $i = $previous_cloud_number; $i <= $#$cores+1; $i++) {
		if ( $count >= $cores->[$i-1]) {					#looking at index in array, not shell number, index of array is one less than shell number
			return ($i);												
		}
	}
}

sub addKmer2cloud { my ($kinfo, $cloud_number) = @_;
	if ($kinfo->{count} >= 2) { 
		if ( !$kinfo->{cloud} || $kinfo->{cloud} > $cloud_number ) {
			# print STDERR "\t\t$kinfo->{id} placed in cloud $cloud_number\n";
			$kinfo->{cloud} = $cloud_number;
			my $rc_kmer = revComp ($kinfo->{id});
			$cloud_info->[$cloud_number-1]++; 
			if ( $kmers->{$rc_kmer} ) { 												# this might already be done, need to test if necessary
				my $rc_kinfo = $kmers->{$rc_kmer};
				addKmer2cloud ($rc_kinfo, $cloud_number); 
			}
		} else { }
	}
}

sub checkifeligible { my ($kinfo, $cores) = @_;
	if ( $kinfo->{count} >= $cores->[-1] ) { return 1; }
	else { return 0; }
}

sub exportPClouds { my $fp = shift; my $kmers = shift; my $motif = shift;
    foreach my $kmer (sortBySubKey2($kmers, 'count')) { my $kinfo = $kmers->{$kmer};
        if ($kinfo->{cloud}) { print $fp "$kmer\t$motif\t$kinfo->{cloud}\t$kinfo->{count}\n"; }
    }
}

sub revComp { my $newseq = reverse scalar shift; $newseq =~ tr/ACGTacgt/TGCAtgca/; return $newseq; }

## Internal Subs ##
sub countKmersInSeqs { my ($seqs, $kmer_len, $kmers ) = @_;
	while (my ($name, $seq) = each (%$seqs)) {
		for (my $i = 0; ($i + $kmer_len) < length $seq; $i++) {
			my $kmer = substr($seq, $i, $kmer_len); 
			if ($kmer =~ /[ACGTacgt]{$kmer_len}/) { $kmers->{$kmer}->{count}++; $kmers->{$kmer}->{id} = $kmer; }
		}
	}
}

# separate on whitespace (/\s+/) but prevents empty record on leading whitespace
sub goodSplit { my @chunks = split(' ', shift); return @chunks; }
sub peek { my $fp = shift; my $pos = tell($fp); my $line = <$fp>; seek($fp, $pos, 0); return $line; }

# pass array pointer, check existence, if not, create anonymous
sub getKeyLen { my $hash = shift; my @keys = keys %{$hash}; return length($keys[0]); }
sub getKeyArray { my ($hash, $key) = @_; if (!$hash->{$key}) { $hash->{$key} = array_ref(); } return $hash->{$key}; }
sub getKeyHash { my ($hash, $key) = @_; if (!$hash->{$key}) { $hash->{$key} = hash_ref(); } return $hash->{$key}; }
sub getIndexHash { my ($array, $index) = @_; if (!$array->[$index]) { $array->[$index] = hash_ref(); return $array->[$index]; } }
sub sortBySubKey { my $hash = shift; my $key = shift; return sort { $hash->{$b}->{$key} <=> $hash->{$a}->{$key} } keys(%$hash); }
sub sortBySubKey2 { my $hash = shift; my $key = shift; return sort { ($hash->{$b}->{$key} <=> $hash->{$a}->{$key}) || ($a cmp $b) } keys(%$hash); }

