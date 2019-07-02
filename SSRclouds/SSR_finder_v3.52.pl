#Perl -Sx "{0}"; Exit {Status}
#!perl

use strict;
use warnings;

require "MultithreadCommands.pl";

use File::Glob ':glob';
use POSIX;
# use Math::Random;
use Data::Dumper;

use globals qw(ProgSetUp finish);
use kClouds qw(buildPClouds importClouds exportPClouds getCloudKmers revCompKmers SSR2Cloud);
use seqTools qw(readFasta makeSeqHash);
use SSRTools qw(makeMotifs makeOligos getAlts concatenate);
use Getopt::Long qw(GetOptions);
use basics qw(readx printx openup hashifempty scalar_ref array_ref hash_ref skip talk record);


#globals
my %program =(
	     id => 1,
	     name => "SSR_finder.pl",
         version => "3.5",
	     nickname => "SSR_finder",
	     authors => "Jonathan Shortt, David D. Pollock, Corey Cox",
	     began => "05/01/2014",
	     modified => "08/6/2016", 	#updateable
	     uses => "Find locations of SSRs and SSR-derived sequences",
	     runrecord => "RunRecord_SSR_finder.txt",
	     computer => "Talisker",	#updateable
	    );

#Settings for using control/default/factory system
my $globals = hash_ref();
my $factoryfile = "factory"; 	# factory file setting is meant to be an example, or fixed if *

#main memory
my $silent = 0; my $quiet = 1; my $murmur = 2; my $talk = 3; my $loud = 5; 
my $read = "<"; my $write = ">"; my $append = ">>"; my $read_append = "+>>"; my $no_record = 1;

### Your variables that will be used through the program (main memory) should go here.
my $clouds = {}; my $kmers = {}; my $ptrs = {};
my $search_kmers = {}; #my $locs = [];  # $locs might not be used ever as a global

my $motifs =  hash_ref(); my $IDs = array_ref(); my $genome_coords = hash_ref();
my $FPdist = hash_ref();            # freqs at which we can expect different length fp
my $FPdata = hash_ref();

#############################
##                         ##
##       Begin Main        ##
##                         ##
#############################

#system("purge");	# clear out inactive memory - on my system each purge takees about 8 seconds to run
ProgSetUp($factoryfile, $globals, \%program);	# read control files, import all values in $globals hash
my $flags = {
    kmerbits => 1,
    recordpos => 1,  ### currently positions are needed to build clouds!!!
    cloudbits => 0,
    cloudid => 1,    ### currently required for proper export and position annotation
    consensus => 1,
    old_himp_clouds => 0,
    rev_com_clouds => 0	## used to add reverse complement of given kmer to same cloud as kmer
};
if ($globals->{'rev_com'} eq 1) { $flags->{'rev_com_clouds'} = 1 };

GetOptions('runmode=s' => \$globals->{'runmode'},
	'dir=s' => \$globals->{'dir'},
	'fasta:s' => \$globals->{'fasta'},
	'out:s' => \$globals->{'out'},
	'merge:s' => \$globals->{'merge'}, 
	'genome:s' => \$globals->{'genome'},
	'assign:s' => \$globals->{'assign'},
	'mode:s' => \$globals->{'mode'},
	'FPTable:s' => \$globals->{'FPTable'},
	'bed:s' => \$globals->{'bed'}, 
	'thread_count:i' => \$globals->{'thread_count'} ) or die "Usage: $0 --runmode [runmode]\n";
if (!$globals->{'runmode'} || !$globals->{'dir'}) { die "Usage: $0 [OPTIONS] --runmode [runmode] --dir [directory]\n"; }

### Create Microsats
talk ("Creating motifs up to $globals->{motifmax}bp in length\n", $murmur, $globals);
makeMotifs ( $motifs, $IDs, $globals->{motifmax}, $globals->{rev_com} );
talk ("Creating oligos of $globals->{oligolength}bp from motifs\n", $murmur, $globals);
makeOligos ( $motifs, $globals->{oligolength});
my $motifs_of_interest = findMotifsOfInterest ();

if ( $globals->{'runmode'} eq "annotatePerfects" ) {
	talk ("Running in mode $globals->{'runmode'}. Beginning annotation of perfect repeats\n", $talk, $globals);
	checkOptions ($globals, 'fasta', 'dir');
	checkInputFiles ($globals->{'fasta'});
	writeLog ($globals);                      
	my $assign_folder = $globals->{'dir'}."/".$globals->{'perfect_assign_dir'}; 
	mkdir $assign_folder if not -d $assign_folder; 
	my $bed_folder = $globals->{'dir'}."/".$globals->{'bed_folder'};
	mkdir $bed_folder if not -d $bed_folder;
	makeAssignFiles($motifs_of_interest, $assign_folder); 
	annotatePerfectAssignFiles($motifs_of_interest, $assign_folder, $bed_folder);	
		
} elsif ($globals->{'runmode'} eq "combinePerfects") {
	talk ("Running in mode $globals->{'runmode'}. Combining perfect annotations.\n", $talk, $globals);
	checkOptions ($globals, 'dir', 'out');
	writeLog ($globals);
	my $beds = array_ref();
	foreach my $motif (@$motifs_of_interest) {
		my $bed = get_file_name ($motif, $globals->{'dir'}, ".bed");
		if (-s $bed) { push (@$beds, $bed); }
		# else { warn "No bed file for $motif found.\n"; }
	}
	if ($#$beds == -1) { die "No files found to merge\n"; }
	talk ("Found $#$beds of $#$motifs_of_interest desired files. Combining now\n", $talk, $globals);
	my $cmd;
	if ($globals->{'merge'}) { 
		$cmd = "cat @$beds | sortBed | bedtools merge -i stdin -d $globals->{'merge'} > $globals->{'out'}"; 
		print STDERR "executing: concatenating, sorting, and merging files of interest\n";
		}
	else { 
		$cmd = "cat @$beds | sortBed  > $globals->{'out'}"; 
		print STDERR "executing: concatenating and sorting files of interest\n"; 
	}
	system $cmd;
	
} elsif ( $globals->{'runmode'} eq "buildClouds" ) {
	talk ("Running in mode $globals->{'runmode'}. Create Clouds from Sequence.\n", $talk, $globals);
	checkOptions ($globals, 'dir', 'genome', 'fasta'); #send a blurb to check options with running instructions?
	checkInputFiles ($globals->{'fasta'}, $globals->{'genome'} );
	writeLog ($globals);
	my $cloud_info_file = startCloudInfoFile($globals->{'dir'});
	makeClouds($motifs_of_interest, $globals->{'dir'}, $cloud_info_file);	
	
} elsif ($globals->{'runmode'} eq "combineClouds" ) {
	talk ("Running in mode $globals->{'runmode'}. Combining clouds.\n", $talk, $globals);
	checkOptions ($globals, 'dir', 'out' );
	writeLog ($globals);
	my $assigns = array_ref();
	foreach my $motif (@$motifs_of_interest) {
		my $assign = get_file_name ($motif, $globals->{'dir'}."/".$globals->{'assign_folder'}, ".assign");
		if (-s $assign) { push (@$assigns, $assign); }
		# else { warn "No bed file for $motif found.\n"; }
	}
	talk ("Found $#$assigns of $#$motifs_of_interest desired files. Combining now\n", $talk, $globals);
	if ($#$assigns == -1) { die "No files found to merge\n"; }
	my $cmd = "cat @$assigns > $globals->{'out'}";
	print STDERR "executing: $cmd\n";
	system $cmd;

} elsif ( $globals->{'runmode'} eq "annotateClouds" ) {
	talk ("Running in mode $globals->{'runmode'}. Annotating clouds.\n", $talk, $globals);
	checkOptions ($globals, 'dir', 'out', 'fasta', 'assign', 'mode', 'FPTable' );
	checkInputFiles ( $globals->{'fasta'}, $globals->{'assign'} );
	writeLog ($globals);
	multiAnnotateAssignFiles( );
	
} elsif ($globals->{'runmode'} eq "exportFPR") {
	talk ("Running in mode $globals->{'runmode'}. Recording FPR from supplied annotation.\n", $talk, $globals);
	checkOptions ($globals, 'dir', 'out', 'bed', 'fasta' );
	checkInputFiles ( $globals->{'bed'}, $globals->{'fasta'} );
	writeLog ($globals);
	getFPFromSim ($FPdist, $FPdata, $globals->{'bed'}, $globals->{'fasta'}, $globals->{'out'}); 
	
} else { talk ("Unrecognized runmode $globals->{'runmode'}\n", $talk, $globals); }

#system("purge"); # clear out inactive memory
finish($globals->{starttime}); # finish program and close out

#############################
##                         ##
##        End Main         ##
##                         ##
#############################

##					##
##   subroutines 	##
##					##

sub writeLog { my $globals = shift;
	my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($globals->{'starttime'});
	my $time = join (":", $hour, $min, $sec);
	my $date = join("_",$mday, $months[$mon], 1900+$year);
	my $combo = $date."_".$time.".log";
	mkdir $globals->{'dir'}."/LOGS" if not -d $globals->{'dir'}."/LOGS";
	open (my $fp, '>', $globals->{'dir'}."/LOGS/".$combo) or die "Couldn't open log for writing\n";
	foreach my $param (sort (keys %{$globals}) ) {
		my $val = $globals->{$param};
		if ($val) { print $fp "$param = $val\n"; }
		# print STDERR "$param = $val\n";
	}
	close ($fp);
}

sub checkOptions { my $globals = shift;
	foreach my $option (@_) {
		if (!$globals->{$option}) {die "Missing option --$option\n"; }
	}
}

sub checkInputFiles { my @files = @_;
	foreach my $file (@files) {
		if ( -s $file ) { print STDERR "Found $file\n"; }
		else { die "$file not found\n"; }
	}	
}

sub findMotifsOfInterest {
	my $motifs_of_interest = [];
	if ( $globals->{motifs_of_interest} ) {
		if (ref($globals->{motifs_of_interest}) eq "ARRAY") { #if (ref($globals->{ind_select}) eq "ARRAY") {
			foreach my $motif ( @{$globals->{motifs_of_interest}} ) {
				if ($motifs->{$motif}) { push (@$motifs_of_interest, $motifs->{$motif}->{id}); }
				else {talk("Unrecognized motif $motif in control file 'motifs of interest'\n", $murmur, $globals); die; }
			}
		} else {
			my $motif = $globals->{motifs_of_interest};
			if ($motifs->{$motif}) { push (@$motifs_of_interest, $motifs->{$motif}->{id}); }
			else {talk("Unrecognized motif $motif in control file 'motifs of interest'\n", $murmur, $globals); die; }
		}
	} else {
		talk("No motifs of interest found in control file. Analysis will use all motifs.\n", $murmur, $globals);
		foreach my $motif (@$IDs) {
			push (@$motifs_of_interest, $motif);
		}
	}
	return $motifs_of_interest;	
}

sub makeAssignFiles { my ($motifs_of_interest, $run_folder) = @_;	
	mkdir $run_folder if not -d $run_folder;
	foreach my $motif ( @$motifs_of_interest ) {
		my $filename = get_file_name($motif, $run_folder, ".assign") ;
		create_assign_file ($motif, $filename );  
	} 
}

sub create_assign_file { my ($motif, $filename ) = @_;
	my $fp = openup ($filename, $write, "\r\e[Kprinting perfect oligos", $globals);
	foreach my $oligo ( keys %{$motifs->{$motif}->{oligofams}} ) {
		print $fp "$oligo\t$motif\t900\n";
	}
}

sub get_file_name { my ($motif, $run_folder, $suffix) = @_;
	my $filename = $run_folder."/fam_".$motif."_"."rc".$globals->{'rev_com'}."_".$globals->{'oligolength'}."-mers".$suffix;
	return $filename;
}

sub annotatePerfectAssignFiles { my ($motifs_of_interest, $assign_folder, $bed_folder)= @_; my $info_array= array_ref();
	foreach my $motif (@$motifs_of_interest) {
		my $assign_file = get_file_name($motif, $assign_folder, ".assign");
		my $bed_file = get_file_name ($motif, $bed_folder, ".bed");
		push (@$info_array, { assign => $assign_file, bed => $bed_file, motif =>$motif } );
	}
	ghettoThreadAnnotateAssignFiles($info_array, $globals->{'fasta'});
}

sub ghettoThreadAnnotateAssignFiles { my ($file_info_array, $fasta)= @_;
	my $thread_count = $globals->{'thread_count'};
	for ( my $index = 1; $index <= $#$file_info_array+1; $index++) { #start at index 1 so the modulus doesn't stop the ghettoThread on the first round
		my $file_info = $file_info_array->[$index-1];
		my $assign_file = $file_info->{'assign'}; my $bed_file = $file_info->{'bed'};
		my $motif = $file_info->{'motif'};
		if ( $index % $globals->{'thread_count'} ) {
			# talk ("Beginning annotation of $motif and moving on\n", $murmur, $globals);
			my $cmd = "./cloudPos -c $assign_file -f $fasta -m 0 -p 0 -r 0 > $bed_file &";
			system $cmd;
		} else {
			talk ("\r\e[KBeginning annotation of $motif and waiting to finish. $index of $#$file_info_array", $murmur, $globals);
			my $cmd = "./cloudPos -c $assign_file -f $fasta -m 0 -p 0 -r 0 > $bed_file "; 
			system $cmd;
		}
	}
	talk ("\nAnnotation complete\n", $murmur, $globals);
}

# prints header for the cloud info file
sub startCloudInfoFile { my $run_folder = shift;
	my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($globals->{'starttime'});
	my $time = join (":", $hour, $min, $sec);
	my $date = join("_",$mday, $months[$mon], 1900+$year);
	my $cloud_info_file = $run_folder."/".$globals->{'cloud_info_file'}."_rc".$globals->{'rev_com'}."_$globals->{'cloud_kmer_length'}"."-mers_".$date."_".$time.".cloudinfo";
	my $fp_out = openup ($cloud_info_file, $write, "recording cloud info to $cloud_info_file.\n", $globals);
	print $fp_out "File\tshell_reduction\t";
	print $fp_out "cloud1_threshold\tcloud2_threshold\tcloud3_threshold\tcloud4_threshold\tcloud5_threshold\t";
	print $fp_out "cloud1_kmercount\tcloud2_kmercount\tcloud3_kmercount\tcloud4_kmercount\tcloud5_kmercount\n";
	return $cloud_info_file;
}

# prepares folders and files to store clouds
sub makeClouds { my ($motifs_of_interest, $run_folder, $cloud_info_file ) = @_; my $cloudinfo = [];
	my $bed_folder = $run_folder."/".$globals->{'bed_folder'}; die "Cannot locate $bed_folder\n" if not -d $bed_folder;
	create_training_and_test_sets ( $bed_folder );
	my $assign_folder = $run_folder."/".$globals->{'assign_folder'}; 
	my $fasta_folder = $run_folder."/".$globals->{'fasta_folder'};
	mkdir $assign_folder if not -d $assign_folder; mkdir $fasta_folder if not -d $fasta_folder;
	my $builds = 0; 
	foreach my $motif (@$motifs_of_interest) {
		# talk("Beginning cloud formation for $motif\n", $murmur, $globals);
		my $cloud_file = get_file_name ($motif, $assign_folder, ".assign");
		my $bed_file = get_file_name ($motif, $bed_folder."_TRAIN", ".bed");
		if (-s $bed_file) { } #talk ("Found $bed_file\n", $murmur, $globals); }
		else { 
			# talk ("No bed for $motif could be found. Shucks. Will not make clouds for $motif.\n", $talk, $globals);
			makePerfectClouds ($motif, $cloud_file);  #not sure if I should do this here or not- this else loop is entered when the bed file doesn't exist so that means no perfects were found so no point in creating the oligos for clouds
			next;
		}
		my $max_core = getMostStringentCoreThreshold ( $bed_file, $motif );
		if ($max_core == 0) {
			# talk("Motif $motif doesn't have enough loci for cloud formation. Skipping\n\n", $murmur, $globals);
			makePerfectClouds ($motif, $cloud_file);
			next;
		}
		$builds++;
		my $fasta = makeCloudFasta ($motif, $fasta_folder, $bed_file, $globals->{'slop_len'}, $globals->{'genome'});
		push (@$cloudinfo, { max_core=> $max_core, fasta=> $fasta, cloud_file=> $cloud_file, motif=> $motif } );
		# print STDERR "\n";
	}
	talk ("Clouds will be built for $builds families of $#$motifs_of_interest total families.\n", $murmur, $globals);
	multiThreadCloudBuilding( $cloudinfo, $cloud_info_file);
}

sub makePerfectClouds { my $motif = shift; my $file = shift;
	open (my $fp, '>', $file) or die "Couldn't open $file for printing\n";
	my $short_perfs = getKLenMotifs ( $motif, $flags->{'rev_com_clouds'}, $globals->{'oligolength'} ); 
	foreach my $short_perf (@$short_perfs) {
		my $count = 999999;
		print $fp "$short_perf\t$motif\t1\t$count\n";		
	}
	close $fp;
}

sub getKLenMotifs { my ($motif, $rev_com_clouds, $kmer_len) = @_; my $perf_oligos = [];
	foreach my $alt_motif ( keys %{getAlts($motif, length $motif, $rev_com_clouds)} ) {
		my $klen_motif = concatenate ($alt_motif, $kmer_len);
		push ( @$perf_oligos, $klen_motif );
	}
	return $perf_oligos;
}

# genome, seq
sub makeCloudFasta { my ($motif, $run_folder, $bed_file, $slop_len, $seq) = @_;
	my $fasta = get_file_name ($motif, $run_folder, ".fasta");
	# need to figure out the genome file- don't want to make one for each
	my $command = "sort -k1,1 -k2,2n $bed_file | bedtools slop -b $slop_len -i stdin -g $globals->{genome} | bedtools merge -i stdin | bedtools getfasta -fi $globals->{'fasta'} -bed stdin -fo $fasta";
	system $command;
	return $fasta;
}

# returns a threshold count for the most stringent cloud based on the number of loci in a bedfile used to form clouds
sub getMostStringentCoreThreshold { my ($bed_file, $motif) = @_;
	my $rows = `wc -l $bed_file`; my $num_rows = (split(' ', $rows))[0]; 
	if ($num_rows < 100) { 									#need to come up with something better here
		# talk("Motif $motif doesn't have a large enough training set. Skipping\n", $murmur, $globals);
		return 0; 
	} 
	my $most_stringent_core_threshold = $num_rows/$globals->{core_set};
	return $most_stringent_core_threshold;
}

#cloud_kmer_len, shell_factor, 
sub multiThreadCloudBuilding { my ($cloudinfo, $cloud_info_file)= @_;
	my $commands = [];
	my $cloud_kmer_len = $globals->{cloud_kmer_length}; my $shell = $globals->{shell_factor}; my $rev_com_flag = $flags->{rev_com_clouds};
	foreach my $motif_hash ( @$cloudinfo ) {
		my $motif = $motif_hash->{motif};
		my $fasta = $motif_hash->{fasta}; my $cloud_file = $motif_hash->{cloud_file};
		my $max_core = $motif_hash->{max_core}; 
		my @args = ($fasta, $cloud_file, $max_core, $cloud_kmer_len, $shell, $rev_com_flag, $cloud_info_file, $motif);
		# my @args = ($training_fasta, $cloud_file, $num_rows, $cloud_kmer_len, $shell, $rev_com_flag, $cloud_info_file, $motif);
		push (@$commands, "perl SSRCloudBuilder.pl @args" );
	}
	MultithreadCommands($commands); # Chooses a reasonable number of threads
}

sub multiAnnotateAssignFiles { 
	readInKmersFromAssign ( $globals->{'assign'}, $kmers);
	# readx ( $globals->{'assign'}, \&readInKmersFromAssign, "extracting kmers" , 1, $kmers );
	multiAnnotate ( $kmers, $globals->{'fasta'}, $globals->{'out'}, $globals->{'mode'}, $globals->{'merge'} );		
}

# opens fasta for reading and bedfile for writing then passes to sub for annotation
sub multiAnnotate { my ($kmers, $genome, $results_bed, $mode, $merge) = @_;
	my $fp_flag = 1;
	if ($mode eq 'sim') { $fp_flag = 0; }
	else { getFPFromTable ($globals->{'FPTable'}); }
	my $fp_in = openup($genome, $read, "annotating SSRs\n", $globals);
	my $fp_out = openup($results_bed, $write, "recording SSR loci\n", $globals);
	getCloudPos_noMem ($fp_in, $kmers, $fp_out, $fp_flag, $merge );
	close ($fp_in);
	close ($fp_out);
}

sub getFPFromTable { my $infile = shift;
	open (my $fp, '<', $infile) or die "Coudn't open $infile to read in false positive rates.\n";
	while (<$fp>) { chomp;
		my ($strin, $len, $fpr) = split (' ', $_);
		if ( !$FPdist->{$strin}->{'cfd'} ) { $FPdist->{$strin}->{'cfd'} = array_ref(); }
		$FPdist->{$strin}->{'cfd'}->[$len] = $fpr;
	}
	close ($fp);
}

# annotates a sequence
sub getCloudPos_noMem { my ($fp_in, $kmers, $fp_out, $fp_flag, $merge) = @_;
    if (! defined $merge) { $merge = 0; }										#move to globals or default but don't leave here
    my $nextLoc = makeLocFinder($fp_in, $kmers);
    my $prevLoc = $nextLoc->();
    while ($prevLoc->{'begin'} >= 0) {
        my $loc = $nextLoc->();
		# select STDOUT; print Dumper ($prevLoc);
        if ($prevLoc->{'contig'} eq $loc->{'contig'}) {
			if ($prevLoc->{'end'} + $merge >= $loc->{'begin'}) {
                push (@{$prevLoc->{'kmers'}}, $loc->{'kinfo'}->{'fams'}); $loc->{'kmers'} = $prevLoc->{'kmers'}; 
				$loc->{'kinfo'} = $prevLoc->{'kinfo'}; 
				push (@{$prevLoc->{'begin_pos'}}, $loc->{'begin'} ); $loc->{'begin_pos'} = $prevLoc->{'begin_pos'};
				$loc->{'begin'} = $prevLoc->{'begin'};
				push (@{$prevLoc->{'end_pos'}}, $loc->{'end'}); $loc->{'end_pos'} = $prevLoc->{'end_pos'};
				$loc->{'pos'} = $prevLoc->{'pos'};
				if ( $prevLoc->{'end'} > $loc->{'end'} ) { $loc->{'end'} = $prevLoc->{'end'}; }
			} else { printLocInfo ($fp_out, $prevLoc, $fp_flag); }
        } else { printLocInfo ($fp_out, $prevLoc, $fp_flag); }
        $prevLoc = $loc;
    }
}

# returns the next kmer from a sequence
sub makeLocFinder { my ($fp_in, $kmers) = @_;
    my $pos = 0; my $contig = ""; my $buffer = ""; my $i = 0;
	my $len = $globals->{'cloud_kmer_length'}; my $perf_len = $globals->{'oligolength'}; 
	my $ilen = 0; my $line = ""; 
    return sub {
        while (1) {
            for (; $i < $ilen; $i++) { my $kmer = substr ($buffer, $i, $len); my $perf_kmer = substr($buffer, $i, $perf_len);
				if ($kmers->{$kmer} || $kmers->{$perf_kmer}) {
					my $loc;
					if ($kmers->{$perf_kmer}) {
						# print STDERR "Found perf $perf_kmer\t$pos\t$kmers->{$perf_kmer}\n";
						$loc = store_loc_info($loc, $pos, $pos + $perf_len, $kmers->{$perf_kmer}, $contig);
					} elsif ($kmers->{$kmer}) { 
						# print STDERR "Found $kmer\t$pos\t$kmers->{$perf_kmer}\n";
						$loc = store_loc_info($loc, $pos, $pos + $len, $kmers->{$kmer}, $contig);
					}
					$pos++; $i++; 
					return $loc;
				} $pos++; 
            }
            $i = 0; $buffer = substr ($buffer, $ilen);
            if ($line = <$fp_in>) { chomp $line; } else { return { 'begin' => -1, 'end' => -1, 'contig' => "" }; }
            if ($line =~ s/^>//) { $contig = $line; $pos = 0; $buffer = "";
                # talk("Found contig $contig; resetting position\n", $murmur, $globals);
                $line = <$fp_in>; chomp $line;
            }
            $buffer .= uc($line); $ilen = length($buffer) - $len;
        }
    }
}

sub store_loc_info { my ($loc, $begin, $end, $kmer_info, $contig ) = @_;
	$loc = { 'begin' => $begin, 'end' => $end, 'kinfo' => $kmer_info, 'contig' => $contig };
	push (@{$loc->{'begin_pos'}}, $begin); push (@{$loc->{'end_pos'}}, $end);
	push (@{$loc->{'kmers'}}, $kmer_info->{'fams'});
	return $loc;
}

# prints a locus
sub printLocInfo { my $fp_out = shift; my $loc = shift; my $fp_flag = shift;
	# select STDOUT; print Dumper ($loc);
	# my $bests = hash_ref();
	my $loc_strins = $loc->{'kmers'};
	my $locFams = getLocFams ( $loc->{'kmers'}, $loc->{'end'} - $loc->{'begin'}, $globals->{'max_fams'} );
	my $strin_lens = getLongestByStringency ($loc);
	print $fp_out "$loc->{'contig'}\t$loc->{'begin'}\t$loc->{'end'}\t$locFams\t+/-\t";
	if ( $fp_flag == 1 ) {
		my $fps = getLocFP ($strin_lens);
		print $fp_out "$fps\t";
	} else { 
		print $fp_out "N/A\t";  
	}
	print_strin_lens ($fp_out, $strin_lens);
	# exit;
}

# prints the length of the longest consecutive region of each stringency in a locus
sub print_strin_lens { my ($fp_out, $strin_lens) = @_;
	my $max_strin = $globals->{'max_strin'};									
	my @strins;
	foreach my $strin (1 .. $max_strin) {
		push (@strins, $strin_lens->{$strin});
	}
	my $strin_len = join("::", @strins);
	print $fp_out "$strin_len\n";
	# print STDERR "$strin_len\n";
}

# returns the false positive rates for every stringency
sub getLocFP { my ($strin_lens) = @_; my $locFP;
	my $locFPs = array_ref();
	foreach my $strin (sort keys %$strin_lens) {
		my $len_fp = getStrinLenFP ($strin, $strin_lens->{$strin});
		push (@$locFPs, $len_fp);
	}
	return ( join("::", @$locFPs) );
}

# returns the longest consecutive length annotated by each stringency
sub getLongestByStringency { my $loc = shift;
	my $max_strin = $globals->{'max_strin'};										
	my $begins = $loc->{'begin_pos'}; my $ends = $loc->{'end_pos'}; my $locKmers = $loc->{'kmers'};
	my $bests = hash_ref(); #this hash_ref holds on to all the longest consecutive locus lengths for each stringency
	for my $strin (1 .. $max_strin) { $bests->{$strin} = 0; }
	my $temp_locs = hash_ref(); for my $strin (1 .. $max_strin) { $temp_locs->{$strin} = hash_ref(); }
	for (my $i = 0; $i <= $#$begins; $i++) {
		my $begin = $begins->[$i]; my $end = $ends->[$i];
		my $kmerFams = $locKmers->[$i];
		my $best_strin = $max_strin;		# need to set this to highest possible cloud stringency.
		foreach my $fam (keys %$kmerFams) {	# get the best stringency at the position
			if ($kmerFams->{$fam} < $best_strin) { $best_strin = $kmerFams->{$fam};} 
		}
		for (my $strin = $best_strin; $strin <= $max_strin; $strin++) {    # save best stringency at position begin and end
			push(@{$temp_locs->{$strin}->{'begins'}}, $begin); 
			push(@{$temp_locs->{$strin}->{'ends'}}, $end);
		}
	} 
	foreach my $strin ( 1.. $max_strin-1) {								# $max_strin-1 because at $max_strin, length will just be length of locus
		if ($temp_locs->{$strin}->{'begins'}) {
			my $strin_begins =  $temp_locs->{$strin}->{'begins'}; my $strin_ends = $temp_locs->{$strin}->{'ends'};
			my $prev_begin = $strin_begins->[0]; my $prev_end = $strin_ends->[0]; # initializing these values for the loop below
			for (my $i = 1; $i <= $#$strin_begins; $i++) {
				my $begin = $strin_begins->[$i]; my $end = $strin_ends->[$i];
				if ($begin <= $prev_end) { $prev_end = $end; }
				else {
					if ($prev_end-$prev_begin > $bests->{$strin}) { $bests->{$strin} = $prev_end-$prev_begin; }
					$prev_begin = $begin; $prev_end = $end; 
				}
			}
			if ($prev_end-$prev_begin > $bests->{$strin}) { $bests->{$strin} = $prev_end-$prev_begin; }
		}
	} $bests->{$max_strin} = $loc->{'end'} - $loc->{'begin'}; # I made it this way so that merged loci would have strin 5 length as the length of the whole locus
	return $bests;											  # what I find instead is that $loc->{'end'} is actually the end position of the last kmer I found
}

# returns a list of colon separated families from the supplied locus
sub getLocFams { my $loc = shift; my $loc_len = shift; my $max_fams = shift;
	my $fam_counts = {};
	for (my $pos = 0; $pos <= $#$loc; $pos++ ) {
		foreach my $fam (keys %{$loc->[$pos]} ) {
			$fam_counts->{$fam}->{'strin_sum'} += $loc->[$pos]->{$fam}; 
			$fam_counts->{$fam}->{'total'}++;
		}
	}
	foreach my $fam ( sort keys %{$fam_counts} ) {
		my $avg_score = $fam_counts->{$fam}->{'strin_sum'}/$fam_counts->{$fam}->{'total'};
		my $loc_perc = $fam_counts->{$fam}->{'total'}/($loc_len - $globals->{'oligolength'} + 1);
		$fam_counts->{$fam}->{'fam'} = $fam;
		$fam_counts->{$fam}->{'score'} = sprintf( "%.2f", $avg_score/$loc_perc);
		$fam_counts->{$fam}->{'mean_cloud'} = sprintf( "%.2f", $avg_score);
		$fam_counts->{$fam}->{'loc_perc'} = sprintf( "%.2f", $loc_perc);
	}
	my @sort_fams = sort {$fam_counts->{$a}->{'score'} <=> $fam_counts->{$b}->{'score'}} keys %$fam_counts;
	my @fam_info;
	for (my $i = 0; $i < $max_fams && $i <= $#sort_fams; $i++) {
		my $fam = $fam_counts->{$sort_fams[$i]};
		push (@fam_info, $fam->{'fam'}.",".$fam->{'score'}.",".$fam->{'mean_cloud'}.",".$fam->{'loc_perc'});
	}
	return ( join("::", @fam_info) );
}

# get pmf and cmf for false positive lengths
sub getFPFromSim { my ( $FPdist, $FPdata, $simfile, $sim_genome, $output ) = @_;
	my $simBP = getGenomeSize ($sim_genome);
	getSimLengths( $FPdata, 12, $simfile );
	getCumulativeStats ($FPdist, $FPdata->{LDcounts}, $simBP, 12 );							# will need to put the 16 in globals if this works
	open ( my $fp, '>', $output) or die "Couldn't open $output for print false positive rates\n";
	foreach my $strin (sort keys %$FPdist) {
		for ( my $len=0; $len<= $#{$FPdist->{$strin}->{'cfd'}}; $len++ ) {
			my $fpr = getStrinLenFP ( $strin, $len );
			print $fp "$strin\t$len\t$fpr\n";
		}
	}
	close ($fp);
	print STDERR "False positive rates have been printed to $output\n";
}

# returns a false positive rate based on the length and stringency supplied
sub getStrinLenFP { my $strin = shift; my $length = shift;
	no warnings 'recursion';
	my $fp;
	if ($length == 0) { $fp = "N/A"; return $fp; }	
	if ( !$FPdist->{$strin}->{cfd}->[$length]) {
		$fp = getStrinLenFP ( $strin, $length-1);
	} else {  
		$fp = $FPdist->{$strin}->{cfd}->[$length];
	}
	return $fp;
}


sub getGenomeSize { my $genome = shift;
	my $capture = `grep ">" $genome | wc -c`;
	my @contig_names_characters = split(' ', $capture );
	$capture = `grep ">" $genome | wc -l`;
	my @contig_count = split(' ', $capture );
	my $contig_characters = $contig_names_characters[0] - $contig_count[0];
	$capture = `wc -c $genome`;
	my @total_character_count = split(' ', $capture);
	$capture = `wc -l $genome`;
	my @newlines_count = split(' ', $capture);
	my $actual_characters = $total_character_count[0] - $newlines_count[0];
	my $size = $actual_characters - $contig_characters;
	print STDERR "Rates based off genome size of $size\n";
	return $size;
}

sub getSimLengths{ my ( $FPdata, $minLength, $simfile ) = @_;
	open (my $fp, '<', $simfile) or die "Couldn't open $simfile to get false positive rates\n";
	getBedLengthsByStringency( $fp, $FPdata, $minLength );
	close ($fp);
}

sub getBedLengthsByStringency{ my ($fpIN, $FPdata, $minLength) = @_;
	# print STDERR "minLength is $minLength\n";
	my $numContaminants = 0;
	$FPdata->{LDcounts} = hash_ref();
	while (<$fpIN>) { chomp;
        my @data = split('\t', $_);
		my @strin_lens = split("::", $data[6]);
		foreach (my $i = 0; $i <= $#strin_lens; $i++) {
			my $strin = $i + 1;
			my $len = $strin_lens[$i];
			if ($i == $#strin_lens) {
				$len = $data[2] - $data[1]; # for last stringency, use entire length of locus (this includes regions that were merged into the annotation)
			}
			if ( $len >= $minLength || $len == 0) {
				$FPdata->{LDcounts}->{$strin}->{ $len }++;
			} else { $numContaminants++; }
		}
    }
    if ($numContaminants == 0){ $numContaminants = "no"; }
	talk("There were $numContaminants contaminants based on length in simulated bed\n", $murmur, $globals);
}

# problem: this treats lengths that do not appear in the simulated genome and lengths that appear fewer than $min_count (defined herein)
# as if they appeared $min_count times. This does not change much for shorter length annotations but does add in extra false positive bases
# for longer length annotations. As written, this is not the most conservative of false positive rate calculators
sub getCumulativeStats{ my ( $FPdist, $strinLenCounts, $genomeBP, $kmer_len ) = @_;
	foreach my $strin ( keys %{$strinLenCounts} ) {
		$FPdist->{$strin} = hash_ref(); 
		foreach my $len ( keys %{$strinLenCounts->{$strin}} ) { 
			my $Dist = $FPdist->{$strin}; 
			my $Counts = $strinLenCounts->{$strin};
			$Dist->{frac} = array_ref(); $Dist->{cfd} = array_ref(); my ( $sumBP, $CuFreq )= 0;
			my $countArray = array_ref();
			foreach my $cnt ( sort {$b<=>$a} keys %{$Counts} ) { $countArray->[$cnt] = $Counts->{$cnt}; }
			my $longestSim = $#$countArray;
			for (my $oLeng = $longestSim; $oLeng >= $kmer_len;  $oLeng--){ # for each oligo length class
				if ( $Counts->{$oLeng} ) {
					$sumBP = ($Counts->{$oLeng} * $oLeng);#Get
					$Dist->{frac}->[$oLeng] = $sumBP/$genomeBP; # 
					$CuFreq  += $Dist->{frac}->[$oLeng];
					$Dist->{cfd}->[$oLeng] = $CuFreq;
				} 
			}
		}
	}
}

sub readInKmersFromAssign { my ($file, $kmers) = @_;
	open (my $fp, '<', $file) or die "Couldn't open $file\n";
	talk("Opened $file to extract kmers\n", $murmur, $globals);
	while (<$fp>) { chomp;
		my ($kmer, $fam, $stringency, $count) = split(' ', $_);
		if ($stringency <= $globals->{'cloud_strin_cutoff'}) {
			$kmers->{$kmer} = hashifempty ($kmers, $kmer);
			$kmers->{$kmer}->{fams}->{$fam} = $stringency;
		}
	}
	close ($fp);
}

# divides loci into training and test sets used in cloud formation and testing
sub create_training_and_test_sets { my $bed_folder = shift;
	my $test_folder = $bed_folder."_TEST"; my $train_folder = $bed_folder."_TRAIN";
	mkdir $test_folder if not -d $test_folder; mkdir $train_folder if not -d $train_folder;
	my $percent = $globals->{'training_perc'};
	my $bed_count = 0;
	foreach my $motif (@$motifs_of_interest) {
		my $bed_file = get_file_name ($motif, $bed_folder, ".bed");
		if (-s $bed_file) {
			my $train_file = get_file_name ($motif, $train_folder, ".bed");
			my $test_file = get_file_name ($motif, $test_folder, ".bed");
			my ($locs, $loc_count) = get_locs ($bed_file);
			my $num_locs_to_sample = int($loc_count*$percent);
			if ($num_locs_to_sample < 1) { next;}
			my $sampled_locs = sample_wout_replacement ($#$locs + 1, $num_locs_to_sample);
			print_loc_files ($locs, $sampled_locs, $train_file, $test_file);
		} 
	}
}

# memorize all locs in a bed file
sub get_locs {my $infile = shift;
	open (my $fp, '<', $infile) or die "Couldn't open $infile to read\n";
	my $locs = []; my $loc_count = 0;
	while (<$fp>) { chomp;
		my ($chrom, $begin, $end) = split(' ', $_);
		push (@$locs, {'chrom'=>$chrom, 'begin'=>$begin, 'end'=>$end} );
		$loc_count++;
	}
	close ($fp);
	return ($locs, $loc_count);
}

# return a hash with x keys that are randomly sampled indexes 0 to n
sub sample_wout_replacement { my $n = shift; my $x = shift;
	my $samples = {};
	my $num_left = $x;
	while ($num_left > 0) {
		my $number = int(rand($n));
		if (!$samples->{$number}) {
			$samples->{$number} = 1;
			$num_left--;
		}
	}
	return $samples;
}

# prints to two separate files based on index of locus
sub print_loc_files { my ($locs, $samples, $outfile1, $outfile2) = @_;
	open (my $out1, '>', $outfile1) or die "Couldn't open $outfile1 for printing\n";
	open (my $out2, '>', $outfile2) or die "Couldn't open $outfile2 for printing\n";
	for (my $i = 0; $i<= $#$locs; $i++) {
		my $loc = $locs->[$i];
		if ($samples->{$i}) { print_loc ($loc, $$out1); }
		else { print_loc ($loc, $$out2); }
	}
	close ($out1); close ($out2);
}

# prints a locus to file
sub print_loc { my $loc = shift; my $fp = shift;
	print $fp "$loc->{'chrom'}\t$loc->{'begin'}\t$loc->{'end'}\n";
}