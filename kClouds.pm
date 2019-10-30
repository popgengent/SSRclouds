#!/usr/local/bin/perl
#
#
#	kClouds.pm
#
# David Pollock, Jaime Merlano & Corey Cox, copyright 30-Apr-2012
#
# kClouds.pm is meant to handle cloud structures for other program processing
#

package kClouds;
use strict;
use warnings;

require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '1.11';
@ISA         = qw(Exporter);
@EXPORT_OK   = qw(buildPClouds importClouds exportPClouds getCloudKmers revCompKmers SSR2Cloud);

use constant KEYNUCS => ('A', 'C', 'G', 'T');

## Exported Subs ##
sub buildPClouds { my ($seqs, $kseqs, $clouds, $kmers, $core, $outer, $klen, $flags) = @_;
    getSeqKmers($seqs, $kseqs, $kmers, $klen, $flags);
    my $train_size = 0;  ### for now we are using everything to build clouds... pass $train_size later.
    my ($train_set, $test_set) = getTrainSeqs($kseqs, $train_size);
    getKmerInfo($kseqs, $kmers, $train_set);
    makePClouds($kmers, $clouds, $core, $outer, $klen, $flags);
}

sub importClouds { my ($fp, @the_rest) = @_;
    my $importSub = getImportSub(peek($fp)); $importSub->(@_);
}

sub exportPClouds { my $fp = shift; my $kmers = shift;
    foreach my $kmer (sortBySubKey($kmers, 'count')) { my $kinfo = $kmers->{$kmer};
        if ($kinfo->{cloud}) { print $fp "$kmer\t$kinfo->{cloud}->{id}\t$kinfo->{count}\n"; }
    }
}

sub getCloudKmers { my $clouds = shift; my $cloud_kmers = {};
    while (my ($seed, $cloudinfo) = each %$clouds) {
        foreach my $kinfo (@{$cloudinfo->{members}}) { $cloud_kmers->{$kinfo->{seq}} = $kinfo; }
    } return $cloud_kmers;
}

sub revCompKmers { my $kmers = shift; my @keys = keys %{$kmers};
    foreach my $kmer (@keys) { my $rev_kmer = revComp($kmer);
        if (!$kmers->{$rev_kmer}) { $kmers->{$rev_kmer} = $kmers->{$kmer}; }
    }
}

sub SSR2Cloud { my ($motifs, $IDs, $clouds, $kmers, $flags,  $motifmin,  $motifmax) = @_; my $seeds = {};
	foreach my $ID (@{$IDs}) {
		my $checklen = length $ID;
		if ($checklen >= $motifmin && $checklen <= $motifmax) {
			foreach my $kmer (keys %{$motifs->{$ID}->{oligofams}}) {	   
				my $count = 900; #just to give a value to $count; may be able to delete it
				my $kinfo = addKmer($kmers, $kmer, $flags, $count);
				if (!$seeds->{$ID}) { $seeds->{$ID} = $kmer;
				seedCloud($clouds, $kmer, $kinfo, $flags, $ID);
				}
				else { addToCloud($clouds, $seeds->{$ID}, $kinfo); }
			}
		}
	}
}

## Internal Subs ##
sub getSeqKmers{  my ($seqs, $kseqs, $kmers, $klen, $flags, @loci) = @_;
	foreach my $seq (keys %{$seqs}){
        for (my $i = 0; ($i + $klen) < length($seqs->{$seq}); $i++){
            my $kmer = substr($seqs->{$seq}, $i, $klen);
            if ($kmer =~ /[ACGTacgt]{$klen}/) {  ## Filter out kmers that have non-canonical bases
                addKmer($kmers, $kmer, $flags, 0, $kseqs, $seq, $i, @loci);
            }
        }
    }
}

sub addKmer{ my ($kmers, $kmer, $flags, $count, $kseqs, $seqid, $pos, @loci) = @_;  ### note order of arguments drastically changed!!!
    my $kinfo = getKeyHash($kmers, $kmer); $kinfo->{count} = $count;
    if (!$kinfo->{seq}) { $kinfo->{seq} = $kmer; }
    if (checkFlag($flags, 'kmerbits')) { splitifempty($kinfo, "bits", $kmer); }
    if (checkFlag($flags, 'recordpos')) { if ($kseqs) { addPosToKseq ($kseqs, $seqid, $kinfo, $pos); } }
#    if (@loci) { push(@{getKeyArray($kinfo, 'loci')}, @loci); }  ### keeping this for if we want to import/export sequence loci
    return $kinfo;
}

sub getTrainSeqs { my $seqs = shift; my $train_size = shift; my $training = []; my $testset = [];
    if (!$train_size){ @$training = keys %$seqs; @$testset = my @empty; } # default: use them all
    else { my $count = 0; my @source = keys %$seqs; my $numseqs = scalar @source;
        if ($numseqs < $train_size) { die "Training set size ($train_size) is greater than number of sequences ($numseqs).\n"; }
        for (my $i = $train_size; $i > 0; $i--) { my $targetID = rand($numseqs); $numseqs--;
            my $target = splice(@source, $targetID, 1); push (@{$training},$target);
        } @$testset = @source;
    } return ($training, $testset)
}

sub getKmerInfo {  my ($kseqs, $kmers, $train) = @_;
    foreach my $seqid (@$train){
        for (my $pos = 0 ; $kseqs->{$seqid}->[$pos]; $pos++){
            my $kinfo = $kseqs->{$seqid}->[$pos]; $kinfo->{count}++;
            addLocus($kinfo, $seqid, $pos);  ### locus information should also be optional...
        }
    }
}

### Subs for building pclouds   ###
sub makePClouds{ my ($kmers, $clouds, $core, $outer, $kmerlen, $flags) = @_;
	foreach my $kmer (sortBySubKey($kmers, 'count')){  #sort by kmer count
        my $kinfo = getKeyHash($kmers, $kmer);
        if (eligible($kinfo, $core)){ my $cloud = $kmer;
            seedCloud($clouds, $cloud, $kinfo, $flags);
            if (checkFlag($flags, 'old_himp_clouds')) { expandCloud_himp_old($kmers, $clouds, $kmer, $core, $outer, $kmerlen); }
            else { expandCloud($kmers, $clouds, $kmer, $core, $outer, $kmerlen, $flags); }
            $clouds->{$cloud}->{consensus} = getConsensus($kmerlen, $clouds->{$cloud}->{members});
        }
    }
}

sub seedCloud { my ($clouds, $seed, $kinfo, $flags, $cloud_id) = @_;
    my $cloudinfo = getKeyHash($clouds, $seed); $cloudinfo->{kmer} =  $seed;
    if (checkFlag($flags, 'cloudid')) { if ($cloud_id) { $cloudinfo->{id} = $cloud_id; } else { $cloudinfo->{id} = scalar keys %$clouds; } }
    if (checkFlag($flags, 'cloudbits')) { splitifempty($cloudinfo, "bits", $seed); }
    addToCloud($clouds, $seed, $kinfo);
}

sub addToCloud { my ($clouds, $cloud, $kinfo) = @_;
    my $cloudinfo = $kinfo->{cloud} = getKeyHash($clouds, $cloud);
    $cloudinfo->{count}++;   ### number of kmers in cloud (optional?,relabel???)  --- count is poor description of this...
    push (@{getKeyArray($cloudinfo, "members")}, $kinfo);
}

sub expandCloud { my ($kmers, $clouds, $seedmer, $core, $outer, $kmerlen, $flags) = @_;
	my $cloudinfo = $kmers->{$seedmer}->{cloud};
	for (my $i = 0; $i < $kmerlen; $i++) {
		foreach my $nuc (KEYNUCS) {
			my $testmer = $seedmer; substr($testmer, $i, 1) = $nuc; 
			if ($kmers->{$testmer}) { my $kinfo = $kmers->{$testmer};
				if ( eligible($kinfo, $outer) ) {
					if ( eligible($kinfo, $core) ) { # addtoclouds needs to stay nested here for eligible
					addToCloud($clouds, $cloudinfo->{kmer}, $kinfo);  # makes current kmer no longer eligible (prevents loop)
					expandCloud($kmers, $clouds, $testmer, $core, $outer, $kmerlen, $flags);
					} else { addToCloud($clouds, $cloudinfo->{kmer}, $kinfo); }
				}
			}
			if (checkFlag($flags, 'rev_com_clouds')) {my $rc_testmer = revComp($testmer);
				if ($kmers->{$rc_testmer}) { my $kinfo = $kmers->{$rc_testmer};
					if ( eligible($kinfo, $outer) ) {
						if ( eligible($kinfo, $core) ) { # addtoclouds needs to stay nested here for eligible
						addToCloud($clouds, $cloudinfo->{kmer}, $kinfo);  # makes current kmer no longer eligible (prevents loop)
						expandCloud($kmers, $clouds, $rc_testmer, $core, $outer, $kmerlen, $flags);
						} else { addToCloud($clouds, $cloudinfo->{kmer}, $kinfo); }
					}
				}
			}
		}
	}
}

sub eligible{ my ($testinfo, $criterion) = @_;
    if ($testinfo && ! $testinfo->{cloud} && $testinfo->{count} >= $criterion){ return 1; } else { return 0; }
}

sub expandCloud_himp_old { my ($kmers, $clouds, $seedmer, $core, $outer, $kmerlen) = @_;
    my $cloudinfo = $kmers->{$seedmer}->{cloud};
    for(my $i=0; $i < $kmerlen; $i++){
        foreach my $nuc (KEYNUCS) {
            my $testmer = $seedmer; substr($testmer, $i, 1) = $nuc;
            if ($kmers->{$testmer}) { my $kinfo = $kmers->{$testmer};
				if ( eligible($kinfo, $outer) ) { addToCloud($clouds, $cloudinfo->{kmer}, $kinfo);
					if ( eligible($kinfo, $core) ) { expandCloud_himp_old($kmers, $clouds, $testmer, $core, $outer, $kmerlen); }
				}
			}
        }
    }
}

### Subs for importing pclouds from file to hash structure.   ###
sub getImportSub { my $num_fields = scalar goodSplit(shift); my $sub;
    if (3 == $num_fields) { $sub = \&importPClouds; }
	if (4 == $num_fields) { $sub = \&importSSRPClouds; }
    elsif ($num_fields > 5) { $sub = \&importHimPClouds; }
    else { print STDERR "\n####\n\n Cloud Format not recognized\n\n####\n"; die; }
    return $sub;
}

sub importPClouds { my $fp = shift; my $clouds = shift; my $kmers = shift; my $flags = shift; my $seeds = {};
    while (<$fp>) { chomp;
        my ($kmer, $cloud_id, $count) = goodSplit($_);
        my $kinfo = addKmer($kmers, $kmer, $flags, $count);
        if (!$seeds->{$cloud_id}) { $seeds->{$cloud_id} = $kmer;
            seedCloud($clouds, $kmer, $kinfo, $flags, $cloud_id);
        }
        else { addToCloud($clouds, $seeds->{$cloud_id}, $kinfo); }
    }
}

sub importSSRPClouds { my $fp = shift; my $clouds = shift; my $kmers = shift; my $flags = shift; my $seeds = {};
	while (<$fp>) { chomp;
		my ($kmer, $fam, $stringency, $count) = goodSplit($_);
	}
}

sub importHimPClouds { my $fp = shift; my $clouds = shift; my $kmers = shift;
    while (<$fp>) { chomp; if ($_ eq "") { next; }
        my ($seed, $num_members, $oligos, $anchor, $supercloud, @members) = split(/\t/, $_);
        ### need to rewrite to use current seedCloud and addToCloud code...
    }
}

### subs for optional cloud features   ###
sub checkFlag { my ($flags, $flag) = @_; if ($flags && $flags->{$flag}) { return 1; } else { return 0; } }

sub getConsensus{ my ($klen, $members) = @_;
    my $nuccounts = hash_ref(); my $consensus; # removed array_ref()
    for (my $i = 0; $i < $klen; $i++) { my $maxcount = 0; my $maxnuc = "";
        foreach my $kinfo (@{$members}){ $nuccounts->{$kinfo->{bits}->[$i]}++; }
        foreach my $nuc (keys %$nuccounts){
            if ($nuccounts->{$nuc} > $maxcount){ $maxcount = $nuccounts->{$nuc}; $maxnuc = $nuc;}
            $nuccounts->{$nuc} = 0;
        }  $consensus .= $maxnuc;
    }
    return $consensus;
}

sub addLocus { my $info = shift; my $seqid = shift; my $pos = shift;
    push (@{getKeyArray($info, "loci")}, { 'seqid' => $seqid, 'position' => $pos});
}

sub addPosToKseq { my $kseqs = shift; my $seqid = shift; my $kinfo = shift; my $pos = shift;
    my $kseqinfo = getKeyArray($kseqs, $seqid); $kseqinfo->[$pos] = $kinfo;
}

sub checkBits { my ($cloudinfo, $cloud) = @_;  ### this seems like it could be unnecessary in general... use splitifempty instead.
    if (!$cloudinfo->{kmer}) { $cloudinfo->{kmer} = $cloud; }
    splitifempty($cloudinfo, "bits", $cloud);
}

sub splitifempty{ my ($inhash, $inkey, $kmer) = @_;
    if (!$inhash->{$inkey}) { $inhash->{$inkey} = array_ref();
        for(my $i = 0; $i < length $kmer; $i++){ $inhash->{$inkey}->[$i] = substr($kmer, $i, 1); }
    }
}

# seperate on whitespace (/\s+/) but prevents empty record on leading whitespace
sub goodSplit { my @chunks = split(' ', shift); return @chunks; }
sub peek { my $fp = shift; my $pos = tell($fp); my $line = <$fp>; seek($fp, $pos, 0); return $line; }
sub revComp { my $newseq = reverse scalar shift; $newseq =~ tr/ACGTacgt/TGCAtgca/; return $newseq; }

# pass array pointer, check existence, if not, create anonymous
sub getKeyLen { my $hash = shift; my @keys = keys %{$hash}; return length($keys[0]); }
sub getKeyArray { my ($hash, $key) = @_; if (!$hash->{$key}) { $hash->{$key} = array_ref(); } return $hash->{$key}; }
sub getKeyHash { my ($hash, $key) = @_; if (!$hash->{$key}) { $hash->{$key} = hash_ref(); } return $hash->{$key}; }
sub getIndexHash { my ($array, $index) = @_; if (!$array->[$index]) { $array->[$index] = hash_ref(); return $array->[$index]; } }
sub sortBySubKey { my $hash = shift; my $key = shift; return sort { $hash->{$b}->{$key} <=> $hash->{$a}->{$key} } keys(%$hash); }

sub hash_ref{ my %hash; return \%hash; }
sub array_ref{ my @array; return \@array; }

### sub getHashAtKey { my ($hash, $key) = @_; if (!$hash->{$key}) { $hash->{$key} = hash_ref(); } return $hash->{$key}; }
### sub getArrayAtKey { my ($hash, $key) = @_; if (!$hash->{$key}) { $hash->{$key} = array_ref(); } return $hash->{$key}; }
### sub getHashAtIndex { my ($array, $index) = @_; if (!$array->[$index]) { $array->[$index] = hash_ref(); return $array->[$index]; } }

sub print_value { my $out = shift; my $value = shift; my $indent = shift; my $ptrs = shift;
    if (!$ptrs){ $ptrs = {}; } elsif ($ptrs->{$value}) { return; }
    my $val_ref = ref $value;
    if ($val_ref eq "") { print $out "$indent$value\n"; }
    elsif ($val_ref eq  "ARRAY") { $ptrs->{$value} = 1; print_array($out, $value, $indent." ", $ptrs); }
    elsif ($val_ref eq  "HASH") { $ptrs->{$value} = 1; print_hash($out, $value, $indent." ", $ptrs); }
}

sub print_hash { my $out = shift; my $hash = shift; my $indent = shift; my $ptrs = shift;
    while (my ($key, $value) = each %{$hash}) { print $out "$indent$key:\n"; print_value($out, $value, $indent." ", $ptrs); }
}

sub print_array { my $out = shift; my $array = shift; my $indent = shift; my $ptrs = shift;
    for (my $i = 0; $i <= $#{$array}; $i++) { print_value($out, $array->[$i], $indent, $ptrs); }
}


1;
