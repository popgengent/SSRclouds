#Perl -Sx "{0}"; Exit {Status}
#!perl

#This script runs a post-SSR annotation analysis to determine the FDR of each annotated locus.

#5 will have FPRs after splitting. 6 will have lengths of annotations at different stringencies.

#Need to calculate fdrs on cumulative counts, not fractional- fpr is already cumulative
# Also need SSR-finder to give fprs for loci under 16 bp
use strict;
use warnings;

use File::Glob ':glob';
use POSIX;
#use Math::Random;
use Data::Dumper;

my $fprs = {};
my $strin_len_counts = {};
my $fdrs = {};

if ($#ARGV != 1) {die "Usage:\nperl get_fdr.pl [bed_file] [genome]\n"; }

my $infile = $ARGV[0];
my $genome = $ARGV[1];

my $genome_size = getGenomeSize($genome);

# read in annotations- store stringency and length dependent false positive rates
# also add up total length of annotations at each stringency and length
# also creates and expected number of false positves based on genome size and false positive rate for each stringency and length
open (my $infp, '<', $infile) or die "Couldn't open $infile for reading\n";
my $count = 0;
while (<$infp>) { chomp;
	my ($chrom, $begin, $end, $fam, $strand, $combined_fps, $combined_lens) = split (' ', $_);
	my @fps = split('::', $combined_fps);
	my @lens = split('::', $combined_lens);
	for (my $i = 0; $i <= $#lens; $i++) {
		my $strin = $i +1;
		my $len = $lens[$i];
		my $fpr = $fps[$i];
		if ($len >0) {
			if (!$fprs->{$strin}->{$len}->{'fpr'}) {
				$fprs->{$strin}->{$len}->{'fpr'} = $fpr;
				$fprs->{$strin}->{$len}->{'exp'} = $fpr*$genome_size;
			} elsif ($fprs->{$strin}->{$len}->{'fpr'} eq 'N/A') {
			} elsif ($fpr != $fprs->{$strin}->{$len}->{'fpr'}){
				die "Error! Fpr at strin $strin and len $len should be $fprs->{$strin}->{$len} not $fpr\n";
			}
			$strin_len_counts->{$strin}->{$len}->{'frac'} += $len;
			# print STDERR "len: $len\tmx: $strin_len_counts->{$strin}->{'max'}\tmn:$strin_len_counts->{$strin}->{'min'}\n";
			# if (!$len || !$strin_len_counts->{$strin}->{'min'} || !$strin_len_counts->{$strin}->{'min'}) {
				# $count++;
				# print STDERR "len: $len\tmx: $strin_len_counts->{$strin}->{'max'}\tmn:$strin_len_counts->{$strin}->{'min'}\n";
				# if ($count > 15) { exit; }
			# }
			if (!$strin_len_counts->{$strin}->{'max'} || $len > $strin_len_counts->{$strin}->{'max'} ) {
				$strin_len_counts->{$strin}->{'max'} = $len;
			}
			if (!$strin_len_counts->{$strin}->{'min'} || $len < $strin_len_counts->{$strin}->{'min'} ) {
				$strin_len_counts->{$strin}->{'min'} = $len;
			}
		}
	}
}
close($infp);

# calculates a false discovery rate using expected number of false positive annotations and cumulative length of annotations for each stringency and length
for my $strin (1..5) {
	my $min = $strin_len_counts->{$strin}->{'min'};
	my $max = $strin_len_counts->{$strin}->{'max'};
	my $total = 0;
	for (my $i = $max; $i >= $min; $i--) {
		my $len = $i;
		if (!$fprs->{$strin}->{$len}) { next; }
		if ($strin_len_counts->{$strin}->{$len}->{'frac'}) {
			$total += $strin_len_counts->{$strin}->{$len}->{'frac'};
		}
		# $strin_len_counts->{$strin}->{$len}->{'cumul'} = $total;
		$fdrs->{$strin}->{$len} = $fprs->{$strin}->{$len}->{'exp'}/$total;
		# print STDERR "$strin\t$len\t$total\t$fprs->{$strin}->{$len}->{'fpr'}\t$fprs->{$strin}->{$len}->{'exp'}\t$fdrs->{$strin}->{$len}\n";
	}
}
# print Dumper($fdrs); exit;

# finds the lowest fdr from the available stringencies and lengths of each annotation and prints to file
open (my $outfp, '>', $infile."_cumul_fdr.bed") or die "Couldn't open outfile for printing\n";
open ($infp, '<', $infile) or die "Couldn't open $infile for fdr calculations and re-printing\n";
while (<$infp>) { chomp;
	my ($chrom, $begin, $end, $fam, $strand, $combined_fps, $combined_lens) = split (' ', $_);
	my @lens = split('::', $combined_lens);
	my $loc_fdrs = [];
	for (my $i = 0; $i <= $#lens; $i++) {
		my $strin = $i +1;
		my $len = $lens[$i];
		if ($len >0) {
			# print STDERR "$strin\t$len\t$chrom\t$begin\t$end\t$strand\t$combined_lens\n";
			# print Dumper ($fdrs->{$strin});
			# print STDERR "$fdrs->{$strin}->{$len}\n";
			# my $fdr = $fprs->{$strin}->{$len}->{'exp'}/($strin_len_counts->{$strin}->{$len}->{'cumul'});
			my $fdr = $fdrs->{$strin}->{$len};
			push (@$loc_fdrs, $fdr);
		}
	}
	my $fdr = getLowest($loc_fdrs);
	print $outfp "$chrom\t$begin\t$end\t$fdr\t$fam\t$strand\t$combined_lens\n";
}
close($infp);
close ($outfp);

exit;

sub getGenomeSize { my $genome = shift;
	my $size_capture = `wc -c $genome`;
	my @fields = split(' ', $size_capture);
	my $size = $fields[0];
	print STDERR "Genome size is $size\n";
	return $size;
}

sub getLowest { my $rates = shift;
	my $lowest = 1;
	foreach my $rate (@$rates) {
		if( $rate < $lowest) {$lowest = $rate; }
	}
	return $lowest;
}
