#!/usr/local/bin/perl
#
#
#	Globals.pm
#
# Part of PerlTrees, under construction 
# David Pollock & Corey Cox, copyright 12/1/05
#
# globals.pm is meant to read in factory and settings files
#

package globals;

use strict;

require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '0.01';
@ISA         = qw(Exporter);
@EXPORT_OK   = qw(ProgSetUp finish);

# main memory
my $append = "append";
my $invariant = "*"; my $variant = "+"; my $numsym = "#";
my $silent = 0; my $quiet = 1; my $murmur = 2; my $talk = 3; my $loud = 5;

my $globals = {}; my $program = {}; my $presets = {};

#### globals.pm could use a significant rewrite.
#### 1) Remove all bareword filehandles
####   - these are essentially depricated in new Perl because the don't play well in namespaces/scope
#### 2) Use openup, talk and record subs (from skel)
####   - unify the way we report to the user and record to file
#### 3) Change to an object to better handle repeated uses of variables
####   - this will significantly simplify the code for openup, talk and record subs in (2)

## Exported Subs ##
sub ProgSetUp{
	my $factoryfile = shift; $globals = shift; $program = shift;
	$globals->{starttime} = time();
	ReadFactory($factoryfile, $presets, $globals);
	ReadControl($globals->{default}, $globals, $presets);
	ReadControl($globals->{control}, $globals, $presets);
	
	srand($globals->{seed});
	SetRecord($globals->{starttime});
}

sub finish { my $starttime = shift;
	my $now = time(); my $delta = $now - $starttime ;
	warn "I'm leaving you now...$delta seconds elapsed.\n";
	exit(0);
}

## Internal Subs ##
sub ReadFactory{
	my $infile = shift; my $settings = shift; my $globals = shift;
	my %sethash; my $setptr; my $rangeptr;
	
	open(FACTORY,"<$infile");
	while(<FACTORY>) {
		chomp($_);
		my ($info, $comment) = split(/$numsym+/, $_);
		my ($key, $type, $constant, $example, @range) = split(/[\s=]+/, $info);
		$settings->{$key} = {};
		if ($key ne " " | $key) {
			$setptr = $settings->{$key};
			$setptr->{type} = $type;
			$setptr->{constant} = $constant;
			$setptr->{example} = $example;
			$setptr->{range} = [];
			my $rangeptr = $setptr->{range};
			push(@{$rangeptr},@range);
			if($constant eq $invariant){
				$globals->{$key} = $example;
			} else { $globals->{$key} = $example; } # note: currently sets globals to factory
		}
	}
	close (FACTORY);
	print "Factory File Read.\n";
}

sub ReadControl{
	my $infile = shift; my $controls = shift; my $presets = shift;
	
	print STDERR "\nReading control file ** $infile **\n";
	
	open(CONTROL,"<$infile");
	while(<CONTROL>) {
		chomp($_);
		my ($info, $comments) = split(/[#]+/,$_);
		$info =~ s/\t//g;
		# if($info){ printf (STDERR "\tread control item: %-20s",($info));if($comments){print STDERR "$comments";}print STDERR "\n";}
		my ($key, $value, @vallist) = split(/[\s:=]+/, $info);
		if (@vallist) {
			my $tempval = $value;
			$value = [];
			print STDERR "\t\tLIST! Key: $key; values: 0=>$tempval";
			$value->[0]=$tempval; my $i=1; foreach my $val (@vallist) {$value->[$i]=$val;print STDERR ", $i=>$val"; $i++;} print STDERR ".\n"; 
		}
		my $thisvar;
		if ($presets && $presets->{$key}){ $thisvar = $presets->{$key}; }
		if ($presets && $presets->{$key} && ($thisvar->{constant} eq $invariant)){
			print "Cannot change $key from control file; reverting to factory setting ($thisvar->{example}).\n";
		} else {
			$controls->{$key} = $value;
		}
	}
	print "Control File $infile Has Been Read.\n";
}

sub SetRecord{
	my $starttime = shift;
    openoutselect($globals->{record_run}, "Recording program run", $append); 
	updateprogram(); displayprogram();
	print "Run at $starttime seconds\n";
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($globals->{starttime});
	my @mos = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my @days = qw( Sun Mon Tue Wed Thu Fri Sat );
	my $month = @mos[$mon]; $year += 1900; my $wkday = @days[$wday];
	print "$wkday, $month $mday, $year. $hour:$min:$sec.\n";
	FactoryOut($presets);
	GlobalsOut($globals);
    close();
}

sub printcontrolhash{
	my $inhashptr = shift; my $blurb = shift;
	my ($key, $value); my $range = "range";
	my $invariant = "*";
	my $variant = "+";
	
	print "$blurb\n";
	foreach $key (keys %{$inhashptr}){
		my $thisvar;
		if ($presets && $presets->{$key}){
			$thisvar = $presets->{$key};
		}
		if ($presets && $presets->{$key} && ($thisvar->{constant} eq $invariant)){
			;
		} else {
			#print "Variable $thisvar->{constant}.\n";
			$value = $inhashptr->{$key};
			if ($key eq $range || $key eq "rates" || $key eq "freqs" || $key eq "mcounts") {
				print "\t$key = ",join (', ', @{$value}), "\n";
			}
			else {  print "\t$key = $value\n"; }
		}
	}
}

sub openoutselect{
	my $outfile=shift; my $blurb=shift; my $do_append = shift;
	my $defaultblurb = "unknown items";
	if (! $blurb) { $blurb = $defaultblurb; }
	if ($do_append eq $append){
		open (fp, ">>$outfile") or warn "Couldn't open $outfile to append.\n";
	} else {
		open (fp, ">$outfile") or warn "Couldn't open $outfile for output.\n";
	}
	print STDERR "printing $blurb to file $outfile.\n";
	select fp;
	return \*fp;
}

sub updateprogram{
	foreach my $key (keys %{$program}) {
		if ($globals->{$key}) { $program->{$key} = $globals->{$key}; }
	}
}

sub displayprogram{
	my @keys =  sort keys(%{$program});
	print "\Program $program->{name}\n";
	foreach my $key (@keys) {
		if ($program->{$key})
		{print "\t$key => $program->{$key}\n";}
	}
}

# Do any extra global initialization necessary, display global variables
# later, this should include reading in the trees
sub InitGlobals {
	my $globs = shift; my $test = shift;
	my $outfile = $test->{testfile};
	open (OUTFILE, ">$outfile") or warn "Couldn't open $outfile for output.\n";
	select OUTFILE;
	display($globs);
	print $globs->{infile},"\n";
}

# Output globals; later this should be nicer output, more specific. Might descriminate
# types of globals (what is default, what was changed, what is factory)
# also, will want output dependent on the "loudness" variable
sub GlobalsOut{
	my $globs = shift; my $blurb=shift;
	if ($globals->{loudness} > $globals->{silent}){  
		if ($blurb) { printcontrolhash($globs, $blurb); } 
		else {printcontrolhash($globs, "\nGlobal Variable settings:"); }
	}
}

sub FactoryOut{
	my $sets = shift; my $key;
	if ($globals->{loudness} > $globals->{quiet} ){  
		print "Factory Set Variables: ", join(', ', (keys %{$sets})), "\n";
		foreach $key (keys %{$sets}){
			if ($key) { printcontrolhash($sets->{$key}, "Information for global constant $key:")};
		}
	} else { print "Factory settings not being printed.\n"; }
}

1;