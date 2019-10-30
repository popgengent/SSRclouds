
#!/usr/local/bin/perl
#
#
#

package basics;
use strict;
use warnings;
use Data::Dumper;

require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '0.01';
@ISA         = qw(Exporter);
@EXPORT_OK   = qw(readx printx openup hashifempty scalar_ref array_ref hash_ref skip talk record);

#main memory
my $silent = 0; my $quiet = 1; my $murmur = 2; my $talk = 3; my $loud = 5; 
my $read = "<"; my $write = ">"; my $append = ">>"; my $read_append = "+>>"; my $no_record = 1;

### sub-routines that call other subroutines should go here ###
sub readx { my ($file,  $sub, $blurb, $no_record, @args) = @_;
    if (my $fp = openup($file, $read, $blurb, $no_record)) { $sub->($fp, @args); close $fp; }
    else { warn "Unable to read input file with readx!\n\n Aborting!!\n\n"; return -1; }
}

sub printx {
	my $flag = shift; my $sub_ref = shift; my $target=shift; my $blurb = shift; my @args = @_;
	talk("Printing Generic $blurb ($flag, $sub_ref, $target, $blurb, @_, etc.)\n\n", $murmur);
	my $file_handle; my $handle;
	if ($flag > 0) {
		if (!$target) { $handle = \*STDERR; }
        if (!$blurb) { $blurb = "unknown information"; }
		if ($target eq \*OUTFILE || $target eq \*STDERR || $target eq \*STDOUT) {
            talk("standard output ID: file pointer is $target\n", $murmur); $handle = $target;
        }
		else {
            $file_handle = openup($target, $write, "for printing via printx");
            talk("output in form of file name:  $target\n\n", $murmur);  $handle = $file_handle;
        }
		if ($handle) { $sub_ref->($handle, @args);  return; }# call print function and exit sub
	}
    if ($file_handle) { close($file_handle); }
    talk("\nNothing printed for $blurb, flag was off ($flag) or can't open file ($target).\n", $murmur);
}

##
##		Basic subroutines needed for many apps
##

sub openup { my $file = shift; my $mode = shift; my $blurb = shift; my $globals = shift; my $no_record = shift; 
	if (! $blurb) { $blurb = "unknown"; }
    open (my $file_handle, $mode, $file) or warn "Couldn't open file: $file. \n";
    if (!$no_record) { talk("Opened $file as $blurb", $murmur, $globals); }
	return $file_handle;
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

sub talk { my $phrase = shift; my $loudness = shift; my $globals = shift;
	if ($loudness <= $globals->{record_loudness}) { record($phrase, $globals); }
	if ($loudness >= $globals->{loudness}) {  print STDOUT $phrase; }
}

sub record { my $phrase = shift; my $globals = shift;
    my $recordfile = openup($globals->{'record_run'}, $append, "", $globals, $no_record);
    my $now = time(); my $delta = $now - $globals->{starttime};
	print $recordfile "$delta: $phrase"; close($recordfile);
}

1;