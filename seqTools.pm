#Perl -Sx "{0}"; Exit {Status}
#!perl

package seqTools;
use strict;
use warnings;

require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '0.01';
@ISA         = qw(Exporter);
@EXPORT_OK   = qw(readFasta makeSeqHash);

my $silent = 0; my $quiet = 1; my $murmur = 2; my $talk = 3; my $loud = 5;
my $read = "<"; my $write = ">"; my $append = ">>"; my $read_append = "+>>"; my $no_record = 1;

sub readFasta { my $target = shift; my $sub = shift; my $blurb = shift; my @args = @_; my $name = ""; my $seq = ""; my $fp;
    if (eval {fileno $target}) { $fp = $target; }
    else { $fp = main::openup($target, $read, "fasta file for ".$blurb); }
    while (<$fp>) { if ($_ =~ s/>//) { chomp; $name = $_; }; last; } # get first sequence name
    while (<$fp>) { chomp;
        if ($_ =~ s/>//) { $sub->($name, $seq, @args); $name = $_; $seq = ""; }
        else { $seq .= uc($_); }
    }
    $sub->($name, $seq, @args); $name = $_; $seq = ""; # run sub on last sequence
    close $fp;
}

# This seems like it has some memory liabilities but it is needed for how we currently build clouds.
sub makeSeqHash { my ($name, $seq, $seqs, $ptrs) = @_;
    my $seqinfo = {}; $seqinfo->{name} = $name;
    $seqs->{$seqinfo.""} = $seq; $ptrs->{$seqinfo.""} = $seqinfo;
}

sub revComp { my $newseq = reverse scalar shift; $newseq =~ tr/ACGTacgt/TGCAtgca/; return $newseq; }

1;