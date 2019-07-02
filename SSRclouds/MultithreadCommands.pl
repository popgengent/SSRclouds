# A more complex example of multithreading multiple system calls
# Now there is a single interface and you can simply pass all the commands 
# you want to call and it will do everything for you.  
#
# The MultithreadCommands subroutine will block until all commands have 
# finished running. 
# 
# Only tested on Windows and Ubuntu 12.10 
#
# http://perldoc.perl.org/perlthrtut.html
# http://perldoc.perl.org/threads.html
#
# Stephen Pollard 12/1/2014
# 
# Usage: 
# require "MultithreadCommands.pl"
# 
# # $commands is a ref to an array of command strings 
#
# MultithreadCommands($commands); # Chooses a reasonable number of threads
# OR 
# MultithreadCommands($commands, $number_of_worker_threads);

use strict;
use warnings;

use threads;
use Thread::Queue;


my $command_queue = Thread::Queue->new();
my $failed_command_queue = Thread::Queue->new();


unless(caller) {
	# First generate all the commands to be run
	# This for loop could loop over a parameter changing
	my $commands = [];
	push @$commands, "intentional_fail_command";
	for (0 .. 10) {
		push @$commands, "echo $_";
	}
	push @$commands, "intentional_fail_command";
	
	MultithreadCommands($commands);
	
}

# Total number of threads used = number of worker threads + boss thread
sub MultithreadCommands {
	my ($commands, $number_of_worker_threads) = @_;
	
	$command_queue->enqueue(@$commands);
	
	# Figure out how many threads to start.
	if (not defined $number_of_worker_threads) {
		if (NumberOfProcessors() >= 4) {
			# I used to subtract one more from this number to account for the 
			# current thread, but the current thread is always sleeping while
			# waiting for the other threads, and so it does not take a 
			# processor
			$number_of_worker_threads = NumberOfProcessors() - 2;
		}
		else {
			$number_of_worker_threads = 2;
		}
	}
	print STDERR "$number_of_worker_threads worker threads will be used\n";
	
	my $threads = [];
	# Create all the threads
	for (1 .. $number_of_worker_threads) {
		push @$threads, threads->create(\&worker_thread);
	}
	
	# print STDERR "Done creating threads. They should be running.\n";
	
	# print STDERR "Now I'll wait until all the threads are done using join().\n";
	
	for my $thread (@$threads) {
		# print STDERR "Waiting for thread " . $thread->tid() . " to finish\n"; 
		$thread->join();
	}
	
	if ($failed_command_queue->pending()) {
		warn "These commands failed:\n";
		while (my $failed_command = $failed_command_queue->dequeue_nb()){
			warn $failed_command . "\n";
		}
	}
	
	# print STDERR "You can continue with analysis knowing that all commands are complete\n";
}

# This is only tested for Windows and Ubuntu 12.10
sub NumberOfProcessors {
	if ($^O eq "MSWin32") {
		return $ENV{NUMBER_OF_PROCESSORS}
	}
	else {
		return int(`grep -c processor /proc/cpuinfo`); 
	}
}

sub worker_thread {
	# Keep shifting commands off the queue and running them until there are
	# none left
	# print STDERR "Starting thread " . threads->tid() . "\n"; 
	while (my $command = $command_queue->dequeue_nb()) {
		# print STDERR "Running '$command' on thread " . threads->tid() . "\n";
		my $failed = system $command;
		if ($failed) {
			# Record any failed commands
			$failed_command_queue->enqueue($command);
		}
	}
	# print STDERR "Thread " . threads->tid() . " finished processing commands\n"
}