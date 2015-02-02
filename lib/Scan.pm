package Scan;
use strict;

##################################################################
#
# Scan
#
# NOTES:
#
# SYNOPSIS:
#
# TODO:
# -
# -
#
# CHANGELOG:
# 2008.02.19
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $scan = Scan->new( $scan_string );
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ($scan_string) = @_;

	if ( defined($scan_string) ) {
		if ( $scan_string =~ /^.+\.\d+\.\d+\.\d+$/ ) {
			$self->{scan_string} = $scan_string;
		}
		else {
			print "  scan string: '$scan_string' format not allowed!\n";
		}
	}
	return $self;
}

# Object Method
# Title	    :
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_scan {
	my $self = shift();

	if ( exists( $self->{scan_string} ) ) {
		if ( defined( $self->{scan_string} ) ) {
			return $self->{scan_string};
		}
	}
	return '';
}

sub get_run {
	my $self            = shift;
	my ($optional_scan) = @_;
	my $scan            = '';
	if ( defined($optional_scan) ) {
		$scan = $optional_scan;
	}
	else {
		$scan = $self->get_scan();
	}

	# 2006-10-14_03_s4.0465.0465.2
	if ( $scan =~ /^(.+)\.\d+\.\d+\.\d+$/ ) {
		$scan =~ /^(.+)\.\d+\.\d+\.\d+$/;
		return $1;
	}
	else {
		n_print("run not defined!\n");
		return '-';
	}
}

sub get_scan_nr1 {
	my $self            = shift;
	my ($optional_scan) = @_;
	my $scan            = '';
	if ( defined($optional_scan) ) {
		$scan = $optional_scan;
	}
	else {
		$scan = $self->get_scan();
	}

	# 2006-10-14_03_s4.0465.0465.2
	if ( $scan =~ /^.+\.(\d+)\.\d+\.\d+$/ ) {
		$scan =~ /^.+\.(\d+)\.\d+\.\d+$/;
		return $1;
	}
	else {
		n_print("scan nr1 not defined!\n");
		return 0;
	}
}

sub get_charge {
	my $self            = shift;
	my ($optional_scan) = @_;
	my $scan            = '';
	if ( defined($optional_scan) ) {
		$scan = $optional_scan;
	}
	else {
		$scan = $self->get_scan();
	}

	# 2006-10-14_03_s4.0465.0465.2
	if ( $scan =~ /^.+\.(\d+)\.\d+\.\d+$/ ) {
		$scan =~ /^.+\.\d+\.\d+\.(\d+)$/;
		return $1;
	}
	else {
		n_print("charge not defined!\n");
		return 0;
	}
}

sub n_print {
	print "  " . $_[0];
}

sub v_print {
	my $self = shift();
	print "  " . $_[0] if $self->{v};
}

sub print_object {    # print the object to the console with dumpvar
	my $self = shift();
	if ( defined( $_[0] ) ) {

		package main;
		require "dumpvar.pl";
		dumpValue( $_[0] );
	}
	else {

		package main;
		require "dumpvar.pl";
		dumpValue($self);
	}
}

1;

__END__

Lukas Reiter




