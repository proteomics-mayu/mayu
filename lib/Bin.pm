package Bin;
use strict;

##################################################################
#
# Bin
#
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
# 2007.06.22
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $b = Bin->new( $verbose, $status );
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $verbose, $status ) = @_;

	if ( defined($verbose) ) {
		$self->{v} = $verbose;
	}
	else {
		$self->{v} = 0;
	}
	if ( defined($status) ) {
		$self->{s} = $status;
	}
	else {
		$self->{s} = 0;
	}

	return $self;
}

# Object Method
# Title	    :  add_testing_type()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_testing_type {
	my $self = shift();
	my ( $one, $two ) = @_;

	unless ( defined($two) ) {
		$self->{testing_type} = 'x=one';
		$self->{one}          = $one;
	}
	else {
		$self->{testing_type} = 'one>x>=two';
		$self->{one}          = $one;
		$self->{two}          = $two;
	}
}

# Object Method
# Title	    :  test()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub test {
	my $self = shift();
	my ($value) = @_;

	if ( $self->{testing_type} eq 'x=one' ) {
		if ( $value == $self->{one} ) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		if ( $value > $self->{one} && $value <= $self->{two} ) {
			return 1;
		}
		else {
			return 0;
		}
	}
}

# Object Method
# Title	    :  get_test_string()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_test_string {
	my $self = shift();

	if ( $self->{testing_type} eq 'x=one' ) {
		my $s = 'x = ' . $self->{one};
		return $s;
	}
	else {
		my $s = $self->{one} . ' < x <= ' . $self->{two};
		return $s;
	}
}

# Object Method
# Title	    :  get_nice_test_string()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nice_test_string {
	my $self = shift();

	if ( $self->{testing_type} eq 'x=one' ) {
		my $s = 'x = ' . $self->{one};
		return $s;
	}
	else {
		my $s = int( $self->{one} ) . ' < x <= ' . int( $self->{two} );
		return $s;
	}
}

# Object Method
# Title	    :  add_att()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_att {
	my $self = shift();
	my ( $name, $value ) = @_;

	$self->{att}{$name} = $value;
}

# Object Method
# Title	    :  get_att_by_name()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_att_by_name {
	my $self = shift();
	my ($name) = @_;

	if ( exists( $self->{att}{$name} ) ) {
		return $self->{att}{$name};
	}
	else {
		$self->n_print("  attribute '$name' does not exists!\n");
		return undef;
	}
}

sub n_print {
	print "  " . $_[0];
}

sub v_print {
	my $self = shift();
	print "  " . $_[0] if $self->{v};
}

1;

__END__

Lukas Reiter




