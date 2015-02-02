package HypergProbDist;
use strict;

##################################################################
#
# HypergProbDist
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
# Usage     :  my $hpd = HypergProbDist->new( $verbose, $status );
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
	
	# TODO choose
	$self->{factorial_switch} = 1000;

	return $self;
}

# Object Method
# Title	    :  complete_hyperg()
# Usage     :
# Function	:
# Returns   :
# Args      :  N: total number of balls (white and black)
#              w: white balls
#              d: draws without replacement
sub complete_hyperg {
	my $self = shift;
	my ( $N, $w, $d ) = @_;

	my @x = ();
	for ( my $x = 0 ; $x <= $d ; $x++ ) {
		push @x, $self->hyperg( $x, $N, $w, $d );
	}

	return \@x;
}

# Object Method
# Title	    :  check_n_hyperg()
# Usage     :
# Function	:  first checks the input and then
#              calculates the probability using the logs
#              of the binomial coefficients.
# Returns   :
# Args      :  x: probability to draw x w balls
#              N: total number of balls (white and black)
#              w: white balls
#              d: draws without replacement
sub check_n_hyperg {
	my $self = shift;
	my ( $x, $N, $w, $d ) = @_;

	# check the input for consistency and correct
	( $x, $N, $w, $d ) = $self->check_hyperg_input( $x, $N, $w, $d );
	
	# hypergeometric probability
	return $self->hyperg( $x, $N, $w, $d );
}

# Title	    :  check_hyperg_input()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub check_hyperg_input {
	my $self = shift;
	my ( $fp, $N, $w, $cf ) = @_;

	$w  = $N  if $w > $N;
	$cf = $N  if $cf > $N;
	$fp = $cf if $fp > $cf;

	return ( $fp, $N, $w, $cf );
}

# Object Method
# Title	    :  hyperg()
# Usage     :
# Function	:  calculates the probability using the logs
#              of the binomial coefficients.
# Returns   :
# Args      :  x: probability to draw x w balls
#              N: total number of balls (white and black)
#              w: white balls
#              d: draws without replacement
sub hyperg {
	my $self = shift;
	my ( $x, $N, $w, $d ) = @_;

	return 0 if $d == 0;

	# natural logarithm of the probability
	my $log_prob =
	  $self->log_binomial( $w,      $x ) +
	  $self->log_binomial( $N - $w, $d - $x ) -
	  $self->log_binomial( $N,      $d );
	
	return exp($log_prob);
}

# Object Method
# Title	    :  log_binomial()
# Usage     :
# Function	:  calculates the log of the binomial
#              coefficient
# Returns   :
# Args      :
sub log_binomial {
	my $self = shift;
	my ( $n, $k ) = @_;

	my $log_bnok = $self->log_factorial($n) -
	  $self->log_factorial($k) -
	  $self->log_factorial( $n - $k );
	  
	return $log_bnok;
}

# Object Method
# Title	    :  log_factorial()
# Usage     :
# Function	:  uses stirling approximation for high numbers
#              to calculate the log of the factorials.
#              n! =~ sqrt(2*PI*n)*(n/e)^n
#              If the numbers are decent the factorials are
#              calculated regularly.
# Returns   :
# Args      :
sub log_factorial {
	my $self = shift;
	my ( $n ) = @_;

	my $switch = $self->{factorial_switch};
	
	if ( $n < $switch ) {
		my $log_fac = $self->exact_log_factorial( $n );
		return $log_fac;
	}
	else {
		my $log_fac = $self->stirling_log_factorial( $n );
		return $log_fac;
	}
}

# Object Method
# Title	    :  exact_log_factorial()
# Usage     :
# Function	:  
# Returns   :
# Args      :
sub exact_log_factorial {
	my $self = shift;
	my ( $n ) = @_;
	
	# log 1
	my $log_fac = 0;
	
	for ( my $i = 2; $i <= $n; $i++ ) {
		$log_fac += log($i);
	}
	return $log_fac;
}

# Object Method
# Title	    :  exact_factorial()
# Usage     :
# Function	:  
# Returns   :
# Args      :
sub exact_factorial {
	my $self = shift;
	my ( $n ) = @_;
	
	my $fac = 1;
	for ( my $i = 2; $i <= $n; $i++ ) {
		$fac *= $i;
	}
	return $fac;
}

# Object Method
# Title	    :  stirling_log_factorial()
# Usage     :
# Function	:  natural logarithm of n! approximated
#              by the stirling approximation.
#              log(n!) =~ log(sqrt(2*PI*n)) + $n*log(n) - n
# Returns   :
# Args      :
sub stirling_log_factorial {
	my $self = shift;
	my ( $n ) = @_;
	
	my $PI = 3.141592653589793;
	
	my $log_stir_fac = log( sqrt( 2*$PI*$n ) ) + $n*log($n) - $n;
	return $log_stir_fac;
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




