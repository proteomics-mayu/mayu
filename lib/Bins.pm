package Bins;
use strict;

use Bin;

##################################################################
#
# Bins
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
# Usage     :  my $b = Bins->new( $verbose, $status );
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

	$self->{bins} = [];

	return $self;
}

# Object Method
# Title	    :  set_equidistant_bins()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_equidistant_bins {
	my $self = shift();
	my ( $ra, $nr_bins, $special_zero_bin ) = @_;

	# 0: equidistant
	# 1: equicount
	my $binning_type = 0;

	$self->set_bins( $ra, $nr_bins, $special_zero_bin, $binning_type );
}

# Object Method
# Title	    :  set_equicount_bins()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_equicount_bins {
	my $self = shift();
	my ( $ra, $nr_bins, $special_zero_bin ) = @_;

	# 0: equidistant
	# 1: equicount
	my $binning_type = 1;

	$self->set_bins( $ra, $nr_bins, $special_zero_bin, $binning_type );
}

# Object Method
# Title	    :  set_bins()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_bins {
	my $self = shift();
	my ( $ra, $total_nr_bins, $special_zero_bin, $binning_type ) = @_;

	# do not make a special zero bin if the total number of bins is one
	$special_zero_bin = 0 if $total_nr_bins == 1;
	
	# create a special bin for only the zero values
	my $normal_bins = $total_nr_bins;
	if ($special_zero_bin) {
		my $ra_filtered_values = [];
		foreach (@$ra) {
			push @$ra_filtered_values, $_;
		}
		$ra = $ra_filtered_values;

		$self->add_zero_bin();

		$normal_bins = $total_nr_bins - 1;
	}

	my $ra_breaks = [];
	if ( $binning_type == 0 ) {
		$ra_breaks = $self->get_equidistant_breaks( $ra, $normal_bins );
	}
	else {
		$ra_breaks = $self->get_equicount_breaks( $ra, $normal_bins );
	}

	my $bin_type = 1;
	for ( my $i = 0 ; $i < $normal_bins ; $i++ ) {
		my $lower_bound = $ra_breaks->[$i];
		my $upper_bound = $ra_breaks->[ $i + 1 ];
		$self->add_bin( $lower_bound, $upper_bound );
	}
}

# Object Method
# Title	    :  get_equidistant_breaks()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_equidistant_breaks {
	my $self = shift();
	my ( $ra, $nr_bins ) = @_;

	my $ra_breaks = [];

	# factor used to calculate the lowest and upper most breaks
	my $f = 10**8;

	# lower_break > x <= upper break
	# -> lowest break has to be lower than the lowest value
	# -> highest break because of numeric problems as well
	my @s = sort { $a <=> $b } @$ra;

	my $min  = $s[0];
	my $max  = $s[-1];
	my $span = $max - $min;
	my $tiny = $span / $f;
	my $part = $span / $nr_bins;

	# lowest break a little smaller then the smallest value
	push @$ra_breaks, $min - $tiny;

	# all but the highest break
	for ( my $i = 1 ; $i < $nr_bins ; $i++ ) {
		push @$ra_breaks, $min + $i * $part;
	}

	# add the highest break
	push @$ra_breaks, $max + $tiny;

	return $ra_breaks;
}

# Object Method
# Title	    :  get_equicount_breaks()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_equicount_breaks {
	my $self = shift();
	my ( $ra, $nr_bins ) = @_;
	
	# sort the entities according to value
	my @sorted = sort { $a <=> $b } @$ra;
	
	# add one element of the same value to the end
	# to circumvent boundary problems
	push @sorted, $sorted[-1];

	my $entries    = @$ra;
	my $nr_per_bin = $entries / $nr_bins;
	my $ra_index   = [];
	my $ra_ebr     = [];
	
	# loop over all breaks (nr_bins + 1)
	for ( my $i = 0 ; $i <= $nr_bins ; $i++ ) {

		# round to integers
		my $index = sprintf( "%.0f", $nr_per_bin * $i );
		my $value = $sorted[$index];
		push @$ra_index, $index;
		push @$ra_ebr,   $value;
	}

	# for each bin calculate the number of proteins rounded
	my $ra_N = [];
	for ( my $i = 0 ; $i < $nr_bins ; $i++ ) {
		my $N = $ra_index->[ $i + 1 ] - $ra_index->[$i] + 1;
		push @$ra_N, $N;
	}

	# correct the last index
	$ra_index->[-1]--;

	# correct the extremes
	my $f    = 1000;
	my $min  = $sorted[0];
	my $max  = $sorted[-1];
	my $span = abs( $max - $min );
	my $av   = $span / ( $nr_bins );
	my $tiny = $av / $f;
	$ra_ebr->[0] -= $tiny;
	$ra_ebr->[-1] += $tiny;

	# ra_index and ra_N not needed
	return $ra_ebr;
}

# Object Method
# Title	    :  add_zero_bin()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_zero_bin {
	my $self = shift();

	my $bin = Bin->new();
	$bin->add_testing_type(0);

	push @{ $self->{bins} }, $bin;
}

# Object Method
# Title	    :  add_bin()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_bin {
	my $self = shift();
	my ( $lower_bound, $upper_bound ) = @_;
	my $bin = Bin->new();
	$bin->add_testing_type( $lower_bound, $upper_bound );

	push @{ $self->{bins} }, $bin;
}

# Object Method
# Title	    :  add_att_to_bin()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_att_to_bin {
	my $self = shift();
	my ( $index, $att_name, $att_value ) = @_;
	my $ra_breaks = [];

	if ( exists( $self->{bins}[$index] ) ) {
		$self->{bins}[$index]->add_att( $att_name, $att_value );
	}
	else {
		print "  bin '$index' does not exist!\n";
	}
}

# Object Method
# Title	    :  get_equidistant_breaks()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_bin_nr {
	my $self = shift();
	my ($value) = @_;

	for ( my $i = 0 ; $i < @{ $self->{bins} } ; $i++ ) {
		if ( $self->{bins}[$i]->test($value) ) {
			return $i;
		}
	}

	# no bin could be found
	return -1;
}

# Object Method
# Title	    :  get_breaks()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_breaks {
	my $self = shift();

	my $ra_breaks = [];

	for ( my $i = 0 ; $i < @{ $self->{bins} } ; $i++ ) {
		push @$ra_breaks, $self->{bins}[$i]->get_test_string();
	}

	return $ra_breaks;
}

# Object Method
# Title	    :  get_bin_by_index()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_bin_by_index {
	my $self = shift();
	my ($index) = @_;

	if ( defined( $self->{bins}[$index] ) ) {
		return $self->{bins}[$index];
	}
	else {
		$self->n_print("bin $index does not exist!\n");
		return undef;
	}
}

# Object Method
# Title	    :  get_nr_bins()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_bins {
	my $self = shift();

	my $nr_bins = @{ $self->{bins} };

	return $nr_bins;
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




