package MayuTools;
use strict;

##################################################################
#
# NOTES:
# A set of tools
#
#
# CHANGELOG:
# 2008.02.25
#            - created
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $tools = MayuTools->new( $verbose, $status );
# Function  :
# Returns   :
# Args      :  $verbose:   print out some stuff or not (1 or 0)
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
# Title	    :  get_file_name_base()
# Usage     :  $t->get_file_name_base();
# Function	:  creates a timestamp file name base
# Returns   :
# Args      :
sub get_file_name_base {
	my $self = shift;
	my (
		$Sekunden, $Minuten,   $Stunden,   $Monatstag, $Monat,
		$Jahr,     $Wochentag, $Jahrestag, $Sommerzeit
	  )
	  = localtime(time);
	$Jahr += 1900;
	$Monat++;
	my ( $Mo, $Ta, $St, $Mi, $Se ) =
	  $self->add_leading_zero(
		[ $Monat, $Monatstag, $Stunden, $Minuten, $Sekunden ], 2 );
	return $Jahr . '-' . $Mo . '-' . $Ta . '_' . $St . '.' . $Mi . '.'
	  . $Se;
}

# my ( $nr ) = add_leading_zero( [ $file_count ], 2 );
#
# use:     my ( $nr_with_leading_zero1, $nr_with_leading_zero2 ) =
#               add_leading_zero( $ra_nr , $format );
# $ra_nr:  a reference to an array with numbers
# $format: number of final digits
# returns: adds leading zero to a number such that the final
#          number has $format digits
sub add_leading_zero {
	my $self = shift;
	my ( $ra, $f ) = @_;
	my @new;
	foreach (@$ra) {

		# number of digits of the current number
		my $l  = length($_);
		my $os = '';
		for ( my $i = 0 ; $i < $f - $l ; $i++ ) {
			$os .= '0';
		}

		push @new, $os . $_;
	}
	return @new;
}

# my ( $nr ) = add_trailing_zero( [ $number1, $number2 ], 2 );
#
# use:     my ( $nr_with_trailing_zero1, $nr_with_trailing_zero2 ) =
#               add_trailing_zero( $ra_nr , $format );
# $ra_nr:  a reference to an array with numbers
# $format: number of final digits
# returns: adds trailing zero to a number such that the final
#          number has $format digits
sub add_trailing_zero {
	my $self = shift;
	my ( $ra, $f ) = @_;
	my @new;
	foreach my $nr (@$ra) {

		# if there is no point
		unless ( $nr =~ /\./ ) {
			$nr .= '.';
		}

		# number of digits of the current number
		my $l  = length($nr);
		my $os = '';
		for ( my $i = 0 ; $i < $f - $l ; $i++ ) {
			$os .= '0';
		}

		push @new, $nr . $os;
	}
	return @new;
}

#-------------------------------------------------------------
# get files with a file suffix
# if no suffix is specified all files are returned
#-------------------------------------------------------------
# $in:     directory with wanted files or a file
# $fs:     file suffix, including the point
# returns: the file names including paths of the files
#          with the file suffix $fs in the $in
sub get_files {
	my $self = shift;
	my ( $in, $fs ) = @_;
	my @return_files = ();

	# return an empty array if the dir input is not defined
	unless ( defined($in) ) {
		return @return_files;
	}

	# input is a file
	elsif ( -f $in ) {

		# file has a file ending
		if ( defined($fs) ) {
			if ( $in =~ /\.[^\.]+$/ ) {
				$in =~ /(\.[^\.]+)$/;
				my $curr_fs = $1;
				if ( $curr_fs eq $fs ) {
					return ($in);
				}
				else {
					return ();
				}
			}
			else {
				return ();
			}
		}
		else {
			return ($in);
		}
	}
	elsif ( -d $in ) {    # input is a directory
		opendir( IN, $in ) or die "get_files(): Couldn't open $in: $!\n";
		my @curr_files = readdir(IN);
		closedir(IN) or warn $!;
		if ( defined($fs) ) {
			foreach my $curr_file (@curr_files) {
				my $conc = $in . $curr_file;
				if ( $curr_file =~ /\.[^\.]+$/ && -f $conc ) {
					$curr_file =~ /(\.[^\.]+)$/;
					my $curr_fs = $1;
					if ( $curr_fs eq $fs ) {
						push @return_files, $conc;
					}
				}
			}
			return @return_files;
		}
		else {
			foreach (@curr_files) {
				my $conc = $in . $_;
				if ( -f $conc ) {    # if it's not . or ..
					push @return_files, $conc;
				}
			}
			return @return_files;
		}
	}
	else {
		return ();
	}
}

# Object Method
# Title	    :  parse_table_file()
# Usage     :  
# Function	:  parse a table file whose columns are separated with
#              $split
# Returns   :
# Args      :
sub parse_table_file {
	my $self = shift;
	my ( $in, $split ) = @_;
	my $ra_ra = [];

	# get_files just returns all the files if no
	my @f = $self->get_files($in);
	foreach my $f (@f) {
		$self->v_print("parsing table file $f\n");
		open( IN, "$f" ) or die $!;
		while (<IN>) {
			my $line = $_;
			chomp($line);
			unless ( $line =~ /^\s*#|^\s*$/ ) {
				my @col = split( /$split/, $line );
				push @$ra_ra, \@col;
			}
		}
		close(IN) or warn $!;
	}
	return $ra_ra;
}

# Object Method
# Title	    :  print_data_in_table_style()
# Usage     :  
# Function	:  prints data in table style with non shifting columns.
#              If the file already exists the data is appended
# Returns   :
# Args      :
sub print_data_in_table_style {
	my $self = shift;
	my ( $ra, $out ) = @_;

	# create the formats
	my @f;
	my $ra_m = $self->get_max_of_each_column($ra);
	for ( my $i = 0 ; $i < @$ra_m ; $i++ ) {
		$f[$i] = "%-" . $ra_m->[$i] . "s";
	}
	my $s1 = '';
	my $s2 = '';
	my $c  = 0;
	foreach (@f) {
		$s1 .= '$f[' . $c . '] ';
		$s2 .= '$_->[' . $c . '],';
		$c++;
	}
	chop($s1);
	chop($s2);

	my $ev = '
	foreach ( @$ra ) {
		printf P "'
	  . $s1 . '\n", ' . $s2 . ';
	}';

	open( P, ">>$out" ) or die $!;
	eval($ev);
	close(P) or warn $!;
}

# Object Method
# Title	    :  print_data_to_stdo_in_table_style()
# Usage     :  
# Function	:  prints data in table style with non shifting columns
#              to standard out
# Returns   :
# Args      :
sub print_data_to_stdo_in_table_style {
	my $self = shift;
	my ($ra) = @_;

	# create the formats
	my @f;
	my $ra_m = $self->get_max_of_each_column($ra);
	for ( my $i = 0 ; $i < @$ra_m ; $i++ ) {
		$f[$i] = "%-" . $ra_m->[$i] . "s";
	}
	my $s1 = '';
	my $s2 = '';
	my $c  = 0;
	foreach (@f) {
		$s1 .= '$f[' . $c . '] ';
		$s2 .= '$_->[' . $c . '],';
		$c++;
	}
	chop($s1);
	chop($s2);

	my $ev = '
	foreach ( @$ra ) {
		printf "'
	  . $s1 . '\n", ' . $s2 . ';
	}';

	eval($ev);
}

# Object Method
# Title	    :  get_max_of_each_column()
# Usage     :  
# Function	:  get the maximum width of each column from
#              a table like structure
# Returns   :
# Args      :
sub get_max_of_each_column {
	my $self = shift;
	my $ra   = $_[0];
	my @max;
	for ( my $i = 0 ; $i < @{ $ra->[0] } ; $i++ ) {
		my @sorted =
		  sort { length( $a->[$i] ) <=> length( $b->[$i] ) } @{$ra};
		push @max, length( $sorted[-1]->[$i] );
	}
	return \@max;
}

# Object Method
# Title	    :  div()
# Usage     :  
# Function	:  
# Returns   :
# Args      :
sub div {
	my $self = shift;
	my ( $x, $y ) = @_;
	return $x if $y == 0;
	return $x / $y;
}

# Object Method
# Title	    :  get_var()
# Usage     :  
# Function	:  
# Returns   :
# Args      :  $ra:  an array of probabilities. the array index
#                    corresponds to the value associated with the
#                    probability.
#              $exp: the expectation value of the probability 
#                    distribution
sub get_var {
	my $self = shift;
	my ( $ra, $exp ) = @_;

	my $var = 0;
	for ( my $i = 0 ; $i < @$ra ; $i++ ) {
		$var += ( ( $i - $exp )**2 ) * $ra->[$i];
	}

	return $var;
}

# Object Method
# Title	    :  get_exp_val()
# Usage     :  
# Function	:  
# Returns   :
# Args      :  $ra:  an array of probabilities. the array index
#                    corresponds to the value associated with the
#                    probability.
sub get_exp_val {
	my $self = shift;
	my ($ra) = @_;

	my $exp = 0;
	for ( my $i = 0 ; $i < @$ra ; $i++ ) {
		$exp += $ra->[$i] * $i;
	}

	return $exp;
}

# Object Method
# Title	    :  get_min()
# Usage     :  
# Function	:  return the minimum of an array of numbers
# Returns   :
# Args      :
sub get_min {
	my $self = shift;
	my ($ra) = @_;
	if ( @$ra > 0 ) {
		my @s = sort { $a <=> $b } @$ra;
		return $s[0];
	}
	else {
		return undef;
	}
}

# Object Method
# Title	    :  get_max()
# Usage     :  
# Function	:  get the maximum of an array of numbers
# Returns   :
# Args      :
sub get_max {
	my $self = shift;
	my ($ra) = @_;
	if ( @$ra > 0 ) {
		my @s = sort { $a <=> $b } @$ra;
		return $s[-1];
	}
	else {
		return undef;
	}
}

# Object Method
# Title	    :  get_sum()
# Usage     :  
# Function	:  get the sum of an array of numbers
# Returns   :
# Args      :
sub get_sum {
	my $self = shift;
	my ($ra) = @_;
	my $sum = 0;
	foreach (@$ra) {
		if ( defined($_) ) {
			$sum += $_;
		}
	}
	return $sum;
}

# Object Method
# Title	    :  mean()
# Usage     :  
# Function	:  calculate the mean of an array of numbers
# Returns   :
# Args      :
sub mean {
	my $self = shift;
	my $ra = $_[0];
	
	return undef unless defined($ra);
	my $nr = @$ra;
	return undef if $nr == 0;
	
	my $c = 0;
	my $s = 0;
	foreach ( @$ra ) {
		if ( $_ =~ /[\d\.]+/ ) {
			$s += $_;
			$c++;
		}
	}
	return undef if $c == 0;
	return $s / $c;
}

# Object Method
# Title	    :  rms_error()
# Usage     :  
# Function	:  calculate the standard deviation of an array
#              of numbers
# Returns   :
# Args      :
sub rms_error {
	my $self = shift;
	my $ra   = $_[0];
	
	my $mean = $self->mean($ra);

	my $c = 0;
	my $s = 0;
	foreach ( @{$ra} ) {
		$s += ( $_ - $mean )**2;
		$c++;
	}
	return undef if $c == 0;
	return 0 if $c == 1;
	return sqrt( $s / ( $c - 1 ) );
}

# Object Method
# Title	    :  shuffle_array()
# Usage     :  
# Function	:  
# Returns   :
# Args      :
sub shuffle_array {
	my $self = shift;
	my ( $ra_in ) = @_;
	
	my @ra = @$ra_in;
	my $ra = \@ra;
    my $i = scalar(@$ra);
    my $j;
    foreach my $item (@$ra )
    {
        --$i;
        $j = int rand ($i+1);
        next if $i == $j;
        @$ra [$i,$j] = @$ra[$j,$i];
    }
    return $ra;
}

# Object Method
# Title     :  linear_inter_or_extrapolate
# Usage     :
# Function  :  v1 a value the corresponds to the first column in the
#              two dimensional array $ra_ra_v1_vn, the linearly
#              interpolated value of the second column is returned
#              min and max describe the minimal and maximal values
#              for the value to be returned.
#              - add a column to the 2D array delta(v1)
#              - sort according to delta(v1)
#              - find the two lowest delta(v1) values
#              - v2 = a*v1 + b
#                a  = (vn2 - vn1) / (v12 - v11)
#                     (yn - y1) / (xn - x1)
#                b  = vn1 - a*v11
#                     y1 - a*x1
#              - check whether x2 - x1 is not 0
#              - if yes find another data pair
#
#              this function is not so performant if the $ra_ra
#              is very big (~10000)
# Returns   :  the inter- or extraplated vn value(s)
# Args      :  $v1:          value corresponding to the first column
#              $ra_ra_v1_vn: 2D array, with n being 1 to n
sub linear_inter_or_extrapolate {
	my $self = shift;
	my ( $v1, $ra_ra_v1_vn ) = @_;

	unless ( defined($v1) && defined($ra_ra_v1_vn) ) {
		print "  undefined input!\n";
		return 0;
	}

	# add a delta v1 column
	# absolut distance of input v1 to the value
	# in the first column of the table
	my @delta_indices = ();
	for ( my $i = 0; $i < @$ra_ra_v1_vn; $i++ ) {
		my $delta = abs( $v1 - $ra_ra_v1_vn->[$i][0] );
		
		# input value is the same as an x value in the table
		if ( $delta == 0 ) {
			my @interpolated_ys = @{$ra_ra_v1_vn->[$i]};
			
			# remove the x value
			shift @interpolated_ys;
			return @interpolated_ys;
		}
		push @delta_indices, [ $delta, $i ];
		
	}

	# number of values depending on the x value
	# substract 1 because of x
	my $nr_ys = @{ $ra_ra_v1_vn->[0] };
	$nr_ys -= 1;
	
	# sort according to the delta (input - value in first column)
	my @sorted_delta = sort { $a->[0] <=> $b->[0] } @delta_indices;

	# find the nearest data point pair that works
	my $try_again      = 1;
	my $ind1           = 0;
	my $ind2           = 1;
	my $nr_data_points = @sorted_delta;
	my $ind1_iterator  = 0;
	
	my $count = 0;
	while ($try_again) {
		
		# guarantee that the loop ends
		$count++;
		$try_again = 0 if $count > $nr_data_points;
		
		# $sorted_delta[$ind1][1]: indices in the input array
		# $sorted_delta[$ind1][0]
		my $delta_x = $ra_ra_v1_vn->[$sorted_delta[$ind1][1]][0] 
			- $ra_ra_v1_vn->[$sorted_delta[$ind2][1]][0];
		if ( $delta_x > 0.000001 ) {
			$try_again = 0;
		}
		else {
			if ( $ind2 == ( $nr_data_points - 1  ) ) {
				if ( $ind1 == $nr_data_points - 1 ) {
					print "  no good data point pair found!\n";
					return 0;
				}
				else {
					$ind1++;
					$ind2 = $ind1 + 1;
				}
			}
			else {
				$ind2++;
			}
		}
	}

	# calculate the inter- extrapolated value
	my @interpolated_ys = ();
	for ( my $i = 0 ; $i <= $nr_ys ; $i++ ) {
		my $delta_x = $ra_ra_v1_vn->[$sorted_delta[$ind1][1]][0] 
			- $ra_ra_v1_vn->[$sorted_delta[$ind2][1]][0];
		my $delta_y = $ra_ra_v1_vn->[$sorted_delta[$ind1][1]][ $i + 1 ] 
			- $ra_ra_v1_vn->[$sorted_delta[$ind2][1]][ $i + 1 ];
		my $a       = $delta_y / $delta_x;
		my $b       = $ra_ra_v1_vn->[$sorted_delta[$ind1][1]][ $i + 1 ] 
			- $a * $ra_ra_v1_vn->[$sorted_delta[$ind1][1]][0];

		my $vn = $a * $v1 + $b;
		push @interpolated_ys, $vn;
	}

	return @interpolated_ys;
}

# Object Method
# Title	    :  set_time()
# Usage     :  
# Function	:  set starttime
# Returns   :
# Args      :
sub set_time {
	my $self = shift;
	my $t = $_[0];
	if ( defined($t) ) {
		$self->{timestamp} = $t;
	}
	else {
		$self->{timestamp} = time();
	}
}

# Object Method
# Title	    :  get_time()
# Usage     :  
# Function	:  time since initialization
# Returns   :
# Args      :
sub get_time {
	my $self = shift;
	if ( exists($self->{timestamp}) ) {
		if (defined($self->{timestamp})) {
			return time() - $self->{timestamp};
		}
		else {
			$self->set_time();
			return 0;
		}
	}
	else {
		$self->set_time();
		return 0;
	}
}

# Object Method
# Title	    :  get_starttime()
# Usage     :  
# Function	:  get the starttime in seconds
# Returns   :
# Args      :
sub get_starttime {
	my $self = shift;
	if ( exists($self->{timestamp}) ) {
		if (defined($self->{timestamp})) {
			return $self->{timestamp};
		}
		else {
			$self->set_time();
			return $self->get_starttime();
		}
	}
	else {
		$self->set_time();
		return $self->get_starttime();
	}
}

# Object Method
# Title	    :  print_time()
# Usage     :  
# Function	:  prints the time from initialization
# Returns   :
# Args      :
sub print_time {
	my $self = shift;
	n_print("time: " . $self->get_time() . "\n");
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




