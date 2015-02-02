package CSVParser;
use strict;

##################################################################
#
# NOTES:
# - records separated with ,
# - strings quoted with "
# - , in strings allowed
# - " in strings escaped with "
# - a newline starts a new line of records
# - if a newline is in a quoted record the newline is integrated
#   into the record
#
#
# SYNOPSIS:
#
#
# TODO:
# -
# -
#
# CHANGELOG:
# 2007.07.13
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $p = CSVParser->new( $verbose, $status,... );
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $verbose, $status, $record_separator, $record_quoting_char,
		$record_quoting_char_escaping )
	  = @_;

	# save the options
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
	if ( defined($record_separator) ) {
		$self->{sep} = $record_separator;
	}
	else {
		$self->{sep} = ',';
	}
	if ( defined($record_quoting_char) ) {
		$self->{quote} = $record_quoting_char;
	}
	else {
		$self->{quote} = '"';
	}
	if ( defined($record_quoting_char_escaping) ) {
		$self->{esc} = $record_quoting_char_escaping;
	}
	else {
		$self->{esc} = '"';
	}

	# stops a record unless currently reading a quoted record
	$self->{newline} = "\n";

	return $self;
}

# Object Method
# Title	    :  parse_file()
# Usage     :  $p->parse_file( $csv );
# Function	:  newline does not start a new record
# Returns   :  returns an array of refs to arrays with records
# Args      :
sub parse_file {
	my $self = shift;
	my ( $csv ) = @_;

	if ( -e $csv ) {
		my ( $prev_char, $quoted, $ra_lines, $ra_records, $ra_curr_record );
	
		open( I, $csv ) or warn $!;
		while (<I>) {
			( $prev_char, $quoted, $ra_lines, $ra_records, $ra_curr_record ) =
			  $self->parse_unfinished_string( $_, $prev_char, $quoted,
				$ra_lines, $ra_records, $ra_curr_record );
		}
		close(I);
	
		# finish the last record if not already done
		unless ( $prev_char =~ /$self->{newline}/ ) {
			my $record = join( "", @$ra_curr_record );
			push @$ra_records, $record;
			push @$ra_lines,   $ra_records;
		}
		
		# if there was no input data
		$ra_lines = [ [ '' ] ] unless defined( $ra_lines );
		$ra_lines = [ [ '' ] ] if @$ra_lines == 0;
	
		return @$ra_lines;
	}
	else {
		print "  csv file '$csv' does not exists!\n" if $self->{v};
		return ( [ '' ] );
	}
}

# Object Method
# Title	    :  parse_string()
# Usage     :  $p->parse_string( $string );
# Function	:  newline does not start a new record
# Returns   :
# Args      :
sub parse_string {
	my $self = shift();
	my ($string) = @_;

	# undefined input
	unless ( defined($string) ) {
		print "  string not defined!\n" if $self->{v};
		return ( [''] );
	}

	# length of input string
	my $L = length($string);

	# input string is empty
	if ( $L == 0 ) {
		print "  empty string!\n" if $self->{v};
		return ( [''] );
	}

	# input string has some length
	else {
		my ( $prev_char, $quoted, $ra_lines, $ra_records, $ra_curr_record )
		  = $self->parse_unfinished_string($string);

		# finish the last record
		my $record = join( "", @$ra_curr_record );
		push @$ra_records, $record;
		push @$ra_lines,   $ra_records;

		return @$ra_lines;
	}
}

# Object Method
# Title	    :  parse_unfinished_string()
# Usage     :  $p->parse_string( $string );
# Function	:  newline does not start a new record
# Returns   :
# Args      :
sub parse_unfinished_string {
	my $self = shift;
	my ( $string, $prev_char, $quoted, $ra_lines, $ra_records,
		$ra_curr_record )
	  = @_;

	# previous character
	my $prev_char = $self->{sep} unless defined($prev_char);

	# currently characters are quoted
	my $quoted = 0 unless defined($quoted);

	# lines of arrays of records
	my $ra_lines = [] unless defined($ra_lines);

	# array of records
	my $ra_records = [] unless defined($ra_records);

	# current record as an array
	my $ra_curr_record = [] unless defined($ra_curr_record);

	# length of input string
	my $L = length($string);

	for ( my $i = 0 ; $i < $L ; $i++ ) {
		my $c = substr( $string, $i, 1 );

		#--------------------------------
		# currently reading unquoted
		#--------------------------------
		unless ($quoted) {

			# separator, a new record starts
			if ( $c =~ /$self->{sep}/ ) {
				my $record = join( "", @$ra_curr_record );
				push @$ra_records, $record;
				$ra_curr_record = [];
			}

			# a new line of records starts
			elsif ( $c =~ /$self->{newline}/ ) {
				my $record = join( "", @$ra_curr_record );
				push @$ra_records, $record;
				my @copy = @$ra_records;
				push @$ra_lines, \@copy;
				$ra_curr_record = [];
				$ra_records     = [];
			}

			# quoting string
			elsif ( $c =~ /$self->{quote}/ ) {

				# a quoted record starts
				if ( $prev_char =~ /$self->{sep}/ ) {
					$quoted = 1;
				}

				# a quote after an escaping quote
				elsif ( $prev_char =~ /$self->{quote}/ ) {
					push @$ra_curr_record, $c;
					$quoted = 1;
				}

				# a quote in the middle of an unquoted record
				else {
					push @$ra_curr_record, $c;
				}
			}

			# a record character
			else {
				push @$ra_curr_record, $c;
			}
		}

		#--------------------------------
		# currently reading quoted record
		#--------------------------------
		else {

			# either an escaping quote or a stopping quote
			if ( $c =~ /$self->{quote}/ ) {
				$quoted = 0;
			}

			# a record character
			else {
				push @$ra_curr_record, $c;
			}
		}

		# set the current character to the previous character
		$prev_char = $c;

	}    # end for loop over string length

	return ( $prev_char, $quoted, $ra_lines, $ra_records,
		$ra_curr_record );
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




