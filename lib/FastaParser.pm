package FastaParser;
use strict;

##################################################################
#
# FastaParser
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
# Usage     :  my $fp = FastaParser->new( $verbose, $status );
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
# Title	    :  parse()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub parse {
	my $self = shift();
	my ($db) = @_;

	my $rh_id_seq;

	if ( defined($db) ) {
		if ( -e $db ) {
			my @seq = ();
			my $id  = '';

			open( DB, $db ) or warn $!;
			while (<DB>) {
				my $line = $_;
				chomp($line);
				if ( substr( $line, 0, 1 ) eq '>' ) {
					$rh_id_seq->{$id} = join( '', @seq )
					  if defined( $seq[0] );
					if ( $line =~ /^>([^\s]+)/ ) {
						$line =~ /^>([^\s]+)/;
						$id = $1;
					}
					else {
						$id = '';
					}
					@seq = ();
				}
				elsif ( $line =~ /[A-Za-z]/ ) {
					push @seq, $line;
				}
			}
			close(DB);
			$rh_id_seq->{$id} = join( '', @seq ) if defined( $seq[0] );

		}
		else {
			print "  Could not find $db!\n";
		}
	}
	else {
		print "  input fasta db not defined!\n";
	}

	return $rh_id_seq;
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




