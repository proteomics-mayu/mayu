package Digest;
use strict;

use PepMass;

##################################################################
#
# Digest
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
# 2008.09.15
#            - added PepMass.pm
# 2007.06.22
#            - created
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $d = Digest->new( $verbose, $status );
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

	# to check during the digest whether peptides are within
	# a mass range
	my $pre_calc_aa_comb = 4;
	my $pm = PepMass->new( $pre_calc_aa_comb, 0, 0 );
	$self->{pep_mass} = $pm;

	return $self;
}

# Object Method
# Title	    :  get_nonred_rh_pep_digest()
# Usage     :
# Function	:  fully tryptic digest of an amino acid sequence
# Returns   :  $rh_pep:  non redundant peptide sequences
#                        from the digest
#              $nr_pep:  total number of non redundant protein
#                        sequences
# Args      :  $aa:      amino acid sequence
#              $min_pep_length: minimal peptide sequence length
sub get_nonred_rh_pep_digest {
	my $self = shift();
	my ( $aa, $min_pep_length, $nr_missed_cleavages, $min_pep_mass,
		$max_pep_mass )
	  = @_;

	# output
	my $rh_pep;
	my $nr_pep = 0;

	my $aa_length = length($aa);
	my $last_pos  = $aa_length - 1;

	# optimized
	my ( $ra_starts, $ra_ends ) =
	  get_trypt_coor( $aa, $aa_length, $last_pos );

	# optimized
	my $nr_ends         = @$ra_ends;
	my $last_ends_index = $nr_ends - 1;
	my $count           = 0;

	# loop through all the starts
	for ( my $i = 0 ; $i < @$ra_starts ; $i++ ) {
		my $highest_end = $i + $nr_missed_cleavages;
		$highest_end = $last_ends_index if $highest_end > $last_ends_index;

		# loop through the ends
		for ( my $j = $i ; $j <= $highest_end ; $j++ ) {

			my $pep_length = $ra_ends->[$j] - $ra_starts->[$i] + 1;
			if ( $pep_length >= $min_pep_length ) {

				my $pep_aa = substr( $aa, $ra_starts->[$i], $pep_length );
				my $mass   = $self->{pep_mass}->get_pep_mass($pep_aa);

				# check the mass constraint for the peptide
				if ( $mass >= $min_pep_mass && $mass <= $max_pep_mass ) {
					$count++;
					$rh_pep->{$pep_aa}++;
				}
			}
		}
	}
	$nr_pep = keys %$rh_pep;

	return ( $rh_pep, $nr_pep );
}

sub get_trypt_coor {
	my ( $aa, $aa_length, $last_pos ) = @_;

	# save coordinates of all possible peptide starts and ends
	my @starts = ();
	my @ends   = ();

	push @starts, 0;

	# use of [KR][^P] not straightforward because of
	# problematic handling of ..KKK..
	while ( $aa =~ /[KR]/g ) {
		my $pos = pos($aa);
		if ( $pos <= $last_pos ) {
			if ( substr( $aa, $pos, 1 ) ne 'P' ) {
				push @starts, $pos;
				push @ends,   $pos - 1;
			}
		}
	}

	push @ends, $last_pos;

	return ( \@starts, \@ends );
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




