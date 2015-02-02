package PepMass;
use strict;

##################################################################
#
# PepMass
#
#
# NOTES:
#
# SYNOPSIS:
#
# TODO:
# - monoisotopic, average
#   as soon as average is added the nlets precalculations for
#   average has to be added to precalculate_nlets( ...,$ra_sets )
# - 
#
# CHANGELOG:
# 2008.09.11
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $d = PepMass->new( $nlets, $v, $s );
# Function  :
# Returns   :
# Args      :  $nlets:  make precalculations for n-lets of
#                       amino acids. 4 seems to be a balanced
#                       value
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $nlets, $v, $s ) = @_;

	# initialize the masses of the building blocks
	# and do some precalculations of building blocks if selected
	# (nlets for amino acids)
	$self->init( $nlets, $v, $s );

	return $self;
}

# Object Method
# Title	    :  init()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub init {
	my $self = shift();
	my ( $nlets, $v, $s ) = @_;

	if ( defined($v) ) {
		$self->{v} = $v;
	}
	else {
		$self->{v} = 0;
	}
	if ( defined($s) ) {
		$self->{s} = $s;
	}
	else {
		$self->{s} = 0;
	}
	
	# mono or average
	$self->{mass_type} = 'mono';

	#-------------------------------------------------------------
	# monoisotopic masses
	#-------------------------------------------------------------
	
	# monoisotopic masses of the amino acids without water
	my $rh_mono_aa_wo_water = {
		'A' => 71.03711,
		'C' => 103.00919,
		'D' => 115.02694,
		'E' => 129.04259,
		'F' => 147.06841,
		'G' => 57.02146,
		'H' => 137.05891,
		'I' => 113.08406,
		'K' => 128.09496,
		'L' => 113.08406,
		'M' => 131.04049,
		'N' => 114.04293,
		'P' => 97.05276,
		'Q' => 128.05858,
		'R' => 156.10111,
		'S' => 87.03203,
		'T' => 101.04768,
		'U' => 149.90419,    # Selenocysteine
		'V' => 99.06841,
		'W' => 186.07931,
		'Y' => 163.06333
	};
	$self->{mono}{aa_wo_water} = $rh_mono_aa_wo_water;

	# monoisotopic masses of elements
	my $rh_mono_el = {
		'H' => 1.007825,
		'O' => 15.994915
	};
	$self->{mono}{el} = $rh_mono_el;

	# monoisotopic masses of organic molecules
	my $rh_mono_org_mol = { 'water' => 18.01056 };
	$self->{mono}{org_mol} = $rh_mono_org_mol;

	#-------------------------------------------------------------
	# precalculate amino acid combinations
	#-------------------------------------------------------------
	
	# the size of the nlets that are precalculated
	# e.g. an nlet is a combination of n amino acids
	$self->{nlet} = 1;

	# determined on T60p 2GB RAM
	$self->{maxnlet} = 5;
	
	# more than seven is too much for any RAM in the
	# near future (20^7 =~ 1.28*10^9 )
	if ( defined($nlets) ) {
		if ( $nlets > 1 ) {

			# the sets for which the precalculations should be done
			my $ra_sets = [ $self->{mono}{aa_wo_water} ];
			$self->precalculate_nlets( $nlets, $ra_sets );
		}
	}
}

# Object Method
# Title	    :  precalculate_nlets()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub precalculate_nlets {
	my $self = shift();
	my ( $nlets, $ra_sets ) = @_;

	if ( $nlets > $self->{maxnlet} ) {
		print "  n-lets is $nlets. set to " . $self->{maxnlet} . "!\n"
		  if $self->{v};
		$nlets = $self->{maxnlet};
	}
	print "  PepMass: precalculating nlets for n = $nlets...\n"
	  if $self->{v};

	# all the sets for which the precalculations should be done
	# e.g.
	# - monoisotopic amino acids
	# - average mass amino acids
	foreach my $rh_set (@$ra_sets) {

		# initialize
		my %nlets    = %$rh_set;
		my $rh_nlets = \%nlets;

		# for nlets - 1 add a variation level
		for ( my $i = 1 ; $i < $nlets ; $i++ ) {

			# the new hash with one level more
			$rh_nlets = $self->add_variation_level( $rh_nlets, $rh_set );
		}

		foreach ( keys %$rh_nlets ) {
			$rh_set->{$_} = $rh_nlets->{$_};
		}
	}

	# save the size of the nlets that are precalculated
	$self->{nlet} = $nlets;
}

# Object Method
# Title	    :  add_variation_level()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_variation_level {
	my $self = shift();
	my ( $rh_nlets, $rh_set ) = @_;

	my $rh_new_nlets;

	foreach my $base ( keys %$rh_nlets ) {

		# add also the smaller nlets
		$rh_new_nlets->{$base} = $rh_nlets->{$base};

		# add one level to the nlets
		foreach my $add ( keys %$rh_set ) {

			# update the key of the new nlet with one more level
			my $new_key = $base . $add;

			# update the mass of the new nlet
			$rh_new_nlets->{$new_key} =
			  $rh_nlets->{$base} + $rh_set->{$add};
		}
	}

	return $rh_new_nlets;
}

# Object Method
# Title	    :  get_mz_by_peps_n_charges()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_mz_by_peps_n_charges {
	my $self = shift();
	my ( $ra_pep, $ra_charge ) = @_;

	my $ra_mz = [];

	foreach my $pep (@$ra_pep) {
		if ( $pep =~ /^[A-Z]+$/ ) {
			foreach my $charge (@$ra_charge) {
				if ( $charge =~ /^\d+$/ ) {
					push @$ra_mz,
					  $self->get_mz_by_pep_n_charge( $pep, $charge );
				}
			}
		}
	}

	return $ra_mz;
}

# Object Method
# Title	    :  get_mz_by_pep_n_charge()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_mz_by_pep_n_charge {
	my $self = shift();
	my ( $pep, $charge ) = @_;

	my $pep_mass = $self->get_pep_mass($pep);
	
	# currently set mass type (default is mono)
	my $mt = $self->{mass_type};

	my $prec_mz =
	  ( $pep_mass + ( $charge * $self->{$mt}{el}{H} ) ) / $charge;

	return $prec_mz;
}

# Object Method
# Title	    :  get_pep_mass()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_pep_mass {
	my $self = shift();
	my ($pep) = @_;

	my $pep_mass = 0;
	
	# currently set mass type (default is mono)
	my $mt = $self->{mass_type};

	if ( defined($pep) ) {
		if ( $pep !~ /^$/ ) {
			
			# defines the sizes of the parts that are added up
			my $nlet = 1;
			if ( exists( $self->{nlet} ) ) {
				$nlet = $self->{nlet} if $self->{nlet} > 1;
			}

			my $pep_len = length($pep);

			my $start = 0;
			while ( $start < $pep_len ) {
				my $end = $start + $nlet - 1;
				$end = $pep_len - 1 if $end > $pep_len - 1;
				my $len = $end - $start + 1;

				my $ss = substr( $pep, $start, $len );

				if ( exists( $self->{$mt}{aa_wo_water}{$ss} ) ) {
					$pep_mass += $self->{$mt}{aa_wo_water}{$ss};
				}

				# update start
				$start += $nlet;
			}

			$pep_mass += $self->{$mt}{org_mol}{water};
		}
	}

	return $pep_mass;
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




