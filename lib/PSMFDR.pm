package PSMFDR;
use strict;

use MayuTools;

##################################################################
#
# PSMFDR
#
# NOTES:
# All the error model classes have to implemement certain methods:
# - get_error_model_identifier
# - create_error_model
#
# calculates PSM fdrs
#
#
# SYNOPSIS:
#
# TODO:
# -
# -
#
# CHANGELOG:
# 2008.03.10
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $em = PSMFDR->new( $verbose, $status, $ratio );
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $verbose, $status, $target_decoy_ratio ) = @_;

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

	# used for error estimation
	if ( defined($target_decoy_ratio) ) {
		$self->{target_decoy_ratio} = $target_decoy_ratio;
	}
	else {
		$self->{target_decoy_ratio} = 1;
	}
	
	my $tools = MayuTools->new();
	$self->{tools} = $tools;

	# used to identify the error model
	$self->{error_model_identifier} = 'PSMFDR';

	# used for processed printing out
	$self->{output}{digits}         = 8;
	$self->{output}{bin_header}     = [];
	$self->{output}{summary_header} =
	  [ 'target_PSM', 'decoy_PSM', 'FP_PSM', 'TP_PSM' ];

	return $self;
}

# Object Method
# Title	    :  is_local_error_model()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub is_local_error_model {
	my $self = shift;

	return 0;
}

# Object Method
# Title	    :  get_error_model_identifier()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_error_model_identifier {
	my $self = shift;

	return $self->{error_model_identifier};
}

# Object Method
# Title	    :  create_error_model()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub create_error_model {
	my $self = shift;

	my $em = PSMFDR->new( $self->{v}, $self->{s},
		$self->{target_decoy_ratio} );

	return $em;
}

# Object Method
# Title	    :  set_psm_set()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_psm_set {
	my $self = shift;
	my ($ra_ra_psm) = @_;

	$self->estimate_error($ra_ra_psm);
}

# Object Method
# Title	    :  estimate_error()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub estimate_error {
	my $self = shift;
	my ( $ra_ra_psm ) = @_;

	my ( $tot_tar, $tot_dec, $tot_fp, $tot_tp, $tot_fp_var, $tot_fdr ) =
	  ( 0, 0, 0, 0, 0, 0 );

	my $nr_target_psm = $self->get_nr_target_psm($ra_ra_psm);
	my $nr_decoy_psm  = $self->get_nr_decoy_psm($ra_ra_psm);

	$self->{summed}{'target_ids'} = $nr_target_psm;
	$self->{summed}{'decoy_ids'}  = $nr_decoy_psm;
	$self->{summed}{'FP'} = $nr_decoy_psm * $self->{target_decoy_ratio};
	$self->{summed}{'TP'} =  $nr_target_psm - $self->{summed}{'FP'};
}

# Object Method
# Title	    :  get_error_bin_table()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_error_bin_table {
	my $self = shift;

	return [];
}

# Object Method
# Title	    :  get_error_summary()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_error_summary {
	my $self = shift;

	my $digits = $self->{output}{digits};
	my $format = "%." . $digits . "f";

	my $tar = $self->{summed}{'target_ids'};
	my $dec = $self->{summed}{'decoy_ids'};
	my $fp  = sprintf( "%.0f", $self->{summed}{'FP'} );
	my $tp  = sprintf( "%.0f", $self->{summed}{'TP'} );

	my $ra_error = [ $tar, $dec, $fp, $tp ];

	return $ra_error;
}

# Object Method
# Title	    :  get_nr_target_psm()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_target_psm {
	my $self = shift;
	my ($ra_ra_psm) = @_;

	my $nr = 0;
	foreach my $ra_psm (@$ra_ra_psm) {
		if ( $ra_psm->[5] == 0 ) {
			$nr++;
		}
	}
	return $nr;
}

# Object Method
# Title	    :  get_nr_decoy_psm()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_decoy_psm {
	my $self = shift;
	my ($ra_ra_psm) = @_;

	my $nr = 0;
	foreach my $ra_psm (@$ra_ra_psm) {
		if ( $ra_psm->[5] == 1 ) {
			$nr++;
		}
	}
	return $nr;
}

# Title	    :  get_summary_header()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_summary_header {
	my $self = shift;
	return $self->{output}{summary_header};
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






