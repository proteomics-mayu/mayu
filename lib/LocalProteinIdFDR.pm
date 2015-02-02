package LocalProteinIdFDR;
use strict;

use MayuTools;

##################################################################
#
# LocalProteinIdFDR
#
# NOTES:
# All the error model classes have to implemement certain methods:
# - get_error_model_identifier
# - create_error_model
# - ...
#
#
# SYNOPSIS:
#
# TODO:
# -
# -
#
# CHANGELOG:
# 2008.03.13
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $em = LocalProteinIdFDR->new( $v, $s,
#                           $target_decoy_ratio );
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
	$self->{error_model_identifier} = 'LocalProteinIdFDR';

	# needs input (global FDR) from
	# to calculate the local FDR using bayes
	$self->{needed_error_model_identifier} = 'ProteinIdFDR';

	# used for processed printing out
	$self->{output}{digits}            = 8;
	$self->{output}{single_hit_header} = [
		'target_protIDs', 'decoy_protIDs',
		'FP_protIDs',     'TP_protIDs',
		'protFDRs'
	];
	$self->{output}{no_single_hit_header} = [
		'target_protIDns', 'decoy_protIDns',
		'FP_protIDns',     'TP_protIDns',
		'protFDRns'
	];

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

	return 1;
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
# Title	    :  get_needed_error_model_identifier()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_needed_error_model_identifier {
	my $self = shift;

	return $self->{needed_error_model_identifier};
}

# Object Method
# Title	    :  create_error_model()
# Usage     :
# Function	:  creates a copy object using the same binning ob
# Returns   :
# Args      :
sub create_error_model {
	my $self = shift;

	my $em =
	  LocalProteinIdFDR->new( $self->{v}, $self->{s},
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
	my ( $fdr, $rh_rh_id_protein_features ) = @_;

	my ( $target_protIDs, $target_protIDns, $target_protID ) = ( 0, 0, 0 );
	my ( $decoy_protIDs, $decoy_protIDns, $decoy_protID ) = ( 0, 0, 0 );

	foreach my $pid ( keys %$rh_rh_id_protein_features ) {
		my $rh = $rh_rh_id_protein_features->{$pid};
		if ( exists( $rh->{decoy} ) && exists( $rh->{NS} ) ) {

			# target protein ids
			if ( $rh->{decoy} == 0 ) {
				$target_protID++;
				if ( $rh->{NS} == 1 ) {
					$target_protIDs++;
				}
				elsif ( $rh->{NS} > 1 ) {
					$target_protIDns++;
				}
			}

			# decoy protein ids
			elsif ( $rh->{decoy} == 1 ) {
				$decoy_protID++;
				if ( $rh->{NS} == 1 ) {
					$decoy_protIDs++;
				}
				elsif ( $rh->{NS} > 1 ) {
					$decoy_protIDns++;
				}
			}
		}
	}
	
	# correct for the target to decoy ratio
	my $decoy_protIDs = $decoy_protIDs * $self->{target_decoy_ratio};
	my $decoy_protIDns = $decoy_protIDns * $self->{target_decoy_ratio};
	my $decoy_protID = $decoy_protID * $self->{target_decoy_ratio};

	# single hits
	my ( $P_s_fp, $P_s, $protFDRs, $FP_protIDs, $TP_protIDs ) =
	  ( 0, 0, 0, 0, 0 );    
	my $P_s_fp = $decoy_protIDs / $decoy_protID unless $decoy_protID == 0;
	my $P_s = $target_protIDs / $target_protID unless $target_protID == 0;

	# bayes
	$protFDRs   = $fdr * $P_s_fp / $P_s unless $P_s == 0;
	$FP_protIDs = $target_protIDs * $protFDRs;
	$TP_protIDs = $target_protIDs - $FP_protIDs;

	# all but single hits
	my ( $P_ns_fp, $P_ns, $protFDRns, $FP_protIDns, $TP_protIDns ) =
	  ( 0, 0, 0, 0, 0 );   
	my $P_ns_fp = $decoy_protIDns / $decoy_protID unless $decoy_protID == 0;
	my $P_ns = $target_protIDns / $target_protID unless $target_protID == 0; 

	# bayes
	$protFDRns   = $fdr * $P_ns_fp / $P_ns unless $P_ns == 0;
	$FP_protIDns = $target_protIDns * $protFDRns;
	$TP_protIDns = $target_protIDns - $FP_protIDns;

	# single hits
	$self->{summed}{single_hit}{'target_ids'} = $target_protIDs;
	$self->{summed}{single_hit}{'decoy_ids'}  = $decoy_protIDs;
	$self->{summed}{single_hit}{'FP'}         = $FP_protIDs;
	$self->{summed}{single_hit}{'TP'}         = $TP_protIDs;
	$self->{summed}{single_hit}{'FDR'}        = $protFDRs;

	# all but single hits
	$self->{summed}{no_single_hit}{'target_ids'} = $target_protIDns;
	$self->{summed}{no_single_hit}{'decoy_ids'}  = $decoy_protIDns;
	$self->{summed}{no_single_hit}{'FP'}         = $FP_protIDns;
	$self->{summed}{no_single_hit}{'TP'}         = $TP_protIDns;
	$self->{summed}{no_single_hit}{'FDR'}        = $protFDRns;
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
	my $self     = shift;
	my $ra_error = [];

	my $digits = $self->{output}{digits};
	my $format = "%." . $digits . "f";

	# single hits
	if ( exists( $self->{summed}{single_hit} ) ) {

		my $target_protIDs = $self->{summed}{single_hit}{'target_ids'};
		my $decoy_protIDs  = $self->{summed}{single_hit}{'decoy_ids'};
		my $FP_protIDs     = $self->{summed}{single_hit}{'FP'};
		my $TP_protIDs     = $self->{summed}{single_hit}{'TP'};
		my $protFDRs       = $self->{summed}{single_hit}{'FDR'};

		push @$ra_error, sprintf( "%.0f", $target_protIDs );
		push @$ra_error, sprintf( "%.0f", $decoy_protIDs );
		push @$ra_error, sprintf( "%.0f", $FP_protIDs );
		push @$ra_error, sprintf( "%.0f", $TP_protIDs );
		push @$ra_error, sprintf( $format, $protFDRs );
	}

	# all but single hits
	if ( exists( $self->{summed}{no_single_hit} ) ) {

		my $target_protIDns = $self->{summed}{no_single_hit}{'target_ids'};
		my $decoy_protIDns  = $self->{summed}{no_single_hit}{'decoy_ids'};
		my $FP_protIDns     = $self->{summed}{no_single_hit}{'FP'};
		my $TP_protIDns     = $self->{summed}{no_single_hit}{'TP'};
		my $protFDRns       = $self->{summed}{no_single_hit}{'FDR'};

		push @$ra_error, sprintf( "%.0f", $target_protIDns );
		push @$ra_error, sprintf( "%.0f", $decoy_protIDns );
		push @$ra_error, sprintf( "%.0f", $FP_protIDns );
		push @$ra_error, sprintf( "%.0f", $TP_protIDns );
		push @$ra_error, sprintf( $format, $protFDRns );
	}

	return $ra_error;
}

# Title	    :  get_summary_header()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_summary_header {
	my $self = shift;

	my $ra_header = [];

	# single hits
	$ra_header = $self->{output}{single_hit_header};

	# all but single hits
	foreach ( @{ $self->{output}{no_single_hit_header} } ) {
		push @$ra_header, $_;
	}

	return $ra_header;
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









