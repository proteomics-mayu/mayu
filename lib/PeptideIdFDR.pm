package PeptideIdFDR;
use strict;

use MayuTools;
use HypergProbDist;

##################################################################
#
# PeptideIdFDR
#
# NOTES:
# All the error model classes have to implemement certain methods:
# - get_error_model_identifier
# - create_error_model
#
#
# SYNOPSIS:
#
# TODO:
# -
# -
#
# CHANGELOG:
# 2008.03.05
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $em = PeptideIdFDR->new( $v, $s, $bin_ob,
#                         $target_decoy_ratio );
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $verbose, $status, $binning_object, $target_decoy_ratio ) = @_;

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

	# binning object
	if ( defined($binning_object) ) {
		$self->{binning_ob} = $binning_object;
	}
	else {
		$self->{binning_ob} = undef;
		print "  need a binning object as input!\n";
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
	$self->{error_model_identifier} = 'PeptideIdFDR';

	# used for processed printing out
	$self->{output}{digits}     = 8;
	$self->{output}{bin_header} = [
		'peptide_bin',    'bin_desc',
		'target_pep',     'target_pepID',
		'decoy_pepID',    'FP_pepID',
		'FP_pepID_stdev', 'TP_pepID',
		'pepFDR'
	];
	$self->{output}{summary_header} = [
		'target_pepID', 'decoy_pepID',
		'FP_pepID',     'FP_pepID_stdev',
		'TP_pepID',     'pepFDR'
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
# Function	:  creates a copy object using the same binning ob
# Returns   :
# Args      :
sub create_error_model {
	my $self = shift;

	my $em =
	  PeptideIdFDR->new( $self->{v}, $self->{s}, $self->{binning_ob},
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

	# the calculation should be done binned
	# total number of peptide ids in db has to be extracted
	# out of the binning object
	if ( defined( $self->{binning_ob} ) ) {
		$self->estimate_error($ra_ra_psm);
	}
	else {
		print "  no binning object found. Cannot calculate errors!\n";
	}
}

# Object Method
# Title	    :  estimate_error()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub estimate_error {
	my $self = shift;
	my ($ra_ra_psm) = @_;
	
	my ( $tot_tar, $tot_dec, $tot_fp, $tot_tp, $tot_fp_var, $tot_fdr ) =
	  ( 0, 0, 0, 0, 0, 0 );
	my $binned_entity = $self->{binning_ob};

	my $ra_target_ids = $self->get_target_peptide_ids($ra_ra_psm);
	my $ra_decoy_ids  = $self->get_decoy_peptide_ids($ra_ra_psm);
	my $nr_bins       = $binned_entity->get_nr_bins();
	my $tot_tar       = @$ra_target_ids;
	my $tot_dec       = @$ra_decoy_ids;

	# loop over the peptide bins
	for ( my $i = 0 ; $i < $nr_bins ; $i++ ) {
		
		my $bin = $binned_entity->get_bin_by_index($i);

		my $bin_desc        = $bin->get_nice_test_string();
		my $target_entities = $bin->get_att_by_name('entities');

		# estimate the number of entities that contain false positives
		my $nr_cf_ids = $tot_dec * $self->{target_decoy_ratio};		
		my ( $fp, $fp_var ) =
		  $self->get_fp_exp_and_var( $target_entities, $tot_tar,
			$nr_cf_ids );
		
		my $tp  = $tot_tar - $fp;
		my $fdr = 0;
		$fdr = $fp / $tot_tar unless $tot_tar == 0;
		$tot_fp     += $fp;
		$tot_tp     += $tp;
		$tot_fp_var += $fp_var;

		$self->{bins}[$i]{'bin_nr'}         = $i;
		$self->{bins}[$i]{'bin_desc'}       = $bin_desc;
		$self->{bins}[$i]{'target_entries'} = $target_entities;
		$self->{bins}[$i]{'target_ids'}     = $tot_tar;
		$self->{bins}[$i]{'decoy_ids'}      = $tot_dec;
		$self->{bins}[$i]{'FP'}             = $fp;
		$self->{bins}[$i]{'FP_stdev'}       = sqrt($fp_var);
		$self->{bins}[$i]{'TP'}             = $tp;
		$self->{bins}[$i]{'FDR'}            = $fdr;
	}
	$tot_fdr = $tot_fp / $tot_tar unless $tot_tar == 0;
	$self->{summed}{'target_ids'} = $tot_tar;
	$self->{summed}{'decoy_ids'}  = $tot_dec;
	$self->{summed}{'FP'}         = $tot_fp;
	$self->{summed}{'FP_stdev'}   = sqrt($tot_fp_var);
	$self->{summed}{'TP'}         = $tot_tp;
	$self->{summed}{'FDR'}        = $tot_fdr;
}

# Object Method
# Title	    :  get_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_fdr {
	my $self = shift;

	if ( exists( $self->{summed}{'FDR'} ) ) {
		if ( defined( $self->{summed}{'FDR'} ) ) {
			my $fdr = $self->{summed}{'FDR'};
			$fdr = 0 if $fdr < 0;
			$fdr = 1 if $fdr > 1;
			return $fdr;
		}
	}
	print "  no FDR found!\n";
	return 1;
}

# Object Method
# Title	    :  get_error_bin_table()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_error_bin_table {
	my $self = shift;

	my $ra_ra_bin_table = [];
	
	# if the errors were never calculated
	return $ra_ra_bin_table unless exists( $self->{bins} );

	my $digits = $self->{output}{digits};
	my $format = "%." . $digits . "f";

	foreach my $bin ( @{ $self->{bins} } ) {
		my $bin_index      = $bin->{'bin_nr'};
		my $bin_desc       = $bin->{'bin_desc'};
		my $target_entries = $bin->{'target_entries'};
		my $target_ids     = $bin->{'target_ids'};
		my $decoy_ids      = $bin->{'decoy_ids'};
		my $fp             = sprintf( "%.0f", $bin->{'FP'} );
		my $fp_stdev       = sprintf( $format, $bin->{'FP_stdev'} );
		my $tp             = sprintf( "%.0f", $bin->{'TP'} );
		my $fdr            = sprintf( $format, $bin->{'FDR'} );
		$fp  = 0 if $fp < 0;
		$tp  = 0 if $tp < 0;
		$fdr = 0 if $fdr < 0;
		$fdr = 1 if $fdr > 1;

		push @$ra_ra_bin_table,
		  [
			$bin_index,  $bin_desc,  $target_entries,
			$target_ids, $decoy_ids, $fp,
			$fp_stdev,   $tp,        $fdr
		  ];
	}

	return $ra_ra_bin_table;
}

# Object Method
# Title	    :  get_error_summary()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_error_summary {
	my $self = shift;
	
	# if the errors were never calculated
	return [] unless exists( $self->{summed} );

	my $digits = $self->{output}{digits};
	my $format = "%." . $digits . "f";

	my $tar      = $self->{summed}{'target_ids'};
	my $dec      = $self->{summed}{'decoy_ids'};
	my $fp       = sprintf( "%.0f", $self->{summed}{'FP'} );
	my $fp_stdev = sprintf( $format, $self->{summed}{'FP_stdev'} );
	my $tp       = sprintf( "%.0f", $self->{summed}{'TP'} );
	my $fdr      = sprintf( $format, $self->{summed}{'FDR'} );

	$fp  = 0 if $fp < 0;
	$tp  = 0 if $tp < 0;
	$fdr = 0 if $fdr < 0;
	$fdr = 1 if $fdr > 1;

	my $ra_error = [ $tar, $dec, $fp, $fp_stdev, $tp, $fdr ];

	return $ra_error;
}

# Object Method
# Title	    :  get_target_peptide_ids()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_target_peptide_ids {
	my $self = shift;
	my ($ra_ra_psm) = @_;

	my $rh_peptide_ids;
	foreach my $ra_psm (@$ra_ra_psm) {
		if ( $ra_psm->[5] == 0 ) {
			$rh_peptide_ids->{ $ra_psm->[1] }++;
		}
	}
	my @peptide_ids = sort keys %$rh_peptide_ids;

	return \@peptide_ids;
}

# Object Method
# Title	    :  get_decoy_peptide_ids()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_decoy_peptide_ids {
	my $self = shift;
	my ($ra_ra_psm) = @_;

	my $rh_peptide_ids;
	foreach my $ra_psm (@$ra_ra_psm) {
		if ( $ra_psm->[5] == 1 ) {
			$rh_peptide_ids->{ $ra_psm->[1] }++;
		}
	}
	my @peptide_ids = sort keys %$rh_peptide_ids;

	return \@peptide_ids;
}

# Object Method
# Title	    :  set_data()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_fp_exp_and_var {
	my $self = shift;
	my ( $N, $t, $cf ) = @_;
	
	# bad estimate of search space
	# e.g. semitryptic search
	if ( $t > $N || $cf > $N ) {
		return( $t, 0 );
	}

	my $hpd     = HypergProbDist->new();
	my @fp_prob = ();
	for ( my $fp = 0 ; $fp <= $cf ; $fp++ ) {
		my $tp = $t - $fp;
		my $w  = $N - $tp;
		my $prob = $hpd->check_n_hyperg( $fp, $N, $w, $cf );
		push @fp_prob, $prob;
	}
	
	my $tools = $self->{tools};

	# normalize such that the total prob will be one
	my $sum = $tools->get_sum( \@fp_prob );
	@fp_prob = map( $tools->div( $_, $sum ), @fp_prob );

	my $fp_exp = $tools->get_exp_val( \@fp_prob );
	my $fp_var = $tools->get_var( \@fp_prob, $fp_exp );
	
	return ( $fp_exp, $fp_var );
}

# Title	    :  get_bin_header()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_bin_header {
	my $self = shift;
	return $self->{output}{bin_header};
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

# TODO remove?
sub print_object {    # print the object to the console with dumpvar
	my $self = shift();
	if ( defined( $_[0] ) ) {

		package main;
		require "dumpvar.pl";
		dumpValue( $_[0] );
	}
	else {

		package main;
		require "dumpvar.pl";
		dumpValue($self);
	}
}

1;

__END__

Lukas Reiter






