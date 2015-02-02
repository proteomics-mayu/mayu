package MayuManager;
use strict;

use TandemMSIdSet;

##################################################################
#
# MayuManager
#
# NOTES:
#
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
# Usage     :  my $pm = MayuManager->new( $verbose, $status );
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

	# Saves error models in an array. These error models are
	# used to create error models of the same type for all
	# data selections delivered with TandemMSIdSelectionScheme
	$self->{error_model_proxys} = [];
	$self->{id_sets}            = [];

	return $self;
}

# Object Method
# Title	    :  set_error_model()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_error_model {
	my $self = shift;
	my ($error_model) = @_;

	if ( defined($error_model) ) {
		push @{ $self->{error_model_proxys} }, $error_model;
	}
	else {
		print "  input errror model not defined!\n";
	}
}

# Object Method
# Title	    :  print_registered_error_models()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub print_registered_error_models {
	my $self = shift;

	print "  registered error models:\n";
	my @em_identifiers = ();
	foreach my $em ( @{ $self->{error_model_proxys} } ) {
		push @em_identifiers, $em->get_error_model_identifier();
	}
	print "  " . join( ", ", @em_identifiers ) . "\n";
}

# Object Method
# Title	    :  set_protein_features()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_protein_features {
	my $self = shift;
	my ($prot_feat) = @_;

	$self->{protein_features} = $prot_feat;
}

# Object Method
# Title	    :  set_protein_features()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_protein_feature_output_base {
	my $self = shift;
	my ($out_base) = @_;

	$self->{protein_feature_file} = $out_base;
}

# Object Method
# Title	    :  set_psm_set()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_psm_set {
	my $self = shift;
	my ( $ra_ra_psm, $rh_attributes ) = @_;

	my $id_set = TandemMSIdSet->new( $self->{v}, $self->{s} );

	# binned_entity object (this might be changed in the future)
	if ( exists( $self->{protein_features} ) ) {
		$id_set->set_protein_features( $self->{protein_features} );
	}
	if ( exists( $self->{protein_feature_file} ) ) {
		$id_set->set_protein_feature_file( $self->{protein_feature_file} );
	}

	# add attributes like to the TandemMSIdSet
	# - nr_files
	# - nr_runs
	# - psm_fdr
	# - ...
	foreach my $name ( keys %$rh_attributes ) {
		my $value = $rh_attributes->{$name};
		$id_set->set_attribute( $name, $value );
	}

	# add the error models to the TandemMSIdSet
	my @em = ();
	foreach my $error_model ( @{ $self->{error_model_proxys} } ) {

		# creates a new instance of this error model
		# all error models have to have implemented this
		# method
		my $em = $error_model->create_error_model();
		push @em, $em;
	}
	$id_set->add_error_models( \@em, $ra_ra_psm );

	# add the TandemMSIdSet
	push @{ $self->{id_sets} }, $id_set;
}

# Object Method
# Title	    :  error_model_exists()
# Usage     :
# Function	:  checks whether the required error model
#              is registered
# Returns   :
# Args      :
sub error_model_exists {
	my $self = shift;
	my ($error_model_identifier) = @_;

	my $exists = 0;

	foreach my $em ( @{ $self->{error_model_proxys} } ) {
		if ( $em->get_error_model_identifier() eq $error_model_identifier )
		{
			$exists = 1;
		}
	}

	return $exists;
}

# Object Method
# Title	    :  get_prot_bin_table()
# Usage     :
# Function	:
# Returns   :  a table with following columns:
#              -
# Args      :
sub get_prot_bin_table {
	my $self = shift;
	my ($header) = @_;

	my $ra_ra_table = [];

	foreach my $id_set ( @{ $self->{id_sets} } ) {

		my $ra_att_values = $id_set->get_attribute_values();

		# get binned errors by error model name
		# this does only make sense for error models that bin
		# the data in order to calculate the error model
		my $ra_ra_data =
		  $id_set->get_error_bin_table_by_em_id('ProteinIdFDR');
		foreach my $ra_row (@$ra_ra_data) {
			foreach my $att_value (@$ra_att_values) {
				unshift @$ra_row, $att_value;
			}
			push @$ra_ra_table, $ra_row;
		}
	}

	if ($header) {
		unshift @$ra_ra_table, $self->get_prot_bin_header();
	}

	return $ra_ra_table;
}

# Object Method
# Title	    :  get_table()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_table {
	my $self = shift;
	my ($header) = @_;

	my $ra_ra_table = [];

	foreach my $id_set ( @{ $self->{id_sets} } ) {

		my $ra_att_values = $id_set->get_attribute_values();

		# table with all condensed information
		my $ra_row = $id_set->get_row();
		foreach my $att_value (@$ra_att_values) {
			unshift @$ra_row, $att_value;
		}
		push @$ra_ra_table, $ra_row;
	}

	if ($header) {
		unshift @$ra_ra_table, $self->get_header();
	}
    
	return $ra_ra_table;
}

# Object Method
# Title	    :  get_prot_bin_header()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_prot_bin_header {
	my $self = shift;

	my $em        = $self->get_error_model('ProteinIdFDR');
	my $ra_header = $em->get_bin_header();

	my $ra_att_names = $self->{id_sets}[0]->get_attribute_names();
	foreach my $att_name (@$ra_att_names) {
		unshift @$ra_header, $att_name;
	}

	return $ra_header;
}

# Object Method
# Title	    :  get_header()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_header {
	my $self = shift;

	my $ra_header = [];

	# id set attributes
	my $ra_att_names = $self->{id_sets}[0]->get_attribute_names();
	foreach my $att_name (@$ra_att_names) {
		unshift @$ra_header, $att_name;
	}

	# columns from the error models
	foreach my $em ( @{ $self->{error_model_proxys} } ) {
		my $ra_single_header = $em->get_summary_header();
		foreach (@$ra_single_header) {
			push @$ra_header, $_;
		}
	}

	return $ra_header;
}

# Object Method
# Title	    :  get_error_model()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_error_model {
	my $self = shift;
	my ($error_model_identifier) = @_;

	foreach my $em ( @{ $self->{error_model_proxys} } ) {
		if ( $em->get_error_model_identifier() eq $error_model_identifier )
		{
			return $em;
		}
	}
	print "  could not find the error model!\n";
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






