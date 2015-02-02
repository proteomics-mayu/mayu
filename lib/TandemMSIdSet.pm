package TandemMSIdSet;
use strict;

use MayuTools;

##################################################################
#
# TandemMSIdSet
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
# 2008.03.06
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $ts = TandemMSIdSet->new( $verbose, $status );
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

	my $tools = MayuTools->new();
	$self->{tools} = $tools;

	$self->{error_models} = [];
	$self->{attributes}   = {};

	return $self;
}

# Object Method
# Title	    :  add_error_models()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_error_models {
	my $self = shift;
	my ( $ra_em, $ra_ra_psm ) = @_;

	# keep the order for the printing out later
	for ( my $i = 0 ; $i < @$ra_em ; $i++ ) {
		my $em = $ra_em->[$i];
		unless ( $em->is_local_error_model() ) {
			$em->set_psm_set($ra_ra_psm);
			$self->{error_models}[$i] = $em;
		}
	}

	# check for local error models that need
	# FDRs from other error models
	# keep the order for the printing out later
	for ( my $i = 0 ; $i < @$ra_em ; $i++ ) {
		my $em = $ra_em->[$i];

		# e.g. protein local fdr
		if ( $em->is_local_error_model() ) {
			my $needed_id = $em->get_needed_error_model_identifier();
			my $id        = $em->get_error_model_identifier();

			# calculate the features for the local protein
			if ( $id eq 'LocalProteinIdFDR' ) {

				# get the protein features (nr_scans,...)
				# prot_id
				# nr_files
				# nr_runs
				# mFDR
				# decoy
				# NS
				# NP
				# PAT
				# PSL
				# acNTP
				my ( $rh_rh_id_protein_features,
					$ra_sorted_protein_id_feature_header )
				  = $self->get_protein_id_features($ra_ra_psm);

				if ( $self->{protein_feature_file} ) {
					$self->print_protein_feature_file(
						$rh_rh_id_protein_features,
						$ra_sorted_protein_id_feature_header );
				}

				# search the corresponding global FDR
				my $fdr = $self->get_global_fdr_by_em_id($needed_id);

				$em->set_psm_set( $fdr, $rh_rh_id_protein_features );
				$self->{error_models}[$i] = $em;
			}
		}
	}
}

# Object Method
# Title	    :  set_attribute()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_attribute {
	my $self = shift;
	my ( $name, $value ) = @_;

	if ( defined($name) && defined($value) ) {
		$self->{attributes}{$name} = $value;
	}
	else {
		print "  name or value not defined!\n";
	}
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
# Title	    :  set_protein_feature_file()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set_protein_feature_file {
	my $self = shift;
	my ($out_base) = @_;

	$self->{protein_feature_file} = $out_base;
}

# Object Method
# Title	    :  get_attributes()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_attributes {
	my $self = shift;
	return $self->{attributes};
}

# Object Method
# Title	    :  get_attribute_values()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_attribute_values {
	my $self = shift;

	my $ra_values = [];

	foreach my $att_name ( sort keys %{ $self->{attributes} } ) {
		push @$ra_values, $self->{attributes}{$att_name};
	}

	return $ra_values;
}

# Object Method
# Title	    :  get_attribute_names()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_attribute_names {
	my $self = shift;

	my $ra_names = [];

	foreach my $att_name ( sort keys %{ $self->{attributes} } ) {
		push @$ra_names, $att_name;
	}

	return $ra_names;
}

# Object Method
# Title	    :  get_global_fdr_by_em_id()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_global_fdr_by_em_id {
	my $self = shift;
	my ($id) = @_;

	my $fdr = 1;

	foreach my $em ( @{ $self->{error_models} } ) {
		if ( $em->get_error_model_identifier() eq $id ) {
			$fdr = $em->get_fdr();
		}
	}

	return $fdr;
}

# Object Method
# Title	    :  get_protein_id_features()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_protein_id_features {
	my $self = shift;
	my ($ra_ra_psm) = @_;

	# prot_id
	# nr_files
	# nr_runs
	# mFDR
	# decoy
	# nr_scans
	# nr_pep
	# sot
	# aa
	# ac_ntp
	my $rh_rh_id_protein_features = {};

	# first pass through the data create a useful data structure
	my $rh_rh_rh_pid_pep_scan_count = {};
	my $rh_pid_decoy                = {};
	foreach my $ra_psm (@$ra_ra_psm) {
		$rh_rh_rh_pid_pep_scan_count->{ $ra_psm->[2] }{ $ra_psm->[1] }
		  { $ra_psm->[0] }++;
		$rh_pid_decoy->{ $ra_psm->[2] } = $ra_psm->[5];
	}

	foreach my $pid ( keys %$rh_rh_rh_pid_pep_scan_count ) {

		# this was already calculated while estimating
		# protein size
		my ( $aa, $ac_ntp ) = ( 0, 0 );
		if ( exists( $self->{protein_features} ) ) {
			my $pfeat = $self->{protein_features};
			$aa     = $pfeat->get_feature_by_id( $pid, 'aa' );
			$ac_ntp = $pfeat->get_feature_by_id( $pid, 'corr_ntp' );
		}

		my ( $nr_scans, $sot, $nr_pep ) =
		  $self->get_prot_scan_counts(
			$rh_rh_rh_pid_pep_scan_count->{$pid} );

		# these keys are used by the LocalProteinIdFDR module!
		# don't change the keys!
		my $rh_features = {
			'decoy' => $rh_pid_decoy->{$pid},
			'NS'    => $nr_scans,
			'NP'    => $nr_pep,
			'PAT'   => $sot,
			'PSL'   => $aa,
			'acNTP' => $ac_ntp
		};
		$rh_rh_id_protein_features->{$pid} = $rh_features;
	}

	my @features_header = ( 'decoy', 'NS', 'NP', 'PAT', 'PSL', 'acNTP' );
	my @sorted_features_header = sort @features_header;

	return ( $rh_rh_id_protein_features, \@sorted_features_header );
}

# Object Method
# Title	    :  get_prot_scan_counts()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_prot_scan_counts {
	my $self = shift;
	my ($rh_rh_pep_scan) = @_;

	my ( $nr_scans, $sot, $nr_pep ) = ( 0, 0, 0 );

	# number of distinct peptides and scans
	my $nr_pep = keys %$rh_rh_pep_scan;
	foreach my $pep ( keys %$rh_rh_pep_scan ) {
		foreach my $scan (keys %{ $rh_rh_pep_scan->{$pep} }) {
			$nr_scans += $rh_rh_pep_scan->{$pep}{$scan};
		}
	}

	# 7 distinct cases of scan occupation
	# 0: one scan
	# 1: two scans, one pep
	# 2: two scans, two pep
	# 3: three scans, one pep
	# 4: three scans, two pep
	# 5: three scans, three pep
	# 6: all the rest
	if ( $nr_scans == 1 ) {
		$sot = 0;
	}
	elsif ( $nr_scans == 2 ) {
		if ( $nr_pep == 1 ) {
			$sot = 1;
		}
		elsif ( $nr_pep == 2 ) {
			$sot = 2;
		}
	}
	elsif ( $nr_scans == 3 ) {
		if ( $nr_pep == 1 ) {
			$sot = 3;
		}
		elsif ( $nr_pep == 2 ) {
			$sot = 4;
		}
		else {
			$sot = 5;
		}
	}

	return ( $nr_scans, $sot, $nr_pep );
}

# Object Method
# Title	    :  get_error_bin_table_by_em_id()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_error_bin_table_by_em_id {
	my $self = shift;
	my ($error_model_identifier) = @_;

	my $ra_ra_bin_table = [];

	foreach my $em ( @{ $self->{error_models} } ) {
		if ( $em->get_error_model_identifier() eq $error_model_identifier )
		{
			$ra_ra_bin_table = $em->get_error_bin_table();
		}
	}

	return $ra_ra_bin_table;
}

# Object Method
# Title	    :  get_row()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_row {
	my $self = shift;

	my $ra_full_row = [];

	foreach my $em ( @{ $self->{error_models} } ) {
		my $ra_row = $em->get_error_summary();
		foreach (@$ra_row) {
			push @$ra_full_row, $_;
		}
	}

	return $ra_full_row;
}

# Object Method
# Title	    :  print_protein_feature_file()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub print_protein_feature_file {
	my $self = shift;
	my ( $rh_rh_id_protein_features, $ra_sorted_protein_id_feature_header )
	  = @_;

	my $header  = 1;
	my $csv_sep = "\t";

	my $out_base      = $self->{protein_feature_file};
	my $prot_feat_csv = $out_base . '.csv';

	# if the file exists already the header was already printed
	# and data should just be appended
	my $print_header = 0 if ( -e $prot_feat_csv );
	my $print_header = 1 unless ( -e $prot_feat_csv );

	my $fh_prot_feat_csv = FileHandle->new();
	$fh_prot_feat_csv->open(">>$prot_feat_csv") or die $!;

	my @pid = sort keys %{$rh_rh_id_protein_features};

	# header
	if ($print_header) {
		my $ra_header = ['id'];
		foreach ( @{ $self->get_attribute_names() } ) {
			push @$ra_header, $_;
		}
		foreach my $key (@$ra_sorted_protein_id_feature_header) {
			push @$ra_header, $key;
		}
		print $fh_prot_feat_csv join( $csv_sep, @$ra_header ) . "\n";
	}

	# data
	foreach my $pid (@pid) {
		my $ra_row = [$pid];
		foreach ( @{ $self->get_attribute_values() } ) {
			push @$ra_row, $_;
		}
		my @sorted_keys =
		  sort keys %{ $rh_rh_id_protein_features->{$pid} };
		foreach my $key (@sorted_keys) {
			push @$ra_row, $rh_rh_id_protein_features->{$pid}{$key};
		}
		print $fh_prot_feat_csv join( $csv_sep, @$ra_row ) . "\n";
	}

	$fh_prot_feat_csv->close();
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







