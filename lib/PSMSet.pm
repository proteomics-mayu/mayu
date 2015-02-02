package PSMSet;
use strict;

use File::Basename;

use Scan;
use FastaParser;

##################################################################
#
# PSMSet
#
# NOTES:
# a PSMSet is a space efficient way to store PSM. Data is stored
# in an array of arrays $self->{data}{psm}:
#
# if the input of the PSMSet is to be extended then a class
# variable for the number of input features can be stored.
# This can then be used to return only the standard set of
# fields for each feature to keep the functionality of the
# already implemented functions the same.
#
#
# 0 spectrum
# 1 peptide
# 2 main id
# 3 modification info
# 4 PeptideProphet probability (PPs)
# 5 1: decoy, 0: target
# 6 PSM FDR estimate
#
# SYNOPSIS:
#
# TODO:
# - do not loose the association data to input file
#   do so with a run to input file hash table
#   use scan to extract the run
# - create an object for the ROC
#
# CHANGELOG:
# 2008.02.20
# - added get_psm_by_fdr()
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $s = PSMSet->new( $verbose, $status );
# Function  :
# Returns   :
# Args      :  $verbose:   print out some stuff or not (1 or 0)
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $verbose, $status, $ra_ra_psm, $file ) = @_;

	# delta mass (Da) to test for modifications
	$self->{mod_test_delta} = 0.1;

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

	# initialize
	# the runs, nr_runs, run_file,...
	# correspond to the unfiltered data!
	$self->{data}{psm}      = [];
	$self->{runs}           = {};
	$self->{nr_runs}        = 0;
	$self->{run_file}       = {};
	$self->{file}           = [];
	$self->{base_file}      = [];
	$self->{nr_input_files} = 0;

	if ( defined($ra_ra_psm) ) {
		my $nr = @$ra_ra_psm;
		if ( $nr > 0 ) {
			$self->{data}{psm} = $ra_ra_psm;

			my $scan = Scan->new();
			my %runs;
			foreach $_ ( @{$ra_ra_psm} ) {
				$runs{ $scan->get_run( $_->[0] ) }++;
			}
			my $nr_runs = keys %runs;
			$self->{runs}    = \%runs;
			$self->{nr_runs} = $nr_runs;

			# run to file association
			if ( defined($file) ) {
				foreach my $run ( keys %runs ) {
					$self->{run_file}{$run} = $file;
				}
			}
		}
	}

	# data is derived from one input file
	if ( defined($file) ) {
		push @{ $self->{file} },      $file;
		push @{ $self->{base_file} }, basename($file);
		$self->{nr_input_files} = 1;
	}

	# output default settings
	# 0  ds cutoff
	# 1  PSM FDR
	# 2  FP PSM
	# 3  TP PSM
	# 4  td PSM FDR
	# 5  td FP PSM
	# 6  td TP PSM
	# 7  target PSM
	# 8  decoy PSM
	$self->{roc_header} = [
		'IP/PPs',     'mFDR',  'FP',    'TP',
		'TD_mFDR', 'TD_FP', 'TD_TP', 'target_PSM',
		'decoy_PSM'
	];
	$self->{digits} = 8;

	return $self;
}

# Object Method
# Title	    :  add_psm()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_psm {
	my $self = shift();
	my ( $ra_ra_psm, $file, $rh_runs, $nr_input_runs ) = @_;

	# the runs, nr_runs, run_file,...
	# correspond to the unfiltered data!
	if ( defined($ra_ra_psm) ) {

		$self->add_ra_ra_psm($ra_ra_psm);

		my @runs              = ();
		my $nr_additionl_runs = 0;

		# handle the run names:
		# if the run names are provided use them
		# else use the $ra_ra_psm
		if ( defined($rh_runs) ) {
			@runs = keys %$rh_runs;
		}
		else {
			my $scan = Scan->new();
			my %runs;
			foreach $_ ( @{$ra_ra_psm} ) {
				$runs{ $scan->get_run( $_->[0] ) }++;
			}
			@runs = keys %runs;
		}

		# number of runs
		# if the number of new runs is passed then
		# the number of runs is updated with this number
		# no matter whether there is already identical runs.
		#
		# This can make sense if not really runs but
		# searched runs have to be counted.
		if ( defined($nr_input_runs) ) {
			$nr_additionl_runs = $nr_input_runs;
		}
		else {
			$nr_additionl_runs = @runs;
		}

		# update the information run information
		$self->{nr_runs} += $nr_additionl_runs;
		foreach (@runs) {
			$self->{runs}{$_}++;
		}

		if ( defined($file) ) {

			# run to file association
			if ( defined($file) ) {
				foreach my $run (@runs) {
					$self->{run_file}{$run} = $file;
				}
			}
			push @{ $self->{file} },      $file;
			push @{ $self->{base_file} }, basename($file);
			$self->{nr_input_files}++;
		}
	}
	else {
		print "  add_psm(): no PSM defined!\n";
	}
}

# Object Method
# Title	    :  add_ra_ra_psm()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_ra_ra_psm {
	my $self = shift();
	my ($ra_ra_psm) = @_;

	foreach (@$ra_ra_psm) {
		push @{ $self->{data}{psm} }, $_;
	}
}

# Object Method
# Title	    :  get_psm()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_psm {
	my $self = shift();

	if ( exists( $self->{data}{psm} ) ) {
		return $self->{data}{psm};
	}
	else {
		return [];
	}
}

# Object Method
# Title	    :  get_target_psm()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_target_psm {
	my $self = shift();
	
	my $ra_psm = $self->get_psm();
	
	my $ra_target_psm = [];
	foreach (@$ra_psm) {
		if ( defined( $_->[5] ) ) {
			push @$ra_target_psm, $_ if $_->[5] == 0;
		}
	}
	
	return $ra_target_psm;
}

# Object Method
# Title	    :  get_decoy_psm()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_decoy_psm {
	my $self = shift();
	
	my $ra_psm = $self->get_psm();
	
	my $ra_decoy_psm = [];
	foreach (@$ra_psm) {
		if ( defined( $_->[5] ) ) {
			push @$ra_decoy_psm, $_ if $_->[5] == 1;
		}
	}
	
	return $ra_decoy_psm;
}

# Object Method
# Title	    :  get_psm_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_psm_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	# if not defined return all the data
	$fdr = 1 unless defined( $fdr );
	
	my @filt_psm = ();

	if ( exists( $self->{data}{psm} ) ) {
		foreach my $ra_psm ( @{ $self->{data}{psm} } )
		{
			if ( defined( $ra_psm->[6] ) ) {
				if ( $ra_psm->[6] <= $fdr ) {
					push @filt_psm, $ra_psm;
				}
			}
		}
	}

	return \@filt_psm;
}

# Object Method
# Title	    :  get_psm_by_files_and_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_psm_by_files_and_fdr {
	my $self = shift();
	my ( $ra_files, $fdr ) = @_;

	my @filt_psm = ();

	# gather the corresponding runs
	my $rh_runs;
	foreach my $run ( keys %{ $self->{run_file} } ) {
		my $file = $self->{run_file}{$run};
		foreach my $wanted_file (@$ra_files) {
			if ( $file eq $wanted_file ) {
				$rh_runs->{$run}++;
			}
		}
	}

	# first assemble a list of runs
	my $scan = Scan->new();

	if ( exists( $self->{data}{psm} ) ) {
		foreach my $ra_psm ( @{ $self->{data}{psm} } ) {
			if ( defined( $ra_psm->[6] )
				&& exists( $rh_runs->{ $scan->get_run( $ra_psm->[0] ) } ) )
			{
				if ( $ra_psm->[6] <= $fdr ) {
					push @filt_psm, $ra_psm;
				}
			}
		}
	}

	return \@filt_psm;
}

# Object Method
# Title	    :  get_target_psm_by_files_and_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_target_psm_by_files_and_fdr {
	my $self = shift();
	my ( $ra_files, $fdr ) = @_;

	my @filt_psm = ();

	# gather the corresponding runs
	my $rh_runs;
	foreach my $run ( keys %{ $self->{run_file} } ) {
		my $file = $self->{run_file}{$run};
		foreach my $wanted_file (@$ra_files) {
			if ( $file eq $wanted_file ) {
				$rh_runs->{$run}++;
			}
		}
	}

	# first assemble a list of runs
	my $scan = Scan->new();

	if ( exists( $self->{data}{psm} ) ) {
		foreach my $ra_psm ( @{ $self->{data}{psm} } ) {
			if ( defined( $ra_psm->[6] )
				&& exists( $rh_runs->{ $scan->get_run( $ra_psm->[0] ) } ) )
			{
				if ( $ra_psm->[6] <= $fdr ) {
					if ( defined( $ra_psm->[5] ) ) {
						push @filt_psm, $ra_psm if $ra_psm->[5] == 0;
					}
				}
			}
		}
	}

	return \@filt_psm;
}

# Object Method
# Title	    :  get_psm_by_run_and_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_psm_by_run_and_fdr {
	my $self = shift();
	my ( $ra_runs, $fdr ) = @_;

	my @filt_psm = ();

	my $rh_runs = {};
	foreach my $run (@$ra_runs) {
		$rh_runs->{$run}++;
	}

	# first assemble a list of runs
	my $scan = Scan->new();

	if ( exists( $self->{data}{psm} ) ) {
		foreach my $ra_psm ( @{ $self->{data}{psm} } ) {
			if ( defined( $ra_psm->[6] )
				&& exists( $rh_runs->{ $scan->get_run( $ra_psm->[0] ) } ) )
			{
				if ( $ra_psm->[6] <= $fdr ) {
					push @filt_psm, $ra_psm;
				}
			}
		}
	}

	return \@filt_psm;
}

# Object Method
# Title	    :  get_target_psm_by_run_and_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_target_psm_by_run_and_fdr {
	my $self = shift();
	my ( $ra_runs, $fdr ) = @_;

	my @filt_psm = ();

	my $rh_runs = {};
	foreach my $run (@$ra_runs) {
		$rh_runs->{$run}++;
	}

	# first assemble a list of runs
	my $scan = Scan->new();

	if ( exists( $self->{data}{psm} ) ) {
		foreach my $ra_psm ( @{ $self->{data}{psm} } ) {
			if ( defined( $ra_psm->[6] )
				&& exists( $rh_runs->{ $scan->get_run( $ra_psm->[0] ) } ) )
			{
				if ( $ra_psm->[6] <= $fdr ) {
					if ( defined( $ra_psm->[5] ) ) {
						push @filt_psm, $ra_psm if $ra_psm->[5] == 0;
					}
				}
			}
		}
	}

	return \@filt_psm;
}

# Object Method
# Title	    :  get_mod_psm_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_mod_psm_by_fdr {
	my $self = shift();
	my ( $fdr, $rh_mod_desc ) = @_;
	my @mod_psm = ();

	if ( exists( $self->{data}{psm} ) ) {
		my $all_nr = @{ $self->{data}{psm} };
		foreach my $ra_psm ( @{ $self->{data}{psm} } ) {
			if ( defined( $ra_psm->[6] ) && defined( $ra_psm->[3] ) ) {
				if (
					$ra_psm->[6] <= $fdr
					&& $self->is_mod( $ra_psm->[1], $ra_psm->[3],
						$rh_mod_desc ) == 1
				  )
				{
					push @mod_psm, $ra_psm;
				}
			}
		}
	}

	return \@mod_psm;
}

# Object Method
# Title	    :  get_target_psm_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_target_psm_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my @target_psm = ();
	my $ra_ra_psm  = $self->get_psm_by_fdr($fdr);

	foreach my $ra_psm (@$ra_ra_psm) {
		if ( defined( $ra_psm->[5] ) ) {
			push @target_psm, $ra_psm if $ra_psm->[5] == 0;
		}
	}

	return \@target_psm;
}

# Object Method
# Title	    :  get_decoy_psm_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_decoy_psm_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my @decoy_psm = ();
	my $ra_ra_psm = $self->get_psm_by_fdr($fdr);

	foreach my $ra_psm (@$ra_ra_psm) {
		if ( defined( $ra_psm->[5] ) ) {
			push @decoy_psm, $ra_psm if $ra_psm->[5] == 1;
		}
	}

	return \@decoy_psm;
}

# Object Method
# Title	    :  get_ds_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_ds_by_fdr {   
    my $self = shift();
	my ($fdr) = @_;

	my $ra_ra_target_psm = $self->get_target_psm_by_fdr($fdr);
    my $smallest = 1;
    foreach my $ra_psm ( @$ra_ra_target_psm )
	{
		if ( defined( $ra_psm->[4] ) && $ra_psm->[4] < $smallest ) {
            $smallest = $ra_psm->[4];
		}
	}

	return $smallest;
}

# Object Method
# Title	    :  get_nr_target_psm_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_target_psm_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my $ra_ra_target_psm = $self->get_target_psm_by_fdr($fdr);
	my $nr_target_psm    = @$ra_ra_target_psm;

	return $nr_target_psm;
}

# Object Method
# Title	    :  get_nr_decoy_psm_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_decoy_psm_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my $ra_ra_decoy_psm = $self->get_decoy_psm_by_fdr($fdr);
	my $nr_decoy_psm    = @$ra_ra_decoy_psm;

	return $nr_decoy_psm;
}

# Object Method
# Title	    :  get_nr_psm_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_psm_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my $ra_ra_psm = $self->get_psm_by_fdr($fdr);
	my $nr_psm    = @$ra_ra_psm;

	return $nr_psm;
}

# Object Method
# Title	    :  get_target_pep_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_target_pep_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my %target_pep;
	my $ra_ra_psm = $self->get_psm_by_fdr($fdr);

	foreach my $ra_psm (@$ra_ra_psm) {
		if ( defined( $ra_psm->[5] ) ) {
			$target_pep{ $ra_psm->[1] }++ if $ra_psm->[5] == 0;
		}
	}
	my @target_pep = keys %target_pep;
	return \@target_pep;
}

# Object Method
# Title	    :  get_decoy_pep_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_decoy_pep_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my %decoy_pep;
	my $ra_ra_psm = $self->get_psm_by_fdr($fdr);

	foreach my $ra_psm (@$ra_ra_psm) {
		if ( defined( $ra_psm->[5] ) ) {
			$decoy_pep{ $ra_psm->[1] }++ if $ra_psm->[5] == 1;
		}
	}
	my @decoy_pep = keys %decoy_pep;
	return \@decoy_pep;
}

# Object Method
# Title	    :  get_nr_target_pep_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_target_pep_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my $ra_ra_target_pep = $self->get_target_pep_by_fdr($fdr);
	my $nr_target_pep    = @$ra_ra_target_pep;

	return $nr_target_pep;
}

# Object Method
# Title	    :  get_nr_decoy_pep_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_decoy_pep_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my $ra_ra_decoy_pep = $self->get_decoy_pep_by_fdr($fdr);
	my $nr_decoy_pep    = @$ra_ra_decoy_pep;

	return $nr_decoy_pep;
}

# Object Method
# Title	    :  get_target_prot_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_target_prot_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my %target_prot;
	my $ra_ra_psm = $self->get_psm_by_fdr($fdr);

	foreach my $ra_psm (@$ra_ra_psm) {
		if ( defined( $ra_psm->[5] ) ) {
			$target_prot{ $ra_psm->[2] }++ if $ra_psm->[5] == 0;
		}
	}
	my @target_prot = keys %target_prot;
	return \@target_prot;
}

# Object Method
# Title	    :  get_decoy_prot_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_decoy_prot_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my %decoy_prot;
	my $ra_ra_psm = $self->get_psm_by_fdr($fdr);

	foreach my $ra_psm (@$ra_ra_psm) {
		if ( defined( $ra_psm->[5] ) ) {
			$decoy_prot{ $ra_psm->[2] }++ if $ra_psm->[5] == 1;
		}
	}
	my @decoy_prot = keys %decoy_prot;
	return \@decoy_prot;
}

# Object Method
# Title	    :  get_nr_target_prot_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_target_prot_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my $ra_ra_target_prot = $self->get_target_prot_by_fdr($fdr);
	my $nr_target_prot    = @$ra_ra_target_prot;

	return $nr_target_prot;
}

# Object Method
# Title	    :  get_nr_decoy_prot_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_decoy_prot_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my $ra_ra_decoy_prot = $self->get_decoy_prot_by_fdr($fdr);
	my $nr_decoy_prot    = @$ra_ra_decoy_prot;

	return $nr_decoy_prot;
}

# Object Method
# Title	    :  print_six_nr_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub print_six_nr_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my ( $nr_tpsm, $nr_dpsm, $nr_tpep, $nr_dpep, $nr_tpr, $nr_dpr ) =
	  $self->get_six_nr_by_fdr($fdr);
	print "  target: $nr_tpsm PSM,\t$nr_tpep peptides,\t$nr_tpr"
	  . " proteins at $fdr FDR\n";
	print "  decoy:  $nr_dpsm PSM,\t$nr_dpep peptides,\t$nr_dpr"
	  . " proteins at $fdr FDR\n";
}

# Object Method
# Title	    :  get_six_nr_by_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_six_nr_by_fdr {
	my $self = shift();
	my ($fdr) = @_;

	my $nr_tpsm = $self->get_nr_target_psm_by_fdr($fdr);
	my $nr_dpsm = $self->get_nr_decoy_psm_by_fdr($fdr);
	my $nr_tpep = $self->get_nr_target_pep_by_fdr($fdr);
	my $nr_dpep = $self->get_nr_decoy_pep_by_fdr($fdr);
	my $nr_tpr  = $self->get_nr_target_prot_by_fdr($fdr);
	my $nr_dpr  = $self->get_nr_decoy_prot_by_fdr($fdr);

	return ( $nr_tpsm, $nr_dpsm, $nr_tpep, $nr_dpep, $nr_tpr, $nr_dpr );
}

# Object Method
# Title	    :  filter_by_peptide_length()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub filter_by_peptide_length {
	my $self = shift();
	my ($min_pep_length) = @_;

	my $prior_nr_psm = @{ $self->{data}{psm} };

	return if $prior_nr_psm == 0;

	# 0 spectrum
	# 1 peptide
	# 2 main id
	# 3 modification info
	# 4 PeptideProphet probability (PPs)
	my @length_filtered_psm = ();
	foreach my $ra_psm ( @{ $self->{data}{psm} } ) {
		if ( length( $ra_psm->[1] ) >= $min_pep_length ) {

			push @length_filtered_psm, $ra_psm;
		}
	}

	$self->{data}{psm} = \@length_filtered_psm;

	my $post_nr_psm = @{ $self->{data}{psm} };

	my $diff            = $prior_nr_psm - $post_nr_psm;
	my $diff_percentage = $diff / $prior_nr_psm * 100;
	$diff_percentage = sprintf( "%.5f", $diff_percentage );
	$self->v_print( "$diff_percentage% of PSM removed "
		  . "(length < $min_pep_length aa)\n" );
}

# Object Method
# Title	    :  filter_ambiguous_peptides()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub filter_ambiguous_peptides {
	my $self = shift();
	my ($db) = @_;

	$self->v_print( "removing ambiguous peptides...\n" );
	
	my $prior_nr_psm = @{ $self->{data}{psm} };

	return if $prior_nr_psm == 0;

	my $p = FastaParser->new();
	
	my $rh_id_seq = $p->parse( $db );
	my @seq = ();
	foreach (keys%$rh_id_seq) {
		push @seq, $rh_id_seq->{$_};
	}
	my $joined_seq = join( '@', @seq );
	
	# 0 spectrum
	# 1 peptide
	# 2 main id
	# 3 modification info
	# 4 PeptideProphet probability (PPs)
	my $rh_pep;
	foreach ( @{ $self->{data}{psm} } ) {
		$rh_pep->{$_->[1]}++;
	}
	my $rh_amb_pep;
	foreach my $pep ( keys%$rh_pep ) {
		my $c = 0;
		while ( $joined_seq =~ /$pep/g ) {
			$c++;
			if ( $c == 2 ) {
				$rh_amb_pep->{$pep}++;
				last;
			}
		}
	}
	my @non_amb_psm = ();
	foreach ( @{ $self->{data}{psm} } ) {
		unless (exists($rh_amb_pep->{$_->[1]})) {
			push @non_amb_psm, $_;
		}
	}

	$self->{data}{psm} = \@non_amb_psm;

	my $post_nr_psm = @{ $self->{data}{psm} };

	my $diff            = $prior_nr_psm - $post_nr_psm;
	my $diff_percentage = $diff / $prior_nr_psm * 100;
	$diff_percentage = sprintf( "%.5f", $diff_percentage );
	$self->v_print( "$diff_percentage% of PSM with ambiguous"
		. " peptides removed\n" );
}

# Object Method
# Title	    :  add_target_decoy_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub add_target_decoy_fdr {
	my $self = shift();
	my ( $decoy_id_prefix, $target_decoy_ratio ) = @_;

	# ROC table
	# 0  ds cutoff
	# 1  PSM FDR
	# 2  FP PSM
	# 3  TP PSM
	# 4  td PSM FDR
	# 5  td FP PSM
	# 6  td TP PSM
	# 7  target PSM
	# 8  decoy PSM
	my $ra_ra_ROC = [];

	# 0 spectrum
	# 1 peptide
	# 2 main id
	# 3 modification info
	# 4 PeptideProphet probability (PPs)
	my $rh_ra_PPs_indices;
	for ( my $i = 0 ; $i < @{ $self->{data}{psm} } ; $i++ ) {
		
		# make the index of this PSM retreavable through its
		# discriminant score
		push @{ $rh_ra_PPs_indices->{ $self->{data}{psm}[$i][4] } }, $i;
	}

	# 0 spectrum
	# 1 peptide
	# 2 main id
	# 3 modification info
	# 4 PeptideProphet probability (PPs)
	# 5 1: decoy, 0: target
	# 6 PSM FDR estimate
	my ( $target_count, $decoy_count ) = ( 0, 0 );
	foreach my $pps ( sort { $b <=> $a } keys %$rh_ra_PPs_indices ) {

		# update the counts for all the PSM that have
		# the same PPs and add the decoy flag
		foreach my $i ( @{ $rh_ra_PPs_indices->{$pps} } ) {
			my $id = $self->{data}{psm}[$i][2];
			if ( $id =~ /^$decoy_id_prefix/ ) {
				push @{ $self->{data}{psm}[$i] }, 1;
				$decoy_count++;
			}
			else {
				push @{ $self->{data}{psm}[$i] }, 0;
				$target_count++;
			}
		}

		# PSM FDR calculation for only target PSM
		my $est_psm_fdr = 1;
		my $est_fp      = $decoy_count * $target_decoy_ratio;
		unless ( $target_count == 0 ) {
			$est_psm_fdr = $est_fp / $target_count;
		}
		my $est_tp = $target_count - $est_fp;
		$est_tp = 0 if $est_tp < 0;

		# PSM FDR calculation for target plus decoy PSM
		my $est_total_psm_fdr = 1;
		my $total_psm         = $target_count + $decoy_count;
		my $est_total_fp      = $decoy_count * ( 1 + $target_decoy_ratio );
		unless ( $total_psm == 0 ) {
			$est_total_psm_fdr = $est_total_fp / $total_psm;
		}
		my $est_total_tp = $total_psm - $est_total_fp;
		$est_total_tp = 0 if $est_total_tp < 0;

		# add the PSM FDR
		foreach my $i ( @{ $rh_ra_PPs_indices->{$pps} } ) {
			push @{ $self->{data}{psm}[$i] }, $est_psm_fdr;
		}

		# add row to table for the cutoff PPs
		# 0  ds cutoff
		# 1  PSM FDR
		# 2  FP PSM
		# 3  TP PSM
		# 4  td PSM FDR
		# 5  td FP PSM
		# 6  td TP PSM
		# 7  target PSM
		# 8  decoy PSM
		push @$ra_ra_ROC,
		  [
			$pps,          $est_psm_fdr,       $est_fp,
			$est_tp,       $est_total_psm_fdr, $est_total_fp,
			$est_total_tp, $target_count,      $decoy_count
		  ];
	}

	$self->{data}{roc}           = $ra_ra_ROC;
	$self->{data}{psm_fdr_added} = 1;
}

# Object Method
# Title	    :  filter_by_psm_fdr()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub filter_by_psm_fdr {
	my $self = shift();
	my ($max_psm_fdr) = @_;

	# if the PSM FDR was already calculated and added
	# to the PSM
	if ( exists( $self->{data}{psm_fdr_added} ) ) {
		if ( $self->{data}{psm_fdr_added} ) {

			my $prior_nr_psm = @{ $self->{data}{psm} };

			return if $prior_nr_psm == 0;

			# 0 spectrum
			# 1 peptide
			# 2 main id
			# 3 modification info
			# 4 PeptideProphet probability (PPs)
			# 5 1: decoy, 0: target
			# 6 PSM FDR estimate
			my @psm_fdr_filtered = ();
			foreach my $ra_psm ( @{ $self->{data}{psm} } ) {
				if ( $ra_psm->[6] <= $max_psm_fdr ) {

					push @psm_fdr_filtered, $ra_psm;
				}
			}

			$self->{data}{psm} = \@psm_fdr_filtered;

			my $post_nr_psm = @{ $self->{data}{psm} };

			my $diff            = $prior_nr_psm - $post_nr_psm;
			my $diff_percentage = $diff / $prior_nr_psm * 100;
			$diff_percentage = sprintf( "%.5f", $diff_percentage );
			$self->v_print( "$diff_percentage% of PSM removed "
				  . "(mFDR <= $max_psm_fdr)\n" );

		}
		else {
			print "  psm fdr needed!\n";
		}
	}
	else {
		print "  psm fdr needed!\n";
	}

}

# Object Method
# Title	    :  filter_by_protein_id()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub filter_by_protein_id {
	my $self = shift();
	my ($rh_ids) = @_;

	my $prior_nr_psm = @{ $self->{data}{psm} };

	return if $prior_nr_psm == 0;

	# 0 spectrum
	# 1 peptide
	# 2 main id
	# 3 modification info
	# 4 PeptideProphet probability (PPs)
	my @id_filtered_psm = ();
	foreach my $ra_psm ( @{ $self->{data}{psm} } ) {
		if ( exists( $rh_ids->{ $ra_psm->[2] } ) ) {
			push @id_filtered_psm, $ra_psm;
		}
	}

	$self->{data}{psm} = \@id_filtered_psm;

	my $post_nr_psm = @{ $self->{data}{psm} };

	my $diff            = $prior_nr_psm - $post_nr_psm;
	my $diff_percentage = $diff / $prior_nr_psm * 100;
	$diff_percentage = sprintf( "%.5f", $diff_percentage );
	$self->v_print( "$diff_percentage% of PSM removed "
		  . "(protein id not existing)\n" );
}

# Object Method
# Title	    :  get_nr_files()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_files {
	my $self = shift();

	my $nr_files = 0;
	if ( exists( $self->{nr_input_files} ) ) {
		$nr_files = $self->{nr_input_files};
	}

	return $nr_files;

}

# Object Method
# Title	    :  get_files()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_files {
	my $self = shift();

	my $ra_files = [];

	if ( exists( $self->{file} ) ) {
		if ( defined( $self->{file} ) ) {
			$ra_files = $self->{file};
		}
	}

	return $ra_files;
}

# Object Method
# Title	    :  get_nr_runs()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_nr_runs {
	my $self = shift();

	my $nr_runs = 0;
	if ( exists( $self->{nr_runs} ) ) {
		return $self->{nr_runs};
	}
	elsif ( exists( $self->{data}{psm} ) ) {
		my $scan = Scan->new();
		my %runs;
		foreach $_ ( @{ $self->{data}{psm} } ) {
			$runs{ $scan->get_run( $_->[0] ) }++;
		}
		$nr_runs = keys %runs;
	}

	return $nr_runs;

}

# Object Method
# Title	    :  get_roc()
# Usage     :
# Function	:  make a copy of the roc and return it
# Returns   :
# Args      :
sub get_roc {
	my $self = shift();
	my ($header) = @_;

	my @roc_copy = ();
	if ($header) {
		my @header = @{ $self->{roc_header} };
		unshift @roc_copy, \@header;
	}
	if ( exists( $self->{data}{roc} ) ) {
		foreach my $ra ( @{ $self->{data}{roc} } ) {
			my @row_copy = @$ra;
			push @roc_copy, \@row_copy;
		}
	}
	
	# 0  ds cutoff
	# 1  PSM FDR
	# 2  FP PSM
	# 3  TP PSM
	# 4  td PSM FDR
	# 5  td FP PSM
	# 6  td TP PSM
	# 7  target PSM
	# 8  decoy PSM
	return \@roc_copy;
}

# Object Method
# Title	    :  get_roc_processed()
# Usage     :
# Function	:  make a copy of the roc and return it
#              round numbers
# Returns   :  
# Args      :
sub get_roc_processed {
	my $self = shift();
	my ($header) = @_;

	my $ra_ra_roc = [];
	my $digits    = $self->{digits};
	my $format    = "%." . $digits . "f";

	my $ra_ra = $self->get_roc(0);

	foreach (@$ra_ra) {
		push @$ra_ra_roc,
		  [
			sprintf( $format, $_->[0] ),
			sprintf( $format, $_->[1] ),
			int( $_->[2] ),
			int( $_->[3] ),
			sprintf( $format, $_->[4] ),
			int( $_->[5] ),
			int( $_->[6] ),
			int( $_->[7] ),
			int( $_->[8] )
		  ];
	}

	if ($header) {
		my @header = @{ $self->{roc_header} };
		unshift @$ra_ra_roc, \@header;
	}

	return $ra_ra_roc;
}

# Object Method
# Title	    :  get_run_file()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_run_file {
	my $self = shift();

	if ( exists( $self->{run_file} ) ) {
		if ( defined( $self->{run_file} ) ) {
			return $self->{run_file};
		}
	}

	return {};
}

# Object Method
# Title	    :  get_runs()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_runs {
	my $self = shift();

	if ( exists( $self->{runs} ) ) {
		if ( defined( $self->{runs} ) ) {
			return $self->{runs};
		}
	}

	return {};
}

# Object Method
# Title	    :  get_file_run()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_file_run {
	my $self = shift();

	my $rh_rh_file_run = {};

	if ( exists( $self->{run_file} ) ) {
		if ( defined( $self->{run_file} ) ) {
			foreach my $run ( keys %{ $self->{run_file} } ) {
				my $file = $self->{run_file}{$run};
				$rh_rh_file_run->{$file}{$run}++;
			}
		}
	}

	return $rh_rh_file_run;
}

# Class Method
# Title	    :  is_mod()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub is_mod {
	my $self = shift;
	my ( $pepseq, $rh_pos_mass, $rh_mod_desc ) = @_;

	my $delta = $self->{mod_test_delta};

	if ( defined($rh_pos_mass) ) {
		foreach my $pos ( keys %$rh_pos_mass ) {

			# mass and amino acid of the modified peptide
			my $mass = $rh_pos_mass->{$pos};
			my $aa   = substr( $pepseq, $pos - 1, 1 );

			# wanted modifications
			foreach my $t_aa ( keys %$rh_mod_desc ) {
				if ( defined($t_aa) && defined($aa) ) {
					if ( $t_aa eq $aa ) {
						my $t_mass = $rh_mod_desc->{$t_aa};
						if ( abs( $t_mass - $mass ) < $delta ) {
							return 1;
						}
					}
				}
			}
		}
	}

	return 0;
}

# Object Method
# Title	    :  is_empty()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub is_empty {
	my $self = shift();

	if ( exists( $self->{data}{psm} ) ) {
		if ( defined( $self->{data}{psm} ) ) {
			my $nr = @{ $self->{data}{psm} };
			if ( $nr == 0 ) {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
	else {
		return 1;
	}
}

# TODO remove
sub print_object {    # print the object to the console with dumpvar
	my $self = shift();
	my ($r) = @_;
	if ( defined($r) ) {
		{

			package main;
			require "dumpvar.pl";
			dumpValue($r);
		}
	}
	else {

		package main;
		require "dumpvar.pl";
		dumpValue($self);
	}
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





