package TandemMSIdSelectionScheme;
use strict;

use MayuTools;
use Scan;

##################################################################
#
# TandemMSIdSelectionScheme
#
# NOTES:
# selects TandemMS derived Ids in different ways and passes them
# in each step to MayuManager. Because the PSMSets are quite
# big MayuManager is passed to TandemMSIdSelectionScheme
# and the selected identification sets are processed one after
# another instead of returning an array of identification sets.
#
# selection type for PSM sets:
#  all_id_mFDR_range()
#      all data                                -> a range of mFDR
#
#  x: only special mod peptides (e.g. phospho) -> a range of mFDR
# -1: sort runs according to orthogonality     -> a range of mFDR
# -2: shuffle runs                             -> a range of mFDR
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
# 2008.02.20
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $s = TandemMSIdSelectionScheme->new( $v, $s );
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

	return $self;
}

# Object Method
# Title	    :  all_id_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub all_id_mFDR_range {
	my $self = shift;
	my ( $manager, $psm_set, $max_psm_fdr, $psm_fdr_steps ) = @_;

	$| = 1 if $self->{s};

	my $nr_files = $psm_set->get_nr_files();
	my $nr_runs  = $psm_set->get_nr_runs();

	print "  $nr_runs total runs, looping through id sets...\n"
	  if $self->{v};
	my $psm_fdr_step_size = $max_psm_fdr / ( $psm_fdr_steps - 1 );
	for ( my $i = 0 ; $i < $psm_fdr_steps ; $i++ ) {
		my $step    = $i + 1;
		my $psm_fdr = $i * $psm_fdr_step_size;
		print "\r                                             \r"
		  if $self->{s};
		print "  step $step\t$psm_fdr mFDR" if $self->{s};

		my $ra_ra_psm = $psm_set->get_psm_by_fdr($psm_fdr);

		# feed this set of selected identifications into
		# the TandemMSIdSet
		# TandemMSIdSet will distribute it among the
		# registered error models
		my $rh_attributes = {
			'nr_files' => $nr_files,
			'nr_runs'  => $nr_runs,
			'mFDR'     => $psm_fdr
		};
		$manager->set_psm_set( $ra_ra_psm, $rh_attributes );
	}
	print "\n" if $self->{s};
}

# Object Method
# Title	    :  cumulative_file_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub cumulative_file_mFDR_range {
	my $self = shift;
	my ( $manager, $psm_set, $max_psm_fdr, $psm_fdr_steps ) = @_;

	my $ra_files       = $psm_set->get_files();
	my $rh_rh_file_run = $psm_set->get_file_run();
	my $nr_runs        = $psm_set->get_nr_runs();

	my @sorted_files = sort @$ra_files;

	print "  $nr_runs total runs, looping through id sets...\n"
	  if $self->{v};
	$self->loop_file_mFDR_range( \@sorted_files, $rh_rh_file_run,
		$psm_fdr_steps, $max_psm_fdr, $psm_set, $manager );
}

# Object Method
# Title	    :  cumulative_shuffled_file_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub cumulative_shuffled_file_mFDR_range {
	my $self = shift;
	my ( $manager, $psm_set, $max_psm_fdr, $psm_fdr_steps ) = @_;

	my $ra_files       = $psm_set->get_files();
	my $rh_rh_file_run = $psm_set->get_file_run();
	my $nr_runs        = $psm_set->get_nr_runs();

	# shuffle the array
	my $ra_shuffled_files = $self->{tools}->shuffle_array($ra_files);

	print "  $nr_runs total runs, looping through id sets...\n"
	  if $self->{v};
	$self->loop_file_mFDR_range( $ra_shuffled_files, $rh_rh_file_run,
		$psm_fdr_steps, $max_psm_fdr, $psm_set, $manager );
}

# Object Method
# Title	    :  loop_file_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub loop_file_mFDR_range {
	my $self = shift;
	my ( $ra_files, $rh_rh_file_run, $psm_fdr_steps, $max_psm_fdr,
		$psm_set, $manager )
	  = @_;

	$| = 1 if $self->{s};

	my $nr_files = @$ra_files;

	my $psm_fdr_step_size = $max_psm_fdr / ( $psm_fdr_steps - 1 );

	my @cumulative_files = ();
	for ( my $j = 0 ; $j < $nr_files ; $j++ ) {

		my $file = $ra_files->[$j];
		push @cumulative_files, $file;
		my $nr_f = $j + 1;

		my $nr_runs = 0;
		foreach my $f (@cumulative_files) {
			if ( exists( $rh_rh_file_run->{$f} ) ) {
				my $r = keys %{ $rh_rh_file_run->{$f} };
				$nr_runs += $r;
			}
		}

		for ( my $i = 0 ; $i < $psm_fdr_steps ; $i++ ) {
			my $step    = $i + 1;
			my $psm_fdr = $i * $psm_fdr_step_size;
			print "\r                                             \r"
			  . "\r                                             \r"
			  if $self->{s};
			print
			  "  $nr_f/$nr_files input files, step $step, $psm_fdr mFDR"
			  if $self->{s};

			my $ra_ra_psm =
			  $psm_set->get_psm_by_files_and_fdr( \@cumulative_files,
				$psm_fdr );

			# feed this set of selected identifications into
			# the TandemMSIdSet
			# TandemMSIdSet will distribute it among the
			# registered error models
			my $rh_attributes = {
				'nr_files' => $nr_f,
				'nr_runs'  => $nr_runs,
				'mFDR'     => $psm_fdr
			};
			$manager->set_psm_set( $ra_ra_psm, $rh_attributes );
		}
	}
	print "\n" if $self->{s};
}

# Object Method
# Title	    :  cumulative_run_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub cumulative_run_mFDR_range {
	my $self = shift;
	my ( $manager, $psm_set, $max_psm_fdr, $psm_fdr_steps, $nr_runs_steps )
	  = @_;

	my $rh_rh_run_file = $psm_set->get_run_file();
	my @runs           = sort keys %$rh_rh_run_file;
	my $nr_runs        = @runs;

	print "  $nr_runs total runs, looping through id sets...\n"
	  if $self->{v};
	$self->loop_run_mFDR_range( \@runs, $nr_runs_steps, $psm_fdr_steps,
		$max_psm_fdr, $psm_set, $manager );
}

# Object Method
# Title	    :  cumulative_shuffled_run_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub cumulative_shuffled_run_mFDR_range {
	my $self = shift;
	my ( $manager, $psm_set, $max_psm_fdr, $psm_fdr_steps, $nr_runs_steps )
	  = @_;

	my $nr_files       = $psm_set->get_nr_files();
	my $ra_files       = $psm_set->get_files();
	my $rh_rh_run_file = $psm_set->get_run_file();
	my @runs           = sort keys %$rh_rh_run_file;
	my $nr_runs        = @runs;

	my $ra_shuffled_run = $self->{tools}->shuffle_array( \@runs );

	print "  $nr_runs total runs, looping through id sets...\n"
	  if $self->{v};
	$self->loop_run_mFDR_range( $ra_shuffled_run, $nr_runs_steps,
		$psm_fdr_steps, $max_psm_fdr, $psm_set, $manager );
}

# Object Method
# Title	    :  ortho_run_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub ortho_run_mFDR_range {
	my $self = shift;
	my ( $manager, $psm_set, $max_psm_fdr, $psm_fdr_steps, $nr_runs_steps,
		$out_base, $version )
	  = @_;

	my $ortho_ana_fdr = 0.01;

	my $nr_files       = $psm_set->get_nr_files();
	my $ra_files       = $psm_set->get_files();
	my $rh_rh_file_run = $psm_set->get_file_run();
	my $rh_rh_run_file = $psm_set->get_run_file();
	my @all_runs       = sort keys %$rh_rh_run_file;
	my $nr_runs        = @all_runs;

	my $ra_ra_psm = $psm_set->get_psm_by_fdr($ortho_ana_fdr);

	my $rh_run_est_ms = $self->get_est_ms($ra_ra_psm);

	print "  estimating the orthogonality of each LC-MS/MS run...\n"
	  if $self->{v};
	my $ra_ortho_runs =
	  $self->get_ortho_runs( $ra_ra_psm, $rh_run_est_ms, $out_base,
		$version, \@all_runs );

	print "  $nr_runs total runs, looping through id sets...\n"
	  if $self->{v};
	$self->loop_run_mFDR_range(
		$ra_ortho_runs, $nr_runs_steps, $psm_fdr_steps,
		$max_psm_fdr,   $psm_set,       $manager
	);
}

# Object Method
# Title	    :  get_ortho_runs()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_ortho_runs {
	my $self = shift;
	my ( $ra_ra_psm, $rh_run_est_ms, $out_base, $version, $ra_all_runs ) =
	  @_;

	$| = 1 if $self->{s};

	my $rh_all_pep     = {};
	my $rh_all_prot    = {};
	my $rh_rh_pep_run  = {};
	my $rh_rh_prot_run = {};
	my $rh_rh_run_pep  = {};
	my $rh_rh_run_prot = {};
	my $scan           = Scan->new();
	foreach my $ra (@$ra_ra_psm) {
		$rh_all_pep->{ $ra->[1] }++;
		$rh_all_prot->{ $ra->[2] }++;
		$rh_rh_pep_run->{ $ra->[1] }{ $scan->get_run( $ra->[0] ) }++;
		$rh_rh_prot_run->{ $ra->[2] }{ $scan->get_run( $ra->[0] ) }++;
		$rh_rh_run_pep->{ $scan->get_run( $ra->[0] ) }{ $ra->[1] }++;
		$rh_rh_run_prot->{ $scan->get_run( $ra->[0] ) }{ $ra->[2] }++;
	}

	print "  peptide orthogonality...\n"
	  if $self->{v};
	my $rh_run_pep_diff = {};
	foreach my $pep ( keys %$rh_rh_pep_run ) {
		print "\r  $pep                                             "
		  if $self->{s};
		my @runs    = keys %{ $rh_rh_pep_run->{$pep} };
		my $nr_runs = @runs;
		if ( $nr_runs == 1 ) {
			$rh_run_pep_diff->{ $runs[0] }++;
		}
	}
	print "\n" if $self->{s};

	# add the runs that have no identifications
	foreach my $run (@$ra_all_runs) {
		unless ( exists( $rh_run_pep_diff->{$run} ) ) {

			# empty hashes
			$rh_run_pep_diff->{$run}  = 0;
		}
	}

	print "  protein orthogonality...\n"
	  if $self->{v};
	my $rh_run_prot_diff = {};
	foreach my $prot ( keys %$rh_rh_prot_run ) {
		print "\r  $prot                                             "
		  if $self->{s};
		my @runs    = keys %{ $rh_rh_prot_run->{$prot} };
		my $nr_runs = @runs;
		if ( $nr_runs == 1 ) {
			$rh_run_prot_diff->{ $runs[0] }++;
		}
	}
	print "\n" if $self->{s};

	# add the runs that have no identifications
	foreach my $run (@$ra_all_runs) {
		unless ( exists( $rh_run_prot_diff->{$run} ) ) {

			# empty hashes
			$rh_run_prot_diff->{$run}  = 0;
		}
	}

	# 0: run
	# 1: pep
	# 2: prot
	# 3: est ms2
	# 4: pep diff
	# 5: prot diff
	# 6: pep diff norm
	# 7: prot diff norm
	my $ra_ra_run_ortho_table = [];
	foreach my $run ( sort keys %{$rh_run_pep_diff} ) {
		my ( $nr_pep, $nr_prot, $est_ms2, $pep_diff, $prot_diff,
			$pep_diff_norm, $prot_diff_norm )
		  = ( 0, 0, 0, 0, 0, 0 );
		$nr_pep   = keys %{ $rh_rh_run_pep->{$run} };
		$nr_prot  = keys %{ $rh_rh_run_prot->{$run} };
		$est_ms2  = $rh_run_est_ms->{$run};
		$pep_diff = $rh_run_pep_diff->{$run}
		  if exists( $rh_run_pep_diff->{$run} );
		$prot_diff = $rh_run_prot_diff->{$run}
		  if exists( $rh_run_prot_diff->{$run} );
		$pep_diff_norm  = $pep_diff / $est_ms2  unless $est_ms2 == 0;
		$prot_diff_norm = $prot_diff / $est_ms2 unless $est_ms2 == 0;

		push @$ra_ra_run_ortho_table,
		  [
			$run,           $nr_pep,   $nr_prot,
			$est_ms2,       $pep_diff, $prot_diff,
			$pep_diff_norm, $prot_diff_norm
		  ];
	}
	print "\n" if $self->{s};

	# sort according to pep diff norm
	my @s = sort { $b->[6] <=> $a->[6] } @$ra_ra_run_ortho_table;

	my $ra_ortho_runs = [];
	foreach (@s) {
		push @$ra_ortho_runs, $_->[0];
	}

	my $ortho_file = $out_base . '_ortho_run_' . $version . '.txt';
	unshift @s,
	  [
		'run',           'nr_pep',
		'nr_prot',       'est_scan',
		'pep_diff',      'prot_diff',
		'pep_diff_norm', 'prot_diff_norm'
	  ];
	$self->{tools}->print_data_in_table_style( \@s, $ortho_file );

	return $ra_ortho_runs;
}

# Object Method
# Title	    :  loop_run_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_est_ms {
	my $self = shift;
	my ($ra_ra_psm) = @_;

	my $rh_run_est_max_ms = {};

	my $scan = Scan->new();
	foreach my $ra_psm (@$ra_ra_psm) {
		my $run     = $scan->get_run( $ra_psm->[0] );
		my $scan_nr = $scan->get_scan_nr1( $ra_psm->[0] );
		if ( exists( $rh_run_est_max_ms->{$run} ) ) {
			if ( $scan_nr > $rh_run_est_max_ms->{$run} ) {
				$rh_run_est_max_ms->{$run} = $scan_nr;
			}
		}
		else {
			$rh_run_est_max_ms->{$run} = $scan_nr;
		}
	}

	return $rh_run_est_max_ms;
}

# Object Method
# Title	    :  loop_run_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub loop_run_mFDR_range {
	my $self = shift;
	my ( $ra_runs, $nr_runs_steps, $psm_fdr_steps, $max_psm_fdr, $psm_set,
		$manager )
	  = @_;

	$| = 1 if $self->{s};

	my $psm_fdr_step_size = $max_psm_fdr / ( $psm_fdr_steps - 1 );

	my $nr_runs          = @$ra_runs;
	my $nr_runs_per_step = int( $nr_runs / $nr_runs_steps ) + 1;

	my @cumulative_runs = ();
	for ( my $j = 0 ; $j < $nr_runs ; $j++ ) {

		my $run = $ra_runs->[$j];
		push @cumulative_runs, $run;
		my $nr_current_runs = @cumulative_runs;

		if ( $nr_current_runs % $nr_runs_per_step == 0
			|| ( $j + 1 ) == $nr_runs )
		{

			my $nr_f = 0;

			for ( my $i = 0 ; $i < $psm_fdr_steps ; $i++ ) {
				my $step    = $i + 1;
				my $psm_fdr = $i * $psm_fdr_step_size;
				print "\r                                             \r"
				  . "\r                                             \r"
				  if $self->{s};
				print "  $nr_current_runs/$nr_runs runs, step $step, "
				  . "$psm_fdr mFDR"
				  if $self->{s};

				my $ra_ra_psm =
				  $psm_set->get_psm_by_run_and_fdr( \@cumulative_runs,
					$psm_fdr );

				# feed this set of selected identifications into
				# the TandemMSIdSet
				# TandemMSIdSet will distribute it among the
				# registered error models
				my $rh_attributes = {
					'nr_files' => $nr_f,
					'nr_runs'  => $nr_current_runs,
					'mFDR'     => $psm_fdr
				};
				$manager->set_psm_set( $ra_ra_psm, $rh_attributes );
			}
		}
	}
	print "\n" if $self->{s};
}

# Object Method
# Title	    :  phospho_id_mFDR_range()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub phospho_id_mFDR_range {
	my $self = shift;
	my ( $manager, $psm_set, $max_psm_fdr, $psm_fdr_steps, $nr_runs_steps )
	  = @_;

	$| = 1 if $self->{s};

	my $rh_phospho_mod = {
		'S' => 166.9984,
		'T' => 181.0140,
		'Y' => 243.0297
	};

	my $nr_files = $psm_set->get_nr_files();
	my $nr_runs  = $psm_set->get_nr_runs();

	my $psm_fdr_step_size = $max_psm_fdr / ( $psm_fdr_steps - 1 );
	for ( my $i = 0 ; $i < $psm_fdr_steps ; $i++ ) {
		my $step    = $i + 1;
		my $psm_fdr = $i * $psm_fdr_step_size;
		print "\r                                             \r"
		  . "\r                                             \r"
		  if $self->{s};
		print "  step $step\t$psm_fdr mFDR" if $self->{s};

		my $ra_ra_psm =
		  $psm_set->get_mod_psm_by_fdr( $psm_fdr, $rh_phospho_mod );

		# feed this set of selected identifications into
		# the error models ( registered in the writer )
		my $rh_attributes = {
			'nr_files' => $nr_files,
			'nr_runs'  => $nr_runs,
			'mFDR'     => $psm_fdr
		};
		$manager->set_psm_set( $ra_ra_psm, $rh_attributes );
	}
	print "\n" if $self->{s};
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












