package MascotCSV;
use strict;

use File::Basename;
use CSVParser;

##################################################################
#
# NOTES:
# initially written for Mayu.
#
#
#		"prot_hit_num","prot_acc","prot_desc","prot_score",
#		"prot_mass","prot_matches","pep_query","pep_exp_mz",
#		"pep_exp_mr","pep_exp_z","pep_calc_mr","pep_delta",
#		"pep_miss","pep_score","pep_expect","pep_rank","pep_res_before",
#		"pep_seq","pep_res_after","pep_var_mod"
#
# SYNOPSIS:
#
#
# TODO:
# -
# -
#
# CHANGELOG:
# 2007.07.14
#            - created
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $m = MascotCSV->new( $verbose, $status );
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $verbose, $status ) = @_;

	# save the options
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

	# start with the default values
	my $csv_parser = CSVParser->new();
	$self->{csv_parser} = $csv_parser;

	return $self;
}

# Object Method
# Title	    :  is_mascot_csv_file()
# Usage     :  $m->is_mascot_csv_file( $csv );
# Function	:  uses a csv parser to parse the data
# Returns   :
# Args      :
sub is_mascot_csv_file {
	my $self = shift;
	my ($csv) = @_;

	my $is_mascot_csv_file = 0;

	if ( -e $csv ) {

		open( I, $csv ) or warn $!;
		while (<I>) {
			my $line = $_;
			chomp($line);

			if (   $line =~ /pep_query/
				&& $line =~ /pep_seq/
				&& $line =~ /prot_hit_num/
				&& $line =~ /prot_acc/
				&& $line =~ /pep_var_mod/
				&& $line =~ /pep_score/ )
			{
				$is_mascot_csv_file = 1;
				close(I);
				last;
			}
		}
		close(I);
	}

	return $is_mascot_csv_file;
}

# Object Method
# Title	    :  get_mayu_psm()
# Usage     :  my $ra_p_psm = $m->get_mayu_psm( $csv );
# Function	:  uses a csv parser to parse the data
# Returns   :
# Args      :
# Var       :  depricated! use get_tpp_psm() instead!
sub get_mayu_psm {
	my $self = shift;
	my ($csv) = @_;

	my $ra_ra_mayu_psm = [];

	my $csv_base = basename($csv);
	$csv_base =~ s/\.[^\.]+$//g;

	my ( $ra_ra_header, $ra_col_names, $ra_ra_psm, $rh_colname_indices ) =
	  $self->get_psm($csv);

	return $ra_ra_mayu_psm unless defined($rh_colname_indices);

	if (   exists( $rh_colname_indices->{pep_query} )
		&& exists( $rh_colname_indices->{pep_exp_z} )
		&& exists( $rh_colname_indices->{pep_seq} )
		&& exists( $rh_colname_indices->{prot_acc} )
		&& exists( $rh_colname_indices->{pep_var_mod} )
		&& exists( $rh_colname_indices->{pep_score} )
		&& defined($ra_ra_psm) )
	{
		my $snri = $rh_colname_indices->{pep_query};
		my $chi  = $rh_colname_indices->{pep_exp_z};
		my $psi  = $rh_colname_indices->{pep_seq};
		my $pai  = $rh_colname_indices->{prot_acc};
		my $mi   = $rh_colname_indices->{pep_var_mod};
		my $si   = $rh_colname_indices->{pep_score};

		foreach my $ra_psm (@$ra_ra_psm) {
			my $spec =
			    $csv_base . '.'
			  . $ra_psm->[$snri] . '.'
			  . $ra_psm->[$snri] . '.'
			  . $ra_psm->[$chi];
			my $pep_seq = $ra_psm->[$psi];
			my $prot    = $ra_psm->[$pai];

			# mod
			my @mod = split( ';', $ra_psm->[$mi] );
			my $rh_mod = {};
			for ( my $j = 0 ; $j < @mod ; $j++ ) {
				$rh_mod->{$j} = $mod[$j];
			}
			my $ds          = $ra_psm->[$si];
			my $ra_mayu_psm = [ $spec, $pep_seq, $prot, $rh_mod, $ds ];

			push @$ra_ra_mayu_psm, $ra_mayu_psm;
		}

	}

	return $ra_ra_mayu_psm;
}

# Object Method
# Title	    :  get_tpp_psm()
# Usage     :  my $ra_p_psm = $m->get_tpp_psm( $csv );
# Function	:  uses a csv parser to parse the data
# Returns   :
# Args      :
sub get_tpp_psm {
	my $self = shift;
	my ($csv) = @_;

	my $ra_ra_tpp_psm = [];

	my $csv_base = basename($csv);
	$csv_base =~ s/\.[^\.]+$//g;

	my $ra_col_names = [
		'pep_query',   'pep_exp_z', 'pep_seq', 'prot_acc',
		'pep_var_mod', 'pep_score'
	];

	my $ra_ra_sel_psm = $self->get_psm_by_col_names( $csv, $ra_col_names );

	if ( defined($ra_ra_sel_psm) ) {
		if ( @$ra_ra_sel_psm > 0 ) {
			if ( @{ $ra_ra_sel_psm->[0] } == 6 ) {
				foreach my $ra_psm (@$ra_ra_sel_psm) {
					my $spec =
					    $csv_base . '.'
					  . $ra_psm->[0] . '.'
					  . $ra_psm->[0] . '.'
					  . $ra_psm->[1];
					my $pep_seq = $ra_psm->[2];
					my $prot    = $ra_psm->[3];

					# mod
					my @mod = split( ';', $ra_psm->[4] );
					my $rh_mod = {};
					for ( my $j = 0 ; $j < @mod ; $j++ ) {
						$rh_mod->{$j} = $mod[$j];
					}
					my $ds          = $ra_psm->[5];
					my $ra_mayu_psm =
					  [ $spec, $pep_seq, $prot, $rh_mod, $ds ];

					push @$ra_ra_tpp_psm, $ra_mayu_psm;
				}
			}
		}
	}

	return $ra_ra_tpp_psm;
}

# Object Method
# Title	    :  get_psm_by_col_names()
# Usage     :  my $ra_ra_psm = $m->get_psm_by_col_names( $csv, $ra_col_names );
# Function	:  uses a csv parser to parse the data
# Returns   :
# Args      :
sub get_psm_by_col_names {
	my $self = shift;
	my ( $csv, $ra_wanted_col_names ) = @_;

	# result with only the selected columns
	my $ra_ra_sel_psm = [];

	my $csv_base = basename($csv);
	$csv_base =~ s/\.[^\.]+$//g;

	my ( $ra_ra_header, $ra_col_names, $ra_ra_psm, $rh_colname_indices ) =
	  $self->get_psm($csv);

	# if the column names could not be extracted from the mascot csv file
	# or the column names were not specified
	if ( !defined($rh_colname_indices) || !defined($ra_col_names) ) {
		return $ra_ra_sel_psm;
	}

	# extract the indices
	my $ra_indices = [];
	foreach my $wanted_col_name (@$ra_wanted_col_names) {
		foreach my $mascot_col_name ( keys %$rh_colname_indices ) {
			if ( $mascot_col_name =~ /^$wanted_col_name/ ) {
				push @$ra_indices, $rh_colname_indices->{$mascot_col_name};
			}
		}
	}

	foreach my $ra_psm (@$ra_ra_psm) {
		my $ra_sel_psm = [];
		foreach (@$ra_indices) {
			if ( @$ra_psm > $_ ) {
				push @$ra_sel_psm, $ra_psm->[$_];
			}
			else {
				push @$ra_sel_psm, undef;
			}
		}
		push @$ra_ra_sel_psm, $ra_sel_psm;
	}

	return $ra_ra_sel_psm;
}

# Object Method
# Title	    :  get_psm()
# Usage     :  my ( $header, $col_names, $psm $rh_col_ind)
#                = $m->get_psm( $csv );
# Function	:  uses a csv parser to parse the data
# Returns   :
# Args      :
sub get_psm {
	my $self = shift;
	my ($csv) = @_;

	# header
	my $ra_ra_header = [];

	# column names
	my $ra_col_names        = [];
	my $column_names_passed = 0;

	# peptide spectrum matches
	my $ra_ra_psm = [];

	# parse the csv file
	my @lines = $self->{csv_parser}->parse_file($csv);

	foreach my $ra_line (@lines) {
		my $line = join( " ", @$ra_line );
		if (   $line =~ /pep_query/
			&& $line =~ /pep_seq/
			&& $line =~ /prot_hit_num/
			&& $line =~ /prot_acc/
			&& $line =~ /pep_var_mod/
			&& $line =~ /pep_score/ )
		{
			$column_names_passed = 1;
			$ra_col_names        = $ra_line;
		}
		elsif ($column_names_passed) {
			push @$ra_ra_psm, $ra_line;
		}
		else {
			push @$ra_ra_header, $ra_line;
		}
	}

	# complete the prot_acc, prot_desc, prot_score, prot_mass,
	# prot_matches
	my $rh_colname_indices;
	for ( my $i = 0 ; $i < @$ra_col_names ; $i++ ) {
		$rh_colname_indices->{ $ra_col_names->[$i] } = $i;
	}
	my $pn = $rh_colname_indices->{prot_hit_num};
	my $pa = $rh_colname_indices->{prot_acc};
	my $pd = $rh_colname_indices->{prot_desc};
	my $ps = $rh_colname_indices->{prot_score};
	my $pm = $rh_colname_indices->{prot_mass};
	my $pt = $rh_colname_indices->{prot_matches};

	# mascot does not store the protein information for each
	# PSM. therefore this information has to be completed
	#
	# save the protein information for each protein hit number
	my $rh_protnr_prot;
	foreach my $ra_psm (@$ra_ra_psm) {
		unless ( $ra_psm->[$pa] =~ /^$/ ) {
			$rh_protnr_prot->{ $ra_psm->[$pn] } = [
				$ra_psm->[$pa], $ra_psm->[$pd], $ra_psm->[$ps],
				$ra_psm->[$pm], $ra_psm->[$pt],
			];
		}
	}
	foreach my $ra_psm (@$ra_ra_psm) {
		if ( $ra_psm->[$pa] =~ /^$/ ) {
			$ra_psm->[$pa] = $rh_protnr_prot->{ $ra_psm->[$pn] }[0];
			$ra_psm->[$pd] = $rh_protnr_prot->{ $ra_psm->[$pn] }[1];
			$ra_psm->[$ps] = $rh_protnr_prot->{ $ra_psm->[$pn] }[2];
			$ra_psm->[$pm] = $rh_protnr_prot->{ $ra_psm->[$pn] }[3];
			$ra_psm->[$pt] = $rh_protnr_prot->{ $ra_psm->[$pn] }[4];
		}
	}

	return ( $ra_ra_header, $ra_col_names, $ra_ra_psm,
		$rh_colname_indices );
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




