package MayuProteinSize;
use strict;

use FastaParser;
use Digest;
use Bins;
use BinnedEntity;

##################################################################
#
# MayuProteinSize
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
# 2007.06.22
# -
# -
#
##################################################################

# Constructor
# Title     :  new
# Usage     :  my $s = MayuProteinSize->new( $verbose, $status );
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

	return $self;
}

# Object Method
# Title	    :  protein_size_analysis()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub protein_size_analysis {
	my $self = shift;
	my (
		$tools,
		$db,
		$min_pep_length,
		$decoy_id_prefix,
		$nr_missed_cleavages,
		$min_pep_mass,
		$max_pep_mass,
		$nr_prot_size_bins,
		$dont_corr_id_seq,
		$mirror_decoy_ids_from_target_ids,
		$gene_group_identification_type,
		$gene_group_regex,
		$general_min_pep_length,
		$special_zero_bin,
		$equinr_bins,
		$filter_ids
	  )
	  = @_;

	# $rh_ra_id_protprop
	# 0: aa length
	# 1: identical sequence corrected aa length
	# 2: number of tryptic peptides
	# 3: gene group corrected number of tryptic peptides
	my $rh_ra_id_protprop;

	my $fp = FastaParser->new( $self->{v}, $self->{s} );
	$self->v_print("parsing $db...\n");
	my $rh_id_seq = $fp->parse($db);

	# check whether there is a target id for each decoy id
	#
	# if yes, delete all the decoy entries in the hash $rh_id_seq
	# because the calculations do not have to be done twice
	# 
	# if the decoys are perfectly mirrored:
	# the returned $rh_id_seq corresponds to target ids and sequences only!
	#
	# if the decoys are not perfectly mirrored:
	# the returned $rh_id_seq corresponds to both, target and decoy ids and sequences
	if ($mirror_decoy_ids_from_target_ids) {
		$self->v_print("checking decoy ids for mirror target ids...");
		( $mirror_decoy_ids_from_target_ids, $rh_id_seq ) =
		  $self->check_decoy_mirror( $rh_id_seq, $decoy_id_prefix );
		if ($mirror_decoy_ids_from_target_ids) {
			$self->v_print("ok\n");
		}
	}

	# correct sequence length of groups of identical sequences
	# add to $rh_ra_id_protprop
	# 0: aa length
	# 1: identical sequence corrected aa length
	unless ($dont_corr_id_seq) {
		( $rh_id_seq, $rh_ra_id_protprop ) =
		  $self->correct_identical_sequences($rh_id_seq);
	}

	# group according to gene to correct for identical tryptic peptides
	my $rh_rh_gene_id_seq =
	  $self->group_according_to_gene( $rh_id_seq,
		$gene_group_identification_type,
		$gene_group_regex );

	# calculate the number of tryptic peptides for each protein
	# add to $rh_ra_id_protprop
	# 2: ntp
	# 3: corrected ntp
	$self->v_print("estimating protein sizes...\n");
	$rh_ra_id_protprop = $self->get_protein_sizes(
		$rh_ra_id_protprop,                $rh_rh_gene_id_seq,
		$min_pep_length,                   $decoy_id_prefix,
		$nr_missed_cleavages,              $min_pep_mass,
		$max_pep_mass,                     $nr_prot_size_bins,
		$mirror_decoy_ids_from_target_ids, $general_min_pep_length
	);

	# fixed 20150128 because of Lorenz Blum
	# moved up here in front of the protein binning
	#
	# copy the target entries to decoy entries
	if ($mirror_decoy_ids_from_target_ids) {
		$self->v_print("copying target entries to decoy entries...\n");
		foreach ( keys %$rh_id_seq ) {
			my $decoy_id = $decoy_id_prefix . $_;
			$rh_ra_id_protprop->{$decoy_id} = $rh_ra_id_protprop->{$_};
		}
	}

	# $ra_bin_nrprot
	# index: bin number
	# value: number of total target proteins in that bin
	#
	# add to $rh_ra_id_protprop
	# 4: bin nr
	$self->v_print("getting bins for the proteins...\n");
	my $ra_bin_nrprot = [];
	my $bins_ob;
	( $ra_bin_nrprot, $rh_ra_id_protprop, $bins_ob ) =
	  $self->bin_proteins( $rh_ra_id_protprop,
		$mirror_decoy_ids_from_target_ids,
		$nr_prot_size_bins, $decoy_id_prefix, $special_zero_bin,
		$equinr_bins );
	my $ra_breaks = $bins_ob->get_breaks();



	# create a BinnedEntity object for proteins
	my $binned_prot_entity =
	  $self->get_binned_prot_entity( $rh_ra_id_protprop, $ra_bin_nrprot,
		$bins_ob, $filter_ids );

	# create a BinnedEntity object for peptides with one bin
	my $binned_pep_entity =
	  $self->get_binned_pep_entity($rh_ra_id_protprop);

	return ( $binned_prot_entity, $binned_pep_entity );
}

# Object Method
# Title	    :  get_binned_entity()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_binned_prot_entity {
	my $self = shift;
	my ( $rh_ra_id_protprop, $ra_bin_nrprot, $bins_ob, $filter_ids ) = @_;
                                                                      
	# $filter_ids = 0 means that for proteins from the PSM file which     
	# are not found in the fasta, Mayu will use an average protein length 
	my $binned_entity      = BinnedEntity->new( $filter_ids );            
	my %index_feature_name = (
		'aa'             => 0,
		'id_seq_corr_aa' => 1,
		'ntp'            => 2,
		'corr_ntp'       => 3,
		'bin'            => 4
	);

	while ( my ( $name, $value ) = each %index_feature_name ) {
		my $ra_ra_id_feature = [];
		foreach ( keys %$rh_ra_id_protprop ) {
			push @$ra_ra_id_feature,
			  [ $_, $rh_ra_id_protprop->{$_}[$value] ];
		}
		$binned_entity->add_feature( $ra_ra_id_feature, $name );
	}

	# add the binsob
	for ( my $i = 0 ; $i < @$ra_bin_nrprot ; $i++ ) {
		$bins_ob->add_att_to_bin( $i, 'entities', $ra_bin_nrprot->[$i] );
	}

	$binned_entity->add_bins_ob($bins_ob);

	return $binned_entity;
}

# Object Method
# Title	    :  get_binned_pep_entity()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_binned_pep_entity {
	my $self = shift;
	my ($rh_ra_id_protprop) = @_;

	my $corr_ntp_index = 3;

	my $binned_entity = BinnedEntity->new();

	my $total_nr_pep = 0;
	foreach my $id ( keys %$rh_ra_id_protprop ) {
		$total_nr_pep += $rh_ra_id_protprop->{$id}[$corr_ntp_index];
	}

	my $bins_ob = Bins->new();
	$bins_ob->add_bin( 0, 1000000 );
	$bins_ob->add_att_to_bin( 0, 'entities', $total_nr_pep );
	$binned_entity->add_bins_ob($bins_ob);

	return $binned_entity;
}

# Object Method
# Title	    :  bin_proteins()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub bin_proteins {
	my $self = shift;
	my ( $rh_ra_id_protprop, $mirror_decoy_ids_from_target_ids,
		$nr_prot_size_bins, $decoy_id_prefix, $special_zero_bin,
		$equinr_bins )
	  = @_;

	# number of total target proteins in each bin
	my $ra_bin_nrprot  = [];
	my $corr_ntp_index = 3;

	# use the corrected number of tryptic peptides as protein sizes
	my @protein_sizes = ();
	foreach ( keys %$rh_ra_id_protprop ) {
		push @protein_sizes, $rh_ra_id_protprop->{$_}[$corr_ntp_index];
	}

	# get a Bins object that contains several Bin objects
	my $bins_ob = Bins->new();
	if ($equinr_bins) {
		$bins_ob->set_equicount_bins( \@protein_sizes, $nr_prot_size_bins,
			$special_zero_bin );
	}
	else {
		$bins_ob->set_equidistant_bins( \@protein_sizes,
			$nr_prot_size_bins, $special_zero_bin );
	}

	# count the number of target proteins in each bin
	for ( my $i = 0 ; $i < $nr_prot_size_bins ; $i++ ) {
		$ra_bin_nrprot->[$i] = 0;
	}
	foreach my $id ( keys %$rh_ra_id_protprop ) {
		my $corr_ntp = $rh_ra_id_protprop->{$id}[$corr_ntp_index];
		my $bin_nr   = $bins_ob->get_bin_nr($corr_ntp);
		push @{ $rh_ra_id_protprop->{$id} }, $bin_nr;
		# only count the non decoy proteins!
		# fixed 20150128 because of Lorenz Blum
		if ( $id !~ /^$decoy_id_prefix/ ) {
			$ra_bin_nrprot->[$bin_nr]++;
		}
		# this was the old code
		#$ra_bin_nrprot->[$bin_nr]++;
	}

	return ( $ra_bin_nrprot, $rh_ra_id_protprop, $bins_ob );
}

# Object Method
# Title	    :  get_protein_sizes()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_protein_sizes {
	my $self = shift;
	my (
		$rh_ra_id_protprop,                $rh_rh_gene_id_seq,
		$min_pep_length,                   $decoy_id_prefix,
		$nr_missed_cleavages,              $min_pep_mass,
		$max_pep_mass,                     $nr_prot_size_bins,
		$mirror_decoy_ids_from_target_ids, $general_min_pep_length
	  )
	  = @_;

	# $rh_ra_id_protprop
	# - aa length
	# - corr aa length

	# peptides smaller then this size will never be identified with
	# high confidence
	$min_pep_length = $general_min_pep_length
	  if $general_min_pep_length > $min_pep_length;
	  
	my $digest = Digest->new();

	# print status
	$| = 1 if $self->{s};
	my $c = 0;
	foreach my $gene ( sort keys %$rh_rh_gene_id_seq ) {
		my @ids = sort keys %{ $rh_rh_gene_id_seq->{$gene} };

		# store the peptides for each protein of a gene group
		my $ra_rh_peptides = [];
		foreach my $id (@ids) {
			$c++;
			print "\r                              \r  $c\t" . $id
			  if $self->{s};
			my $prot_seq = $rh_rh_gene_id_seq->{$gene}{$id};
			my ( $rh_pep, $nr_pep ) =
			  $digest->get_nonred_rh_pep_digest( $prot_seq,
				$min_pep_length, $nr_missed_cleavages, $min_pep_mass,
				$max_pep_mass );
			push @$ra_rh_peptides, $rh_pep;

			# add ntp
			push @{ $rh_ra_id_protprop->{$id} }, $nr_pep;
		}

		# make the correction for gene groups
		for ( my $i = 1 ; $i < @ids ; $i++ ) {
			for ( my $j = 0 ; $j < $i ; $j++ ) {
				my $rh_test_pep = $ra_rh_peptides->[$i];
				my $rh_main_pep = $ra_rh_peptides->[$j];
				foreach my $test_pep ( keys %$rh_test_pep ) {
					if ( exists( $rh_main_pep->{$test_pep} ) ) {
						delete( $rh_test_pep->{$test_pep} );
					}
				}
			}
		}

		# add corr ntp
		for ( my $i = 0 ; $i < @ids ; $i++ ) {
			my $id          = $ids[$i];
			my $nr_corr_pep = keys %{ $ra_rh_peptides->[$i] };

			# add protein length and ntp
			push @{ $rh_ra_id_protprop->{$id} }, $nr_corr_pep;
		}
	}

	return $rh_ra_id_protprop;
}

# Object Method
# Title	    :  group_according_to_gene()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub group_according_to_gene {
	my $self = shift;
	my ( $rh_id_seq, $gene_group_identification_type, $gene_group_regex ) =
	  @_;

	my $rh_rh_gene_id_seq;

	# suffix
	if ($gene_group_identification_type) {

		foreach my $id ( keys %$rh_id_seq ) {
			my $gene = $id;
			$gene =~ s/$gene_group_regex$//;
			$rh_rh_gene_id_seq->{$gene}{$id} = $rh_id_seq->{$id};
		}
	}

	# praefix
	elsif ($gene_group_identification_type) {

		foreach my $id ( keys %$rh_id_seq ) {
			my $gene = $id;
			$gene =~ s/^$gene_group_regex//;
			$rh_rh_gene_id_seq->{$gene}{$id} = $rh_id_seq->{$id};
		}
	}

	# middle
	else {

		foreach my $id ( keys %$rh_id_seq ) {
			my $gene = $id;
			$gene =~ s/$gene_group_regex//;
			$rh_rh_gene_id_seq->{$gene}{$id} = $rh_id_seq->{$id};
		}
	}

	my $nr_groups = keys %$rh_rh_gene_id_seq;
	$self->v_print("$nr_groups protein groups...\n");

	return $rh_rh_gene_id_seq;
}

# Object Method
# Title	    :  correct_identical_sequences()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub correct_identical_sequences {
	my $self = shift;
	my ($rh_id_seq) = @_;

	# properties of the ids
	my $rh_ra_id_protprop;

	# for identical sequences the non alphabetically first
	# sequence is deleted ( = '' )
	my $rh_id_corrseq;

	# add aa length to protein properties
	# and create temporary hash
	my $rh_seq_id;
	while ( my ( $id, $seq ) = each %$rh_id_seq ) {
		push @{ $rh_ra_id_protprop->{$id} }, length($seq);
		push @{ $rh_seq_id->{$seq} },        $id;
	}

	# sort the ids alphabetically
	my $count_total_equal = 0;
	my $count_corrected   = 0;
	while ( my ( $seq, $ra_ids ) = each %$rh_seq_id ) {
		my @sorted_ids = sort @$ra_ids;
		my $nr_ids     = @sorted_ids;
		$count_total_equal += $nr_ids if $nr_ids > 1;

		# loop over all ids for the current sequence
		for ( my $i = 0 ; $i < $nr_ids ; $i++ ) {
			my $corrected_seq = '';

			if ( $i == 0 ) {
				$corrected_seq = $seq;
			}
			else {
				$count_corrected++;
			}

			# correct the protein sequence in the hash and
			# add the corrected aa length to the protein
			$rh_id_corrseq->{ $sorted_ids[$i] } = $corrected_seq;
			push @{ $rh_ra_id_protprop->{ $sorted_ids[$i] } },
			  length($corrected_seq);
		}
	}
	$self->v_print( "$count_total_equal total identical sequences, "
		  . "$count_corrected sequences corrected to ''\n" );

	return ( $rh_id_corrseq, $rh_ra_id_protprop );
}

# Object Method
# Title	    :  check_decoy_mirror()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub check_decoy_mirror {
	my $self = shift;
	my ( $rh_id_seq, $decoy_id_prefix ) = @_;

	my $rh_decoy_ids_without_prefix;
	my $rh_target_ids;

	foreach ( keys %$rh_id_seq ) {
		if ( $_ =~ /^$decoy_id_prefix/ ) {
			$_ =~ /^$decoy_id_prefix(.+)$/;
			$rh_decoy_ids_without_prefix->{$1} = 1;
		}
		else {
			$rh_target_ids->{$_} = 1;
		}
	}

	my $is_mirrored = 1;
	foreach ( keys %{$rh_decoy_ids_without_prefix} ) {
		unless ( exists( $rh_target_ids->{$_} ) ) {
			$is_mirrored = 0;
		}
	}

	# if there is a target id for each decoy id then
	# delete all the decoy sequences and deduce the protein
	# size (and other features) from its target counterpart
	if ($is_mirrored) {
		my $rh_only_target_id_seq;
		foreach ( keys %$rh_target_ids ) {
			$rh_only_target_id_seq->{$_} = $rh_id_seq->{$_};
		}
		return ( 1, $rh_only_target_id_seq );
	}
	else {
		$self->v_print( "some decoy entries do not have "
			  . "corresponding target entries!\n" );
		$self->v_print( "-> protein sizes have to be calculated for "
			  . "target and decoy sequences\n" );
		return ( 0, $rh_id_seq );
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





