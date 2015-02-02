package MayuPepXMLHandler;
use warnings;
use strict;
use File::Basename;
use FileHandle;
use XML::SAX::Base;
use base qw(XML::SAX::Base);    # inherit from XML::SAX::Base

######################################################################
#
# MayuPepXMLHandler
#
# returns:
# --------
# 0 spectrum:
# 1 peptide:                          (raw sequence)
# 2 main id:                          (alphabetically first protein)
# 3 modification info:                $rh_mod pos -> mass
# 4 PeptideProphet probability (PPs):
#
# example:
# --------
# <spectrum_query spectrum="20050415_S2_15_16.12400.12400.3" start_scan="12400"
#  end_scan="12400" precursor_neutral_mass="2587.3920" assumed_charge="3"
#  index="256">
#  <search_result>
#   <search_hit hit_rank="1" peptide="VHPTGVEGGAYFEPAIITGLSDEAR"
#    peptide_prev_aa="R" peptide_next_aa="A" protein="Y69F12A.2a"
#    num_tot_proteins="2" num_matched_ions=" 38" tot_num_ions=" 96"
#    calc_neutral_pep_mass="2586.7920" massdiff="+0.6" num_tol_term="2"
#    num_missed_cleavages="0" is_rejected="0" protein_descr="CE34422
#    WBGene00000118 locus:alh-12     aldehyde dehydrogenase status:Confirmed
#    SW:Q7Z1Q3 protein_id:AAP46268.1">
#    <alternative_protein protein="Y69F12A.2b" protein_descr="CE34423
#     WBGene00000118 locus:alh-12      status:Confirmed SW:Q7Z1Q2
#     protein_id:AAP46269.1" num_tol_term="2"/>
#    <modification_info modified_peptide="IADFGM[147]AKCADNSSKK">
#     <mod_aminoacid_mass position="6" mass="147.192001"/>
#     <mod_aminoacid_mass position="9" mass="330.399109"/>
#    </modification_info>
#    <search_score name="xcorr" value="5.622"/>
#    <search_score name="deltacn" value="0.548"/>
#    <search_score name="deltacnstar" value="0"/>
#    <search_score name="spscore" value="2115.5"/>
#    <search_score name="sprank" value="1"/>
#    <analysis_result analysis="peptideprophet">
#     <peptideprophet_result probability="1.0000"
#      all_ntt_prob="(1.0000,1.0000,1.0000)">
#      <search_score_summary>
#       <parameter name="fval" value="8.2034"/>
#       <parameter name="ntt" value="2"/>
#       <parameter name="nmc" value="0"/>
#       <parameter name="massd" value="0.600"/>
#       <parameter name="icat" value="0"/>
#      </search_score_summary>
#     </peptideprophet_result>
#    </analysis_result>
#   </search_hit>
#  </search_result>
# </spectrum_query>
#
# TODO: - check the add the missing attribute to the sub
#       - 
######################################################################

my @element_stack;      # remembers element names
my %chosen_el_stack;    # remembers chosen element names
my $ra_ra_psm;          # takes results

my $PPs_cutoff;         # input in set()
my $pepxml_file;        # input in set()
my $parser;             # input in set()
my $verbose;            # input in set()
my $status;             # input in set(), print a status while parsing
my $status_prefix = '';
my $counter       = 0;
my $lineMod       = 500;    # every $lineMod print line number

# Object Method
# Title	    :  set()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub set {
	my $self = shift;
	(
		$PPs_cutoff, $pepxml_file, $parser, $verbose, $status,
		$status_prefix
	  )
	  = @_;
	$| = 1 if $status;    # no output buffering for status
}

# Object Method
# Title	    :  set()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub start_document {
	my ( $self, $doc ) = @_;

	# reset the class variables since they exist over the
	# livespan of an object
	@element_stack = ();
	undef %chosen_el_stack;
	@$ra_ra_psm = ();
}

# Object Method
# Title	    :  start_element()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub start_element {
	my ( $self, $properties ) = @_;

	# print progress of parsing
	if ($status) {
		$counter++;
		if ( $counter % $lineMod == 0 ) {
			print "\r                                               \r";
			my @location = $parser->location();
			my $nr       = $location[0]->{'LineNumber'};
			print "  " . $status_prefix . $nr;
		}
	}

	# note: the hash %{$properties} will lose attribute order
	# remember the name by pushing onto the stack
	my $n   = $properties->{'Name'};
	my $att = $properties->{'Attributes'};
	push( @element_stack, $n );

	if ( $n eq 'spectrum_query' ) {

		if ( defined( $att->{spectrum} ) ) {
			$chosen_el_stack{spectrum} = $att->{spectrum};
		}
		else { warnMissingAttribute(); }

	}
	elsif ( $n eq 'search_result' ) {

		# this does not have to be defined in the pepxml file
		if ( exists( $att->{search_id} ) && defined( $att->{search_id} ) )
		{
			$chosen_el_stack{search_id} = $att->{search_id};
		}
		else {
			$chosen_el_stack{search_id} = 1;
		}

	}
	elsif ( $n eq 'search_hit' ) {

		# hit rank
		if ( defined( $att->{hit_rank} ) ) {
			$chosen_el_stack{hit_rank} = $att->{hit_rank};
		}
		else { warnMissingAttribute(); }

		# peptide
		if ( defined( $att->{peptide} ) ) {
			$chosen_el_stack{peptide} = $att->{peptide};
		}
		else { warnMissingAttribute(); }

		# main_protein
		if ( defined( $att->{protein} ) ) {
			$chosen_el_stack{protein} = $att->{protein};
		}
		else { warnMissingAttribute(); }

	}
	elsif ( $n eq 'mod_aminoacid_mass' ) {

		if ( defined( $att->{position} ) && defined( $att->{mass} ) ) {
			$chosen_el_stack{ref_to_hash_with_mod_info}
			  ->{ $att->{position} } = $att->{mass};
		}
		else { warnMissingAttribute(); }

	}
	elsif ( $n eq 'peptideprophet_result'
		&& stack_contains('spectrum_query') )
	{
		$chosen_el_stack{probability} = $att->{probability};
	}

}

# Object Method
# Title	    :  end_element()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub end_element {
	my ( $self, $properties ) = @_;
	my $n = $properties->{'Name'};

	if ( $n eq 'search_hit' ) {
		if (   defined( $chosen_el_stack{spectrum} )
			&& defined( $chosen_el_stack{peptide} )
			&& defined( $chosen_el_stack{protein} )
			&& defined( $chosen_el_stack{ref_to_hash_with_mod_info} )
			&& defined( $chosen_el_stack{probability} ) )
		{
			if (   $chosen_el_stack{probability} >= $PPs_cutoff
				&& $chosen_el_stack{hit_rank} == 1 )
			{

				push @$ra_ra_psm,
				  [
					$chosen_el_stack{spectrum},
					$chosen_el_stack{peptide},
					$chosen_el_stack{protein},
					$chosen_el_stack{ref_to_hash_with_mod_info},
					$chosen_el_stack{probability}
				  ];
			}

		}

		# delete these entries because some tags don't have all of them
		# this would lead to a carry over from an old tag
		delete( $chosen_el_stack{peptide} );
		delete( $chosen_el_stack{protein} );
		delete( $chosen_el_stack{ref_to_hash_with_mod_info} );
		$chosen_el_stack{ref_to_hash_with_mod_info} = {};
		delete( $chosen_el_stack{probability} );
		delete( $chosen_el_stack{hit_rank} );

	}

	pop(@element_stack);
}

# Object Method
# Title	    :  return_type()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub return_type {
	my $self = shift;

	return [ 'scan', 'pepseq', 'main_cds', 'rh_mod', 'ds' ];
}

# Object Method
# Title	    :  name()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub name {
	my $self = shift;
	return 'PepXMLHandler';
}

# Object Method
# Title	    :  entity_reference()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub entity_reference {
	my ( $self, $properties ) = @_;
}

# Object Method
# Title	    :  end_document()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub end_document {
	my ( $self, $properties ) = @_;
	print "\r                                  \r" if $status;
}

# Object Method
# Title	    :  warnMissingAttribute()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub warnMissingAttribute {
	n_print("missing attribute!\n");
	my @location = $parser->location();
	n_print( "line number: " . $location[0]->{'LineNumber'} . "\n" );
	exit;
}

# Object Method
# Title	    :  results()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub results {
	my $self = shift;
	my @copy = @$ra_ra_psm;    # $ra_ra_psm
	return \@copy;
}

# Class Method
# Title	    :  stack_top()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub stack_top {
	my $guess = shift;
	return $element_stack[$#element_stack] eq $guess;
}

# Class Method
# Title	    :  stack_contains()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub stack_contains {
	my $guess = shift;
	foreach (@element_stack) {
		return 1 if ( $_ eq $guess );
	}
	return 0;
}

sub n_print {
	print "  " . $_[0];
}

sub v_print {
	print "  " . $_[0] if $verbose;
}

1;

__END__

Lukas Reiter


