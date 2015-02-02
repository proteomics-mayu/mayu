package MayuPepXMLParser;
use warnings;
use strict;
use File::Basename;
use FileHandle;

######################################################################
#
# MayuPepXMLParser
# A parser that does not need a PepXML parser from cpan but
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
######################################################################

# Constructor
# Title     :  new
# Usage     :  my $p = MayuPepXMLParser->new();
# Function  :
# Returns   :
# Args      :
sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my ( $initial_pps_cutoff, $pepxml, $verbose, $status, $status_prefix )
	  = @_;

	# cutoff
	if ( defined($initial_pps_cutoff) ) {
		$self->{cutoff} = $initial_pps_cutoff;
	}
	else {
		$self->{cutoff} = 0;
	}

	# inputfile
	if ( defined($pepxml) ) {
		$self->{input_pepxml} = $pepxml;
	}
	else {
		$self->{input_pepxml} = '';
		print "  no input file!\n";
	}

	# printing
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
	if ( defined($status_prefix) ) {
		$self->{s_prefix} = $status_prefix;
	}
	else {
		$self->{s_prefix} = '';
	}

	return $self;
}

# Object Method
# Title	    :  parse()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub parse {
	my $self = shift;

	my $ra_ra_psm = [];

	# no output buffering if line number is to be printed
	my $status        = $self->{s};
	my $status_prefix = $self->{s_prefix};
	$| = 1 if $status;
	
	my $foundip = 0;

	if ( -e $self->{input_pepxml} ) {

		# INPUT_RECORD_SEPARATOR
		$/ = '</spectrum_query>';
		open( I, $self->{input_pepxml} ) or warn $!;
		while (<I>) {
			my $line = $_;
			$line =~ s/\n//g;

			my ($spec);
			if ( $line =~ /spectrum="([^"]+)"/ ) {
				$line =~ /spectrum="([^"]+)"/;
				$spec = $1;
			}

			# split into search hits
			my @search_hits = split( "<search_hit ", $line );

			# remove the spectrum
			shift @search_hits;

			# extract the search hit information
			foreach my $search_hit (@search_hits) {
				my ( $pep, $id, $pps );
				my $rh_mod = {};
				if ( $search_hit =~ /[^<>]*\speptide="([^"]+)"/ ) {
					$search_hit =~ /[^<>]*\speptide="([^"]+)"/;
					$pep = $1;
				}
				if ( $search_hit =~ /[^<>]*\sprotein="([^"]+)"/ ) {
					$search_hit =~ /[^<>]*\sprotein="([^"]+)"/;
					$id = $1;
				}
				if ( $search_hit =~
					/<peptideprophet_result probability="([^"]+)"/ )
				{
					$search_hit =~
					  /<peptideprophet_result probability="([^"]+)"/;
					$pps = $1;
				}
				if ( $search_hit =~
					/<interprophet_result probability="([^"]+)"/ )
				{
					if(!$foundip) 
                    { 
                        print "  Found IprophetProbability. Using this instead of PeptideProphetProbability.\n";
                        print "  Caveat: Do not mix iprophet and peptideprophet files, there is no check!\n";
                        $foundip = 1; 
                    }
					$search_hit =~
					/<interprophet_result probability="([^"]+)"/;
                    $pps = $1;
				}
				my @mod = split( "<mod_aminoacid_mass ", $search_hit );
				foreach my $mod (@mod) {
					if (   $mod =~ /position="([^"]+)"/
						&& $mod =~ /mass="([^"]+)"/ )
					{
						$mod =~ /position="([^"]+)"/;
						my $pos = $1;
						$mod =~ /mass="([^"]+)"/;
						my $mass = $1;
						$rh_mod->{$pos} = $mass;
					}
				}				
				if (   defined($pps)
					&& defined($spec)
					&& defined($pep)
					&& defined($id) )
				{
					if ( $pps >= $self->{cutoff} ) {
						push @$ra_ra_psm,
						  [ $spec, $pep, $id, $rh_mod, $pps ];
					}
				}
			}
			print "\r                                                \r"
			  if $status;
			print "  " . $status_prefix . $. if $status;
		}
		close(I);
		print "\r" if $status;
	}

	return $ra_ra_psm;
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

sub n_print {
	print "  " . $_[0];
}

1;

__END__

Lukas Reiter



