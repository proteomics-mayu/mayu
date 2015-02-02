#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Cwd;

##################################################################
#
# Mayu is a software package to determine the protein
# identification false discovery rate (protFDR) additionally to
# the peptide-spectrum match false discovery rate (mFDR).
# For a comprehensive help about the currently implemented
# features type:
#
# perl Mayu.pl -manual
#
#
# Lukas Reiter
# Manfred Claassen
#
##################################################################
#
# SOFTWARE
#
#	Lukas Reiter - Hengartner Laboratory
#	  lukas.reiter@molbio.uzh.ch
#	  Institute of Molecular Biology
#	  Winterthurerstrasse 190
#	  University of Zuerich - Irchel
#	  CH-8057 Zürich
#	+++++++++++++++++++++++++++++++++++++++++
#	Located at:
#	  Institute for Molecular Systems Biology
#	  Aebersold Laboratory
#	  Wolfgang-Pauli-Str. 16
#	  ETH Hönggerberg, HPT C 75
#	  CH-8093 Zürich
#	  Tel: +41 44 633 39 45
#
##################################################################

##################################################################
# INPUT AND OPTIONS
##################################################################

# no output buffering
$| = 1;

my $version = '1.06';

#-------------------------------------------------------------
# protein size estimation
#-------------------------------------------------------------
my $mirror_decoy_ids_from_target_ids = 1;

# gene group identification type
# 0: prefix
# 1: suffix
# 2: middle
my $gene_group_identification_type = 1;
my $gene_group_regex               = '[a-z]';

# use a special bin for zero protein size
my $special_zero_bin = 1;

# peptides smaller than this size will never be identified with
# high confidence and are therefore neglected in the protein
# size estimation (additionally to the mass window)
my $general_min_pep_length = 4;

#-------------------------------------------------------------
# input PSM parsing
#-------------------------------------------------------------
# initial discriminant score cutoff
# (PeptideProphet probability cutoff)
my $initial_ds_cutoff = 0;

# used for printing out the number of identifications
my $status_fdr = 0.01;

# file ending for the table input files
my $table_file_ending = '.csv';

# separator for table input (and output)
my $table_separator = ',';

# 1: a header is assumed to be there
# 0: there is a check for a header
my $csv_header = 0;

# 1: filter the PSM for only the protein ids fround in the
# provided fasta database
# 0: alternatively use an average protein length for binning
# the protein identifications into protein search space
# size bins.
my $filter_ids = 0;

#-------------------------------------------------------------
# id selection
#-------------------------------------------------------------
# 0: all data and mFDR range (default)
# 1: cumulative input files and mFDR range
# 2: cumulative shuffled input files and mFDR range
# 3: cumulative runs and mFDR range
# 4: cumulative shuffled runs and mFDR range
# 5: cumulative runs (orthogonality sorted) and mFDR range
# (10: phospho peptides (orthogonality sorted) and mFDR range)
my $rh_id_selection_helper = {
	'0'  => 'all data and mFDR range (default)',
	'1'  => 'cumulative input files and mFDR range',
	'2'  => 'cumulative shuffled input files and mFDR range',
	'3'  => 'cumulative runs and mFDR range',
	'4'  => 'cumulative shuffled runs and mFDR range',
	'5'  => 'cumulative runs (orthogonality sorted) and mFDR range',
	'10' => 'phospho peptides (orthogonality sorted) and mFDR range',
};

#-------------------------------------------------------------
# PSM filtering
#-------------------------------------------------------------
my $rh_psm_target_decoy = {
	't'  => 'target PSM',
	'td' => 'target and decoy PSM',
	'd'  => 'decoy PSM'
};

#-------------------------------------------------------------
# variables for command line input
#-------------------------------------------------------------
my (
	$pepxml_in,          $csv_in,           $db,
	$min_pep_length,     $decoy_id_prefix,  $tar_dec_ratio,
	$max_psm_fdr,        $psm_fdr_steps,    $nr_missed_cleavages,
	$min_pep_mass,       $max_pep_mass,     $nr_prot_size_bins,
	$file_name_base,     $id_set_selection, $nr_runs_steps,
	$ids_out_tag,        $equidist_bins,    $equicount_bins,
	$dont_corr_id_seq,   $run_r,            $use_xml_parser,
	$p_mfdr,             $p_bin_prot,       $p_prot_feat,
	$print_input_output, $print_cumulative, $remove_ambiguous, $v,
	$s,                  $help,             $manual
);

# options are described in the usage() sub: perl Mayu.pl -h
GetOptions(
	'A=s'       => \$pepxml_in,            # input pepxml.xml file or dir
	'B=s'       => \$csv_in,               # input table.csv, comma sep
	'C=s'       => \$db,                   # target-decoy database
	'D=i'       => \$min_pep_length,       # peptides have to be >= length
	'E=s'       => \$decoy_id_prefix,      # id prefix for decoy hits
	'F=f'       => \$tar_dec_ratio,        # target to decoy ratio (e.g. 1)
	'G=f'       => \$max_psm_fdr,          # maximal PSM FDR
	'H=i'       => \$psm_fdr_steps,        # number of analysis steps
	'I=i'       => \$nr_missed_cleavages,  # number of missed cleavages
	'J=f'       => \$min_pep_mass,         # minimal peptide mass
	'K=f'       => \$max_pep_mass,         # maximal peptide mass
	'L=i'       => \$nr_prot_size_bins,    # number of protein size bins
	'M=s'       => \$file_name_base,       # use this as file name base
	'N=i'       => \$id_set_selection,     # id set selection
	'O=i'       => \$nr_runs_steps,        # cumulative runs
	'P=s'       => \$ids_out_tag,          # output filtered ids
	'equidist'  => \$equidist_bins,        # equicount / equidist
	'dcis'      => \$dont_corr_id_seq,     # correct for identical seq
	'xmlparser' => \$use_xml_parser,       # use a proper xml parser
	'runR'      => \$run_r,                # run R to finish analysis
	'PmFDR'     => \$p_mfdr,               # print the mFDR files
	'PbinProt'  => \$p_bin_prot,           # print bin protFDR files
	'PprotFeat' => \$p_prot_feat,          # print prot feature files
	'Pio'       => \$print_input_output,   # print pepxml input to csv
	'cumul'     => \$print_cumulative,     # print cumulative input
	'remamb'    => \$remove_ambiguous,     # remove ambiguous PSMs
	'verbose'   => \$v,
	'status'    => \$s,
	'help'      => \$help,
	'manual'    => \$manual
);

##################################################################
# SET DEFAULT VALUES OR CHANGE IMPOSSIBLE VALUES
##################################################################
$pepxml_in           = ''     unless defined($pepxml_in);
$csv_in              = ''     unless defined($csv_in);
$db                  = ''     unless defined($db);
$min_pep_length      = 0      unless defined($min_pep_length);
$min_pep_length      = 0      unless $min_pep_length > 0;
$decoy_id_prefix     = 'rev_' unless defined($decoy_id_prefix);
$tar_dec_ratio       = 1      unless defined($tar_dec_ratio);
$max_psm_fdr         = 0.01   unless defined($max_psm_fdr);
$max_psm_fdr         = 1      unless $max_psm_fdr <= 1;
$psm_fdr_steps       = 11     unless defined($psm_fdr_steps);
$psm_fdr_steps       = 11     unless $psm_fdr_steps >= 2;
$nr_missed_cleavages = 0      unless defined($nr_missed_cleavages);
$min_pep_mass        = 400    unless defined($min_pep_mass);
$min_pep_mass        = 0      unless $min_pep_mass >= 0;
$max_pep_mass        = 6000   unless defined($max_pep_mass);
$max_pep_mass        = 0      unless $max_pep_mass >= 0;
$nr_prot_size_bins   = 10     unless defined($nr_prot_size_bins);
$nr_prot_size_bins   = 10     unless $nr_prot_size_bins >= 1;
$id_set_selection    = 0      unless defined($id_set_selection);
$nr_runs_steps       = 2      unless defined($nr_runs_steps);
$nr_runs_steps       = 2      unless $nr_runs_steps >= 2;
$ids_out_tag         = ''     unless defined($ids_out_tag);
$equidist_bins       = 0      unless defined($equidist_bins);
$equicount_bins      = 1      unless $equidist_bins == 1;
$equicount_bins      = 0      unless $equidist_bins == 0;
$dont_corr_id_seq    = 0      unless defined($dont_corr_id_seq);
$use_xml_parser      = 0      unless defined($use_xml_parser);
$run_r               = 0      unless defined($run_r);
$p_mfdr              = 0      unless defined($p_mfdr);
$p_bin_prot          = 0      unless defined($p_bin_prot);
$p_prot_feat         = 0      unless defined($p_prot_feat);
$print_input_output  = 0      unless defined($print_input_output);
$print_cumulative    = 0      unless defined($print_cumulative);
$remove_ambiguous    = 0      unless defined($remove_ambiguous);
$v                   = 0      unless defined($v);
$s                   = 0      unless defined($s);
$help                = 0      unless defined($help);
$manual              = 0      unless defined($manual);

##################################################################
# CHECK INPUT
##################################################################
if ($manual) {
	manual();
	exit;
}
elsif ($help) {
	usage();
	exit;
}
elsif ( ( -e $pepxml_in || -e $csv_in ) && -e $db ) {
	main();
}
else {
	print "\n\n  need existing input for (-A or -B) and (-C)!\n\n";
	print "  -A: '$pepxml_in'\n";
	print "  -B: '$csv_in'\n";
	print "  -C: '$db'\n\n";
	usage();
}

##################################################################
# MAIN
##################################################################

# Title	    :  main()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub main {

	#-------------------------------------------------------------
	# load modules and show used options
	#-------------------------------------------------------------
	# directory where the script is located
	my $mayu_dir = dirname($0);
	$mayu_dir =~ s/\\/\//g;
	
	# current working directory
	my $wd = getcwd();

	# wd needs to be the Mayu directory for proper loading
	use lib "lib";

	# load modules
	use MayuTools;
	use MascotCSV;
	use MayuProteinSize;
	use PSMSet;
	
	# proper xml parser needs to be installed
	if ($use_xml_parser) {
		my $parser = "XML::Parser::PerlSAX";
		eval("use $parser;");
		die("  Could not find $parser!") if $@ ne "";
		my $handler = "MayuPepXMLHandler";    # XML::SAX::Base
		eval("use $handler;");
		die("  Could not find $handler!") if $@ ne "";
	}
	# no installation required
	else {
		use MayuPepXMLParser;
	}
	use TandemMSIdSelectionScheme;
	use MayuManager;
	use PSMFDR;
	use PeptideIdFDR;
	use ProteinIdFDR;
	use LocalProteinIdFDR;
	use RPlots;

	# helper subroutines
	my $tools = MayuTools->new( 0, 0 );
	$tools->set_time();

	# output file name base
	my $out_base  = '';
	my $timestamp = $tools->get_file_name_base();
	if ( defined($file_name_base) ) {
		$out_base = $file_name_base;
	}
	else {
		$out_base = $timestamp;
	}

	# print all input and main options
	print_options($out_base);

	#-------------------------------------------------------------
	# estimate the protein sizes that are accessible to the
	# search engine. Return BinnedEntity object(s)
	#-------------------------------------------------------------
	print "  ------------------------------------\n" if $v;
	print "  protein size\n"                         if $v;
	print "  ------------------------------------\n" if $v;
	my $Mayu_prot_size = MayuProteinSize->new( $v, $s );
	my ( $bin_prot, $bin_pep ) = $Mayu_prot_size->protein_size_analysis(
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
		$equicount_bins,
		$filter_ids
	);
	print "\n" if $v;

	# TODO - add static modifications and types of modifications
	#        to the PSMSet class
	#      - only possible for pepxml, not for other input formats
	# <aminoacid_modification aminoacid="M" massdiff="15.9994"
	#	mass="147.1920" variable="Y" symbol="*"/>
	# <aminoacid_modification aminoacid="C" massdiff="8.9339"
	#	mass="339.3330" variable="Y" symbol="#"/>
	# <aminoacid_modification aminoacid="C" massdiff="227.2603"
	#	mass="330.3991" variable="N"/>
	#-------------------------------------------------------------
	# get input PSM in the form of a PSMSet object
	#-------------------------------------------------------------
	print "  ------------------------------------\n" if $v;
	print "  input PSM (mFDR)\n"                     if $v;
	print "  ------------------------------------\n" if $v;
	my ( $psm_set, $ra_roc_files ) = get_psm_set(
		$tools,             $pepxml_in,         $csv_in,
		$initial_ds_cutoff, $decoy_id_prefix,   $tar_dec_ratio,
		$max_psm_fdr,       $table_file_ending, $table_separator,
		$out_base,          $bin_prot, $db
	);
	print "\n" if $v;
	if ( $psm_set->is_empty() ) {
		print "  no input data! "
		  . "(file not existing or data of low quality)\n";
		exit;
	}

	#-------------------------------------------------------------
	# MayuManager
	# set the error models to the MayuManager and add binning
	# entities if the model requires them.
	# MayuManager will initialize the error models together with
	# the PSMSets created by the TandemMSIdSelectionScheme.
	# All the error models provide similar functionality for
	# - error calculation
	# - local error calculation
	# - data querying
	#
	# some error models:
	# - PSMFDR
	# - PeptideIdFDR
	# - ProteinIdFDR
	# (- NaivePeptideIdFDR)
	# (- NaiveProteinIdFDR)
	#
	# functions provided by all of the error models:
	# - get_error_model_identifier()
	# - create_error_model()
	# - set_psm_set(), ( $ra_ra with 6 columns )
	# - get_error_summary()
	# functions provided by some of the error models:
	# - get_error_bin_table()
	#-------------------------------------------------------------
	# TandemMSIdSelectionScheme
	# selects the identifications and passes them to the
	# MayuManager
	#-------------------------------------------------------------
	print "  ------------------------------------\n" if $v;
	print "  protein identification FDR (protFDR)\n" if $v;
	print "  ------------------------------------\n" if $v;
	my $p_manager = MayuManager->new( $v, $s );
	$p_manager->set_protein_features($bin_prot);
	my $feat_base = $out_base . '_feat_prot_' . $version;
	$p_manager->set_protein_feature_output_base($feat_base)
	  if $p_prot_feat;

	# set the error models
	my $psm_em = PSMFDR->new( $v, $s, $tar_dec_ratio );
	$p_manager->set_error_model($psm_em);
	my $pep_em = PeptideIdFDR->new( $v, $s, $bin_pep, $tar_dec_ratio );
	$p_manager->set_error_model($pep_em);
	my $prot_em = ProteinIdFDR->new( $v, $s, $bin_prot, $tar_dec_ratio );
	$p_manager->set_error_model($prot_em);
	my $lp_em = LocalProteinIdFDR->new( $v, $s, $tar_dec_ratio );
	$p_manager->set_error_model($lp_em);

	# print out the used error models
	$p_manager->print_registered_error_models() if $v;

	# decide on the id selection scheme
	my $sel_scheme = TandemMSIdSelectionScheme->new( $v, $s );

	#  0: all data and mFDR range (default)
	#  1: cumulative input files and mFDR range
	#  2: cumulative shuffled input files and mFDR range
	#  3: cumulative runs and mFDR range
	#  4: cumulative shuffled runs and mFDR range
	#  5: cumulative runs (orthogonality sorted) and mFDR range
	# 10: phospho peptides (orthogonality sorted) and mFDR range
	if ( $id_set_selection == 1 ) {
		print "\n  DATA SELECTION SCHEME:\n" if $v;
		print "  cumulative input files" . " -> set of mFDR\n\n"
		  if $v;
		$sel_scheme->cumulative_file_mFDR_range( $p_manager, $psm_set,
			$max_psm_fdr, $psm_fdr_steps );
	}
	elsif ( $id_set_selection == 2 ) {
		print "\n  DATA SELECTION SCHEME:\n" if $v;
		print "  cumulative shuffled input files" . " -> set of mFDR\n\n"
		  if $v;
		$sel_scheme->cumulative_shuffled_file_mFDR_range( $p_manager,
			$psm_set, $max_psm_fdr, $psm_fdr_steps );
	}
	elsif ( $id_set_selection == 3 ) {
		print "\n  DATA SELECTION SCHEME:\n" if $v;
		print
		  "  cumulative runs in $nr_runs_steps steps (sorted according "
		  . "to files and then runs) -> set of mFDR\n"
		  if $v;
		$sel_scheme->cumulative_run_mFDR_range(
			$p_manager,     $psm_set, $max_psm_fdr,
			$psm_fdr_steps, $nr_runs_steps
		);
	}
	elsif ( $id_set_selection == 4 ) {
		print "\n  DATA SELECTION SCHEME:\n" if $v;
		print "  cumulative shuffled runs in $nr_runs_steps steps "
		  . " -> set of mFDR\n\n"
		  if $v;
		$sel_scheme->cumulative_shuffled_run_mFDR_range(
			$p_manager,     $psm_set, $max_psm_fdr,
			$psm_fdr_steps, $nr_runs_steps
		);
	}
	elsif ( $id_set_selection == 5 ) {
		print "\n  DATA SELECTION SCHEME:\n" if $v;
		print "  cumulative runs in $nr_runs_steps steps (sorted according"
		  . " to orthogonality) -> set of mFDR\n\n"
		  if $v;
		$sel_scheme->ortho_run_mFDR_range( $p_manager, $psm_set,
			$max_psm_fdr, $psm_fdr_steps, $nr_runs_steps, $out_base,
			$version );
	}

	# TODO implement
	elsif ( $id_set_selection == 10 ) {
		print "\n  DATA SELECTION SCHEME:\n"                   if $v;
		print "  complete phospho data set -> set of mFDR\n\n" if $v;
		$sel_scheme->phospho_id_mFDR_range( $p_manager, $psm_set,
			$max_psm_fdr, $psm_fdr_steps );
	}

	# this is the default
	else {
		print "\n  DATA SELECTION SCHEME:\n"           if $v;
		print "  complete data set -> set of mFDR\n\n" if $v;
		$sel_scheme->all_id_mFDR_range( $p_manager, $psm_set, $max_psm_fdr,
			$psm_fdr_steps );
	}
	print "  $feat_base\n" if $p_prot_feat && $v;
	print "\n"             if $v;

	#-------------------------------------------------------------
	# print out
	#-------------------------------------------------------------
	print "  ------------------------------------\n" if $v;
	print "  printing files\n"                       if $v;
	print "  ------------------------------------\n" if $v;

	# print out a file with local FDR of proteins with similar 
	# sizes if this error model was registered
	if ( $p_manager->error_model_exists('ProteinIdFDR') && $p_bin_prot ) {
		my $psize_locFDR_base = $out_base 
			. '_prot_size_local_FDR_' . $version;
		my ( $prot_size_bin_csv, $prot_size_bin_txt ) =
		  print_prot_size_bin_file( $p_manager, $psize_locFDR_base, $tools );
		print "  $prot_size_bin_csv, $prot_size_bin_txt\n" if $v;
	}

	# print out main summary file
	my $main_base = $out_base . '_main_' . $version;
	my ( $fdr_csv, $fdr_txt ) =
	  print_main_file( $p_manager, $main_base, $tools );
	print "  main output files:\n" if $v;
	print "  $fdr_csv, $fdr_txt\n" if $v;
	print "\n"                     if $v;

	# print out a table of filtered psm
	# e.g. mFDR=0.01:t or pepFDR=0.02:td or protFDR=0.05:t
	unless ( $ids_out_tag =~ /^$/ ) {
		my $psm_file_base = $out_base . '_psm_';
		my (
			$id_csv_file, $nr_psm,       $fdr_type,
			$fdr_value,   $target_decoy, $mFDR
		  )
		  = print_psm_file( $psm_set, $ids_out_tag, $p_manager,
			$psm_file_base, $max_psm_fdr, $tools );
		print "  $nr_psm PSM printed out\n" if $v;
		print "\n"                          if $v;
	}

	#-------------------------------------------------------------
	# R plots
	#-------------------------------------------------------------
	if ( $run_r ) {
		print "  ------------------------------------\n" if $v;
		print "  R plots\n"                              if $v;
		print "  ------------------------------------\n" if $v;
		my $r_plots          = RPlots->new( $v, $s );
		my $main_template_in = $mayu_dir . '/templates/main.R';

		# prepare the path of the r library
		my $r_lib = $mayu_dir . '/r_lib/';
		
		#-------------------------------------------------------------
		# PmFDR
		#-------------------------------------------------------------
		if ( $p_mfdr ) {
						
			my $roc_files_r_vector = 'c("' . join('","', @$ra_roc_files) . '")';
		
			# prepare the file name of the r .pdf output file
			my $r_out_base = $out_base . '_ds_vs_mFDR_' . $version;
			my $r_out  = $r_out_base . '.pdf';
			
			# create replace hash
			my $rh_replace = {
				'<REPLACE_WITH_LIB_PATH>' => $r_lib,
				'"<REPLACE_WITH_ROC_FILES>"' => $roc_files_r_vector,
				'<REPLACE_WITH_OUT_PATH>' => getcwd(),
				'<REPLACE_WITH_OUTFILE>' => $r_out
			};
			
			my $template_in = $mayu_dir . '/templates/PmFDR.R';
			
			# if not defined then the new template is written to the cwd
			my $template_out_dir = undef;
	
			$r_plots->run_r_template( $template_in, $rh_replace,
				$template_out_dir );
		}
		
		#-------------------------------------------------------------
		# main
		#-------------------------------------------------------------
		
		# prepare the file name of the r .pdf output file
		my $main_r_out_base = $out_base . '_main_' . $version;
		my $main_r_out  = $main_r_out_base . '.pdf';
		
		# create replace hash
		my $rh_replace = {
			'<REPLACE_WITH_LIB_PATH>' => $r_lib,
			'<REPLACE_WITH_INPUT_FILE>' => $fdr_csv,
			'<REPLACE_WITH_OUT_PATH>' => getcwd(),
			'<REPLACE_WITH_OUTFILE>' => $main_r_out
		};

		# if not defined then the new template is written to the cwd
		my $template_out_dir = undef;

		$r_plots->run_r_template( $main_template_in, $rh_replace,
			$template_out_dir );
	}

	print "\n";
	print "  " . $tools->get_time() . " seconds run time\n";
	print "\n";

}

# Title	    :  print_psm_file()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub print_psm_file {
	my ( $psm_set, $ids_out_tag, $p_manager, $out_base, $max_psm_fdr,
		$tools ) = @_;

	# extract the type of filtering
	my ( $fdr_type, $fdr_value, $target_decoy ) =
	  get_ids_filtering($ids_out_tag);

	# output file name
	my $id_csv_file = $out_base
	  		. $fdr_type
	  		. $fdr_value . '_'
	  		. $target_decoy . '_'
	  		. $version . '.csv';

	# determine the corresponding mFDR cutoff
	my $mFDR = 0.01;    # default
	if ( $fdr_type ne 'mFDR' ) {
		my $header          = 1;
		my $ra_ra_fdr_table = $p_manager->get_table($header);

		# remove the header and get index of columns
		my $ra_header = shift @$ra_ra_fdr_table;
		my ( $mfdr_col, $fdr_col );
		for ( my $i = 0 ; $i < @$ra_header ; $i++ ) {
			if ( $ra_header->[$i] =~ /^mFDR$/ ) {
				$mfdr_col = $i;
			}
			if ( $ra_header->[$i] =~ /^$fdr_type$/ ) {
				$fdr_col = $i;
			}
		}

		my $ra_ra_FDR_mFDR = [];
		foreach my $ra (@$ra_ra_fdr_table) {
			push @$ra_ra_FDR_mFDR, [ $ra->[$fdr_col], $ra->[$mfdr_col] ];
		}

		# get an mFDR cutoff
		my @mFDR =
		  $tools->linear_inter_or_extrapolate( $fdr_value,
			$ra_ra_FDR_mFDR );
		
		$mFDR = $mFDR[0];
	}
	else {
		$mFDR = $fdr_value;
	}
	
	# check whether the selected PSM can be printed out
	# but print out the PSM anyway
	if ( $mFDR > $max_psm_fdr ) {
		print
		  "  your desired mFDR $mFDR for output (-P ...) is beyond the "
		  . "chosen mFDR $max_psm_fdr used for calculations (-G ...) !\n";
		print "  use a higher mFDR with -G for proper output "
		  . "(e.g. -G $mFDR)!\n";
		print "\n";
	}
	
	# print out the selected printing of PSM
	print "  filtered data:\n"                                 if $v;
	print "  printing a set of filtered PSM to $id_csv_file\n" if $v;
	print "  filtering PSM with $fdr_type = $fdr_value, "
	  . "corresponding mFDR = $mFDR\n"
	  if $v;
	print "  keeping " . $rh_psm_target_decoy->{$target_decoy} . "\n"
	  if $v;

	# extract the data from $psm_set
	my $ra_ra_psm = [];
	my $nr_psm    = 0;
	if ( $target_decoy =~ /^t$/ ) {
		$ra_ra_psm = $psm_set->get_target_psm_by_fdr($mFDR);
		$nr_psm    = $psm_set->get_nr_target_psm_by_fdr($mFDR);
	}
	elsif ( $target_decoy =~ /^td$/ ) {
		$ra_ra_psm = $psm_set->get_psm_by_fdr($mFDR);
		$nr_psm    = $psm_set->get_nr_psm_by_fdr($mFDR);
	}
	else {
		$ra_ra_psm = $psm_set->get_decoy_psm_by_fdr($mFDR);
		$nr_psm    = $psm_set->get_nr_decoy_psm_by_fdr($mFDR);
	}

	open( O, ">$id_csv_file" ) or warn $!;
	print O join( $table_separator,
		( 'scan', 'pep', 'prot', 'mod', 'score', 'decoy', 'mFDR' ) )
	  . "\n";
	foreach my $ra_psm (@$ra_ra_psm) {
		$ra_psm->[3] = get_string_from_mod( $ra_psm->[3] );
		print O join( $table_separator, @$ra_psm ) . "\n";
	}
	close(O);

	return (
		$id_csv_file, $nr_psm,       $fdr_type,
		$fdr_value,   $target_decoy, $mFDR
	);
}

# Title	    :  get_ids_filtering()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_ids_filtering {
	my ($ids_out_tag) = @_;

	# default values
	# filter with 0.05 on protFDR and report only target hits
	my ( $fdr_type, $fdr_value, $target_decoy ) = ( 'protFDR', 0.05, 't' );

	my @parts = split( ':', $ids_out_tag );
	if ( @parts == 2 ) {
		my @fdr_parts = split( '=', $parts[0] );
		if ( @fdr_parts == 2 ) {
			if (   $fdr_parts[0] =~ /^mFDR$/
				|| $fdr_parts[0] =~ /^pepFDR$/
				|| $fdr_parts[0] =~ /^protFDR$/ )
			{
				$fdr_type = $fdr_parts[0];
			}
			if ( $fdr_parts[1] =~ /^[\.\d]+$/ ) {
				if ( $fdr_parts[1] >= 0 && $fdr_parts[1] <= 1 ) {
					$fdr_value = $fdr_parts[1];
				}
			}
		}

		if (   $parts[1] =~ /^t$/
			|| $parts[1] =~ /^d$/
			|| $parts[1] =~ /^td$/ )
		{
			$target_decoy = $parts[1];
		}
	}

	return ( $fdr_type, $fdr_value, $target_decoy );
}

# Title	    :  print_prot_size_bin_file()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub print_prot_size_bin_file {
	my ( $p_manager, $out_base, $tools ) = @_;

	my $csv_sep        = "\t";
	my $prot_size_csv  = $out_base . '.csv';
	my $prot_size_txt  = $out_base . '.txt';

	my $fh_prot_size_bin_csv = FileHandle->new();
	$fh_prot_size_bin_csv->open(">$prot_size_csv") or die $!;

	my $header                    = 1;
	my $ra_ra_prot_size_bin_table =
	  $p_manager->get_prot_bin_table($header);
	foreach (@$ra_ra_prot_size_bin_table) {
		print $fh_prot_size_bin_csv join( $csv_sep, @$_ ) . "\n";
	}

	$fh_prot_size_bin_csv->close();

	# print files with equally spaced column
	$tools->print_data_in_table_style( $ra_ra_prot_size_bin_table,
		$prot_size_txt );

	return ( $prot_size_csv, $prot_size_txt );
}

# Title	    :  print_main_file()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub print_main_file {
	my ( $p_manager, $out_base, $tools ) = @_;

	my $csv_sep  = "\t";
	my $fdr_csv  = $out_base . '.csv';
	my $fdr_txt  = $out_base . '.txt';

	my $fh_fdr_csv = FileHandle->new();
	$fh_fdr_csv->open(">$fdr_csv") or die $!;

	my $header          = 1;
	my $ra_ra_fdr_table = $p_manager->get_table($header);
	foreach (@$ra_ra_fdr_table) {
		print $fh_fdr_csv join( $csv_sep, @$_ ) . "\n";
	}

	$fh_fdr_csv->close();

	# print files with equally spaced column
	$tools->print_data_in_table_style( $ra_ra_fdr_table, $fdr_txt );

	return ( $fdr_csv, $fdr_txt );
}

# Title	    :  get_psm_set()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_psm_set {
	my (
		$tools,             $pepxml_in,         $csv_in,
		$initial_ds_cutoff, $decoy_id_prefix,   $tar_dec_ratio,
		$max_psm_fdr,       $table_file_ending, $table_separator,
		$out_base,          $bin_prot, $db
	  )
	  = @_;

	# combine the PSM of all the input files after filtering
	# and processing in this object
	my $complete_psm_set = PSMSet->new( $v, $s );
	my @roc_files = ();

	if ( -e $pepxml_in ) {
		my @pepxml = $tools->get_files( $pepxml_in, '.xml' );
		my $nr_pepxml = @pepxml;
		print "  no xml input found!\n" if $nr_pepxml == 0;
		foreach my $pepxml ( sort @pepxml ) {
			my ( $psm_set, $roc_file ) =
			  get_psm_set_from_pepxml( $pepxml, $initial_ds_cutoff,
				$decoy_id_prefix, $tar_dec_ratio, $max_psm_fdr, $tools,
				$out_base, $db );
			push @roc_files, $roc_file;

			$complete_psm_set->add_psm( $psm_set->get_psm(), $pepxml,
				$psm_set->get_runs() );

			if ($print_cumulative) {
				print "  cumulative:\n";
				$complete_psm_set->print_six_nr_by_fdr($status_fdr);
			}
			print "\n" if $v;
		}
	}
	if ( -e $csv_in ) {
		my @csv = $tools->get_files( $csv_in, $table_file_ending );
		my $nr_csv = @csv;
		print "  no table input found!\n" if $nr_csv == 0;
		foreach my $csv ( sort @csv ) {
			my ( $psm_set, $roc_file ) =
			  get_psm_set_from_table( $csv, $initial_ds_cutoff,
				$decoy_id_prefix, $tar_dec_ratio, $max_psm_fdr, $tools,
				$table_separator, $out_base, $db );
			push @roc_files, $roc_file;

			$complete_psm_set->add_psm( $psm_set->get_psm(), $csv,
				$psm_set->get_runs() );

			if ($print_cumulative) {
				print "  cumulative:\n";
				$complete_psm_set->print_six_nr_by_fdr($status_fdr);
			}
			print "\n" if $v;
		}
	}

	# filter out PSM whose protein ids were not found in the protein
	# fasta database
	if ($filter_ids) {
		print "  filtering out protein ids that were not found in"
		  . " the target decoy protein database\n"
		  if $v;
		my $rh_ids = $bin_prot->get_ids();
		$complete_psm_set->filter_by_protein_id($rh_ids);
	}

	# print an overview of identifications for a given FDR
	print "  cumulative input data:\n" if $v;
	$complete_psm_set->print_six_nr_by_fdr($status_fdr) if $v;

	return ( $complete_psm_set, \@roc_files );
}

# Title	    :  get_psm_set_from_pepxml()
# Usage     :
# Function	:
# Returns   :  a PSM set object
# Args      :
sub get_psm_set_from_pepxml {
	my ( $pepxml, $initial_pps_cutoff, $decoy_id_prefix, $tar_dec_ratio,
		$max_psm_fdr, $tools, $out_base, $db )
	  = @_;

	my $pepxml_file = basename($pepxml);

	# extracted from pepxml
	# 0 spectrum
	# 1 peptide
	# 2 main id (protein)
	# 3 modification info
	# 4 PeptideProphet probability (PPs)
	my $ra_ra_psm = [];

	my $status_prefix = "parsing $pepxml_file... ";
	if ($use_xml_parser) {
		my $Mayu_handler = MayuPepXMLHandler->new();

		# create parser and set the handler
		my $parser = XML::Parser::PerlSAX->new(
			Handler      => $Mayu_handler,
			ErrorContext => 2,
		);
		my $status_prefix = "parsing $pepxml_file... ";
		$Mayu_handler->set( $initial_pps_cutoff, $pepxml, $parser, $v,
			$s, $status_prefix );
		my %parser_args = ( Source => { SystemId => $pepxml } );
		my $result = $parser->parse(%parser_args);
		$ra_ra_psm = $Mayu_handler->results();
	}
	else {
		my $Mayu_parser =
		  MayuPepXMLParser->new( $initial_pps_cutoff, $pepxml, $v, $s,
			$status_prefix );
		$ra_ra_psm = $Mayu_parser->parse();
	}
	print "\r  $pepxml_file parsed                        \n" if $s;
	print "  $pepxml_file parsed\n"                           if $v && !$s;

	# print the PSM to a csv file
	if ($print_input_output) {
		my $csv_file = $pepxml;
		$csv_file =~ s/\.xml$/_mayu.csv/g;
		open( O, ">$csv_file" ) or warn $!;
		foreach my $ra_psm (@$ra_ra_psm) {
			print O join( $table_separator, 
				( $ra_psm->[0], $ra_psm->[1], $ra_psm->[2], 
				get_string_from_mod( $ra_psm->[3] ), $ra_psm->[4] ) ) 
				. "\n";
		}
		close(O);
	}

	# calculate the PSM FDR based on the target decoy approach,
	# use the PeptideProphet score as discriminant score
	my $psm_set = PSMSet->new( $v, $s, $ra_ra_psm, $pepxml );
	my $roc_file = '';
	
	( $psm_set, $roc_file ) = preprocess_psm_set(
		$psm_set,       $min_pep_length, $decoy_id_prefix,
		$tar_dec_ratio, $status_fdr,     $max_psm_fdr,
		$pepxml_file,   $out_base,       $tools, $db
	);

	return ( $psm_set, $roc_file );
}

# Title	    :  get_psm_set_from_table()
# Usage     :
# Function	:
# Returns   :  a PSM set object
# Args      :
sub get_psm_set_from_table {
	my ( $csv, $initial_ds_cutoff, $decoy_id_prefix, $tar_dec_ratio,
		$max_psm_fdr, $tools, $table_separator, $out_base, $db )
	  = @_;

	my $status_fdr = 0.01;

	my $csv_file = basename($csv);

	# extracted from table
	# 0 spectrum
	# 1 peptide
	# 2 main id (protein)
	# 3 modification info
	# 4 PeptideProphet probability (PPs)
	my $ra_ra_psm =
	  parse_csv( $csv, $initial_ds_cutoff, $table_separator, $csv_header );
	print "\r  $csv_file parsed                        \n" if $s;
	print "  $csv_file parsed\n"                           if $v && !$s;
	
	# print the PSM to a csv file
	if ($print_input_output) {
		my $csv_file = $csv;
		$csv_file =~ s/\.csv$/_mayu.csv/g;
		open( O, ">$csv_file" ) or warn $!;
		foreach my $ra_psm (@$ra_ra_psm) {
			print O join( $table_separator, 
				( $ra_psm->[0], $ra_psm->[1], $ra_psm->[2], 
				get_string_from_mod( $ra_psm->[3] ), $ra_psm->[4] ) ) 
				. "\n";
		}
		close(O);
	}
	
	my $psm_set = PSMSet->new( $v, $s, $ra_ra_psm, $csv );
	
	my $roc_file = '';

	( $psm_set, $roc_file ) = preprocess_psm_set(
		$psm_set,       $min_pep_length, $decoy_id_prefix,
		$tar_dec_ratio, $status_fdr,     $max_psm_fdr,
		$csv_file,      $out_base,       $tools, $db
	);

	return ( $psm_set, $roc_file );
}

# Title	    :  parse_csv()
# Usage     :
# Function	:  parses a table file with PSM info.
#              the file format has to be a mascot output file
#              or in the following way:
#              optional header
#              spectrum, peptide, protein, modification, discr score
#              ...
#
#              $header turns the header on or off
#              , is defined by $table_separator
#              modification: pos1=mass1:pos2=mass2
# Returns   :  an array of references to arrays
#              0 spectrum
#              1 peptide
#              2 main id (protein)
#              3 modification info
#              4 discriminant score (ds)
# Args      :
sub parse_csv {
	my ( $csv, $initial_ds_cutoff, $table_separator, $header ) = @_;

	my $ra_ra_psm = [];

	if ( -e $csv ) {

		# check for mascot csv file
		my $mascot_csv = MascotCSV->new( 0, 0 );

		if ( $mascot_csv->is_mascot_csv_file($csv) ) {
			print "  $csv is mascot file!\n";
			$ra_ra_psm = $mascot_csv->get_mayu_psm($csv);
		}
		else {
			$ra_ra_psm = parse_Mayu_csv( $csv, $header );
		}
	}
	else {
		print "  $csv does not exist!\n";
		return $ra_ra_psm;
	}

	return $ra_ra_psm;
}

# Title	    :  parse_Mayu_csv()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub parse_Mayu_csv {
	my ( $csv, $header ) = @_;

	my $ra_ra_psm = [];

	my $spectrum_col     = 0;
	my $peptide_col      = 1;
	my $protein_col      = 2;
	my $modification_col = 3;
	my $ds_col           = 4;
	my $td_col           = 5;
	my $mfdr_col         = 6;

	my $header_was_removed = 0;
	
	open( CSV, $csv ) or warn $!;
	while (<CSV>) {
		my $line = $_;
		chomp($line);
		
		# check the format
		my @col = split( /$table_separator/, $line );
		my $nr = @col;
		if ( $nr == 7 ) {
			print "  input is " . basename( $0 ) . " output file. "
				. "Using (1 - mFDR score) as discriminant score!\n";
		}
		
		# remove the header 
		if ( $header && !$header_was_removed ) {
			$header_was_removed = 1;
			next;
		}
		else {
			my @col = split( /$table_separator/, $line );
			my $nr = @col;
			if ( defined( $col[0] ) ) {
				if ( $col[0] =~ /^scan$/ && !$header_was_removed ) {
					$header_was_removed = 1;
					next;
				}
				else {
					if ( $nr == 5 ) {
						push @$ra_ra_psm,
						  [
							$col[$spectrum_col],
							$col[$peptide_col],
							$col[$protein_col],
							get_mod_from_string( $col[$modification_col] ),
							$col[$ds_col]
						  ];
					}
					# use the mFDR score
					elsif ( $nr >= 7 ) {
						push @$ra_ra_psm,
						  [
							$col[$spectrum_col],
							$col[$peptide_col],
							$col[$protein_col],
							get_mod_from_string( $col[$modification_col] ),
							1 - $col[$mfdr_col]
						  ];
					}
				}
			}
		}
	}
	close(CSV);

	return $ra_ra_psm;
}

# Title	    :  get_mod_from_string()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_mod_from_string {
	my ( $mod_string, $mod_sep, $mod_split ) = @_;

	$mod_sep   = ':' unless defined($mod_sep);
	$mod_split = '=' unless defined($mod_split);

	my $rh = {};

	if ( defined($mod_string) ) {
		my @mod = split( /$mod_sep/, $mod_string );
		foreach my $mod (@mod) {
			my @pos_mass = split( /$mod_split/, $mod );
			if ( @pos_mass == 2 ) {
				$rh->{ $pos_mass[0] } = $pos_mass[1];
			}
			else {
				print "  something is wrong with modification in the "
				  . "csv input: '$mod'\n";
			}
		}
	}

	return $rh;
}

# Title	    :  get_string_from_mod()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub get_string_from_mod {
	my ( $rh_mod, $mod_sep, $mod_split ) = @_;

	$mod_sep   = ':' unless defined($mod_sep);
	$mod_split = '=' unless defined($mod_split);
	
	my $s = '';
	
	if ( defined($rh_mod) ) {
		my @mod_pos = keys %$rh_mod;
	
		my @mod = ();
		foreach my $pos (@mod_pos) {
			push @mod, $pos . $mod_split . $rh_mod->{$pos};
		}
		$s = join( $mod_sep, @mod )
	}

	return $s;
}

# Title	    :  preprocess_psm_set()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub preprocess_psm_set {
	my (
		$psm_set,       $min_pep_length, $decoy_id_prefix,
		$tar_dec_ratio, $status_fdr,     $max_psm_fdr,
		$input_file,   $out_base,       $tools, $db
	  )
	  = @_;

	$psm_set->filter_by_peptide_length($min_pep_length);
	
	if ( $remove_ambiguous ) {
		$psm_set->filter_ambiguous_peptides($db);
	}

	$psm_set->add_target_decoy_fdr( $decoy_id_prefix, $tar_dec_ratio );

	# get the numbers for all psm
	my (
		$nr_all_tpsm, $nr_all_dpsm, $nr_all_tpep,
		$nr_all_dpep, $nr_all_tpr,  $nr_all_dpr
	  )
	  = $psm_set->get_six_nr_by_fdr(1);
	print "  non FDR filtered target: $nr_all_tpsm PSM,\t$nr_all_tpep "
	  . "peptides,\t$nr_all_tpr proteins\n"
	  if $v;
	print "  non FDR filtered decoy:  $nr_all_dpsm PSM,\t$nr_all_dpep "
	  . "peptides,\t$nr_all_dpr proteins\n"
	  if $v;

	# give overview of size of the parsed data
	my ( $nr_tpsm, $nr_dpsm, $nr_tpep, $nr_dpep, $nr_tpr, $nr_dpr ) =
	  $psm_set->get_six_nr_by_fdr($status_fdr);

	# filter the data set with FDR
	$psm_set->filter_by_psm_fdr($max_psm_fdr);

	# print an overview of identifications for a given FDR
	print "  target: $nr_tpsm PSM,\t$nr_tpep peptides,\t$nr_tpr"
	  . " proteins at $status_fdr FDR\n"
	  if $v;
	print "  decoy:  $nr_dpsm PSM,\t$nr_dpep peptides,\t$nr_dpr"
	  . " proteins at $status_fdr FDR\n"
	  if $v;

	# print out the ROC data to a file
	$input_file =~ s/\.[^\.]+$//;
	my $roc_file =
	  $out_base . '_' . $input_file . '_ds_vs_mFDR_' . $version . '.txt';
	my $header    = 1;
	my $ra_ra_roc = $psm_set->get_roc_processed($header);
	$tools->print_data_in_table_style( $ra_ra_roc, $roc_file ) if $p_mfdr;

	return ( $psm_set, $roc_file );
}

# Title	    :  print_options()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub print_options {
	my ($out_base) = @_;
	my $st = time();
	print "\n  $0 $version\n";
	print "\n"                                                      if $v;
	print "  initial DS cutoff: $initial_ds_cutoff\n"               if $v;
	print "  status fdr: $status_fdr\n"                             if $v;
	print "\n  INPUT:\n"                                            if $v;
	print "  using xml parser: $use_xml_parser\n"                   if $v;
	print "  -A pepxml input: $pepxml_in\n"                         if $v;
	print "  -B table input: $csv_in\n"                             if $v;
	print "  -C search database: $db\n"                             if $v;
	print "  -D minimal peptide length: $min_pep_length\n"          if $v;
	print "\n  TARGET-DECOY OPTIONS:\n"                             if $v;
	print "  -E decoy id prefix: $decoy_id_prefix\n"                if $v;
	print "  -F target to decoy ratio: $tar_dec_ratio\n"            if $v;
	print "  -G maximal PSM FDR for analysis: $max_psm_fdr\n"       if $v;
	print "  -H PSM FDR steps: $psm_fdr_steps\n"                    if $v;
	print "\n  PROTEIN FDR CALCULATION:\n"                          if $v;
	print "  -I number of missed cleavages: $nr_missed_cleavages\n" if $v;
	print "  -J minimal peptide mass: $min_pep_mass\n"              if $v;
	print "  -K maximal peptide mass: $max_pep_mass\n"              if $v;
	print "  -L number of protein size bins: $nr_prot_size_bins\n"  if $v;
	print "  -N id set selection type: $id_set_selection\n"         if $v;

	if ( exists( $rh_id_selection_helper->{$id_set_selection} ) ) {
		print "     " . $rh_id_selection_helper->{$id_set_selection} . "\n"
		  if $v;
	}
	print "  -O cumulative runs steps: $nr_runs_steps\n"
	  if $v
	  && ( $id_set_selection == 3
		|| $id_set_selection == 4
		|| $id_set_selection == 5 );
	print "  don't correct identical sequences: $dont_corr_id_seq\n" if $v;
	print "  use equidist protein binning: $equidist_bins\n"         if $v;
	print "\n  POST ANALYSIS:\n"                                     if $v;
	print "  run R analysis: $run_r\n"                               if $v;
	print "\n  OUTPUT:\n"                                            if $v;
	print "  -P $ids_out_tag\n"           if $ids_out_tag !~ /^$/ && $v;
	print "  file name base: $out_base\n" if $v;
	print "  print mFDR file: $p_mfdr\n"  if $v;
	print "  print protein size bin file: $p_bin_prot\n"            if $v;
	print "  print protein feature file: $p_prot_feat\n"            if $v;
	print "  print input pepxml to csv: $print_input_output\n"      if $v;
	print "  print cumulative input to stdout: $print_cumulative\n" if $v;
	print "  print status: $s\n"                                    if $v;
	print "\n"                                                      if $v;
	print "\n"                                                      if $v;
}

# Title	    :  usage()
# Usage     :
# Function	:
# Returns   :
# Args      :
sub usage {
	print "\n\n\t" . basename($0) . " $version";
	print "
	------------\n";

	print '
	For more information and command lines for the example data set
	type "perl ' . basename($0) 
	. ' -manual" or consult the manual.
	
	GENERAL:
	required: (-A or -B) and -C
	there are default values for the rest of input options
	
	SETTINGS:
	-A :  file or directory with pepXML search result files (.xml)
	-B :  file or directory with Mayu or Mascot csv search result
	      files (.csv)
	-C :  target decoy search database in fasta format
	-E :  decoy id prefix (decoy ids in the faste database
	      have to start with this prefix, e.g. "rev_")
	-G :  defines upper limit of error analysis (maximal mFDR)
	-H :  defines the resolution of error analysis (mFDR steps)
	-P :  select a list of filtered PSM to be printed out, e.g.:
	      mFDR=0.01:t =    filter data with PSM FDR 0.01 and
	                       print out the corresponding target PSM
	      pepFDR=0.02:td = filter data with peptide identification
	                       FDR 0.02 and print out the corresponding 
	                       target and decoy PSM
	      protFDR=0.05:t = filter data with protein identification
	                       FDR 0.05 and print out the corresponding
	                       target PSM
	      
	OPTIONS:
	  -PmFDR         :  print out the mFDR table files separately 
	                    for each input file
	  -PbinProt      :  print out size binned protein file
	  -PprotFeat     :  print out a protein feature file
	                    (can be very large files)
	  -runR          :  generate plots from the table output
	                    (R package needs to be installed and R
	                    recognized as a command on the command line)
	  -verbose       :  print out more
	  -status        :  print out status
	  -help          :  show this help
	  -manual        :  detailed description of options
	';
	print "\n\n";
}

# Title	    :  manual()
# Usage     :
# Function	:  more detailed description of the options
# Returns   :
# Args      :
sub manual {
	print "\n\n\t" . basename($0) . " $version";
	print "
	------------\n";

	print '
	GENERAL:
	Mayu is a software package to determine protein
	identification false discovery rates (protFDR) and
	peptide identification false discovery rates (pepFDR)
	additionally to the peptide-spectrum match false discovery 
	rate (mFDR).
	
	Lukas Reiter
	Manfred Claassen

	SOFTWARE:
	Lukas Reiter - Hengartner Laboratory
	lukas.reiter@molbio.uzh.ch
	Institute of Molecular Biology
	Winterthurerstrasse 190
	University of Zurich - Irchel
	CH-8057 Zurich
	+++++++++++++++++++++++++++++++++++++++
	Located at:
	Institute for Molecular Systems Biology
	Aebersold Laboratory
	Wolfgang-Pauli-Str. 16
	ETH Hoenggerberg
	CH-8093 Zurich
	
	INSTALLATION:
	Perl needs to be installed on the system. E.g. use activestate perl 
	distribution (www.activestate.com). R statistical package needs to 
	be installed for graphical output (-runR). R command needs to be 
	recognized on the command line.
	
	PREREQUISITE:
	- data searched against a single concatenated target decoy database
	- search results as pepxml, mascot csv or Mayu table
	
	RECOMMENDED:
	- all data searched with similar search options
	- protein databases preferably low in redundant sequences
	
	EXAMPLES:
	1. standard analysis, main analysis table printed out
	   "perl Mayu.pl -B example.csv -C tardecdb.fa -v -s"
	2. plot graphics using the R statistical package
	   "perl Mayu.pl -B example.csv -C tardecdb.fa -v -s -runR"
	3. remove peptides smaller than 10 amino acids from target and decoy PSM
	   "perl Mayu.pl -B example.csv -C tardecdb.fa -D 10 -v -s"
	4. do calculations of error rates in 51 steps between 0 and 5% PSM FDR
	   "perl Mayu.pl -B example.csv -C tardecdb.fa -G 0.05 -H 51 -v -s"
	5. print out more result tables in separate files
	   "perl Mayu.pl -B example.csv -C tardecdb.fa -PmFDR -PbinProt -PprotFeat"
	6. start a long run on a unix system and log the standard output
	   "nohup perl Mayu.pl -B example.csv -C tardecdb.fa ... -v > log.txt &"
	7. print out target and decoy PSM, target PSM with a PSM FDR of 0.01
	   "perl Mayu.pl -B example.csv -C tardecdb.fa -P mFDR=0.01:td"
	8. print out target PSM whose protein ids correspond to a protFDR of 0.05
	   "perl Mayu.pl -B example.csv -C tardecdb.fa -P protFDR=0.05:t"
	9. use pepxml as input
	   "perl Mayu.pl -A sequest_pepxml.xml -C tardecdb.fa"
	10.pepxml as input, print out a .csv file of input for faster reanalysis
	   "perl Mayu.pl -A sequest_pepxml.xml -C tardecdb.fa -Pio"
	11.sort the LC-MS/MS runs by orthogonality (run is recognized by its
	   scan base) and perform the analysis on cumulative data sets in 11 steps
	   "perl Mayu.pl -B example.csv -C tardecdb.fa -N 5 -O 11"
	
	SETTINGS:
	required: (-A or -B) and -C
	-A :  pepXML file or directory with pepXML files (.xml)
	      from a search against the target decoy protein
	      database.
	-B :  Mayu or Mascot csv file or directory (.csv, comma separated)
	      with PSM from a search against the target decoy
	      protein database. Check the pdf manual or example
	      file for the Mayu format.
	-C :  target decoy search database in fasta format.
	      Protein decoy entries have to start with a 
	      prefix that marks them unambiguously (e.g. rev_).
	-D :  minimal peptide length for analysis.
	      Remove peptides that are shorter then this length
	-E :  decoy id prefix.
	      Decoy ids in the faste database have to start with this 
	      prefix (e.g. rev_).
	-F :  ratio of target to decoy entries (usually 1)
	-G :  maximal mFDR for Mayu analysis
	      Start analysis from zero mFDR and go with -H steps
	      to -G.
	-H :  mFDR steps (resolution of the analysis).
	      Make the protein identification FDR analysis for
	      -H mFDR cutoffs.
	-I :  number of missed cleavages used for database search.
	      Use the same number of missed cleavages for tryptic
	      peptides as for the database search.
	-J :  minimal peptide mass.
	      Usually determined by the settings of the mass 
	      spectrometer.
	-K :  maximal peptide mass.
	      Usually determined by the settings of the mass 
	      spectrometer.
	-L :  number of protein size bins.
	      Bin the proteins into -L protein bins and perform
	      the protein identification FDR estimation in these
	      bins separately.
	-M :  file name base.
	      Use this output file name base.
	-N :  id set selection.
	      This describes the way in which identifications are
	      selected from the input data, e.g. according to 
	      input file or orthogonality sorted runs additionally
	      to an mFDR cutoff.
	      For the options 3, 4 and 5 consider that roughly 
	      -H * -O * -L calculation entities have to be performed! 
	      Therefore start with decent step numbers.
	      0: all data and mFDR range (default)
	         take all data and apply a range of mFDR cutoffs
	      1: cumulative input files and mFDR range
	         take 1 to n cumulative data sets corresponding
	         to the n files (alphabetically sorted)
	      2: cumulative shuffled input files and mFDR range
	         take 1 to n cumulative data sets corresponding
	         to the n files (shuffled files)
	      3: cumulative runs and mFDR range
	         take 1 to -O cumulative data sets corresponding
	         to the -O*x runs (alphabetically sorted according
	         to files and then runs)
	      4: cumulative shuffled runs and mFDR range
	         take 1 to -O cumulative data sets corresponding
	         to the -O*x runs (shuffled runs)
	      5: cumulative runs (orthogonality sorted) and mFDR range
	         take 1 to -O cumulative data sets corresponding
	         to the -O*x runs. The runs are sorted according to
	         an orthogonality measure describing the uniquely
	         contributed peptides to the data set at 1% mFDR
	         normalized to the data amount in that run.
	         (orthogonality data stored in ..._ortho_run.txt)
	-O :  number of cumulative runs steps.
	      Used with option -N 3, 4 or 5
	-P :  select a filtered list of PSM to be printed out, e.g.:
	      mFDR=0.01:t =    filter data with PSM FDR 0.01 and
	                       print out the corresponding target PSM
	      pepFDR=0.02:td = filter data with peptide identification
	                       FDR 0.02 and print out the corresponding 
	                       target and decoy PSM
	      protFDR=0.05:t = filter data with protein identification
	                       FDR 0.05 and print out the corresponding
	                       target PSM
	      
	OPTIONS:
	For output file formats see also the pdf manual.
	  -equidist      :  equidist binning (according to estimated
	                    number of tryptic peptides) of proteins 
	                    instead of equicount (convergence is faster 
	                    with equicount binning).
	  -dcis          :  do not correct identical sequences in
	                    their protein size.
	                    If not selected (recommended) proteins
	                    that have identical sequences will be
	                    corrected in size to zero.
	  -xmlparser     :  use a proper xml parser for pepXML parsing.
	                    This is slower but much more robust 
	                    and the parser needs to be installed 
	                    first (libxml). It is strongly recommended 
	                    to use this option! Read the manual and
	                    ask your administrator for help.
	  -PmFDR         :  print out the mFDR table files for 
	                    each input file which relates the 
	                    discriminant score (e.g. PeptideProphet
	                    probability) to an FDR estimated with 
	                    the target decoy strategy.
	  -PbinProt      :  print out size binned protein file.
	                    protFDR and additional information on
	                    the protein size partition.
	  -PprotFeat     :  print out a protein feature file
	                    (can be very large files)
	                    For each protein identification a set of
	                    features that can be used for local FDR
	                    calculation. features:
	                    - decoy (1 if decoy, 0 if target)
	                    - NS (number of PSM)
	                    - NP (number of distinct peptides)
	                    - PAT (PSM alignment type)
	                    - PSL (corrected protein sequence length in aa)
	                    - acNTP (corrected number of tryptic peptides)
	  -Pio           :  print input pepxml to a csv table
	                    file in the same directory. This file
	                    contains the PSM and can be used as
	                    input in the next run (-B instead of -A) with 
	                    the advantage of faster parsing.
	  -cumul         :  print cumulative input to stdout.
	                    Number of PSM, peptides and proteins with each
	                    cumulative input file
	  -remamb        :  remove ambiguous peptides prior to analysis
	  -runR          :  generate plots from the table output
	                    (R package needs to be installed and R
	                    recognized as a command on the command line)
	  -verbose       :  print out more
	  -status        :  print out even more
	                    do not use this option if the standard
	                    out is written/logged to a file
	  -help          :  show help
	  -manual        :  show this manual
	';
	print "\n\n";
}

__END__


Lukas Reiter


