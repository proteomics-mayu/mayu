# main.R
# 
# R template to generate graphics from the main Mayu output file.
# generated graphics:
# PSM:
# - # PSM vs. mFDR
# - # TP PSM vs. mFDR
# - # FP PSM vs. mFDR
# peptide:
# - # peptides vs. pepFDR
# - # TP peptides vs. pepFDR
# - # FP peptides vs. pepFDR
# protein:
# - # proteins vs. protFDR
# - # TP proteins vs. protFDR
# - # FP proteins vs. protFDR
#
# TODO 
# - sensitivity
# - comparison plot no sh <-> all 
#   use min of all curves
#
# Author: Lukas Reiter
###############################################################################

############################################
# INIT                                     #
############################################

# clean working environment
rm(list = ls())

############################################
# SUBS                                     #
############################################

#-------------------------------------------
# path of the r library
#-------------------------------------------
t.lib_path <- "<REPLACE_WITH_LIB_PATH>"

t.v_lib <- c( 
		"lib_3d_plot_v0.2.R" # 3d plotting
)
for ( t.i in 1:length(t.v_lib)) {
	source( paste( t.lib_path, t.v_lib[t.i], sep="" ) )
}

############################################
# INPUT/OPTIONS                            #
############################################

# data identifier
t.date <- date()
t.date <- gsub( ' ', '_', t.date, perl=TRUE )
t.date <- gsub( ':', '.', t.date, perl=TRUE )

# choice of printing to stdout or to pdf
plot_to_pdf <- 1

# describes the underlying filtered data entity
t.col_nr_runs <- "nr_runs"
t.col_nr_files <- "nr_files"
t.col_mFDR <- "mFDR"

# PSM information
t.col_target_PSM <- "target_PSM"
t.col_decoy_PSM <- "decoy_PSM"
t.col_FP_PSM <- "FP_PSM"
t.col_TP_PSM <- "TP_PSM"

# peptide information
t.col_target_pepID <- "target_pepID"
t.col_decoy_pepID <- "decoy_pepID"
t.col_FP_pepID <- "FP_pepID"
t.col_TP_pepID <- "TP_pepID"
t.col_pepFDR <- "pepFDR"

# protein information
t.col_target_protID <- "target_protID"
t.col_decoy_protID <- "decoy_protID"
t.col_FP_protID <- "FP_protID"
t.col_TP_protID <- "TP_protID"
t.col_protFDR <- "protFDR"

# single PSM protID information
t.col_target_protIDs <- "target_protIDs"
t.col_decoy_protIDs <- "decoy_protIDs"
t.col_FP_protIDs <- "FP_protIDs"
t.col_TP_protIDs <- "TP_protIDs"
t.col_protFDRs <- "protFDRs"

# non single PSM protID information
t.col_target_protIDns <- "target_protIDns"
t.col_decoy_protIDns <- "decoy_protIDns"
t.col_FP_protIDns <- "FP_protIDns"
t.col_TP_protIDns <- "TP_protIDns"
t.col_protFDRns <- "protFDRns"

# plotting options
t.legend_cex <- 0.6

#-------------------------------------------
# select input .csv main Mayu file
#-------------------------------------------
t.main_file <- "<REPLACE_WITH_INPUT_FILE>"
#t.main_file <- file.choose()

#-------------------------------------------
# path for output
#-------------------------------------------
t.out_path <- "<REPLACE_WITH_OUT_PATH>"

#-------------------------------------------
# output pdf file with r plots
#-------------------------------------------
pdf_plots_file_name <- "<REPLACE_WITH_OUTFILE>"
#pdf_plots_file_name <- paste( dirname( t.main_file ), "/", 
#	t.date, "_PANDORA_main.pdf", sep="")

############################################
# MAIN                                     #
############################################

# open the pdf file after you choose the first file
# such the pdf file will be written to this directory
if (plot_to_pdf) {
	pdf( paste( t.out_path, "/", pdf_plots_file_name, sep="" ) )
}

#-------------------------------------------
# plots mFDR vs various other properties
#-------------------------------------------
# TODO plot sensitivity

t.d <- read.table(
	file=t.main_file,
	sep="\t",
	header=T
	)

if ( t.d[,1] > 0 ) {
	
	#-------------------------------------------
	# make the 3d plots
	# - color code = nr_runs
	#-------------------------------------------
	
	# a vector of the data entities
	t.unique_nr_runs <- as.vector( unique( t.d[,t.col_nr_runs] ) )
	
	t.v_title <- c(
			"Peptide-Spectrum Matches",
			"Proteins",
			"Proteins",
			"Proteins",
			"Proteins",			
			"Peptides",
			"Peptides",
			"Peptides",
			"Single Hits \n(Single PSM Protein Identifications)",
			"Single Hits \n(Single PSM Protein Identifications)",
			"Single Hits \n(Single PSM Protein Identifications)",
			"All But Single Hits \n(Without Single PSM Protein Identifications)",
			"All But Single Hits \n(Without Single PSM Protein Identifications)",
			"All But Single Hits \n(Without Single PSM Protein Identifications)"
	)
	t.v_xlab <- c(
			"peptide-spectrum match FDR",
			"peptide-spectrum match FDR",
			"peptide-spectrum match FDR",
			"protein identification FDR",
			"protein identification FDR",
			"peptide-spectrum match FDR",
			"peptide-spectrum match FDR",
			"peptide identification FDR",
			"peptide-spectrum match FDR",
			"peptide-spectrum match FDR",
			"single hit FDR",
			"peptide-spectrum match FDR",
			"peptide-spectrum match FDR",
			"all but single hit FDR"
		)
	t.v_ylab <- c(
			"true positive peptide-spectrum matches",
			"protein identification FDR",
			"true positive protein identifications",
			"true positive protein identifications",
			"target protein identifications",
			"peptide identification FDR",
			"true positive peptide identifications",
			"true positive peptide identifications",
			"single hit FDR",
			"true positive single hits",
			"true positive single hits",
			"all but single hit FDR",
			"true positive all but single hits",
			"true positive all but single hits"
		)
	t.x <- c(
			t.col_mFDR,
			t.col_mFDR,
			t.col_mFDR,
			t.col_protFDR,
			t.col_protFDR,
			t.col_mFDR,
			t.col_mFDR,
			t.col_pepFDR,
			t.col_mFDR,
			t.col_mFDR,
			t.col_protFDRs,
			t.col_mFDR,
			t.col_mFDR,
			t.col_protFDRns
	)
	t.y <- c( 
			t.col_TP_PSM,
			t.col_protFDR,
			t.col_TP_protID,
			t.col_TP_protID,
			t.col_target_protID,
			t.col_pepFDR, 
			t.col_TP_pepID,
			t.col_TP_pepID,
			t.col_protFDRs, 
			t.col_TP_protIDs,
			t.col_TP_protIDs,
			t.col_protFDRns, 
			t.col_TP_protIDns,
			t.col_TP_protIDns
	)
	t.v_legend_title <- rep( "data set size", length(t.y) )
	
	for ( t.i in 1:length(t.y) ) {
		t.v_x <- c()
		t.v_y <- c()
		t.curve_index <- c()
		t.v_legend <- c()
		for ( t.j in 1:length(t.unique_nr_runs) ) {
			t.curr_nr_runs <- t.unique_nr_runs[t.j]
			t.v_legend <- c( t.v_legend, paste(t.curr_nr_runs, "LC-MS/MS runs") )
			t.sub <- subset( t.d, t.d[,t.col_nr_runs] == t.curr_nr_runs )
			t.v_x <- c( t.v_x, t.sub[,t.x[t.i]] )
			t.v_y <- c( t.v_y, t.sub[,t.y[t.i]] )
			t.curve_index <- c( t.curve_index, rep( t.curr_nr_runs, length( t.sub[,t.x[t.i]] ) ) )
		}
	
		plot_curves_from_3_vectors_inc_legend(
			t.v_x=t.v_x,
			t.v_y=t.v_y,
			t.curve_index=t.curve_index,
			t.plot_y_min=0, 
			t.plot_x_min=0,
			t.title=t.v_title[t.i],
			t.v_legend=t.v_legend,
			t.xlab=t.v_xlab[t.i],
			t.ylab=t.v_ylab[t.i],
			t.legend_title = t.v_legend_title[t.i],
			t.legend_cex=t.legend_cex,
			t.legend_box_lty=1,
			t.legend_position = "bottomright"
		)
	}
	
	#-------------------------------------------
	# compare mFDR with mFDR and only multi hits
	#-------------------------------------------
	t.max_nr_runs <- max( as.vector( unique( t.d[,t.col_nr_runs] ) ) )
	t.sub <- subset( t.d, t.d[,t.col_nr_runs] == t.max_nr_runs )
	
	t.x <- c( 
		t.col_protFDR,
		t.col_protFDRns
		)
	t.y <- c( 
		t.col_TP_protID,
		t.col_TP_protIDns
		)
	t.v_l <- c(
		"mFDR",
		"mFDR and only multi hits"
		)
	
	t.v_x <- c()
	t.v_y <- c()
	t.curve_index <- c()
	t.v_legend <- c()
	for ( t.i in 1:length(t.x) ) {
		t.v_x <- c( t.v_x, t.sub[,t.x[t.i]] )
		t.v_y <- c( t.v_y, t.sub[,t.y[t.i]] )
		t.curve_index <- c( t.curve_index, rep( t.i, length( t.sub[,t.x[t.i]] ) ) )
		t.v_legend <- c( t.v_legend, paste( t.max_nr_runs, "runs,", t.v_l[t.i] ) )
	}
	plot_curves_from_3_vectors_inc_legend(
		t.v_x=t.v_x,
		t.v_y=t.v_y,
		t.curve_index=t.curve_index,
		t.v_lty=c(0,0),
		t.plot_y_min=0, 
		t.plot_x_min=0,
		t.title="Effect of Removing Single Hits\non Estimated True Positive Protein Identifications",
		t.v_legend=t.v_legend,
		t.xlab="protFDR",
		t.ylab="TP_protID",
		t.legend_cex=t.legend_cex,
		t.legend_box_lty=1,
		t.legend_position="bottomright",
		t.legend_title = "with resp. without single hits"
	)
}

# do not keep the pdf file blocked
if (plot_to_pdf) {
	dev.off()
}



