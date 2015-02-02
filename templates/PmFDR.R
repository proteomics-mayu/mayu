# PmFDR.R
# 
# R template to generate graphics for the input files. run the program with
# the option -PmFDR for generating this graphics.
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

# column names of the mFDR files
t.col_ds <- "PPs"
t.col_mFDR <- "mFDR"

#plotting options
t.legend_cex <- 0.6

#-------------------------------------------
# select input .txt roc files
#-------------------------------------------
t.v_ds_vs_mFDR_files <- "<REPLACE_WITH_ROC_FILES>"
#t.v_ds_vs_mFDR_files <- choose.files())

#-------------------------------------------
# path for output
#-------------------------------------------
t.out_path <- "<REPLACE_WITH_OUT_PATH>"

#-------------------------------------------
# output pdf file with r plots
#-------------------------------------------
pdf_plots_file_name <- "<REPLACE_WITH_OUTFILE>"
#pdf_plots_file_name <- paste( dirname( t.v_ds_vs_mFDR_files[1] ), "/", 
#	t.date, "_PmFDR.pdf", sep="" )

############################################
# MAIN                                     #
############################################

# open the pdf file after you choose the first file
# such the pdf file will be written to this directory
if (plot_to_pdf) {
	pdf( paste( t.out_path, "/", pdf_plots_file_name, sep="" ) )
}

#-------------------------------------------
# mFDR vs. discriminant score cutoff for
# all input files
#-------------------------------------------
if ( length(t.v_ds_vs_mFDR_files) > 0 ) {
	t.v_ds <- c()
	t.v_mFDR <- c()
	t.v_roc_file_group <- c()
	for ( t.i in 1:length(t.v_ds_vs_mFDR_files)) {
		t.curr_roc_file <- t.v_ds_vs_mFDR_files[t.i]
		t.df_roc <- read.table(
				file=t.curr_roc_file,
				header=T
		)
		t.v_ds <- c( t.v_ds, t.df_roc[,t.col_ds] )
		t.v_mFDR <- c( t.v_mFDR, t.df_roc[,t.col_mFDR] )
		t.v_roc_file_group <- c( t.v_roc_file_group, 
				rep( basename( t.curr_roc_file ), length(t.df_roc[,1]) ) )
	}
	
	plot_curves_from_3_vectors_inc_legend(
			t.v_x=t.v_mFDR,
			t.v_y=t.v_ds,
			t.curve_index=t.v_roc_file_group,
			t.v_pch=rep("",length(t.v_ds_vs_mFDR_files)),
			t.plot_y_min=0, 
			t.plot_x_min=0,
			t.title="Search Result(s) - Mayu Input File(s)",
			t.v_legend=basename(t.v_ds_vs_mFDR_files),
			t.xlab="peptide-spectrum match FDR",
			t.ylab="discriminant score cutoff",
			t.legend_cex=t.legend_cex,
			t.legend_box_lty=1,
			t.legend_position="topright",
			t.legend_title="input search result files",
			t.v_legend_pch=rep("-",length(t.v_ds_vs_mFDR_files))
	)
	
	plot_curves_from_3_vectors_inc_legend(
			t.v_x=t.v_mFDR,
			t.v_y=t.v_ds,
			t.curve_index=t.v_roc_file_group,
			t.v_pch=rep("",length(t.v_ds_vs_mFDR_files)),
			t.plot_y_min=0, 
			t.plot_x_min=0,
			t.plot_x_max=0.05,
			t.title="zoom: discriminant score in dependence\nof the peptide-spectrum match FDR",
			t.v_legend=basename(t.v_ds_vs_mFDR_files),
			t.xlab="mFDR",
			t.ylab="discriminant score cutoff",
			t.legend_cex=t.legend_cex,
			t.legend_box_lty=1,
			t.legend_position="topright",
			t.legend_title="input search result files",
			t.v_legend_pch=rep("-",length(t.v_ds_vs_mFDR_files))
	)
}

# do not keep the pdf file blocked
if (plot_to_pdf) {
	dev.off()
}



