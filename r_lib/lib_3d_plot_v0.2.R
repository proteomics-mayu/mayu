# 3d_plot_v0.2.R
#
# Author: Lukas Reiter
###############################################################################

# Function
# @title:    plot_curves_from_4_vectors_3_versions
# @param:    
# @usage:    
# @function: 
# @returns:  
plot_curves_from_4_vectors_3_versions <- function( 
	t.v_x=c(1), t.v_y=c(1), t.curve_index=c(1), t.sd_index=c(1),
	t.plot_y_min=NULL, t.plot_y_max=NULL, 
	t.plot_x_min=NULL, t.plot_x_max=NULL, 
	t.v_colors=c(1), t.v_pch=c(20), t.v_lty=c(1), t.v_legend=c(""), 
	t.title="", t.xlab="", t.ylab="", t.tck=NA, 
	cex.main=1,
	t.draw_arrows=1, t.arrow_code=3, t.arrow_angle=70, t.arrow_length=0.03, t.arrow_lwd=1,
	t.legend_border=10, t.legend_factor=1.5, t.legend_box_lty=0, t.legend_cex=1,
	t.intelligent_range=0,
	t.legend_position=NULL, t.legend_title=NULL, t.v_legend_pch=NULL,
	...
	) {
	
#	cat("plot_curves_from_4_vectors_3_versions()\n");
	
	# calculate the means and sds for each sd index
	t.v_unique_sd <- sort( as.vector( unique(t.sd_index) ) )
	
	# new vectors
	t.v_new_x <- c()
	t.v_mean <- c()
	t.v_sd <- c()
	t.v_new_curve_index <- c()
	
	# TODO one could also simply use the t.curve_index
	for ( t.i in 1:length(t.v_unique_sd) ) {
		t.curr_sd <- t.v_unique_sd[t.i]
		t.curr_sd_index <- which( t.sd_index == t.curr_sd )
		
		t.v_unique_x <- sort( as.vector( unique( t.v_x[t.curr_sd_index] ) ) )
		
		# unique values on the x-axes
		for ( t.j in 1:length(t.v_unique_x) ) {
			t.curr_x <- t.v_unique_x[t.j]
		
			t.x_index <- which( t.v_x[t.curr_sd_index] == t.curr_x )
			
			t.v_new_x <- c( t.v_new_x, t.curr_x )
			
			t.v_mean <- c( t.v_mean, mean( t.v_y[t.curr_sd_index][t.x_index], na.rm=T ) )
			
			t.sd <- sd( t.v_y[t.curr_sd_index][t.x_index], na.rm=T )
			if ( is.na(t.sd) ) {
				t.v_sd <- c( t.v_sd, 0 )
			} else {
				t.v_sd <- c( t.v_sd, t.sd )
			}
			
			t.v_new_curve_index <- c( t.v_new_curve_index, t.curr_sd )
		}
	}
	
	# make a plot of the means and the standard deviations
	plot_curves_from_3_vectors_3_versions(
		t.v_new_x,           # x values
		t.v_mean,            # y values
		t.v_sd,              # standard deviations in y direction
		t.v_new_curve_index, # curve index
		t.plot_y_min,
		t.plot_y_max,
		t.plot_x_min,
		t.plot_x_max, 
		t.v_colors,          # colors of the curves
		t.v_pch,             # point types of the curves
		t.v_lty,             # line styles of the curves
		t.v_legend,
		t.title,
		t.xlab, 
		t.ylab, 
		t.tck,
		cex.main,
		t.draw_arrows,
		t.arrow_code,
		t.arrow_angle, 
		t.arrow_length, 
		t.arrow_lwd,
		t.legend_border,
		t.legend_factor,
		t.legend_box_lty,
		t.legend_cex,
		t.intelligent_range,
		t.legend_position,
		t.legend_title,
		t.v_legend_pch,
		...
		)
}

# used by var_ana_v0.1.R
# Function
# @title:    plot_curves_from_3_vectors_3_versions
# @param:    
# @usage:    
# @function: 
# @returns:  
plot_curves_from_3_vectors_3_versions <- function( 
	t.v_x=c(1), t.v_y=c(1), t.v_sd=c(1), t.curve_index=c(1),
	t.plot_y_min=NULL, t.plot_y_max=NULL, 
	t.plot_x_min=NULL, t.plot_x_max=NULL, 
	t.v_colors=c(1), t.v_pch=c(20), t.v_lty=c(1), t.v_legend=c(""), 
	t.title="", t.xlab="", t.ylab="", t.tck=NA, 
	cex.main=1,
	t.draw_arrows=1, t.arrow_code=3, t.arrow_angle=70, t.arrow_length=0.03, t.arrow_lwd=1,
	t.legend_border=10, t.legend_factor=1.5, t.legend_box_lty=0, t.legend_cex=1,
	t.intelligent_range=0,
	t.legend_position=NULL, t.legend_title=NULL, t.v_legend_pch=NULL,
	...
	) {
	
#	cat("plot_curves_from_3_vectors_3_versions()\n");
	
	# input with #points values
	if ( length(t.v_y) != length(t.v_x) ) {
		t.v_y <- rep( 1, length(t.v_x) )
	}
	if ( length(t.v_sd) != length(t.v_x) ) {
		t.v_sd <- rep( 0, length(t.v_x) )
		t.draw_arrows <- 0
	}
	if ( length(t.curve_index) != length(t.v_x) ) {
		t.curve_index <- rep( 1, length(t.v_x) )
	}
	# input with #curves values
	t.nr_curves <- length( unique(t.curve_index) )
	if ( length(t.v_colors) != t.nr_curves ) {
		library(fields)
#		t.v_colors <- tim.colors( t.nr_curves )
		t.v_colors <- rainbow(t.nr_curves, start=0, end=0.6)
	}
	if ( length(t.v_pch) != t.nr_curves ) {
		# 20 as default
		t.v_pch <- rep(20,t.nr_curves)
	}
	if ( length(t.v_lty) != t.nr_curves ) {
		# 1 as default
		t.v_lty <- rep(1,t.nr_curves)
	}
	if ( length(t.v_legend) != t.nr_curves ) {
		t.v_legend <- paste( "curve", 1:t.nr_curves )
	}
	
	plot_curves_from_3_vectors(
		t.v_x, t.v_y, t.v_sd, t.curve_index, 
		t.plot_y_min, t.plot_y_max, t.plot_x_min, t.plot_x_max, 
		t.v_colors, t.v_pch, t.v_lty,
		t.title, t.xlab, t.ylab, t.tck, 
		cex.main,
		t.draw_arrows, t.arrow_code, t.arrow_angle, t.arrow_length, t.arrow_lwd,
		t.intelligent_range,
		...
		)
	plot_curves_from_3_vectors_inc_legend_and_legend_only(
		t.v_x, t.v_y, t.v_sd, t.curve_index, 
		t.plot_y_min, t.plot_y_max, t.plot_x_min, t.plot_x_max, 
		t.v_colors, t.v_pch, t.v_lty, t.v_legend,
		t.title, t.xlab, t.ylab, t.tck, 
		cex.main,
		t.draw_arrows, t.arrow_code, t.arrow_angle, t.arrow_length, t.arrow_lwd,
		t.legend_border, t.legend_factor, t.legend_box_lty, t.legend_cex,
		t.intelligent_range,
		t.legend_position, t.legend_title, t.v_legend_pch,
		...
		)
}

# used by var_ana_v0.1.R
# Function
# @title:    plot_curves_from_3_vectors
# @param:    
# @usage:    
# @function: 
# @returns:  
plot_curves_from_3_vectors <- function(
	t.v_x=c(1), t.v_y=c(1), t.v_sd=c(1), t.curve_index=c(1),
	t.plot_y_min=NULL, t.plot_y_max=NULL, 
	t.plot_x_min=NULL, t.plot_x_max=NULL, 
	t.v_colors=c(1), t.v_pch=c(20), t.v_lty=c(1),
	t.title="", t.xlab="", t.ylab="", t.tck=NA, 
	cex.main=1,
	t.draw_arrows=1, t.arrow_code=3, t.arrow_angle=70, t.arrow_length=0.03, t.arrow_lwd=1,
	t.intelligent_range=0,
	...
	) {
	
#	cat("plot_curves_from_3_vectors()\n");
	
	# input with #points values
	if ( length(t.v_y) != length(t.v_x) ) {
		t.v_y <- rep( 1, length(t.v_x) )
	}
	if ( length(t.v_sd) != length(t.v_x) ) {
		t.v_sd <- rep( 0, length(t.v_x) )
		t.draw_arrows <- 0
	}
	if ( length(t.curve_index) != length(t.v_x) ) {
		t.curve_index <- rep( 1, length(t.v_x) )
	}
	# input with #curves values
	t.nr_curves <- length( as.vector( unique(t.curve_index) ) )
	if ( length(t.v_colors) != t.nr_curves ) {
		library(fields)
#		t.v_colors <- tim.colors( t.nr_curves )
		t.v_colors <- rainbow(t.nr_curves, start=0, end=0.6)
	}
	if ( length(t.v_pch) != t.nr_curves ) {
		# 20 as default
		t.v_pch <- rep(20,t.nr_curves)
	}
	if ( length(t.v_lty) != t.nr_curves ) {
		# 1 as default
		t.v_lty <- rep(1,t.nr_curves)
	}
	
	# adjust the drawing range
	t.x_min <- min(t.v_x, na.rm=T)
	t.x_max <- max(t.v_x, na.rm=T)
	t.y_min <- min(t.v_y, na.rm=T)
	t.y_max <- max(t.v_y, na.rm=T)
	if ( !is.null(t.plot_y_min) == 1 ) {
		t.y_min <- t.plot_y_min
	}
	if ( !is.null(t.plot_y_max) == 1 ) {
		t.y_max <- t.plot_y_max
	}
	if ( !is.null(t.plot_x_min) == 1 ) {
		t.x_min <- t.plot_x_min
	}
	if ( !is.null(t.plot_x_max) == 1 ) {
		t.x_max <- t.plot_x_max
	}
	if ( t.draw_arrows == 1 & t.intelligent_range ) {
		t.tmp_y_min <- min(t.v_y - t.v_sd, na.rm=T)
		t.tmp_y_max <- max(t.v_y + t.v_sd, na.rm=T)
		if ( t.tmp_y_min < t.y_min ) {
			t.y_min <- t.tmp_y_min
		}
		if ( t.tmp_y_max > t.y_max ) {
			t.y_max <- t.tmp_y_max
		}
	}
	
	# TODO for debugging purposes
#	print(t.v_x)
#	print(t.v_y)
#	print(t.v_sd)
#	print(t.curve_index)
#	print(t.plot_y_min)
#	print(t.plot_y_max)
#	print(t.plot_x_min)
#	print(t.plot_x_max)
#	print(t.v_colors)
#	print(t.v_pch)
#	print(t.v_lty)
#	print(t.title)
#	print(t.xlab)
#	print(t.ylab)
#	print(t.tck)
#	t.draw_arrows,
#	t.arrow_code,
#	t.arrow_angle, 
#	t.arrow_length, 
#	t.arrow_lwd,
#	t.legend_border,
#	t.legend_factor,
#	t.legend_box_lty,
#	t.legend_cex,
#	t.intelligent_range,
	
	# set up the plot with writings and size
	plot( 1 ~ 1, 
		main=t.title,
		cex.main=cex.main,
		xlim=c( t.x_min, t.x_max ), ylim=c( t.y_min, t.y_max ), 
		xlab=t.xlab, ylab=t.ylab, tck=t.tck, type="n",
		...
		)
	
	# draw all the curves into the plot
	t.v_indices <- as.vector( unique( t.curve_index ) )
	t.nr_curves <- length( t.v_indices )
#	print("t.v_indices")
#	print(t.v_indices)
	for ( t.curr_index_of_index in 1:t.nr_curves ) {
	
#		print( paste( "index", t.curr_index_of_index ) )
		# old, changed 31.3.2008
		# t.curr_curve_index <- which( t.curve_index == t.curr_index_of_index )

		# new
		t.curr_index_value <- t.v_indices[t.curr_index_of_index]
		t.curr_curve_index <- which( t.curve_index == t.curr_index_value )
#		print("current indices")
#		print(t.curr_curve_index)
		
#		print("values")
#		print(t.v_x[t.curr_curve_index])
#		print(t.v_y[t.curr_curve_index])
		
		points( t.v_y[t.curr_curve_index] ~ t.v_x[t.curr_curve_index], 
			col=t.v_colors[t.curr_index_of_index],
			pch=t.v_pch[t.curr_index_of_index]
			)
		lines( t.v_y[t.curr_curve_index] ~ t.v_x[t.curr_curve_index], 
			col=t.v_colors[t.curr_index_of_index],
			lty=t.v_lty[t.curr_index_of_index]
			)
		if ( t.draw_arrows == 1 ) {
			arrows( t.v_x[t.curr_curve_index], t.v_y[t.curr_curve_index] - t.v_sd[t.curr_curve_index], 
				t.v_x[t.curr_curve_index], t.v_y[t.curr_curve_index] + t.v_sd[t.curr_curve_index],
				code = t.arrow_code, col=t.v_colors[t.curr_index_of_index], angle = t.arrow_angle,
				length = t.arrow_length, lwd=t.arrow_lwd )
		}
	}

}

# used by var_ana_v0.1.R
# Function
# @title:    plot_curves_from_3_vectors_inc_legend
# @param:    
# @usage:    
# @function: 
# @returns:  
plot_curves_from_3_vectors_inc_legend <- function(	
	t.v_x=c(1), t.v_y=c(1), t.v_sd=c(1), t.curve_index=c(1),
	t.plot_y_min=NULL, t.plot_y_max=NULL, 
	t.plot_x_min=NULL, t.plot_x_max=NULL, 
	t.v_colors=c(1), t.v_pch=c(20), t.v_lty=c(1), t.v_legend=c(""), 
	t.title="", t.xlab="", t.ylab="", t.tck=NA, 
	cex.main=1,
	t.draw_arrows=1, t.arrow_code=3, t.arrow_angle=70, t.arrow_length=0.03, t.arrow_lwd=1,
	t.legend_border=10, t.legend_factor=1.5, t.legend_box_lty=0, t.legend_cex=1,
	t.intelligent_range=0,
	t.legend_position=NULL, t.legend_title=NULL, t.v_legend_pch=NULL,
	...
	) {
	
#	cat("plot_curves_from_3_vectors_inc_legend()\n");
	
	plot_curves_from_3_vectors(
		t.v_x, t.v_y, t.v_sd, t.curve_index, 
		t.plot_y_min, t.plot_y_max, t.plot_x_min, t.plot_x_max, 
		t.v_colors, t.v_pch, t.v_lty,
		t.title, t.xlab, t.ylab, t.tck, 
		cex.main,
		t.arrow_code, t.arrow_angle, t.arrow_length, t.arrow_lwd,
		t.intelligent_range,
		...
		)
	
	t.nr_curves <- length( unique(t.curve_index) )
	t.usr <- par("usr")
	add_legend( 
		t.nr_curves, t.v_legend, t.v_colors, t.v_pch, 
		t.usr, t.legend_box_lty, t.legend_border, t.legend_factor, t.legend_cex,
		t.legend_position, t.legend_title, t.v_legend_pch
		)
}

# used by var_ana_v0.1.R
# Function
# @title:    plot_curves_from_3_vectors_inc_legend_and_legend_only
# @param:    
# @usage:    
# @function: 
# @returns:  
plot_curves_from_3_vectors_inc_legend_and_legend_only <- function(	
	t.v_x=c(1), t.v_y=c(1), t.v_sd=c(1), t.curve_index=c(1),
	t.plot_y_min=NULL, t.plot_y_max=NULL, 
	t.plot_x_min=NULL, t.plot_x_max=NULL, 
	t.v_colors=c(1), t.v_pch=c(20), t.v_lty=c(1), t.v_legend=c(""), 
	t.title="", t.xlab="", t.ylab="", t.tck=NA, 
	cex.main=1,
	t.draw_arrows=1, t.arrow_code=3, t.arrow_angle=70, t.arrow_length=0.03, t.arrow_lwd=1,
	t.legend_border=10, t.legend_factor=1.5, t.legend_box_lty=0, t.legend_cex=1,
	t.intelligent_range=0,
	t.legend_position=NULL, t.legend_title=NULL, t.v_legend_pch=NULL,
	...
	) {
	
#	cat("plot_curves_from_3_vectors_inc_legend_and_legend_only()\n");
	
	plot_curves_from_3_vectors(
		t.v_x, t.v_y, t.v_sd, t.curve_index, 
		t.plot_y_min, t.plot_y_max, t.plot_x_min, t.plot_x_max, 
		t.v_colors, t.v_pch, t.v_lty,
		t.title, t.xlab, t.ylab, t.tck, 
		cex.main,
		t.arrow_code, t.arrow_angle, t.arrow_length, t.arrow_lwd,
		t.intelligent_range,
		...
		)
	
	t.nr_curves <- length( unique(t.curve_index) )
	t.usr <- par("usr")
	add_legend( 
		t.nr_curves, t.v_legend, t.v_colors, t.v_pch, 
		t.usr, t.legend_box_lty, t.legend_border, t.legend_factor, t.legend_cex,
		t.legend_position, t.legend_title, t.v_legend_pch
		)
	plot_legend(
		t.nr_curves, t.v_legend, t.v_colors, t.v_pch, 
		t.legend_box_lty, t.legend_border, t.legend_factor, t.legend_cex,
		t.title, t.xlab, t.ylab, t.tck, 
		cex.main,
		t.legend_position, t.legend_title, t.v_legend_pch
		)
}

# used by var_ana_v0.1.R
# Function
# @title:    add_legend
# @param:    
# @usage:    
# @function: 
# @returns:  
add_legend <- function(
	t.nr_curves=1, t.v_legend=1, t.v_colors=1, t.v_pch=20, 
	t.usr, t.legend_box_lty=0, t.legend_border=10, t.legend_factor=1.5, t.legend_cex=1,
	t.legend_position=NULL, t.legend_title=NULL, t.v_legend_pch=NULL
	) {
	
	# some better defaults than 1
	if ( length(t.v_legend) != t.nr_curves ) {
		t.v_legend <- paste( "curve", 1:t.nr_curves )
	}
	if ( length(t.v_colors) != t.nr_curves ) {
		library(fields)
#		t.v_colors <- tim.colors( t.nr_curves )
		t.v_colors <- rainbow(t.nr_curves, start=0, end=0.6)
	}
	
	# pch
	t.v_this_pch <- rep(20,t.nr_curves)
	if ( is.null( t.v_legend_pch ) ) {
		if ( length(t.v_pch) == t.nr_curves ) {
			t.v_this_pch <- t.v_pch
		}
	} else {
		if ( length(t.v_legend_pch) == t.nr_curves ) {
			t.v_this_pch <- t.v_legend_pch
		}
	}
	
	# added 8.8.2008
	if ( !is.null(t.legend_position) == 1 ) {
		legend(
				t.legend_position, 
				legend = t.v_legend,
				col = t.v_colors, 
				title = t.legend_title,
				box.lty = t.legend_box_lty,
				pch = t.v_this_pch,
				cex = t.legend_cex
		)
	} else { 
		
		t.span_x <- abs( t.usr[2] - t.usr[1] )
		t.span_y <- abs( t.usr[3] - t.usr[4] )
		
		# for each legend entry
		for ( t.i in 1:t.nr_curves ) {
			
			# x, y, txt,
			t.v_x <- t.usr[1] + t.span_x / t.legend_border
			t.v_y <- t.usr[4] - t.span_y / t.legend_border - 
					( t.i - 1 ) * ( ( t.span_y - t.span_y / t.legend_border ) / 
						( t.nr_curves * t.legend_factor ) )
			
			legend( t.v_x, t.v_y, 
					t.v_legend[t.i],
					col = t.v_colors[t.i], 
					box.lty = t.legend_box_lty,
					pch=t.v_this_pch[t.i],
					cex = t.legend_cex
			)
		}
		
	}
	
}

# used by var_ana_v0.1.R
# Function
# @title:    plot_legend
# @param:    
# @usage:    
# @function: 
# @returns:  
plot_legend <- function(
	t.nr_curves=1, t.v_legend=1, t.v_colors=1, t.v_pch=20, 
	t.legend_box_lty=0, t.legend_border=10, t.legend_factor=1.5, t.legend_cex=1,
	t.title="", t.xlab="", t.ylab="", t.tck=NA, 
	cex.main=1,
	t.legend_position=NULL, t.legend_title=NULL, t.v_legend_pch=NULL
	) {
	
	# some better defaults than 1
	if ( length(t.v_legend) != t.nr_curves ) {
		# 20 as default
		t.v_legend <- 1:t.nr_curves
	}
	if ( length(t.v_colors) != t.nr_curves ) {
		library(fields)
#		t.v_colors <- tim.colors( t.nr_curves )
		t.v_colors <- rainbow(t.nr_curves, start=0, end=0.6)
	}
	if ( length(t.v_pch) != t.nr_curves ) {
		# 1 as default
		t.v_pch <- rep(20,t.nr_curves)
	}
	
	plot(1,type="n",
		main=t.title,
		cex.main=cex.main,
		xlab=t.xlab, ylab=t.ylab, tck=t.tck
		)
	t.usr <- par("usr")
	add_legend(
		t.nr_curves, t.v_legend, t.v_colors, t.v_pch, 
		t.usr, t.legend_box_lty, t.legend_border, t.legend_factor, t.legend_cex,
		t.legend_position, t.legend_title, t.v_legend_pch
		)
}


