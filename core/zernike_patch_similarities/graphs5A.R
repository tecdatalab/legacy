plot_aux <- function(datax, datay, maintitle, filename_base, filename_prefix,
										 xtitle, ytitle) {
	imagename = paste(filename_prefix, "-", filename_base, ".png", sep="");

	png(filename=imagename, width=500, height=500);

	plot(datax, datay, cex=1, pch=19, xlab=xtitle, ylab=ytitle, main=maintitle);
	dev.off();
}

all_plots <- function(data_filename, maintitle, filename_base) {
	headers = c("invfile1","invfile2","cc","patchintfraction1","totalintfraction1",
	  					"patchintfraction2", "totalintfraction2","centroiddistance",
							"spheredistance");
	data = read.table(data_filename, strip.white=TRUE, fill=TRUE,
										col.names=headers);

	data$avg_patch_interface_fraction <-
		(data$patchintfraction1 + data$patchintfraction2) / 2;
	data$avg_total_interface_fraction <-
		(data$totalintfraction1 + data$totalintfraction2) / 2;

	plot_aux(data$centroiddistance, data$cc, maintitle, filename_base,
					 "centroid-distances", "Distance between patch centroids", "3DZD CC");
	plot_aux(data$spheredistance, data$cc, maintitle, filename_base,
					 "sphere-center-distances", "Distance between sphere centers", "3DZD CC");
	plot_aux(data$avg_patch_interface_fraction, data$cc, maintitle, filename_base,
					 "interface-atoms-over-patch-total",
					 "Avg fraction of atoms in patch that belong to interface", "3DZD CC");
	plot_aux(data$avg_total_interface_fraction, data$cc, maintitle, filename_base,
					 "interface-atoms-over-total-interfaces",
					 "Avg fraction of the interface contained in patch", "3DZD CC");
}

all_plots("5Aint-1KXP-5/all.txt", "1KXP, 5A patches", "1KXP-5");
all_plots("5Aint-1KXP-10/all.txt", "1KXP, 10A patches", "1KXP-10");
all_plots("5Aint-1KXP-15/all.txt", "1KXP, 15A patches", "1KXP-15");
all_plots("5Aint-1KXP-20/all.txt", "1KXP, 20A patches", "1KXP-20");
all_plots("5Aint-1IB1-5/all.txt", "1IB1, 5A patches", "1IB1-5");
all_plots("5Aint-1IB1-10/all.txt", "1IB1, 10A patches", "1IB1-10");
all_plots("5Aint-1IB1-15/all.txt", "1IB1, 15A patches", "1IB1-15");
all_plots("5Aint-1IB1-20/all.txt", "1IB1, 20A patches", "1IB1-20");
all_plots("5Aint-1WQ1-5/all.txt", "1WQ1, 5A patches", "1WQ1-5");
all_plots("5Aint-1WQ1-10/all.txt", "1WQ1, 10A patches", "1WQ1-10");
all_plots("5Aint-1WQ1-15/all.txt", "1WQ1, 15A patches", "1WQ1-15");
all_plots("5Aint-1WQ1-20/all.txt", "1WQ1, 20A patches", "1WQ1-20");
