plot_aux <- function(datax, datay, maintitle, filename_base, filename_suffix,
										 xtitle, ytitle) {
	imagename = paste("plot-", filename_base, "-", filename_suffix, ".png", sep="");

	png(filename=imagename, width=500, height=500);

	plot(datax, datay, pch=".", xlab=xtitle, ylab=ytitle, main=maintitle);
	dev.off();
}

all_plots <- function(data_filename, maintitle, filename_base) {
	headers = c("invfile1","invfile2","cc","patchintfraction1","totalintfraction1",
	  					"patchintfraction2", "totalintfraction2","centroiddistance");
	data = read.table(data_filename, strip.white=TRUE, fill=TRUE,
										col.names=headers);

	data$avg_patch_interface_fraction <-
		(data$patchintfraction1 + data$patchintfraction2) / 2;
	data$avg_total_interface_fraction <-
		(data$totalintfraction1 + data$totalintfraction2) / 2;

	plot_aux(data$centroiddistance, data$cc, maintitle, filename_base,
					 "distances", "Distance between patch centroids", "3DZD CC");
	plot_aux(data$avg_patch_interface_fraction, data$cc, maintitle, filename_base,
					 "interface-atoms-over-patch-total",
					 "Avg fraction of atoms in patch that belong to interface", "3DZD CC");
	plot_aux(data$avg_total_interface_fraction, data$cc, maintitle, filename_base,
					 "interface-atoms-over-total-interfaces",
					 "Avg fraction of the interface contained in patch", "3DZD CC");
}

all_plots("1AY7-15/all.txt", "1AY7, 15A patches", "1AY7-15");
all_plots("1AY7-20/all.txt", "1AY7, 20A patches", "1AY7-20");
all_plots("1KXP-15/all.txt", "1KXP, 15A patches", "1KXP-15");
all_plots("1KXP-20/all.txt", "1KXP, 20A patches", "1KXP-20");
