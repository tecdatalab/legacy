create_graphs <- function(filename,output_prefix)
{
	data = read.table(filename, strip.white=TRUE, fill=TRUE,
			col.names=c("rank","rmsd","pdb_euc","em_euc"));

	pdb_cor = cor(data$rmsd, data$pdb_euc);
	em_cor = cor(data$rmsd, data$em_euc);
	both_cor = cor(data$pdb_euc, data$em_euc);
	
	# RMSD VS Both 3DZD Euclidean distances
	output = paste(output_prefix, "-rmsdVSeuc.corr.png", sep="");
	png(filename = output, width = 700, height = 700);
	maintitle = paste(output_prefix, "RMSD/3DZD Euclidean Distance Correlation");
	subtitle = paste("PDB-PDB CC(blue):", sprintf("%.2f",pdb_cor), "EM-EM CC(red):", sprintf("%.2f",em_cor), "Both CC:", sprintf("%.2f",both_cor));
	plot(data$rmsd, data$pdb_euc, xlab="RMSD to Native",
		ylab="3DZD Euclidean distance", main=maintitle, sub=subtitle,
		pch=20, col="blue", ylim=c(0,max(max(data$pdb_euc),max(data$em_euc))));
	points(data$rmsd, data$em_euc, pch=15, col="red");
	dev.off();

	# Both Euclidean distances plotted against each other
	output = paste(output_prefix, "-eucVSeuc.corr.png", sep="");
	png(filename = output, width = 700, height = 700);
	maintitle = paste(output_prefix, "PDB/EM 3DZD Euclidean Distance Correlation");
	subtitle = paste("Correlation:", sprintf("%.2f",both_cor));
	plot(data$em_euc, data$pdb_euc, xlab="Euclidean distance from EM 3DZD",
		ylab="Euclidean distance from PDB 3DZD", main=maintitle, sub=subtitle,
		pch=20, col="blue");
	dev.off();

	# Print the correlation values for further processing
	cat(output_prefix,both_cor, pdb_cor, em_cor, em_cor - pdb_cor, "\n");
}
params <- commandArgs(trailingOnly=TRUE);
if(length(params) == 2) {
	for(param in params)
	{
		eval(parse(text=param));
	}
	create_graphs(input, output_prefix);
} else {
	cat("Usage: R --slave --vanilla --file=em_rmsd_3dzd_correlation_plots.R --args \"input='1A0R.eucstats'\" \"output_prefix='1A0R'\"\n");
}
