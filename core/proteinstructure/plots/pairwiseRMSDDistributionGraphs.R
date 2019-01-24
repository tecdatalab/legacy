rmsd_distribution_graph <- function(pdb, description, data)
{
	# hardcode fnat to 0.1 to see if we have some hope
	fnatthreshold = 0.1;
	# threshold bounds
	lowerbound = 0;
	upperbound = 10;
	# assign default to total, in case there is a read error (like when we don't have info in the file)
	total = 0;
	#x-axis: each of the RMSD threshold values
	datax = c(lowerbound:upperbound);
	#y-axis: the number of predictions that have a RMSD lower than the associated threshold
	datay = array(0,upperbound - lowerbound + 1);
	if(is(data$irmsd[1],"numeric")) { # we don't do any formal check but at least see if the first value is indeed a number
		updateindex = 1;
		# column 3 represents RMSD
		total = length(data$irmsd);
		for(threshold in datax)
		{
			counter = 0;
			for(prediction in 1:total)
			{
				if(data$irmsd[prediction] < threshold && (!is.nan(data$fnat[prediction]))
					&& data$fnat[prediction] > fnatthreshold)
				{
					counter = counter + 1;
				}
			}
			datay[updateindex] = counter;
			updateindex = updateindex + 1;
		}
	}
	# to make the analysis easier see how many predictions have RMSD <4 (hardcoded to position 5)
	good = datay[5];

	print(good);

	# and now create the image file
	xtitle = "iRMSD threshold";
	ytitle = "Number of hits";
	maintitle = paste(pdb, description);
	subtitle = paste("Total predictions =", total, " Below 4A=", good);
	imagename = paste(maintitle, ".png", sep="");
	png(filename = imagename, width = 500, height = 500);

	# first create an empy plot but use 2 anchor points to set the graph boundaries
	anchorx = c(lowerbound, upperbound);
	anchory = c(0,  max(datay));

	# type="n" avoids the actual display of the two points
	plot(anchorx,anchory, type="n", xlab=xtitle, ylab=ytitle, main=maintitle, sub=subtitle);

	# and now draw the line on that plot
	lines(datax, datay, type="l");
	# close the png output
	dev.off();
}

pairwise_rmsd_graphs <- function(pdb, chains)
{
	# generate one graph per chain combination, go through each of them
	for(rec in 1:(length(chains)-1))
	{
		for(lig in (rec + 1):length(chains))
		{
			fileprefix = paste(chains[rec], "-", chains[lig], sep="");
			filename = paste(fileprefix, ".stats", sep="");
			# print(paste("Processing ", filename));
			data = read.table(filename, strip.white=TRUE, fill=TRUE, col.names=c("id","irmsd","rmsd","fnat"));
			rmsd_distribution_graph(pdb, fileprefix, data);
		}
	}
}
params <- commandArgs(trailingOnly=TRUE);
if(length(params) == 2) {
	for(param in params)
	{
		eval(parse(text=param));
	}
	pairwise_rmsd_graphs(pdb, chains);
} else {
	cat("Usage: R --slave --vanilla --file=pairwiseRMSDDistributionGraphs.R --args \"pdb='1A0R'\" chains=\"c('B','G','P')\"\n");
}
