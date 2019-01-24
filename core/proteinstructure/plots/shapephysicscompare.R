rmsd_distribution_graph <- function(pdb, description, data, physicsdata)
{
	# threshold set to 4A
	threshold = 4;
	# rank bounds
	lowerbound = 1;
	upperbound = 2000;
	# assign default to total, in case there is a read error (like when we don't have info in the file)
	total = 0;
	#x-axis: each of the RMSD threshold values
	datax = seq(lowerbound,upperbound);
	#y-axis: the counts that represent the number of hits up to that rank
	datay = array(0,length(datax));
	physicsdatay = array(0,length(datax));
	counter = 0;
	physicscounter = 0;
	if(is(data$rmsd[1],"numeric")) { # we don't do any formal check but at least see if the first value is indeed a number
		updateindex = 1;
		for(rank in datax)
		{
			if(data$rmsd[rank] < threshold)
			{
				counter = counter + 1;
			}
			if(physicsdata$rmsd[rank] < threshold)
			{
				physicscounter = physicscounter + 1;
			}
			datay[updateindex] = counter;
			physicsdatay[updateindex] = physicscounter;
			updateindex = updateindex + 1;
		}
	}
	#cat(datay[1000]);
	#cat("\n");
	#cat(physicsdatay[1000]);
	# and now create the image file
	xtitle = "Rank";
	ytitle = "Number of hits";
	maintitle = paste(pdb, description);
	imagename = paste(maintitle, ".png", sep="");
	png(filename = imagename, width = 600, height = 500);

	# first create an empy plot but use 2 anchor points to set the graph boundaries
	anchorx = c(lowerbound, upperbound);
	anchory = c(0,  max(max(datay),max(physicsdatay)));

	# type="n" avoids the actual display of the two points
	plot(anchorx,anchory, type="n", xlab=xtitle, ylab=ytitle, cex.lab=1.65, cex.axis=1.65);#, main=maintitle);

	# and now draw the line on that plot
	#par(pch=24,col="red"); #triangle
	lines(datax, datay, type="l", lty="solid");
	#par(pch=21,col="blue"); #circle
	#lines(datax, physicsdatay, type="l", lty="4F");
	lines(datax, physicsdatay, type="l", lty="solid", lwd=8);

	# add line legend
	#legend("bottomright", c("Shape", "Physics"), lty=1:2, cex=1.65);
	legend("bottomright", c("Shape", "Physics"), lwd=c(1,8), cex=1.65);

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
			physicsfilename = paste(fileprefix, ".pstats", sep="");
			print(paste("Processing ", filename));
			data = read.table(filename, strip.white=TRUE, fill=TRUE, col.names=c("id","irmsd","rmsd","fnat"));
			physicsdata = read.table(physicsfilename, strip.white=TRUE, fill=TRUE, col.names=c("id","rmsd","irmsd","fnat", "score", "notused"));
			# sort it according to the score
			sortedphysicsdata = physicsdata[sort.list(physicsdata$score), ];
			rmsd_distribution_graph(pdb, fileprefix, data, sortedphysicsdata);
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
	cat("Usage: R --slave --vanilla --file=shapephysicscompare.R --args \"pdb='1A0R'\" chains=\"c('B','G','P')\"\n");
}
