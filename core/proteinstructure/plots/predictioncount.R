# This function creates a line plot for each of the rows contained
# in the dataset from "start" to "end". Only rows that contain at least
# one value above "threshold" will be drawn
# Each of the rows represents the number of times a specific prediction
# is used in the GA process
predictioncountgraph <- function(data, start, end, threshold, title)
{
	# first create an empy plot but use 2 anchor points to set the graph boundaries
	anchorx = c(0, dim(data)[2]-3);
	anchory = c(0,  max(data[start:end,4:dim(data)[2]]));

	# generate the name of the image that will contain the plot
	imagename = paste(title, ".png", sep="");
	# initialize the png output
	png(filename = imagename, width = 500, height = 500);

	# type="n" avoids the actual display of the two points
	plot(anchorx,anchory, type="n", xlab="Generations", ylab="Count", main=title);

	# Positions 1,2 and 3 correspond to the Receptor, Ligand and Prediction Number
	totallabels = dim(data)[2]-3;
	# create values from 1 to N to identify each of the ranges
	labels = c(1:totallabels)

	for(rowindex in start:end)
	{
		# 4:lastcolumn because the first 3 don't represent counts
		values = data[rowindex,4:dim(data)[2]];

		#only lines that have a value over the threshold are drawn
		drawline = FALSE;
		for(countindex in 1:totallabels)
		{
			if(values[countindex] > threshold)
			{
				drawline = TRUE;
				break;
			}
		}
		if(drawline)
		{
			lines(labels, values, type="l");
		}
	}
	# close the png output
	dev.off();
}

# This function loads the prediction count data from "filename"
# and calls "predictioncountgraph" once for each pair wise combination (A-B, B-C, etc)
# Since the first two columns determine the receptor and ligand involved, and the file
# is sorted, then we can determine the ranges used for each pair wise combination
generateallgraphs <- function(filename, threshold)
{
	# load the data first
	countdata = read.csv(filename);

	# as expected we start from 1 and look for the first row
	# that has a different receptor/ligand combination
	currentstart = 1;
	currentend = 1;

	totalrows = dim(countdata)[1];
	
	# generate graphs until we reach the end of the file, i.e. when currentstart gets to the end
	while(currentstart <= totalrows)
	{
		# first look for the end of the current range
		# if we get to the end of the file then the range ends there
		# currentend+1 in both cases because we want currentend to point to the last element in this range
		while((currentend + 1) <= totalrows &&
			countdata[currentstart,1] == countdata[currentend + 1,1] &&
			countdata[currentstart,2] == countdata[currentend + 1,2])
		{
			# keep incrementing until they are different
			currentend = currentend + 1;
		}
		# at this point start and end point to the beginning and end of a range
		# generate the graph providing a title of the form "Chain1-Chain2"
		title = paste(countdata[currentstart,1], "-", countdata[currentstart,2]);
		predictioncountgraph(countdata, currentstart, currentend, threshold, title);

		#after generating the graph, look for the next range starting at currentend + 1;
		currentstart = currentend + 1;
		currentend = currentstart;
	}
}
#generateallgraphs("test2.csv",50);