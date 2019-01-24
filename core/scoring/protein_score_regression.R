# Assume the data contain the 12 terms plus a RMSD column
# rmsd
TrainWeights <-	function(data, type, rmsd_threshold)
{
	if(type == "linear") {
		mod <- lm(RMSD~vdw+vdw_attr+vdw_rep+elec+elec_sr_attr+elec_lr_attr+elec_sr_rep+elec_lr_rep+hbp_ss+solv+sasa+acp,data);
	} else { #logistic
		data$hit <- (data$RMSD <= rmsd_threshold)*1;
		mod <- glm(as.factor(hit)~vdw+vdw_attr+vdw_rep+elec+elec_sr_attr+elec_lr_attr+elec_sr_rep+elec_lr_rep+hbp_ss+solv+sasa+acp,family="binomial"(link="logit"),data);
	}
	smod <- step(mod,trace=F);

	cf <- coefficients(smod);

	feature_names <- c('(Intercept)', 'vdw', 'vdw_attr', 'vdw_rep',
		'elec', 'elec_sr_attr', 'elec_lr_attr', 'elec_sr_rep', 'elec_lr_rep',
		'hbp_ss', 'solv', 'sasa', 'acp');
	weights <- rep(0,13);

	# it is possible that some of the original features are not used as coefficients
	# update those that are actually returned
	weights[feature_names %in% names(cf)] <- cf;

	weights;

#	cf <- cf[names(cf) %in% vnames];
#	kk[vnames[2:13] %in% names(cf)]	<-	-cf*1000;

#	print(kk);

}

# Program entry point

params <- commandArgs(trailingOnly=TRUE);
if(length(params) < 3) {
	cat("Usage: R --slave --vanilla --file=protein_score_regression.R --args \"training='datafile'\" \"test='datafile'\" \"type=<'linear'/'logistic'>\" [\"hit_cutoff=4.0\"] [\"logit_cutoff=0\"]\n");
} else {
	# by default hits are 4A structures...
	hit_cutoff = 4.0;
	# and classification is set at >0 or <0 for logistic regression
	logit_cutoff = 0;

	for(param in params) {
		eval(parse(text=param));
	}

	training_set <- read.table(training, sep=",", header=T);
	test_set <- read.table(test, sep=",", header=T);
	
	weights <- TrainWeights(training_set, type, hit_cutoff);
	
	cat("Intercept/Weights: ", weights, "\n");

	# Validation 
	score_terms_only <- cbind(1, test_set$vdw, test_set$vdw_attr, test_set$vdw_rep, test_set$elec,
				test_set$elec_sr_attr, test_set$elec_lr_attr, test_set$elec_sr_rep,
				test_set$elec_lr_rep, test_set$hbp_ss, test_set$solv, test_set$sasa,
				test_set$acp);
	
	predicted_values <- score_terms_only %*% weights;


	if(type == "linear") {
		predicted_class <- (predicted_values <= hit_cutoff) * 1;
	} else {
		predicted_class <- (predicted_values > logit_cutoff) * 1;
	}
	actual_class <- (test_set$RMSD <= hit_cutoff) * 1;

	positives <- sum(actual_class);
	negatives <- length(actual_class) - positives;

	combined <- predicted_class + actual_class;

	true_positives <- sum((combined == 2) * 1);
	true_negatives <- sum((combined == 0) * 1);
	false_positives <- sum(((combined == 1) * 1) * predicted_class);
	false_negatives <- sum(((combined == 1) * 1) * actual_class);

	cat("TN: ", true_negatives, " FN: ", false_negatives, "\n");

	cat("Negative Predictive Value [TN / (TN + FN)]: ", true_negatives / (true_negatives + false_negatives), "\n");
	cat("False Discovery Rate [FP / (FP + TP)]: ", false_positives / (false_positives + true_positives), "\n");
	cat("Accuracy: ", (true_positives + true_negatives) / (positives + negatives), "\n");

	cat("TP: ", true_positives, "/", positives, "(", true_positives / positives, ")",
		" FP: ", false_positives, "/", negatives, "(", false_positives / negatives, ")\n");

}
