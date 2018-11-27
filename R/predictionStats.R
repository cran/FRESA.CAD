predictionStats_ordinal <-  function(predictions,plotname="",...)
{
    cat(plotname,"\n")
	dpoints <- nrow(predictions)
	tint <- qt(0.975,dpoints - 1)/sqrt(dpoints)
    ScoresOutcome <- predictions[,1]
    if (nchar(plotname) > 0) 
    {
      boxplot(predictions[,2] ~ ScoresOutcome,main = plotname,xlab = "Outcome" ,ylab = "Prediction",...)
    }
    ct <-  DescTools::KendallTauB(ScoresOutcome,as.integer(predictions[,2]+0.5),conf.level = 0.95)
    bias <- mean(predictions[,2] - ScoresOutcome)
    rstd <- sqrt(mean((predictions[,2] - ScoresOutcome)^2) - bias^2)
    Bias <- c(bias,bias - tint*rstd,bias + tint*rstd)
	theScores <- as.numeric(names(table(ScoresOutcome)))
	BMAE <- NULL;
	for (s in theScores)
	{
		BMAE <- rbind(BMAE,MAE95ci(predictions[ScoresOutcome==s,])); 
	}
	BMAE <- colMeans(BMAE);
	class95ci <- ClassMetric95ci(predictions);
    kp <- irr::kappa2(cbind(as.integer(predictions[,2] + 0.5),ScoresOutcome),"squared")
    zdis <- 2*kp$value/kp$statistic;
    Kapp <- c(kp$value, kp$value - zdis, kp$value + zdis)
    results <- list(Kendall = ct,
					Bias = Bias,
					BMAE = BMAE,
					Kapp = Kapp,
					class95ci = class95ci,
					KendallTauB = ct,
					Kappa.analysis = kp
					);
    return(results);
}

predictionStats_binary <-  function(predictions, plotname="", center=FALSE,...)
{
    cat(plotname,"\n")
	if (any(is.na(predictions[,2])))
	{	
		predictions <- predictions[!is.na(predictions[,2]),]
	}	
    if (center) predictions[,2] <- predictions[,2] - 0.5;
	if (min(predictions[,2]) >= 0)
	{
		predictions[,2] <- predictions[,2] - 0.5;
	}
	
    pm <- NULL;
	citest <- NULL;
    if (nchar(plotname) > 1)
    {
      pm <- plotModels.ROC(predictions,main = plotname,...);
      cis <- ci.auc(pm$roc.predictor)
    }
    else
    {
      pm <- pROC::roc(as.vector(predictions[,1]),predictions[,2]);
      cis <- ci.auc(pm);
      pm$predictionTable <- table(predictions[,2] < 0,1 - predictions[,1]);
		if (nrow(pm$predictionTable) == 1)
		{
			if ((rownames(pm$predictionTable) == "0") || (rownames(pm$predictionTable) == "FALSE"))
			{
				pm$predictionTable <- rbind(c(0,0),pm$predictionTable);
			}
			else
			{
				pm$predictionTable <- rbind(pm$predictionTable,c(0,0));
			}
			rownames(pm$predictionTable) <- c("0","1")
		}
    }
	class95ci <- ClassMetric95ci(cbind(predictions[,1],predictions[,2] >= 0));

#    print(pm$predictionTable)
    if (length(pm$predictionTable) > 2 )
    {
      ci <- epiR::epi.tests(pm$predictionTable);
      accc <- ci$elements$diag.acc;
      berror <- class95ci$berci;
      sensitivity <- ci$elements$sensitivity;
      specificity <- ci$elements$specificity;
    }
    else
    {
      accc <- c(0.5,0.5,0.5);
      berror <- c(0.5,0.5,0.5);
      sensitivity <- c(0.0,0.0,0.0);
      specificity <- c(0.0,0.0,0.0);
      names(accc) <- c("est","lower","upper")
      names(berror) <- c("est","lower","upper")
    }
    aucs <- cis[c(2,1,3)];
    results <- list(accc = accc,berror = berror,
					aucs = aucs,sensitivity = sensitivity,
					specificity = specificity,
					ROC.analysis = pm,
					CM.analysis = ci,
					ClassMetrics = class95ci
					);
    return(results);
}

predictionStats_regression <-  function(predictions, plotname="",...)
{
      cat(plotname,"\n")
	  dpoints <- nrow(predictions)
	  chsqup <- sqrt(dpoints/qchisq(0.025, df = dpoints))
	  chsqdown <- sqrt(dpoints/qchisq(0.975, df = dpoints))
	  tint <- qt(0.975,dpoints - 1)/sqrt(dpoints)

	  if (nchar(plotname) > 0) 
      {
         plot(predictions[,2] ~ predictions[,1],main = plotname,xlab = "Outcome", ylab = "Prediction",...)
      }
	  ct <- cor.test(predictions[,1],predictions[,2],method = "pearson");
	  corci <- c(ct$estimate,ct$conf.int);
	  bias <- mean(predictions[,2] - predictions[,1]);
	  rmse <- sqrt(mean((predictions[,2] - predictions[,1])^2));
	  rstd <- sqrt(rmse^2 - bias^2);
	  biasci <- c(bias,bias - tint*rstd,bias + tint*rstd);
	  RMSEci <- c(rmse,chsqdown*rmse,chsqup*rmse);
	  spearmanci <- sperman95ci(predictions);
	  MAEci <- MAE95ci(predictions);
	  results <- list(corci = corci, 
						biasci= biasci,
						RMSEci=RMSEci,
						spearmanci=spearmanci,
						MAEci=MAEci,
						pearson=ct
						);
	  return(results);
}
