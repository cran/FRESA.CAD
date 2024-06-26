predictionStats_survival <-  function(predictions, plotname="", atriskthr=1.0, ...)
{
	if(!sum(predictions[,4] == 0) == length(predictions[,4]))
	{
		if (!requireNamespace("survminer", quietly = TRUE)) {
			install.packages("survminer", dependencies = TRUE)
		}

		data1 <- data.frame(times=predictions[,1],status=predictions[,2],preds= 1.0/predictions[,4])
		CIRisk <- concordance95ci(datatest = data1)
		data2 <- data.frame(times=numeric(nrow(predictions)),status=predictions[,2],preds= predictions[,3])
		CILp <- concordance95ci(datatest = data2)

		onlycases <- predictions[,2]>0
		datatest <- data.frame(times=predictions[onlycases,1],preds= 1.0/predictions[onlycases,4])
		spearmanCI <- sperman95ci(datatest)
		
		if ( !inherits(atriskthr,"numeric") ) 
		{
			atriskthr <- median(predictions[,4]);
		}
		groups = predictions[,4] >= atriskthr
		labelsplot <- c("Other",sprintf("At Risk > %5.2f",atriskthr));
		paletteplot <- c("#00bbff", "#ff0000")
		newData <- data.frame(times=predictions[,1],status=predictions[,2],preds=predictions[,4],groups = groups);
		Curves <- survival::survfit(survival::Surv(times, status) ~ groups,newData)

        LogRankE <- EmpiricalSurvDiff(times=newData$times,
                  status=newData$status,
                  groups=newData$groups,
                  plots=plotname!="",main=plotname)
		if(plotname!="")
		{
			 graph <- survminer::ggsurvplot(Curves, data=newData, conf.int = TRUE, legend.labs=labelsplot,
                        palette = paletteplot,
                        ggtheme = ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20)),
                        title = plotname,
                        risk.table = TRUE,
                        tables.height = 0.2,
                        tables.theme = survminer::theme_cleantable())
      		print(graph)
		}
		
		LogRank <- survival::survdiff(Surv(predictions[,1], predictions[,2]) ~ predictions[,4] >= atriskthr )
		LogRank <- cbind(LogRank$chisq,  1 - pchisq(LogRank$chisq, length(LogRank$n) - 1));
		colnames(LogRank) <- cbind("chisq","pvalue");
		return( list(CIRisk = CIRisk,CILp=CILp,spearmanCI=spearmanCI, LogRank = LogRank, Curves = Curves,LogRankE = LogRankE,groups=groups) );
	}
	else{
		return( list(CIRisk = rep(0,nrow(predictions)), LogRank = rep(0,nrow(predictions)), Curves = NULL) );
	}
}

concordance95ci <- function(datatest,nss=1000)
{
  sz <- nrow(datatest)
  sesci <- c(0,0,0);
  isROCAUC <- sum(datatest$times)==0;
  if (sz>2)
  {
    ses <- numeric(nss);
    for (i in 1:nss)
    {
      bootsample <- datatest[sample(sz,sz,replace=TRUE),];
      if(isROCAUC)
      {
        ses[i] <- rcorr.cens(bootsample$preds, bootsample$status)[1]
      }
      else
      {
        ses[i] <- rcorr.cens(bootsample$preds, survival::Surv(bootsample$times,bootsample$status))[1]
      }
    }
    sesci <- quantile(ses, probs = c(0.5,0.025,0.975),na.rm = TRUE);
    sesci[1]<-mean(ses)
    names(sesci)<-c("median","lower","upper");
  }
  return (sesci);
}

predictionStats_ordinal <-  function(predictions,plotname="",...)
{
    cat(plotname,"\n")
	dpoints <- nrow(predictions)
	tint <- qt(0.975,dpoints - 1)/sqrt(dpoints)
    ScoresOutcome <- predictions[,1]
    if (nchar(plotname) > 0) 
    {
      boxplot(predictions[,2] ~ ScoresOutcome,main = plotname,...)
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

metric95ci <- function(metric,nss=1000,ssize=0)
{
	sz <- length(metric);
	if (ssize == 0)
	{
		ssize <- sz;
	}
	ssize <- min(sz,ssize);
	metricci <- c(0,0,0);
	if (sz>1)
	{
		meanMetric <- numeric(nss);
		for (i in 1:nss)
		{
		  bootsample <- metric[sample(sz,ssize,replace=TRUE)];
		  meanMetric[i] <- mean(bootsample);
		}
		metricci <- quantile(meanMetric, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	}
	return (metricci);
}

corcen95ci <- function(dataTable,nss=1000)
{
	sz <- nrow(dataTable);
	metricci <- c(0.5,0.5,0.5);
	if (sz>1)
	{
		meanMetric <- numeric(nss);
		for (i in 1:nss)
		{
		  bootsample <- dataTable[sample(sz,sz,replace=TRUE),];
		  meanMetric[i] <- rcorr.cens(bootsample[,2],bootsample[,1], outx=FALSE)[1];
		}
		metricci <- quantile(meanMetric, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	}
	return (metricci);
}

predictionStats_binary <-  function(predictions, plotname="", center=FALSE,...)
{
#    cat(plotname,"\n")
	if (!requireNamespace("epiR", quietly = TRUE)) {
			   install.packages("epiR", dependencies = TRUE)
	}
	cstat <- NULL;
	cstatCI <- c(0.5,0.5,0.5);
	medianTest <- NULL;
	parameters <- list(...);
	thrval <- 0.5;

	if (!is.null(parameters$thr))
	{
		thrval <- parameters$thr;
	}
	else
	{
		if ((min(predictions[,2]) < -0.01) | (max(predictions[,2]) > 1.01))
		{
			thrval <- 0.0;
		}
	}

	
	if (ncol(predictions)>2)
	{
		numberOfModels <- table(predictions[,2]);
		numberOfModels <- as.integer(names(numberOfModels));
		cstat <- numeric(length(numberOfModels))
		modSize <- 0;
		for (mi in numberOfModels)
		{
			  mtest <- predictions[,2] == mi;
			  modSize <- modSize+sum(mtest);
			  cstat[mi] <- rcorr.cens(predictions[mtest,3],predictions[mtest,1], outx=FALSE)[1];
		}
		modSize <- modSize/length(numberOfModels);
		boxstaTest <- try(boxplot(as.numeric(as.character(predictions[,3]))~rownames(predictions),plot = FALSE));
		if (!inherits(boxstaTest, "try-error"))
		{
			medianTest <- cbind(predictions[boxstaTest$names,1],boxstaTest$stats[3,]);
			rownames(medianTest) <- boxstaTest$names;
		}
		else
		{
			warning("boxplot test failed");
			medianTest <- cbind(predictions[,1],rep(0,nrow(predictions)));
			rownames(medianTest) <- rownames(predictions);
		}
		colnames(medianTest) <- c("Outcome","Median");
		smpCI <- as.integer(nrow(medianTest)/modSize+0.5);
#		cat("Avg size:", modSize,"Samp size:",smpCI,"\n");
		if (length(cstat) > 1) 
		{
			cstatCI <- metric95ci(cstat,ssize=smpCI);
		}
		else 
		{
			cstatCI <- c(cstat[1],0.0,1.0);
		}
		predictions	<- medianTest;
	}
	else
	{
		cstatCI <- corcen95ci(predictions,200 + 800*(nrow(predictions) < 1000) );
	}
	if (any(is.na(predictions[,2])))
	{	
		predictions <- predictions[!is.na(predictions[,2]),]
	}	
    if (center) predictions[,2] <- predictions[,2] - 0.5;
	
    pm <- NULL;
	citest <- NULL;
    if (nchar(plotname) > 1)
    {
      pm <- plotModels.ROC(predictions,main = plotname,...);
      cis <- ci.auc(pm$roc.predictor)
    }
    else
    {
      pm <- pROC::roc(as.vector(predictions[,1]),predictions[,2],quiet = TRUE);
      cis <- ci.auc(pm);
      pm$predictionTable <- table(predictions[,2] < thrval,1 - predictions[,1]);
		if (nrow(pm$predictionTable) == 1)
		{
			if ((rownames(pm$predictionTable) == "0") || (rownames(pm$predictionTable) == "TRUE"))
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
	class95ci <- ClassMetric95ci(cbind(predictions[,1],predictions[,2] >= thrval),1000 + 1000*(nrow(predictions) < 1000) );

#    print(pm$predictionTable)
    if (length(pm$predictionTable) > 2 )
    {
      ci <- epiR::epi.tests(pm$predictionTable);
      accc <- ci$detail[5,c(2:4)];
      berror <- class95ci$berci;
      sensitivity <- ci$detail[3,c(2:4)];
      specificity <- ci$detail[4,c(2:4)];
	  rownames(accc) <- NULL
	  rownames(sensitivity) <- NULL
	  rownames(specificity) <- NULL
    }
    else
    {
      accc <- c(0.5,0.5,0.5);
      cIndexSet <- c(0.5,0.5,0.5);
      cstatCI <- c(0.5,0.5,0.5);
      berror <- c(0.5,0.5,0.5);
      sensitivity <- c(0.0,0.0,0.0);
      specificity <- c(0.0,0.0,0.0);
      names(accc) <- c("est","lower","upper")
      names(berror) <- c("est","lower","upper")
    }
    aucs <- cis[c(2,1,3)];
	names(aucs) <- c("est","lower","upper")
    results <- list(accc = accc,berror = berror,
					aucs = aucs,sensitivity = sensitivity,
					specificity = specificity,
					ROC.analysis = pm,
					CM.analysis = ci,
					ClassMetrics = class95ci,
					cIndexSet = cstat,
					cIndexCI = cstatCI,
					medianTest = medianTest
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
