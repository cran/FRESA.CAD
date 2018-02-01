getSignature <-
function(data,varlist=NULL,Outcome=NULL,target=c("All","Control","Case"),CVFolds=3,repeats=9,distanceFunction=signatureDistance,...)
{
#It will use a back elimination procedure to return the list of features that maximize the CV ROC AUC based on a nearest centroid scheme.
#varlist a list of candidate features
#data the data frame where a binary outcome is defined 
#CVfolds the number of folds to be used
#the number of repetitions of each CV folds
#distanceFunction is funtion that returns the distance betweeen the signature template and each row of a dataframe
#... parameters to the distanceFuntion

trimvalue = 0.05; #trim 5% of the tails
theProbs <- c(0.025,0.16, 0.5,0.84, 0.975);
sortVarDistance <- function(data,varlist,Outcome,distanceFunction=signatureDistance,...)
{
	outDistance <- function (x,template,method) {md <- abs(cor(x,template,method=method,use="pairwise.complete.obs")); return (md)}

	indx <- NULL;
	starsig <- colMeans(data[,varlist],na.rm=TRUE);
	plot(starsig);
	matt <- data[,varlist];
	colnames(matt) <- c(1:length(varlist));
	outc <- as.vector(data[,Outcome]);
	names(varlist) <- c(1:length(varlist));
  	metric <- as.numeric(apply(matt,2,outDistance,outc,method="spearman"));
	for (i in 2:(length(varlist)-1))
	{
		metric[is.na(metric)] <- 0;
		maxc <- which.max(metric);
		indx <- c(indx,colnames(matt)[maxc]);
		outc <- as.vector(matt[,maxc]);
		matt <- matt[,-maxc];
		metric <- as.numeric(apply(matt,2,outDistance,outc,method="spearman"));
	}
	indx <- c(indx,colnames(matt));
	print(indx);
	varlist <- varlist[indx];
	starsig <- colMeans(data[,varlist],na.rm=TRUE);
	plot(starsig);
	return (varlist);
}

CVDistance <- function (casesample,controlsample,CVFolds,totRepeats,target)
{
	controlDistance <- numeric();
	caseDistance <- numeric();
	outcomes <-  numeric();
	IDs <- character();
	stdg <- apply(rbind(casesample,controlsample),2,sd,na.rm = TRUE);
	casesamplesize <- nrow(casesample);
	controlsamplesize <- nrow(controlsample);
	for (i in 1:totRepeats) 
	{
		j <- 1 + ((i-1) %% CVFolds)
		if ( j == 1)
		{
			casefolds <- cvTools::cvFolds(casesamplesize, CVFolds,1,  "random");
			controlfolds <- cvTools::cvFolds(controlsamplesize, CVFolds,1,  "random");
		}
		CaseTrainSet <- casesample[casefolds$subsets[casefolds$which != j,],];
		CaseBlindSet <- casesample[casefolds$subsets[casefolds$which == j,],];
		ControlTrainSet <- controlsample[controlfolds$subsets[controlfolds$which != j,],];
		ControlBlindSet <- controlsample[controlfolds$subsets[controlfolds$which == j,],];

		testData <- rbind(CaseBlindSet,ControlBlindSet);
		controlTemplate <- apply(ControlTrainSet,2,quantile,probs = theProbs,na.rm = TRUE);
		caseTemplate <- apply(CaseTrainSet,2,quantile,probs = theProbs,na.rm = TRUE);
		
		controlDistance <- append(controlDistance,distanceFunction(controlTemplate,testData,...));
		caseDistance <- append(caseDistance,distanceFunction(caseTemplate,testData,...));

		outcomes <- append(outcomes,c(rep(1,nrow(CaseBlindSet)),rep(0,nrow(ControlBlindSet))));
		IDs <- append(IDs,rownames(testData));
	}

	if (target=="Control") 
	{
		distance = controlDistance;
		cAUC <- pROC::roc(outcomes,distance,auc=TRUE,direction="<")$auc
		cont <- distance[outcomes==0];
		case <- distance[outcomes==1];
		zdis <- (mean(case,trim=trimvalue,na.rm=TRUE)-mean(cont,trim=trimvalue,na.rm=TRUE))/(1.0e-10+sqrt((var(case,na.rm=TRUE)+var(cont,na.rm=TRUE))/2));
	}
	else
	{
		if (target=="Case") 
		{
			distance = caseDistance;
			cAUC <- pROC::roc(outcomes,distance,auc=TRUE,direction=">")$auc
			cont <- distance[outcomes==0];
			case <- distance[outcomes==1];
			zdis <- (mean(cont,trim=trimvalue,na.rm=TRUE)-mean(case,trim=trimvalue,na.rm=TRUE))/(1.0e-10+sqrt((var(case,na.rm=TRUE)+var(cont,na.rm=TRUE))/2));
		}
		else
		{
			distance = controlDistance-caseDistance;
			sen <- sum((distance>=0)*outcomes)/sum(outcomes);
			spe <- sum((distance<0)*(outcomes==0))/sum(outcomes==0);
			cAUC <- (0.90*(sen+spe)/2+0.1*pROC::roc(outcomes,distance,auc=TRUE,direction="<")$auc);
#			cat("Sen",sen,"Spe",spe,"AUC",pROC::roc(outcomes,distance,auc=TRUE,direction="<")$auc,"(",sum(outcomes),":",sum(outcomes==0),")\n")
			cont <- distance[outcomes==0];
			case <- distance[outcomes==1];
			zdis <- (mean(case,trim=trimvalue,na.rm=TRUE)-mean(cont,trim=trimvalue,na.rm=TRUE))/(1.0e-10+sqrt((var(case,na.rm=TRUE)+var(cont,na.rm=TRUE))/2));
		}
	}

	result <- list(	controlDistance=controlDistance,
					caseDistance=caseDistance,
					outcomes=outcomes,
					IDs=IDs,
					cAUC=cAUC,
					zdis=zdis
					);
	return (result);
}
	
	if (is.null(varlist)) varlist=colnames(data)[-1];
	target <- match.arg(target);
	if (is.null(Outcome))
	{
		Outcome=varlist[1];
		varlist <- varlist[-1];
	}
#	varlist <- sortVarDistance(data,varlist,Outcome,signatureDistance,...);

	casesample = subset(data[,c(Outcome,varlist)],get(Outcome)  == 1);
	controlsample = subset(data[,c(Outcome,varlist)],get(Outcome) == 0);
	casesamplesize <- nrow(casesample);
	controlsamplesize <- nrow(controlsample);
	featuresize <- length(varlist);
	stsize <- min(5,featuresize-1);
	lastFeatureSize <- featuresize;
	totRepeats <- CVFolds*repeats;
	AUCevolution <- numeric();
	Zevolution <- numeric();
	featureSizeEvolution <- numeric();
	ES=0;
	maxES <- 0;
	maxAUC <- 0;
	mmasES <- 0;
	oldES <-0;
	keepfeature <- character();
	bestkeep <- character();
	CVOutput <- NULL;
		
	featuresize <- stsize+1;
	lastFeatureSize <- featuresize;
	snames <- varlist[1:featuresize];
	bestkeep <- snames;
	minkeep <- bestkeep;
	cvdis <- CVDistance(casesample[,snames],controlsample[,snames],CVFolds,totRepeats,target);
	maxAUC <- cvdis$cAUC;
	maxES <- (0.01*cvdis$zdis + 0.99*cvdis$cAUC)*(featuresize-1)/featuresize;

	if (lastFeatureSize<length(varlist))
	{
		fowardkeep <- bestkeep;
		keepfeature <- character();
		featureID=lastFeatureSize+1;
		rdeltaES <- 0.1;
		oldES <- maxES;
		sizef <- length(bestkeep);
		stop <- FALSE
		while (!stop)
		{	
			snames <- c(fowardkeep,keepfeature,varlist[featureID]);
			cvdis <- CVDistance(casesample[,snames],controlsample[,snames],CVFolds,totRepeats,target);
			sz <- length(snames);
			ES <- (0.01*cvdis$zdis + 0.99*cvdis$cAUC)*(sz-1)/sz;
			if (ES >= 0.99*maxES) 
			{
				keepfeature <- append(keepfeature,varlist[featureID]);
				bestkeep <- snames;
				sizef <- length(snames);
				CVOutput <- data.frame(cvdis$IDs,cvdis$outcomes,cvdis$caseDistance,cvdis$controlDistance);
				rdeltaES <- 0.9*rdeltaES+0.01;
				if (cvdis$cAUC >= maxAUC)
				{
					minkeep <- bestkeep;
					maxAUC <- cvdis$cAUC
				}
				if (ES > maxES) 
				{
					maxES <- ES;
				}
			}
			else
			{
				rdeltaES <- 0.8*rdeltaES+0.2*(maxES-oldES)/maxES;
			}
			cat(sprintf("%4d %s %4d %s %8.3f %s %8.3f %s %8.3f %s %8.5f \n",featureID,"Number of features:",sizef,"Max AUC:",maxAUC,"AUC:",cvdis$cAUC,"Z:",cvdis$zdis,"Rdelta:",rdeltaES));
			AUCevolution <- append(AUCevolution,cvdis$cAUC);
			Zevolution <- append(Zevolution,cvdis$zdis);
			featureSizeEvolution <- append(featureSizeEvolution,featureID);
			featureID <- featureID+1;
			if ((rdeltaES < 1.0e-4) || (featureID>length(varlist)))
			{
				stop <- TRUE;
			}
			oldES <- maxES;
		}
	}
	

	colnames(CVOutput) <- c("ID","Outcome","Case.Distance","Control.Distance");

	Ntemplate <- apply(controlsample[,minkeep],2,quantile,probs = theProbs,na.rm = TRUE);
	Ptemplate <- apply(casesample[,minkeep],2,quantile,probs = theProbs,na.rm = TRUE);

	LNtemplate <- apply(controlsample[,bestkeep],2,quantile,probs = theProbs,na.rm = TRUE);
	LPtemplate <- apply(casesample[,bestkeep],2,quantile,probs = theProbs,na.rm = TRUE);

	
	result <- list(controlTemplate=LNtemplate,
				   caseTamplate=LPtemplate,
				   OptControlTemplate=Ntemplate,
				   OptCaseTamplate=Ptemplate,
				   AUCevolution=AUCevolution,
				   Zevolution=Zevolution,
				   featureSizeEvolution=featureSizeEvolution,
				   featureList=bestkeep,
				   CVOutput=CVOutput,
				   maxES=maxES
				   );
  
	return (result);
}
