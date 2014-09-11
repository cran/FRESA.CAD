#
plotModels.ROC <-
function(modelPredictions,number.of.models=0,specificities=c(0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),...) 
{
	par(mfrow=c(1,1))
	rocadded = 0;
	if (number.of.models == 0)
	{
		number.of.models=max(modelPredictions[,"Model"])
	}
	auclist <- vector()
	sumSen <- NULL;
	blindSen <- NULL;
	auc1 <- roc(modelPredictions[,"Outcome"],modelPredictions[,"Blind_Prediction"],col="darkblue",auc=TRUE,plot=TRUE,smooth=FALSE,...)$auc
	par(new=TRUE)
	ley.names <- c(paste("Blind: Coherence (",sprintf("%.3f",auc1),")"))
	ley.colors <- c("darkblue")
	ley.lty <- c(1)
	for (n in 1:number.of.models)
	{
		blindmodel <- modelPredictions[which(modelPredictions[,3] == n),];
		if ( (sum(blindmodel[,"Outcome"]==1) > 3) && (sum(blindmodel[,"Outcome"]==0) > 3))
		{
			auclist <- append(auclist,roc(blindmodel[,"Outcome"],blindmodel[,"Blind_Prediction"],auc=TRUE,plot=TRUE,col="lightgray",lty=4,lwd=1)$auc)
			par(new=TRUE)
			sen <- roc(blindmodel[,"Outcome"],blindmodel[,"Blind_Prediction"],auc=TRUE,plot=FALSE,ci=TRUE,of='se',specificities=specificities,boot.n=100,smooth=FALSE,lty=3,lwd=1)$ci[,2]
			if (n == 1) 
			{
				blindSen <- sen;
			}
			else
			{
				blindSen <- rbind(sen,blindSen);
			}
			rocadded = rocadded +1;
		}
	}
	auc = 0;
	if (rocadded>0)
	{
		boxplot(blindSen,add=TRUE, axes = FALSE,boxwex=0.04,at=specificities);
		sumSen <- colMeans(blindSen,na.rm = TRUE);
		sennames <- names(sumSen);
		sumSen <- append(0,sumSen);
		sumSen <- append(sumSen,1);
		sennames <- append("1",sennames);
		sennames <- append(sennames,"0");
		names(sumSen) <- sennames;
		spevalues <- as.numeric(names(sumSen));
		lines(spevalues,sumSen,col="red",lwd=2.0);
		for (i in 2:length(spevalues))
		{
			auc = auc + (spevalues[i-1]-spevalues[i])*(sumSen[i-1]+(sumSen[i]-sumSen[i-1])/2)
		}
		ley.names <- append(ley.names,c("Blind: ROCs",paste("Blind: Mean Sensitivities(",sprintf("%.3f",auc),")")));
		ley.colors <- append(ley.colors,c("lightgray","red"));
		ley.lty <- append(ley.lty,c(4,1));
	}

	legend(0.6,0.30, ley.names,col = ley.colors, lty = ley.lty,bty="n")

	result <- list(ROC.AUCs=auclist,
	mean.sensitivities=sumSen,
	model.sensitivities=blindSen,
	specificities=specificities,
	senAUC=auc)
	return (result)
}
