#' @method update fitFRESA
summary.fitFRESA <- function(object,type=c("Improvement","Residual"),ci=c(0.025,0.975),data=NULL,...)
{
	coefficients <- NULL;
	type <- match.arg(type);
	coffset <- 1.0*(object$type == "COX") 
	outcome <- as.character(all.vars(object$formula)[1+coffset]);
	FRESAsummary <- NULL;
	bv <- NULL;
	if (length(object$coefficients)>1)
	{
		if (is.null(data)) 
		{
			data <- object$model;
		}
		if ((object$type=="LM") || (type=="Residual"))
		{
			coefficients <- object$coefficients[-1];
			ncoef <- length(coefficients);
			cis <- matrix(rep(NA,ncoef*3),ncoef,3);
			cinames <- c("lower","mean","upper");
			if (hasName(object,"baggingAnalysis"))
			{
				vres <- list();
				vres$unitrainMSE <- object$baggingAnalysis$uMS_values;
				vres$redtrainMSE <- object$baggingAnalysis$rMS_values;
				vres$NeRIs <- object$baggingAnalysis$NeRI_values;
				vres$tP.value <- object$baggingAnalysis$pt_values;
				vres$WilcoxP.value <- object$baggingAnalysis$pWilcox_values;
				vres$FP.value <- object$baggingAnalysis$pF_values;
				vres$BinP.value <- object$baggingAnalysis$pBin_values;
				vres$FullTrainMSE <- object$baggingAnalysis$mMSE_values;
				if (object$baggingAnalysis$n_bootstrap>2)
				{
					tvalue1 <- qt(ci[1],object$baggingAnalysis$n_bootstrap-1);
					tvalue2 <- qt(ci[2],object$baggingAnalysis$n_bootstrap-1);
					cis = cbind(object$baggingAnalysis$coefficients+tvalue1*object$baggingAnalysis$coefstd,object$baggingAnalysis$coefficients,object$baggingAnalysis$coefficients+tvalue2*object$baggingAnalysis$coefstd);
					cinames <- c("lower","mean","upper");
				}
			}
			else
			{
				vres <- getVar.Res(object,data=data,Outcome=outcome,type=object$type);
				bv <- bootstrapValidation_Res(model.formula=object$formula,Outcome=outcome,data=data,type=object$type,...)
				coffset <- 1.0*(object$type != "COX") 
				for (i in 1:ncoef)
				{
					cis[i,] = as.vector(quantile(bv$s.coef[,i+coffset], probs = c(ci[1],0.5,ci[2]), na.rm = TRUE,names = FALSE, type = 7));
				}
				cinames <- c("lower","median","upper");
			}
			coefficients=cbind(coefficients,cis,vres$unitrainMSE,vres$redtrainMSE,vres$FullTrainMSE,vres$NeRIs,vres$FP.value,vres$tP.value,vres$BinP.value,vres$WilcoxP.value);
			colnames(coefficients) <- c("Estimate",cinames,"u.MSE","r.MSE","model.MSE","NeRI","F.pvalue","t.pvalue","Sign.pvalue","Wilcox.pvalue");
			residaulsMSE <- mean((object$response[,1]-predict(object))^2);
			Rsquare <- var(object$response[,1]);
			Rsquare <- (Rsquare-residaulsMSE)/Rsquare;
			FRESAsummary <- list(coefficients=coefficients,MSE=residaulsMSE,R2=Rsquare,bootstrap=bv);
		}
		else
		{
			coefficients <- object$coefficients[-1];
			ncoef <- length(coefficients);
			cis <- matrix(rep(NA,ncoef*3),ncoef,3);
			cinames <- c("lower","mean","upper");
			if (hasName(object,"baggingAnalysis"))
			{
				vres <- list();
				vres$uniTrainAccuracy <- object$baggingAnalysis$uAcc_values;
				vres$redtrainAccuracy <- object$baggingAnalysis$rAcc_values;
				vres$uniTrainAUC <- object$baggingAnalysis$uAUC_values;
				vres$redtrainAUC <- object$baggingAnalysis$rAUC_values;
				vres$IDIs <- object$baggingAnalysis$idi_values;
				vres$NRIs <- object$baggingAnalysis$nri_values;
				vres$z.IDIs <- object$baggingAnalysis$zidi_values;
				vres$z.NRIs <- object$baggingAnalysis$znri_values;
				vres$fullTrainAccuracy <- object$baggingAnalysis$mACC_values;
				vres$fullTrainAUC <- object$baggingAnalysis$mAUC_values;
				if (object$baggingAnalysis$n_bootstrap>2)
				{
					tvalue1 <- qt(ci[1],object$baggingAnalysis$n_bootstrap-1);
					tvalue2 <- qt(ci[2],object$baggingAnalysis$n_bootstrap-1);
					cis = cbind(object$baggingAnalysis$coefficients+tvalue1*object$baggingAnalysis$coefstd,object$baggingAnalysis$coefficients,object$baggingAnalysis$coefficients+tvalue2*object$baggingAnalysis$coefstd)
					cis <- exp(cis)
					cinames <- c("lower","OR","upper");
				}
			}
			else
			{
				vres <- getVar.Bin(object,data=data,Outcome=outcome,type=object$type);
				bv <- bootstrapValidation_Bin(model.formula=object$formula,Outcome=outcome,data=data,type=object$type,...)
				coffset <- 1.0*(object$type != "COX") 
				for (i in 1:ncoef)
				{
					cis[i,] = as.vector(quantile(bv$s.coef[,i+coffset], probs = c(ci[1],0.5,ci[2]), na.rm = TRUE,names = FALSE, type = 7));
				}
				cinames <- c("lower","median","upper");
			}
			coefficients=cbind(coefficients,cis,vres$uniTrainAccuracy,vres$redtrainAccuracy,vres$fullTrainAccuracy,vres$uniTrainAUC,vres$redtrainAUC,vres$fullTrainAUC,vres$IDIs,vres$NRIs,vres$z.IDIs,vres$z.NRIs);
			colnames(coefficients) <- c("Estimate",cinames,"u.Accuracy","r.Accuracy","full.Accuracy","u.AUC","r.AUC","full.AUC","IDI","NRI","z.IDI","z.NRI");
			pred <- 1*(predict(object)>0);
			Accuracy <- mean(1.0*(object$response[,1]==pred));
			tAUC <- sum((object$response[,1]==pred)*object$response[,1])/sum(object$response[,1]);
			tAUC <- 0.5*(tAUC+sum((object$response[,1]==pred)*(object$response[,1]==0))/sum(object$response[,1]==0));
			FRESAsummary <- list(coefficients=coefficients,Accuracy=Accuracy,tAUC=tAUC,bootstrap=bv);
		}
	}
	return (FRESAsummary);
}
