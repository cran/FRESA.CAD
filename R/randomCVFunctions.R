#Filter Methods


FRESAcacheEnv <- new.env();

correlated_Remove <- function(data= NULL,fnames= NULL,thr=0.999)
{
#	cat(thr,":",length(fnames)," ->");
	if (length(fnames)>1)
	{
#		print(fnames);
		colsd <- apply(data[,fnames],2,sd,na.rm = TRUE);
		fnames <- fnames[colsd > 0];
		
		corm <- abs(cor(data[,fnames],method="spearman"));
		corm[diag(length(fnames)) == 1] <- 0;
		keep <- numeric(length(fnames));
		for (i in length(fnames):1)
		{
			keep[i] <- (sum(corm[,i] > thr) == 0);
			if (!keep[i])  
			{
				corm[i,] <- 0;
			}
		}
		attributes(fnames) <- list(removed=fnames[keep == 0]);
		fnames <- fnames[keep == 1];
	}
#	cat(length(fnames),"\n");

	return (fnames);
}

correlated_RemoveToLimit <- function(data,unitPvalues,limit=0,thr=0.975,maxLoops=50,minCorr=0.80)
{
	if ((limit >= 0) && (thr < 1.0))
	{
		ntest <- 0;
		cthr <- thr;
		if (limit > 0)
		{
			slimit <- limit;
			if (limit <= 1)
			{
				slimit <- limit*nrow(data);
			}
			if (slimit < 2) 
			{
				slimit <- 2;
			}
			if (slimit >= nrow(data))
			{
				slimit = nrow(data) - 1;
			}
			if (length(unitPvalues) > (10*slimit))
			{
				unitPvalues <- unitPvalues[1:(10*slimit)];
			}
			while ( (length(unitPvalues) > slimit) && (ntest < maxLoops) && (cthr > minCorr) )
			{
				unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues),cthr)];
				ntest = ntest + 1;
				cthr = cthr*thr;
			}
			if (length(unitPvalues) > slimit)
			{
				unitPvalues <- unitPvalues[1:slimit];
			}
		}
		else
		{
			unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues),cthr)];
		}
	}
	return (unitPvalues);

}


univariate_Logit <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH", uniTest=c("zIDI","zNRI"),limit=0,...,n = 0)
{
	varlist <-colnames(data);
	varlist <- varlist[Outcome != varlist];
	varlist <- cbind(varlist,varlist);
	uniTest <- match.arg(uniTest);
	univ <- ForwardSelection.Model.Bin(nrow(varlist),1.0,0.0,1,"1",Outcome,varlist,data,1,type="LOGIT",selectionType=uniTest);
	unitPvalues <- (1.0 - pnorm(univ$base.Zvalues));
	names(unitPvalues) <-  varlist[,1];
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	unitPvalues <- unitPvalues[order(unitPvalues)];
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues < pvalue];
#	print(unitPvalues)
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
   return(unitPvalues);
}

univariate_residual <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",uniTest=c("Ftest","Binomial","Wilcox","tStudent"),type=c("LM","LOGIT"),limit=0,...,n = 0)
{
	varlist <- colnames(data);
	varlist <- varlist[Outcome != varlist];
	varlist <- cbind(varlist,varlist)
	uniTest <- match.arg(uniTest);
	type <- match.arg(type);

	univ <- ForwardSelection.Model.Res(nrow(varlist),1.0,0.0,1,"1",Outcome,varlist,data,1,type=type,testType=uniTest);
	unitPvalues <- (1.0 - pnorm(univ$base.Zvalues));
#	print(unitPvalues);
	names(unitPvalues) <- varlist[,1];
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	unitPvalues <- unitPvalues[order(unitPvalues)];
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues < pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	return(unitPvalues);
}

univariate_Wilcoxon <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{

	varlist <-colnames(data);
	case <- subset(data,get(Outcome) == 1);
	control <- subset(data,get(Outcome) == 0);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- numeric(length(varlist));
	names(unitPvalues) <- varlist;

	for (j in varlist) 
	{
		 pval <- wilcox.test(control[,j],case[,j],na.action = na.exclude)$p.value; 
		 if (inherits(pval, "try-error")) {pval <- 1.0;}
		 if (is.null(pval)) { pval <- 1.0; }
		 if (is.na(pval)) {pval <- 1.0;}
		 unitPvalues[j] <- pval;
	}
  
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	unitPvalues <- unitPvalues[order(unitPvalues)];
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues < pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	
   return(unitPvalues);
}

univariate_tstudent <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{

	varlist <-colnames(data);
	case <- subset(data,get(Outcome) == 1);
	control <- subset(data,get(Outcome) == 0);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- numeric(length(varlist));
	names(unitPvalues) <- varlist;
	for (j in varlist) 
	{
		pval <- try(t.test(control[,j],case[,j],na.action = na.exclude)$p.value);
		if (inherits(pval, "try-error")) {pval <- 1.0;}
		if (is.null(pval)) { pval <- 1.0; }
		if (is.na(pval)) {pval <- 1.0;}
		unitPvalues[j] <-  pval;	
	}
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	unitPvalues <- unitPvalues[order(unitPvalues)];
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues < pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
   return(unitPvalues);
}

univariate_correlation <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH", method = "kendall",limit=0,...,n = 0)
{

	varlist <-colnames(data);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- numeric(length(varlist));
	names(unitPvalues) <- varlist;

	for (j in varlist) 
	{
		pval <- try(cor.test(data[,Outcome],data[,j],na.action = na.exclude,method = method)$p.value);
		if (inherits(pval, "try-error")) {pval <- 1.0;}
		if (is.null(pval)) { pval <- 1.0; }
		if (is.na(pval)) {pval <- 1.0;}
		unitPvalues[j] <-  pval;
	}
  
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	unitPvalues <- unitPvalues[order(unitPvalues)];
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues < pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
   return(unitPvalues);
}

mRMR.classic_FRESA <- function(data=NULL, Outcome=NULL,feature_count=0,...)
{
	if (!requireNamespace("mRMRe", quietly = TRUE)) {
	   install.packages("mRMRe", dependencies = TRUE)
	} 

	if (feature_count == 0)
	{
		feature_count <- min(c(nrow(data)-1,ncol(data)-1));
	}
	fs <- try(mRMRe::mRMR.classic(data = mRMRe::mRMR.data(data), target_indices = grep(Outcome, colnames(data)), feature_count = feature_count,...));
	if ( !inherits(fs, "try-error"))
	{
		result <- as.vector(mRMRe::scores(fs)[[1]]);
		names(result) <- colnames(data)[mRMRe::solutions(fs)[[1]]];
		result <- result[!is.na(result)];
		gz <- (result >= 0);
		if (sum(1*gz) > 5)	result <- result[gz];
	}
	else
	{
		warning("mRMR.classic Error. I'll make all columns numeric and try again\n")
		data[,1:ncol(data)] <- sapply(data,as.numeric)
		fs <- try(mRMRe::mRMR.classic(data = mRMRe::mRMR.data(data), target_indices = grep(Outcome, colnames(data)), feature_count = feature_count,...));
		if ( !inherits(fs, "try-error"))
		{
			result <- as.vector(mRMRe::scores(fs)[[1]]);
			names(result) <- colnames(data)[mRMRe::solutions(fs)[[1]]];
			result <- result[!is.na(result)];
			gz <- (result >= 0);
			if (sum(1*gz) > 5)	result <- result[gz];
		}
		else
		{
			warning("mRMR.classic Error. Returning the first elements of the data-frame\n")
			result <- 1:feature_count;
			varlist <- colnames(data);
			varlist <- varlist[Outcome != varlist];
			names(result) <- varlist[result];
			result <- result/feature_count;
		}
	}
	return(result);
}

sperman95ci <- function(datatest,nss=4000)
{
	sz <- nrow(datatest)
	sesci <- c(0,0,0);
	if (sz>2)
	{
		ses <- numeric(nss);
		for (i in 1:nss)
		{
		  bootsample <- datatest[sample(sz,sz,replace=TRUE),];
		  ses[i] <- cor(bootsample[,1],bootsample[,2],method="spearman");
		}
		sesci <- quantile(ses, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	}
	return (sesci);
}

MAE95ci <- function(datatest,nss=4000)
{
	sz <- nrow(datatest)
	sesci <- c(0,0,0);
	if (sz>1)
	{
		ses <- numeric(nss);
		for (i in 1:nss)
		{
		  bootsample <- datatest[sample(sz,sz,replace=TRUE),];
		  ses[i] <- mean(abs(bootsample[,1]-bootsample[,2]));
		}
		sesci <- quantile(ses, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	}
	return (sesci);
}
	
ClassMetric95ci <- function(datatest,nss=4000)
{
	  sz <- nrow(datatest)
	  accci <- c(0,0,0);
	  senci <- c(0,0,0);
	  aucci <- c(0,0,0);
	  berci <- c(0,0,0);
	  preci <- c(0,0,0);
	  F1ci <- c(0,0,0);
	  tbstat <- as.matrix(table(datatest[,1]));
	  nscores <- nrow(tbstat);
	  scores <- as.numeric(rownames(tbstat));
	  if (sz > 1)
	  {
		acc <- numeric(nss);
		sen <- numeric(nss);
		csen <- numeric(nscores);
		cspe <- numeric(nscores);
		auc <- numeric(nss);
		ber <- numeric(nss);
		pre <- numeric(nss);
		F1 <- numeric(nss);
		cpre <- numeric(nscores);
		cF1 <- numeric(nscores);
		for (i in 1:nss)
		{
		  bootsample <- datatest[sample(sz,sz,replace = TRUE),];
		  acc[i] <- mean(bootsample[,1] == bootsample[,2]);
		  for (n in 1:nscores)
		  {
			scoresamples <- sum(bootsample[,1] == scores[n])
			if (scoresamples > 0)
			{
				csen[n] <- sum( (bootsample[,1] == scores[n]) & (bootsample[,2] == scores[n]) )/scoresamples;
			}
			else
			{
				csen[n] <- 0;
			}
			if (scoresamples < sz)
			{
				cspe[n] <- sum( (bootsample[,1] != scores[n]) & (bootsample[,2] != scores[n]) )/(sz - scoresamples);
			}
			else
			{
				cspe[n] <- 0;
			}
			tp <- sum( (bootsample[,2] == scores[n]) );
			if (tp > 0)
			{			
				cpre[n] <- sum( (bootsample[,1] == scores[n]) & (bootsample[,2] == scores[n]) )/tp ;
			}
			else
			{
				cpre[n] <- NA;
			}
			if (is.na(csen) || is.na(cpre))
			{
				cF1[n] <- NA;
			}
			else
			{
				cF1[n] <- 2.0*csen[n]*cpre[n]/(csen[n]+cpre[n]);
			}
			sen[i] <- sen[i] + csen[n];
			auc[i] <- auc[i] + (csen[n] + cspe[n])/2;
			ber[i] <- ber[i] + 1.0-(csen[n] + cspe[n])/2;
			pre[i] <- pre[i] + cpre[n];
			F1[i] <- F1[i] + cF1[n];
		  }
		  sen[i] <- sen[i]/nscores;
		  auc[i] <- auc[i]/nscores;
		  ber[i] <- ber[i]/nscores;
		  pre[i] <- pre[i]/nscores;
		  F1[i] <- F1[i]/nscores;
		}
		accci <- quantile(acc, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		senci <- quantile(sen, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		aucci <- quantile(auc, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		berci <- quantile(ber, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		F1ci <- quantile(F1, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		preci <- quantile(pre, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	  }
	  return(list(accci = accci,senci = senci,aucci = aucci,berci = berci,preci = preci,F1ci = F1ci));
}


