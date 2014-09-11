bootstrapVarElimination <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),bootLoops=250,bootFraction=1.00) 
{
  	seltype <- match.arg(selectionType)



	boot.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",dataframe,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),bootLoops,bootFraction) 
	{
		seltype <- match.arg(selectionType)
		type <- match.arg(type);
	  
		varsList <- as.list(attr(terms(object),"variables"))
#		termList <- as.list(attr(terms(object),'term.labels'))
		
		climpvalue = pvalue;
		cthr = abs(qnorm(pvalue));
		removeID = 0;

		outCome = paste(varsList[2]," ~ 1");
		startlist = 3 ;
		frm1 = outCome;
		for ( i in startlist:length(varsList))
		{
			frm1 <- paste(frm1,paste(" + ",varsList[i]));
		}
		ftmp <- formula(frm1);
		idiCV <- bootstrapValidation(bootFraction,bootLoops,ftmp,Outcome,dataframe,type,plots=FALSE)
		startSearch = startlist + startOffset;
		frm1 = outCome;
		if ((startSearch-1) >= startlist)
		{
			for ( i in startlist:(startSearch-1))
			{
				frm1 <- paste(frm1,paste(" + ",varsList[i]));
			}
		}
		minlcl = cthr;
		idx = 2;
		who = 0;
		idlist=startOffset+1;
		for ( i in startSearch:length(varsList))
		{
			if (any(is.na(idiCV$z.IDIs[,idlist])))
			{
				who = i;
			}
			else
			{
				if (seltype=="zIDI")
				{
					ci <- as.vector(quantile(idiCV$z.IDIs[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
					ci2 <- as.vector(quantile(idiCV$test.z.IDIs[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
				}
				else
				{
					ci <- as.vector(quantile(idiCV$z.NRIs[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
					ci2 <- as.vector(quantile(idiCV$test.z.NRIs[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
				}
				cmin = min(ci[idx],ci2[idx]);
				if (cmin <= minlcl)
				{
					minlcl = cmin;
					who = i;
				}
			}
			idlist=idlist+1;
		}
		if ((length(varsList) == startSearch) && (who == startSearch)) 
		{
			who = -1;
			removeID = -1;
		}
		for ( i in startSearch:length(varsList))
		{
			if (who != i)
			{
				frm1 <- paste(frm1,paste(" + ",varsList[i]));
			}
			else
			{
				removeID=idlist;
			}
		}
		cat ("Formula ",frm1,"\n")
		ftmp <- formula(frm1);
		fullModel <- modelFitting(ftmp,dataframe,type)
		if (inherits(fullModel, "try-error"))
		{
			fullModel <- modelFitting(object,dataframe,type)
		}
		result <- list(Model=fullModel,Removed=removeID,BootModel=idiCV,backfrm=frm1);

		return (result)
	}


	changes=1;
	loops=0;
    model <- object;
	mydataFrame <- data;
	myOutcome <- Outcome;
	while ((changes>0) && (loops<100))
	{
		bk <- boot.var.IDISelection(model,pvalue,Outcome=myOutcome,dataframe=mydataFrame,startOffset,type,seltype,bootLoops,bootFraction);
		changes = as.integer(bk$Removed);
		if (changes>0)
		{
		  loops = loops + 1;      
		}
		model = bk$Model;
	}
	modelReclas <- getVarReclassification(model,dataframe=mydataFrame,Outcome=myOutcome,type);
	if ((bootFraction<1) ||  (changes <0 ))
	{
		idiCV <- bootstrapValidation(1.0,bootLoops,model$formula,myOutcome,mydataFrame,type,plots=FALSE);
	}
	else
	{
		idiCV <- bk$BootModel;
	}
	print(summary(model));

	result <- list(back.model=model,loops=loops,reclas.info=modelReclas,bootCV=idiCV,back.formula=bk$backfrm,lastRemoved=changes);
	return (result);
}