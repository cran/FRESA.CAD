backVarElimination <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI")) 
{
  	seltype <- match.arg(selectionType)

	back.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",dataframe,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI")) 
	{
		seltype <- match.arg(selectionType)
		type <- match.arg(type);
	  
		varsList <- as.list(attr(terms(object),"variables"))
		
		
		cthr = abs(qnorm(pvalue));
		removeID = 0;

		outCome = paste(varsList[2]," ~ ");
		startlist = 3 ;
		frm1 = outCome;
		if (length(varsList)>=startlist)
		{
			for ( i in startlist:length(varsList))
			{
				frm1 <- paste(frm1,paste(" + ",varsList[i]));
			}
		}
		else
		{
			frm1 <- paste(frm1," 1 ")
		}
		
		ftmp <- formula(frm1);
		bckform <- frm1;
		fullModel <- modelFitting(ftmp,dataframe,type)
		if ( !inherits(fullModel, "try-error"))
		{
			fullPredict <- predictForFresa(fullModel,newdata=dataframe,type = 'prob');
			startSearch = startlist + startOffset;
			if (length(varsList)>startSearch)
			{
				for ( i in startSearch:length(varsList))
				{
				
					frm1 = outCome;
					for ( j in startlist:length(varsList))
					{
						if (i!=j)
						{
							frm1 <- paste(frm1,paste(" + ",varsList[j]));
						}
					}
					ftmp <- formula(frm1);
					redModel <- modelFitting(ftmp,dataframe,type)
					if ( !inherits(redModel, "try-error"))
					{
						redPredict <- predictForFresa(redModel,newdata=dataframe,type = 'prob');
						iprob <- improveProb(redPredict,fullPredict,dataframe[,Outcome]);
						if (seltype=="zIDI") 
						{
							ztst = iprob$z.idi;
						}
						else
						{
							ztst = iprob$z.nri;
						}
						if (ztst<cthr)
						{
							cthr = ztst;
							removeID = i;
						}
					}
				}
			}
		}

		if ((length(varsList) == startSearch) && (removeID == startSearch)) 
		{
			removeID = -1;
		}
		
		if (removeID > 0)
		{
			frm1 = outCome;
			for ( i in startlist:length(varsList))
			{
				if (i != removeID)
				{
					frm1 = paste(frm1,paste(" + ",varsList[i]));
				}
			}
			ftmp <- formula(frm1);
			bckform <- frm1;
			fullModel <- modelFitting(ftmp,dataframe,type)
		}
		cat("removed : ",removeID,"Final Model: \n")
		print(summary(fullModel));
		result <- list(Model=fullModel,Removed=removeID,backfrm=bckform);

		return (result)
	}


	changes=1;
	loops=0;
    model <- object;
	mydataFrame <- data;
	myOutcome <- Outcome;
	while ((changes>0) && (loops<100))
	{
		bk <- back.var.IDISelection(model,pvalue,Outcome=myOutcome,dataframe=mydataFrame,startOffset,type,seltype);
		changes = as.integer(bk$Removed);
		if (changes>0)
		{
		  loops = loops + 1;      
		}
		model = bk$Model;
	}
	print(summary(model));
	modelReclas <- getVarReclassification(model,dataframe=mydataFrame,Outcome=myOutcome,type);
	result <- list(back.model=model,loops=loops,reclas.info=modelReclas,back.formula=bk$backfrm,lastRemoved=changes);
	return (result);
}