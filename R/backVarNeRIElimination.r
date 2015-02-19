backVarNeRIElimination <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),setIntersect=1) 
{

	back.var.NeRISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),setIntersect=1) 
	{
		type <- match.arg(type);
	  
		varsList <- as.list(attr(terms(object),"variables"))
		termList <- as.list(attr(terms(object),"term.labels"))
		removeID = 0;

		outCome = paste(varsList[2]," ~ ");
		if (setIntersect==0) 
		{
			outCome = paste(outCome," 0  ");
		}
		else
		{
			outCome = paste(outCome," 1  ");
		}
		startlist = 3 ;
		frm1 = outCome;
		cpv = pvalue; 
		if (length(varsList)>=startlist)
		{
			for ( i in startlist:length(varsList))
			{
				frm1 <- paste(frm1,paste(" + ",varsList[i]));
			}
		}
#		cat ("Len: ",length(termList)," : ",frm1,"\n")
		ftmp <- formula(frm1);
		bckform <- frm1;
		if (length(termList)>1)
		{
			fullModel <- modelFitting(ftmp,data,type)
			fullResiduals <- residualForNeRIs(fullModel,data,Outcome);

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
					redModel <- modelFitting(ftmp,data,type)
					if (inherits(redModel, "try-error"))
					{
						redModel <- fullModel
					}
					

					redResiduals <- residualForNeRIs(redModel,data,Outcome);
					iprob <- improvedResiduals(redResiduals,fullResiduals,testType);
					if (iprob$p.value>cpv)
					{
						cpv = iprob$p.value;
						removeID = i;
					}
				}
			}
			if ((length(varsList) == startSearch) && (removeID == startSearch)) 
			{
				removeID = -1;
			}

			if (removeID>0)
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
				cat("formula :",frm1,"\n")
				fullModel <- modelFitting(ftmp,data,type)
			}
		}
		else
		{
			fullModel <- object;
		}
#		cat("Reduced: \n")
#		print(summary(fullModel));
		result <- list(Model=fullModel,Removed=removeID,backfrm=bckform);

		return (result)
	}


	changes=1;
	loops=0;
    model <- object;
	changes2<-0
	while ((changes>0) && (loops<100))
	{
		bk <- back.var.NeRISelection(model,pvalue,Outcome=Outcome,data=data,startOffset,type,testType,setIntersect);
		changes = as.integer(bk$Removed);
		if (changes>0)
		{
		  loops = loops + 1;
			changes2<- as.character(as.list(attr(terms(model),"variables")))[which(!(as.character(as.list(attr(terms(model),"variables")))%in%as.character(as.list(attr(terms(bk$Model),"variables")))))]
			if (length(changes2)>1)
			{
				changes2<-changes2[2]
			}
		}
		if (changes < 0)
		{
			changes2<- changes
		}
		 model = bk$Model;
	}
	cat(bk$backfrm," : Reduced Model:\n")
#	print(summary(model));

	
	modelReclas <- getVarNeRI(model,data=data,Outcome=Outcome,type=type);

	
	result <- list(back.model= model,loops=loops,reclas.info=modelReclas,back.formula=formula(bk$backfrm),bootCV=NULL,lastRemoved=changes2);
	return (result);
}