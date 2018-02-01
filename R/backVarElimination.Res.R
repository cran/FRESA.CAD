backVarElimination_Res <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),setIntersect=1,adjsize=1) 
{

	back.var.NeRISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),setIntersect=1) 
	{
		type <- match.arg(type);
	  
		varsList <- unlist(as.list(attr(terms(object),"variables")))
		termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")
		
		modsize <- length(termList);
		
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
		frm1 = outCome;
		cpv = pvalue; 
		if (length(termList)>0)
		{
			for ( i in 1:length(termList))
			{
				frm1 <- paste(frm1,paste("+",termList[i]));
			}
		}
#		cat ("Len: ",length(termList)," : ",frm1,"\n")
		ftmp <- formula(frm1);
		bckform <- frm1;
		startSearch = 1 + startOffset;
		if (length(termList)>1)
		{
			FullModel <- modelFitting(ftmp,data,type,TRUE)
			FullResiduals <- residualForFRESA(FullModel,data,Outcome);

			if (length(termList)>startSearch)
			{
				for ( i in startSearch:length(termList))
				{
				
					frm1 = outCome;
					for ( j in 1:length(termList))
					{
						if (i!=j)
						{
							frm1 <- paste(frm1,paste("+",termList[j]));
						}
					}
					ftmp <- formula(frm1);
					redModel <- modelFitting(ftmp,data,type,TRUE)
					if (inherits(redModel, "try-error"))
					{
						redModel <- FullModel
					}
					

					redResiduals <- residualForFRESA(redModel,data,Outcome);
					iprob <- improvedResiduals(redResiduals,FullResiduals,testType);
					if (iprob$p.value>cpv)
					{
						cpv = iprob$p.value;
						removeID = i;
					}
				}
			}
			if ((length(termList) == startSearch) && (removeID == startSearch)) 
			{
				removeID = -1;
			}

			if (removeID>0)
			{
				frm1 = outCome;
				for ( i in 1:length(termList))
				{
					if (i != removeID)
					{
						frm1 = paste(frm1,paste("+",termList[i]));
					}
				}
				ftmp <- formula(frm1);
				bckform <- frm1;
#				cat("formula :",frm1,"\n")
				FullModel <- modelFitting(ftmp,data,type,TRUE)
			}
		}
		else
		{
			FullModel <- object;
		}
#		cat("Reduced: \n")
#		print(summary(FullModel));
		result <- list(Model=FullModel,Removed=removeID,backfrm=bckform);

		return (result)
	}

	bkobj <- NULL;
	beforeFSC.formula <- NULL;
	if (adjsize>1)
	{
		bkobj <- backVarElimination_Res(object,pvalue,Outcome,data,startOffset,type,testType,setIntersect,adjsize=1); # remove features that do not improve residuals
		object <- bkobj$back.model;
		adjsize = floor(adjsize);
		adjsize <- min(adjsize,ncol(data)-1);
		beforeFSC.formula <- bkobj$string.formula;
#		cat("Adjusted Size:",adjsize,"\n");
	}

	changes=1;
	loops=0;
    model <- object;
	beforeFSCmodel <- object;
	changes2<-0
	while ((changes>0) && (loops<100))
	{
		p.elimin <- pvalue;
		if (adjsize>0)
		{		
			modsize <- length(attr(terms(model),"term.labels"));	
			if (modsize<1) modsize=1;
			qvalue <- 4.0*pvalue;
			if (qvalue < 0.05) qvalue=0.05 # lests keep the minimum q-value to 0.1
			p.elimin <- min(pvalue,modsize*qvalue/adjsize) # BH alpha  the elimination p-value
		}

		bk <- back.var.NeRISelection(model,p.elimin,Outcome=Outcome,data=data,startOffset,type,testType,setIntersect);
#		cat("Used p :",p.elimin,"Formula <- ", bk$backfrm,"\n");
		changes = as.integer(bk$Removed);
		if (changes>0)
		{
		  loops = loops + 1;
			changes2<- attr(terms(model),"term.labels")[which(!(attr(terms(model),"term.labels") %in% attr(terms(bk$Model),"term.labels")))]
			model = bk$Model;
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
#	cat("Reduced Model:",bk$backfrm,"\n")
#	print(summary(model));

	modelReclas <- getVar.Res(model,data=data,Outcome=Outcome,type=type);
	
	result <- list(back.model= model,
	loops=loops,
	reclas.info=modelReclas,
	back.formula=formula(bk$backfrm),
	bootCV=NULL,
	lastRemoved=changes2,
	number.of.independent=adjsize,
	at.opt.model=beforeFSCmodel,
	string.formula=bk$backfrm,
	beforeFSC.formula=formula(beforeFSC.formula));
	return (result);
}