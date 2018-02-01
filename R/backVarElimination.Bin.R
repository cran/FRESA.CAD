backVarElimination_Bin <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),adjsize=1) 
{
  	seltype <- match.arg(selectionType)

	back.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI")) 
	{
		seltype <- match.arg(selectionType)
		type <- match.arg(type);
	  
		varsList <- unlist(as.list(attr(terms(object),"variables")))
		termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")
		
		
		cthr = abs(qnorm(pvalue));
		removeID = 0;

		outCome = paste(varsList[2]," ~ ");
		frm1 = outCome;
		if (length(termList)>0)
		{
			for ( i in 1:length(termList))
			{
				frm1 <- paste(frm1,paste("+",termList[i]));
			}
		}
		else
		{
			frm1 <- paste(frm1," 1 ")
		}
		
		ftmp <- formula(frm1);
		bckform <- frm1;
		FullModel <- modelFitting(ftmp,data,type,TRUE)
		startSearch = 1 + startOffset;
		if ( !inherits(FullModel, "try-error"))
		{
			FullPredict <- predict.fitFRESA(FullModel,data,'prob');
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
					if ( !inherits(redModel, "try-error"))
					{
						redPredict <- predict.fitFRESA(redModel,data,'prob');
						iprob <- .Call("improveProbCpp",redPredict,FullPredict,data[,Outcome]);
						if (seltype=="zIDI") 
						{
							ztst = iprob$z.idi;
						}
						else
						{
							ztst = iprob$z.nri;
						}
						if (is.na(ztst)) ztst=0;
						if (ztst<cthr)
						{
							cthr = ztst;
							removeID = i;
						}
					}
				}
			}
		}

		if ((length(termList) == startSearch) && (removeID == startSearch)) 
		{
			removeID = -1;
		}
		
		if (removeID > 0)
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
			FullModel <- modelFitting(ftmp,data,type,TRUE)
		}
#		cat("removed : ",removeID,"Final Model: \n")
#		print(summary(FullModel));
		result <- list(Model=FullModel,Removed=removeID,backfrm=bckform);

		return (result)
	}

	bkobj <- NULL;
	beforeFSC.formula <- NULL;
	if (adjsize>1)
	{
		bkobj <- backVarElimination_Bin(object,pvalue,Outcome,data,startOffset,type,selectionType,adjsize=1); 
		object <- bkobj$back.model;
		adjsize = floor(adjsize);
		adjsize <- min(adjsize,ncol(data)-1);
		beforeFSC.formula <- bkobj$string.formula;
#		cat("Adjusted Size:",adjsize,":",bkobj$beforeFSC.formula,"\n");
	}

	changes=1;
	loops=0;
    model <- object;
	beforeFSCmodel <- object;
	mydataFrame <- data;
	myOutcome <- Outcome;
	changes2 <- 0
	while ((changes>0) && (loops<100))
	{
		p.elimin <- pvalue;
		if (adjsize>1)
		{
			modsize <- length(attr(terms(model),"term.labels"));	
			if (modsize<1) modsize=1;
			qvalue <- 4.0*pvalue;
			if (qvalue < 0.05) qvalue=0.05 # lests keep the minimum q-value to 0.1
			p.elimin <- min(pvalue,modsize*qvalue/adjsize) # BH alpha the elimination p-value
		}

		bk <- back.var.IDISelection(model,p.elimin,Outcome=myOutcome,data=mydataFrame,startOffset,type,seltype);

		changes = as.integer(bk$Removed);
		if (changes>0)
		{
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
		  loops = loops + 1;
	}
#	print(summary(model));
#	cat("Reduced Model:",bk$backfrm,"\n")
	modelReclas <- getVar.Bin(model,data=mydataFrame,Outcome=myOutcome,type);
	result <- list(back.model=model,
	loops=loops,
	reclas.info=modelReclas,
	back.formula=formula(bk$backfrm),
	lastRemoved=changes2,
	at.opt.model=beforeFSCmodel,
	string.formula=bk$backfrm,
	beforeFSC.formula=formula(beforeFSC.formula));
	return (result);
}