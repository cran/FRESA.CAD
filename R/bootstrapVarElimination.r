bootstrapVarElimination <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops=250,fraction=1.00,print=TRUE,plots=TRUE) 
{
  	seltype <- match.arg(selectionType)



	boot.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops,fraction) 
	{
		seltype <- match.arg(selectionType)
		type <- match.arg(type);
	  
		varsList <- as.list(attr(terms(object),"variables"))
		
		climpvalue = pvalue;
		cthr = abs(qnorm(pvalue));
		removeID = 0;

		outCome = paste(varsList[2]," ~ 1");
		startlist = 3 ;
		frm1 = outCome;
		idiCV=NULL;
		if (startlist <= length(varsList))
		{
			for ( i in startlist:length(varsList))
			{
				frm1 <- paste(frm1,paste(" + ",varsList[i]));
			}
			ftmp <- formula(frm1);
		
#			cat ("Formula ",frm1,"\n")
			idiCV <- bootstrapValidation(fraction,loops,ftmp,Outcome,data,type,plots=plots)
	#		summary(idiCV);
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
					if (who != -1)
					{
						frm1 <- paste(frm1,paste(" + ",varsList[i]));
					}
				}
				else
				{
					removeID=idlist;
					cat ("Removed ",paste(" -> ",varsList[i]),"\n");
				}
			}
		}
#		cat ("Rm: Formula ",frm1,"\n")
		ftmp <- formula(frm1);
		fullModel <- modelFitting(ftmp,data,type,TRUE)
#		print(summary(fullModel))
		if (inherits(fullModel, "try-error"))
		{
			cat("Warning: Model fitting error\n"); 
		}
		result <- list(Model=fullModel,Removed=removeID,BootModel=idiCV,backfrm=frm1);

		return (result)
	}


	changes=1;
	loopsAux=0;
    model <- object;
	mydataFrame <- data;
	modelReclas <- NULL;
	myOutcome <- Outcome;
	changes2 <- 0
	while ((changes>0) && (loopsAux<100))
	{
		bk <- boot.var.IDISelection(model,pvalue,Outcome=myOutcome,data=mydataFrame,startOffset,type,seltype,loops,fraction);
		if (!inherits(bk$Model, "try-error"))
		{
			changes = as.integer(bk$Removed);
			if (changes>0)
			{
				loopsAux = loopsAux + 1;
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
		else
		{
			model = bk$Model;
			changes = 1;
			loopsAux = loopsAux + 1
		}
	}
	model <- modelFitting(formula(bk$backfrm),data,type);
#	cat("Last Removed ->",bk$Removed,"\n")
	if (bk$Removed>=0)
	{
		modelReclas <- getVarReclassification(model,data=mydataFrame,Outcome=myOutcome,type);
	}
	if ((fraction<1) ||  (changes == 0 ))
	{
		idiCV <- bootstrapValidation(1.0,loops,formula(bk$backfrm),myOutcome,mydataFrame,type,plots=plots);
	}
	else
	{
		idiCV <- bk$BootModel;  #		Copy Bootstrap model
	}
#	print(summary(model));
	cat("Reduced Formula:",bk$backfrm,"\n")
	result <- list(back.model=model,loops=loopsAux,reclas.info=modelReclas,bootCV=idiCV,back.formula=formula(bk$backfrm),lastRemoved=changes2);
	return (result);
}