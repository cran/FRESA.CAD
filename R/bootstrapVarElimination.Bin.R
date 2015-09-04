bootstrapVarElimination_Bin <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops=250,fraction=1.00,print=TRUE,plots=TRUE,adjsize=1) 
{
  	seltype <- match.arg(selectionType)



	boot.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops,fraction,idiCV=NULL) 
	{
		seltype <- match.arg(selectionType)
		type <- match.arg(type);
	  
		varsList <- as.list(attr(terms(object),"variables"))
		
		climpvalue = max(5*pvalue,0.35); #
		cthr = abs(qnorm(pvalue));
		removeID = 0;

		outCome = paste(varsList[2]," ~ 1");
		startlist = 3 ;
		frm1 = outCome;
		
		if (startlist <= length(varsList))
		{
			for ( i in startlist:length(varsList))
			{
				frm1 <- paste(frm1,paste(" + ",varsList[i]));
			}
			ftmp <- formula(frm1);
		
#			cat ("Formula ",frm1,"\n")
			if (is.null(idiCV)) idiCV <- bootstrapValidation_Bin(fraction,loops,ftmp,Outcome,data,type,plots=plots)
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
					cmin = min(ci[2],ci2[2]);
					if (cmin <= minlcl) 
					{
						minlcl = cmin;
						who = i;
					}
					if ((ci[1] < cthr)&&(minlcl == cthr))
					{
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
#					cat ("Removed ",paste(" -> ",varsList[i]),"\n");
				}
			}
		}
#		cat ("Rm: Formula ",frm1,"\n")
		ftmp <- formula(frm1);
		FullModel <- modelFitting(ftmp,data,type)
#		print(summary(FullModel))
		if (inherits(FullModel, "try-error"))
		{
			cat("Warning: Model fitting error\n"); 
		}
		if (removeID>0) idiCV <- bootstrapValidation_Bin(fraction,loops,ftmp,Outcome,data,type,plots=plots)
		result <- list(Model=FullModel,Removed=removeID,BootModel=idiCV,backfrm=frm1);

		return (result)
	}

	bkobj <- NULL;
	
	maxROCAUC <- 1.0;
	stopROCACU <- 0.975;	# threshold for stop for decrease in AUC
	if (adjsize>1)
	{
		bkobj <- bootstrapVarElimination_Bin(object,pvalue,Outcome,data,startOffset,type,selectionType,loops,fraction,print,plots,adjsize=1); 
		object <- bkobj$back.model;
		adjsize = floor(adjsize);
		adjsize <- min(adjsize,ncol(data)-1);
		if (!is.null(bkobj$bootCV)) 
		{
			maxROCAUC <- bkobj$bootCV$blind.ROCAUC$auc;
		}
		stopROCACU <- 0.80;	# threshold for stop for decrease in AUC
	}
	else
	{
		idiCV <- bootstrapValidation_Bin(fraction,loops,formula(object),Outcome,data,type,plots=plots)
		maxROCAUC <- idiCV$blind.ROCAUC$auc;
	}

#	cat("AUC :",maxROCAUC,"\n");
	startAUC = maxROCAUC;
	
	changes=1;
	loopsAux=0;
    model <- object;
	mydataFrame <- data;
	modelReclas <- NULL;
	myOutcome <- Outcome;
	changes2 <- 0

	varsList <- as.list(attr(terms(object),"variables"))
	outCome = paste(varsList[2]," ~ 1");
	startlist = 3 ;
	frm1 = outCome;
	for ( i in startlist:length(varsList))
	{
		frm1 <- paste(frm1,paste(" + ",varsList[i]));
	}
	beforeFSCmodel.formula <- frm1;
	model.formula <- frm1;
	
	
	beforeFSCmodel <- object;

	wts <- rep(1,length(beforeFSCmodel$coefficients));
	names(wts) <- names(beforeFSCmodel$coefficients);
	lastboot = NULL;
#	cat ("Start AUC :",startAUC,"\n")
	while ((changes>0) && (loopsAux<100))
	{
		p.elimin <- pvalue;
		if (adjsize>1)
		{
			modsize <- length(as.list(attr(terms(model),"term.labels")));	
			if (modsize<1) modsize=1;
			p.elimin <- min(pvalue,2*modsize*pvalue/adjsize) # # BH alpha  the elimination p-value
		}

		bk <- boot.var.IDISelection(model,p.elimin,Outcome=myOutcome,data=mydataFrame,startOffset,type,seltype,loops,fraction,lastboot);
#		cat("Used p :",p.elimin,"Formula <- ", bk$backfrm,"\n");

		if (!inherits(bk$Model, "try-error"))
		{
			lastboot <- bk$BootModel;
			changes = as.integer(bk$Removed);
			weight <- 1.0/2.0;
#			cat ("AUC :",bk$BootModel$blind.ROCAUC$auc,"\n")
			if (!is.null(bk$BootModel))
			{			
				if (bk$BootModel$blind.ROCAUC$auc < stopROCACU*maxROCAUC)
				{
					changes = 0;
				}
				if (bk$BootModel$blind.ROCAUC$auc <= maxROCAUC)
				{
					weight <- 1.0;
				}
				else
				{
					maxROCAUC <- bk$BootModel$blind.ROCAUC$auc;
					weight <- 0.0;
				}
			}
			if (changes>0)
			{
				loopsAux = loopsAux + 1;
				changes2<- as.character(as.list(attr(terms(model),"variables")))[which(!(as.character(as.list(attr(terms(model),"variables")))%in%as.character(as.list(attr(terms(bk$Model),"variables")))))]
				model = bk$Model;
				model.formula <- bk$backfrm;

				if (length(changes2)>1)
				{
					changes2<-changes2[2]
				}
				if (adjsize>1)
				{
					for (i in 1:length(beforeFSCmodel$coefficients))
					{
						notadded = TRUE;
						for (j in 1:length(model$coefficients))
						{
							if (names(beforeFSCmodel$coefficients)[i] == names(model$coefficients)[j])
							{
								beforeFSCmodel$coefficients[i] <- (weight*beforeFSCmodel$coefficients[i] + (1.0-weight)*model$coefficients[j]); # it will average the two
								notadded=FALSE;
							}
						}
						if (notadded)
						{
							beforeFSCmodel$coefficients[i] <- weight*beforeFSCmodel$coefficients[i]/2; # it will average with zero
							wts[i] = wts[i]/2;
						}
					}
				}
			}
			if (changes < 0)
			{
				changes2<- changes
			}
		}
		else
		{
			changes = 1;
			loopsAux = loopsAux + 1
		}
	}
	if (str_count(model.formula,"\\+")>1)
	{
		modelReclas <- getVar.Bin(model,data=mydataFrame,Outcome=myOutcome,type);
	}
	idiCV <- bootstrapValidation_Bin(1.0000,loops,formula(model.formula),myOutcome,mydataFrame,type,plots=plots);
	if ((adjsize>1)&&(print == TRUE))
	{
		cat("Before FSC Mod:",beforeFSCmodel.formula,"\n")
		cat("Reduced Model :",model.formula,"\n")
		cat("Adjust size:",adjsize,"\n");
		cat("Start AUC:",startAUC,"last AUC:",idiCV$blind.ROCAUC$auc,"Max AUC:",maxROCAUC,"\n")
	}
	result <- list(back.model=model,
	loops=loopsAux,
	reclas.info=modelReclas,
	bootCV=idiCV,
	back.formula=formula(model.formula),
	lastRemoved=changes2,
	beforeFSC.model=beforeFSCmodel,
	beforeFSC.formula=formula(beforeFSCmodel.formula));
	return (result);
}