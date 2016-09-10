bootstrapVarElimination_Bin <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops=250,fraction=1.00,print=TRUE,plots=TRUE,adjsize=1) 
{
  	seltype <- match.arg(selectionType)



#	boot.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops,fraction,idiCV=NULL) 
	boot.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops,fraction) 
	{
		seltype <- match.arg(selectionType)
		type <- match.arg(type);
	  
		varsList <- as.list(attr(terms(object),"variables"))
		
		quartileValue = 0.49; #
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
			modelFrame <- model.frame(ftmp,data);		
			mfpos=1;
			if (type=='COX') mfpos=2; 

		
#			cat ("Formula ",frm1,"\n")
			idiCV <- bootstrapValidation_Bin(fraction,loops,ftmp,Outcome,data,type,plots=plots)
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
#				cat(as.character(varsList[i])," Min thr: ",minlcl," : ");
				{
					zcorrel <- cor.test(modelFrame[,mfpos], modelFrame[,as.character(varsList[i])], method = "spearman",na.action=na.omit)$p.value
					zcorrel <- -(qnorm(zcorrel))
#					print(idiCV$z.IDIs[,idlist]);
					if (seltype=="zIDI")
					{
						ci <- as.vector(quantile(idiCV$z.IDIs[,idlist], probs = c(quartileValue, 0.5, 1-quartileValue), na.rm = TRUE,names = FALSE, type = 7));
						ci2 <- median(idiCV$test.z.IDIs[,idlist], na.rm = TRUE);
					}
					else
					{
						ci <- as.vector(quantile(idiCV$z.NRIs[,idlist], probs = c(quartileValue, 0.5, 1-quartileValue), na.rm = TRUE,names = FALSE, type = 7));
						ci2 <- median(idiCV$test.z.NRIs[,idlist], na.rm = TRUE);
					}
#					cat(zcorrel,"Ztrain :",ci[2],"Ztest :",ci2[2],"\n");
					if (any(is.na(c(zcorrel,ci[2],ci2)))) 
					{
						who=i;
					}
					else
					{
						if ((zcorrel>ci[2]) && (ci2>0))
						{
							cmin=zcorrel;
						}
						else
						{
							cmin = min(ci[1],ci2);
						}
						if (cmin <= minlcl) 
						{
							minlcl = cmin;
							who = i;
						}
					}
				}
				idlist=idlist+1;
			}
			if ((length(varsList) == startSearch) && (who == startSearch)) 
			{
#				who = -1;
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

		idiCV <- bootstrapValidation_Bin(fraction,loops,ftmp,Outcome,data,type,plots=plots)
		
		result <- list(Removed=removeID,BootModelAUC=idiCV$blind.ROCAUC$auc,backfrm=frm1);

		return (result)
	}

	bkobj <- NULL;
	
	maxROCAUC <- 1.0;
	stopROCACU <- 0.90;	# threshold for stop for decrease in AUC
	if (adjsize>1)
	{
		bkobj <- bootstrapVarElimination_Bin(object,pvalue,Outcome,data,startOffset,type,selectionType,loops,fraction,print,plots,adjsize=1); 
		object <- modelFitting(bkobj$back.formula,data,type);
		adjsize = floor(adjsize);
		maxROCAUC <- bkobj$bootCV$blind.ROCAUC$auc;
		stopROCACU <- 0.90;	# threshold for stop for decrease in AUC
#		cat("Feature size:",adjsize,"\n");
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
	if (length(varsList) >= startlist)
	{
		for ( i in startlist:length(varsList))
		{
			frm1 <- paste(frm1,paste(" + ",varsList[i]));
		}
	}
	beforeFSCmodel.formula <- frm1;
	model.formula <- frm1;
	
	
	beforeFSCmodel <- object;

	wts <- rep(1,length(beforeFSCmodel$coefficients));
	names(wts) <- names(beforeFSCmodel$coefficients);
#	cat ("Start AUC :",startAUC,"\n")
	changed <- 0;
	while ((changes>0) && (loopsAux<100))
	{
		p.elimin <- pvalue;
		if (adjsize>1)
		{
			modsize <- length(as.list(attr(terms(model),"term.labels")));	
			if (modsize<1) modsize=1;
			qvalue <- 2*pvalue;
			if (qvalue < 0.1) qvalue=0.1 # lests keep the minimum q-value to 0.1
			p.elimin <- min(pvalue,modsize*qvalue/adjsize) # # BH alpha  the elimination p-value
		}

		bk <- boot.var.IDISelection(model,p.elimin,Outcome=myOutcome,data=mydataFrame,startOffset,type,seltype,loops,fraction);
#		cat("Used p :",p.elimin,"Formula <- ", bk$backfrm,"\n");
		nmodel = modelFitting(formula(bk$backfrm),data,type);

		if (!inherits(nmodel, "try-error"))
		{
			changes = as.integer(bk$Removed);
			weight <- 1.0;
#			cat ("AUC :",bk$BootModelAUC,"\n")
			if (!is.null(bk$BootModelAUC))
			{			
				if (bk$BootModelAUC < stopROCACU*maxROCAUC)
				{
					changes = 0;
				}
				if (bk$BootModelAUC <= maxROCAUC)
				{
					weight <- 1.0;
				}
				else
				{
					maxROCAUC <- bk$BootModelAUC;
					weight <- 0.0;
				}
			}
			if (changes>0)
			{
				loopsAux = loopsAux + 1;
				changes2<- as.character(as.list(attr(terms(model),"variables")))[which(!(as.character(as.list(attr(terms(model),"variables")))%in%as.character(as.list(attr(terms(nmodel),"variables")))))]
				model <- nmodel;
				model.formula <- bk$backfrm;

				if (length(changes2)>1)
				{
					changes2<-changes2[2]
				}
				if (adjsize>1)
				{
					changed <- 1;
					if ((length(beforeFSCmodel$coefficients)>0)&&(length(model$coefficients)>0))
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
								beforeFSCmodel$coefficients[i] <- weight*beforeFSCmodel$coefficients[i]; # it will average with zero
								wts[i] = wts[i]*weight;
							}
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