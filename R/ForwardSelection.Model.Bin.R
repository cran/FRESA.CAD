ForwardSelection.Model.Bin <-
function(size=100,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,variableList,data,maxTrainModelSize=20,type=c("LM","LOGIT","COX"),timeOutcome="Time",selectionType=c("zIDI","zNRI"),cores=4,randsize = 0)
{
#	    R_CStackLimit = -1;
	type <- match.arg(type)
	seltype <- match.arg(selectionType)
	
#		cat(covariates," <- Covariates\n");


	Outcome<-as.character(Outcome);

	if (type == "COX")
		timeOutcome<-as.character(timeOutcome)
	else
		timeOutcome=".";
	if (is.na(size))
	{
	  stop("Number of variables to be used is not defined\n")
	}


#		cat(timeOutcome," <- Time Outcome\n");
	
	vnames <- as.vector(variableList[,1]);
	mcnt=0;
	i=1;
	if (length(vnames)<size) size = length(vnames)
	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates,"+",covariates[i])
		}
	}

	while ((mcnt==0)&&(i<=size))
	{
		mcnt = mcnt+str_count(vnames[i],"\\*");
		i = i + 1;
	}
	
	if (mcnt>0)
	{
		if (nrow(variableList)<size) size = nrow(variableList)
		if (timeOutcome == ".") 
		{
			frm <- paste(Outcome,"~",acovariates);
		}
		else
		{
			frm <- paste(Outcome,"~",acovariates,"+",timeOutcome);
		}
		for (i in 1:size)
		{
			frm <- paste(frm,"+",vnames[i])
		}
#			cat(frm,"\n")
		modelFrame <- model.frame(formula(frm),data);
	}
	else
	{
		modelFrame <- data;	
	}

#		cat(frm,"\n")
	
	colNames=colnames(modelFrame);
	if (randsize >= 0)
	{
#			cat("Forward\n")
		output<-.Call("ReclassificationFRESAModelCpp",size, fraction, pvalue, loops, covariates, Outcome,as.vector(variableList[,1]), maxTrainModelSize, type, timeOutcome, seltype,data.matrix(modelFrame),colNames,cores);
	}
	else
	{
#			cat("Random Forward\n")
		if (timeOutcome != ".") modelFrame[,timeOutcome] <- runif(nrow(modelFrame));
		output <-.Call("ReclassificationFRESAModelCpp",size, fraction, pvalue, loops, "1", paste("RANDOM",Outcome,sep="") ,as.vector(variableList[,1]), maxTrainModelSize, type, timeOutcome, seltype,data.matrix(modelFrame),colNames,cores);
	}

	random.fraction <- 1.0
	if (randsize<0)
	{
		covcount <- 1;
		mcount <- 0;
		pfind <- 0;
		randsize <- 0;
#			print(output$formula.list);
		for (i in 1:loops)
		{
			plusc = str_count(output$formula.list[i],"\\+");
			if (plusc>=covcount)
			{
				randsize <- randsize + plusc - covcount;
				mcount <- mcount+1;
				pfind <- pfind + 1*(plusc>covcount);
			}
		}
		randsize = (nrow(variableList)/size)*(randsize/loops);
		random.fraction <- pfind/loops;
		cat ("\n Vars:",nrow(variableList),"Size:",size,sprintf(", Fraction= %6.4f,  Average random size = %6.2f, Size:%6.2f",random.fraction,randsize,randsize/pvalue),"\n");
	}
	else
	{
		if (randsize==0) randsize = pvalue*nrow(variableList);
	}
	
	zthr = abs(qnorm(pvalue)); 
	zthrO = abs(qnorm(pvalue*pvalue));
	zthr2 = zthr;
	if (fraction<1) 
	{
		zthr2 = zthr*sqrt(fraction);
	}
	if (zthr2<abs(qnorm(0.1))) zthr2 = abs(qnorm(0.1));



		 baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
		baseForm = paste("Surv(",timeOutcome,",",Outcome,")");
	}

	baseForm = paste(baseForm,"~",acovariates);


		mynames <- output$mynames + 1 
		topvar <- table(mynames);
#		print(mynames)
		if (length(topvar)>1)
		{
			topvar <- topvar[order(-topvar)];
		}
		topvarID <- as.numeric(rownames(topvar));


		
		frm1 <- baseForm;
		frm1 <- paste(frm1,"+");
		frm1 <- paste(frm1,vnames[topvarID[1]]);
		ftmp <- formula(frm1);
#		cat(frm1," <- Start Formula \n")
		bestmodel <- modelFitting(ftmp,data,type,TRUE)
		
		if ( !inherits(bestmodel, "try-error"))
		{
			bestpredict <- predict.fitFRESA(bestmodel,data,'prob');

			vnames_model <- vector();
			model_zmin <- vector();
			varlist <- vector();
			inserted = 1;
			
			vnames_model <- append(vnames_model,vnames[topvarID[1]]);
			varlist <- append(varlist,topvarID[1]);
			model_zmin <- append(model_zmin,NA);
			
			if (length(topvar)>1)
			{
				for ( i in 2:length(topvar))
				{
					if ((topvar[i] > 0) && (inserted < maxTrainModelSize))
					{
						kinserted = 0
						kins = 0 
						frma <- paste(frm1,"+");
						frma <- paste(frma,vnames[topvarID[i]]);
#							cat(frma,":",topvarID[i],"\n");
						ftmp <- formula(frma);
						newmodel <- modelFitting(ftmp,data,type,TRUE)
						if ( !inherits(newmodel, "try-error"))
						{
							iprob <- .Call("improveProbCpp",bestpredict,predict.fitFRESA(newmodel,data,'prob'),data[,Outcome]);
							if (seltype=="zIDI") 
							{
								zmin = iprob$z.idi;
							}
							else
							{
								zmin = iprob$z.nri;
							}
							if (is.numeric(zmin) && !is.na(zmin))
							{ 
								if (zmin>zthr)
								{
									bestpredict <- predict.fitFRESA(newmodel,data,'prob');
									frm1 <- frma;
									vnames_model <- append(vnames_model,vnames[topvarID[i]]);
									model_zmin <- append(model_zmin,zmin);
									varlist <- append(varlist,topvarID[i]);
									inserted = inserted + 1;
									kins = 1;
								}
							}									
						}
					}
				}
			}
		}
		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,data,type,TRUE)

	
	base.Zvalues <- output$Base.values;
	rownames(base.Zvalues) <- vnames[1:nrow(base.Zvalues)];
	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	ranked.var=topvar,
	z.selectionType=model_zmin,
	formula.list=output$formula.list,
	random.formula.size=randsize,
	random.fraction = random.fraction,
	variableList=variableList,
	base.Zvalues=base.Zvalues
	);
#		cat ("Final :",frm1,"\n")
	return (result);
}
	