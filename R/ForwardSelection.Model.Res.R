ForwardSelection.Model.Res <-
function(size=100,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,variableList,data,maxTrainModelSize=20,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),timeOutcome="Time",cores = 4,randsize = 0)
{
#	R_CStackLimit = -1;

	if (is.na(size))
	{
		stop("Size: Number of variables to be explored is not defined\n")
	}


	type <- match.arg(type)
	testType <- match.arg(testType)
	Outcome<-as.character(Outcome);

	if (type == "COX")
		timeOutcome<-as.character(timeOutcome)
	else
		timeOutcome=".";

	vnames <- as.vector(variableList[,1]);
	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates,"+",covariates[i])
		}
	}
	mcnt=0;
	i=1;
	if (length(vnames)<size) size = length(vnames)
	while ((mcnt==0)&&(i<=size))
	{
		mcnt = mcnt+str_count(vnames[i],"\\*");
		i = i + 1;
	}

	if ( mcnt>0 )
	{
#		cat(Outcome," :",nrow(data),":",nrow(variableList)," Model Frames\n")
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
			frm <- paste(frm,"+",vnames[i]);
		}
#		cat(frm,"\n")
		modelFrame <- model.frame(formula(frm),data);
	}
	else
	{
		modelFrame <- data;
	}

	colNames=colnames(modelFrame);
	if (randsize >= 0)
	{
		output<-.Call("ForwardResidualModelCpp",size, fraction, pvalue, loops, covariates, Outcome,vnames, maxTrainModelSize, type, timeOutcome, testType,data.matrix(modelFrame),colNames,cores);
	}
	else
	{
		if (timeOutcome!=".") modelFrame[,timeOutcome] <- runif(nrow(modelFrame));
		output<-.Call("ForwardResidualModelCpp",size, fraction, pvalue, loops, "1", paste("RANDOM",Outcome,sep=""),vnames, maxTrainModelSize, type, timeOutcome, testType,data.matrix(modelFrame),colNames,cores);
	}

	random.fraction <- 1.0
	if (randsize<0)
	{
		covcount <- 1;
		mcount=0;
		pfind <- 0;
		randsize <- 0;
#		print(output$formula.list);
		for (i in 1:loops)
		{
		    plusc = str_count(output$formula.list[i],"\\+");
			if (plusc>=covcount)
			{
				randsize <- randsize + plusc - covcount;
				pfind <- pfind + 1*(plusc>covcount);
				mcount <- mcount+1;
			}
		}
		random.fraction <- pfind/loops;
		randsize = (nrow(variableList)/size)*(randsize/loops);
		cat ("\n Vars:",nrow(variableList),"Size:",size,sprintf(", Fraction= %6.3f,  Average random size = %6.2f, Size:%6.2f",random.fraction,randsize,randsize/pvalue),"\n");
	}
	else
	{
		if (randsize==0) randsize = pvalue*nrow(variableList);
	}

	mynames <- output$mynames + 1;
	formula.list <- output$formula.list

	pthr = pvalue;
	pthrO = pvalue*pvalue;
	baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome,",",Outcome,")");
	}

	baseForm = paste(baseForm,"~",acovariates);

	

		pthr2 = 1-pnorm(sqrt(fraction)*abs(qnorm(pthr)));
		if (pthr2>0.1) pthr2 = 0.1;


		topvar <- table(mynames);
		
		frm1 <- baseForm;
		vnames_model <- vector();
		model_ziri <- vector();
		if (length(topvar)>1)
		{
			topvar <- topvar[order(-topvar)];
			topvarID <- as.numeric(rownames(topvar));

			frm1 <- paste(frm1,"+");
			frm1 <- paste(frm1,vnames[topvarID[1]]);
			
			ftmp <- formula(frm1);
			bestmodel <- modelFitting(ftmp,data,type,TRUE)
#			cat(frm1,"b \n")

			bestResiduals <- residualForFRESA(bestmodel,data,Outcome);

			vnames_model <- append(vnames_model,vnames[topvarID[1]]);
			model_ziri <- append(model_ziri,1);
			varlist <- vector();
			varlist <- append(varlist,topvarID[1]);
			inserted = 1
			for ( i in 2:length(topvar))
			{
				if(topvar[i] > 0)
				{
					frma <- paste(frm1,"+");
					frma <- paste(frma,vnames[topvarID[i]]);
	#				cat(frma," b \n");

					
					ftmp <- formula(frma);
					newmodel <- modelFitting(ftmp,data,type,TRUE)
					kins = 0
					if ( !inherits(newmodel, "try-error"))
					{
						iprob <- .Call("improvedResidualsCpp",bestResiduals,residualForFRESA(newmodel,data,Outcome),testType,0);
						piri <- iprob$p.value;
						if (piri<pthr)
						{
							bestResiduals <- residualForFRESA(newmodel,data,Outcome);
							frm1 <- paste(frm1,"+");
							frm1 <- paste(frm1,vnames[topvarID[i]]);
							varlist <- append(varlist,topvarID[i]);
							vnames_model <- append(vnames_model,vnames[topvarID[i]]);
							model_ziri <- append(model_ziri,abs(qnorm(piri)));
							inserted = inserted + 1;
							kins = 1
						}	
					}
				}
			}
		}


		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,data,type,TRUE);

	#	cat(frm1," Final \n");
	base.Zvalues <- -1.0*qnorm(as.vector(output$Base.values));
	names(base.Zvalues) <- vnames[1:size];

	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	ranked.var=topvar,
	z.NeRIs=model_ziri,
	formula.list=formula.list,
	random.formula.size=randsize,
	random.fraction = random.fraction,
	variableList=variableList,
	base.Zvalues=base.Zvalues
	);
	
	return (result);
}
