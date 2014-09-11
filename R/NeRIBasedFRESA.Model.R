NeRIBasedFRESA.Model <-
function(size=100,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,variableList,dataframe,maxTrainModelSize=10,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),timeOutcome="Time",loop.threshold=20,interaction = 1)
{

	if (is.na(size))
	{
		stop("Size: Number of variables to be explored is not defined\n")
	}
	  type <- match.arg(type)
	  
	 pthr = pvalue;
	 pthrO = pvalue*pvalue;
	 baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome);
	  baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	}


	baseForm = paste(baseForm,paste(" ~ ",covariates));

	pthr2 = 1-pnorm(sqrt(fraction)*abs(qnorm(pthr)));
	if (pthr2>0.1) pthr2 = 0.1;


	sizecases = nrow(dataframe);
	 
	totSamples = fraction*sizecases;
	 
	mynames <- vector();

	vnames <- as.vector(variableList[,1]);
	vnameslen <- length(vnames)
	if (size > vnameslen) size = vnameslen;

	randomorderNames <- sample(1:size, size, replace=FALSE)
	vnames <- vnames[randomorderNames]; #randomize the order of input
	ovnames <- vnames;


	
	formulas <- vector();

	for ( doOver in 1:loops)
	{ 

		vnames <- ovnames;
		if (loops>1)
		{


			theSamples <- sample(1:sizecases, totSamples, replace=TRUE);
			mysample <- dataframe[theSamples,];

			myTestsample <- dataframe[-theSamples,];
			myTestsample <- myTestsample[sample(1:nrow(myTestsample), totSamples, replace=TRUE),];

		}
		else
		{
			mysample <- dataframe;
			myTestsample <- dataframe;
		}
			
		inname = "Inserted";


		frm1 <- baseForm;
		# frm1 <- paste(frm1," + ");
		# frm1 <- paste(frm1,vnames[size]);
		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,mysample,type)
		bestResiduals <- residualForNeRIs(bestmodel,myTestsample,Outcome);
		bestTrainResiduals <- residualForNeRIs(bestmodel,mysample,Outcome);
#		frm1 <- baseForm;

		changes = 1;
		minpiri = 1;
		jmax = -1;
		inserted = 0;
		varlist <- vector();

		while (changes>0)
		{
			changes = 0;
			minpiri = 1;
			jmax = 0;
			nlistidx = 0;
			iprob <- improvedResiduals(bestResiduals,bestResiduals,testType);
			iprobTrain <- improvedResiduals(bestTrainResiduals,bestTrainResiduals,testType);            
			for (j in 1:size)
			{

				hfrec = 1;
				if (doOver > loop.threshold) 
				{
					frec <- topvar[toString(j)];
					if (!is.na(frec))
					{
					    if ((frec/doOver) < 1.0/(2*loop.threshold+1.0)) 
						{ 
							hfrec = 0;
						}
					}
					else
					{
						hfrec = 0;
					}
				}

				if ((vnames[j] != inname) && (hfrec>0))
				{
					gfrm1 <- paste(frm1," + ",vnames[j]);
					ftmp <- formula(gfrm1);
					newmodel <- modelFitting(ftmp,mysample,type)

					if ( !inherits(newmodel, "try-error"))
					{
						izecoef = length(newmodel$coef);
						if (!is.na(newmodel$coef[izecoef]))
						{
							testResiduals <- residualForNeRIs(newmodel,newdata=myTestsample,Outcome)
							trainResiduals <- residualForNeRIs(newmodel,newdata=mysample,Outcome)
							iprob <- improvedResiduals(bestResiduals,testResiduals,testType);            
							iprobTrain <- improvedResiduals(bestTrainResiduals,trainResiduals,testType);            
							if ( !is.na(iprob$p.value) &&  !is.na(iprobTrain$p.value) )
							{
								piri <- max(iprob$p.value,iprobTrain$p.value);
								if (piri  < minpiri)
								{
									 jmax = j;
									 minpiri = piri;
									 nlistidx = 0;
								}
							}
						}	 
					}
					if ((interaction == 2) && (inserted>0))
					{
						for (nlist in 1:inserted)
						{
							if (jmax != j) 
							{
								gfrm1 <- paste(frm1," + ",ovnames[j]," + I(",ovnames[varlist[nlist]],"*",ovnames[j],")")
							}
							else
							{
								gfrm1 <- paste(frm1," + I(",ovnames[varlist[nlist]],"*",ovnames[j],")")
							}
#							cat("nlist:",vnames[varlist[nlist]]," topVar:",vnames[topvarID[i]]," Form:",gfrm1,"\n");
							ftmp <- formula(gfrm1);
							newmodel <- modelFitting(ftmp,mysample,type)
							
							if ( !inherits(newmodel, "try-error"))
							{
								izecoef = length(newmodel$coef);
								if (!is.na(newmodel$coef[izecoef]))
								{
									if (jmax != j)
									{
										iprob <- improvedResiduals(bestResiduals,residualForNeRIs(newmodel,newdata=myTestsample,Outcome),testType);            
										iprobTrain <- improvedResiduals(bestTrainResiduals,residualForNeRIs(newmodel,newdata=mysample,Outcome),testType);            
									}
									else
									{
										iprob <- improvedResiduals(testResiduals,residualForNeRIs(newmodel,newdata=myTestsample,Outcome),testType);            
										iprobTrain <- improvedResiduals(trainResiduals,residualForNeRIs(newmodel,newdata=mysample,Outcome),testType);            
									}
									if ( !is.na(iprob$p.value) &&  !is.na(iprobTrain$p.value) )
									{
										piri <- max(iprob$p.value,iprobTrain$p.value);
										if (jmax != j)
										{
											piri = sqrt(piri);
										}
										if (piri  < minpiri)
										{
											 minpiri = piri;
											 nlistidx = nlist;
										}
									}								 
								}
							}
						}
						if (nlistidx>0)
						{
							jmax = j;
						}
					}
				}
			}
			if ((jmax > 0) && (minpiri < pthr2) && (vnames[jmax] != inname) && (inserted<=maxTrainModelSize))   
			{
				if ((interaction == 2) && (inserted>0) && (nlistidx>0))
				{
					gfrm1 <- paste(frm1,"+",ovnames[jmax]," + I(",ovnames[varlist[nlistidx]],"*",ovnames[jmax],")")
					mynames <- append(mynames,varlist[nlistidx]);
				}
				else
				{
					gfrm1 <- paste(frm1,"+",ovnames[jmax]);
				}
				ftmp <- formula(gfrm1);
				bestmodel <- modelFitting(ftmp,mysample,type)
				bestResiduals <- residualForNeRIs(bestmodel,newdata=myTestsample,Outcome);
				bestTrainResiduals <- residualForNeRIs(bestmodel,newdata=mysample,Outcome);
				frm1 <- gfrm1;
#				cat (frm1,"\n");

				changes <- changes + 1;
				mynames <- append(mynames,jmax);

				varlist <- append(varlist,jmax);
				vnames[jmax] = inname;
				jmax = 0;
				minpiri = 1;
				inserted = inserted + 1;
			}
		}
		cat (frm1,"\n");
		formulas <- append(formulas,frm1);
		topvar <- table(mynames);
		if (length(topvar)>1)
		{
			topvar <- topvar[order(-topvar)];
			barplot(topvar);
			titname <- paste ( "Var Frequency ",doOver);
			title(main=titname);
		}
	}

	vnames <- as.vector(variableList[,1]);

	topvar <- table(mynames);
	
	frm1 <- baseForm;
	vnames_model <- vector();
	model_ziri <- vector();
	if (length(topvar)>1)
	{
		topvar <- topvar[order(-topvar)];
		topvarID <- as.numeric(rownames(topvar));
		tovarnames <- topvarID;
		for (i in 1:length(topvar)) {tovarnames[i]=randomorderNames[topvarID[i]]; }
		rownames(topvar) <- tovarnames;
		topvarID <- tovarnames;

		barplot(topvar);
		titname <- paste ( "Var Frequency Completed");
		title(main=titname);
		frm1 <- paste(frm1," + ");
		frm1 <- paste(frm1,vnames[topvarID[1]]);
		cat(frm1," \n");

		
		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,dataframe,type)


		bestResiduals <- residualForNeRIs(bestmodel,newdata=dataframe,Outcome);

		vnames_model <- append(vnames_model,vnames[topvarID[1]]);
		model_ziri <- append(model_ziri,1);
		varlist <- vector();
		varlist <- append(varlist,topvarID[1]);
		inserted = 1
		for ( i in 2:length(topvar))
		{
			frma <- paste(frm1," + ");
			frma <- paste(frma,vnames[topvarID[i]]);

			
			ftmp <- formula(frma);
			newmodel <- modelFitting(ftmp,dataframe,type)
			kins = 0
			if ( !inherits(newmodel, "try-error"))
			{
				iprob <- improvedResiduals(bestResiduals,residualForNeRIs(newmodel,newdata=dataframe,Outcome),testType);
				piri <- iprob$p.value;
				if (piri<pthr)
				{
					bestResiduals <- residualForNeRIs(newmodel,newdata=dataframe,Outcome);
					frm1 <- paste(frm1," + ");
					frm1 <- paste(frm1,vnames[topvarID[i]]);
					varlist <- append(varlist,topvarID[i]);
					vnames_model <- append(vnames_model,vnames[topvarID[i]]);
					model_ziri <- append(model_ziri,abs(qnorm(piri)));
#					print(summary(newmodel));
					inserted = inserted + 1;
					kins = 1
				}	
				if (interaction == 2)
				{
					for (nlist in 1:inserted)
					{
						if (kins==1)
						{
							pthrOl=pthr;
							frma <- paste(frm1," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
						}
						else
						{
							frma <- paste(frm1," + ",vnames[topvarID[i]]," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
							pthrOl=pthrO;
						}
						ftmp <- formula(frma);
						newmodel <- modelFitting(ftmp,dataframe,type)
						iprob <- improvedResiduals(bestResiduals,residualForNeRIs(newmodel,newdata=dataframe,Outcome),testType);
						piri <- iprob$p.value;
						if (is.numeric(piri) && !is.na(piri) && (piri<pthrOl))
						{
							bestResiduals <- residualForNeRIs(newmodel,newdata=dataframe,Outcome);
							frm1 <- frma;
							vnames_model <- append(vnames_model,vnames[topvarID[i]]);
							model_ziri <- append(model_ziri,abs(qnorm(piri)));
							if (kins == 0)
							{
								varlist <- append(varlist,topvarID[i]);
								inserted = inserted + 1;
							}
							kins =1
						}	
					}
				}
			}
		}
	}
	ftmp <- formula(frm1);
	bestmodel <- modelFitting(ftmp,dataframe,type)
	print(summary(bestmodel));
	
	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	ranked.var=topvar,
	z.iri=model_ziri,
	formulas.list=formulas);
	return (result);
}
