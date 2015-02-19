updateNeRIModel <-
function(Outcome,covariates="1",pvalue=c(0.05,0.02),VarFrequencyTable,variableList,data,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent"), lastTopVariable= 0,timeOutcome="Time",interaction=1,maxTrainModelSize=0)
{
	type <- match.arg(type)
  
	vnames <- as.vector(variableList[,1]);
	topvarID <- as.numeric(rownames(VarFrequencyTable));
	vnames_model <- vector();
	nsize <- nrow(data)

	if (maxTrainModelSize == 0)
	{
		maxTrainModelSize = nsize/5;
	}

	nsize <- nrow(data)
	
	baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome);
	  baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	}
	varlist <- vector();

	frm1 = paste(baseForm,paste(" ~ ",covariates));
	frm1 <- paste(frm1," + ");
	frm1 <- paste(frm1,vnames[topvarID[1]]);
	
	ftmp <- formula(frm1);
	varlist <- append(varlist,topvarID[1])
	bestmodel <- modelFitting(ftmp,data,type,TRUE)
	startIndex = 2;
	if ( inherits(bestmodel, "try-error"))
	{
		frm1 <- paste(frm1," + ",vnames[topvarID[2]]);
		varlist <- append(varlist,topvarID[2])
		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,data,type,TRUE);
		startIndex = 3;
		topvarID[2]=0;
	}
	else
	{
		topvarID[1]=0;
	}
#	cat("Update Formula: ",frm1,"\n")
#	print(summary(bestmodel))
	
	bestResiduals <- residualForNeRIs(bestmodel,data,Outcome);

	model_ziri <- vector();

	loops = 0;
	changes = 1;
	if (lastTopVariable < 1) lastTopVariable = length(VarFrequencyTable);
	if (lastTopVariable > length(VarFrequencyTable)) lastTopVariable = length(VarFrequencyTable);
	inserted = 1
	kins=1
	cpyformula <- frm1;
	termsinserted = 1;
	for (pval in 1:length(pvalue))
	{
		cthr = pvalue[pval];
		pthrO = cthr*cthr;
		cat ("Update at:",cthr,"\n");
		while ((termsinserted <= maxTrainModelSize)&&((loops<5) ||((changes>0) && (loops<100))))
		{
			changes = 0;

			theTrainSamples <- sample(1:nsize, nsize, replace=TRUE);
			myTrainSample <- data[theTrainSamples,];
			myTestSample <- data[-theTrainSamples,];
			myTestSample <- myTestSample[sample(1:nrow(myTestSample), nsize, replace=TRUE),];
			
#			myTrainSample <- data[sample(1:nsize, nsize, replace=TRUE),]
#			myTestSample <- data[sample(1:nsize, nsize, replace=TRUE),];


			ftmp <- formula(frm1);
			bestmodel <- modelFitting(ftmp,myTrainSample,type,TRUE)
			if ((loops == 0)&&(inherits(bestmodel, "try-error")))
			{
				frm1 <- paste(frm1," + ",vnames[topvarID[startIndex]]);
				cat("Update Formula 1: ",frm1,"\n")
				varlist <- append(varlist,topvarID[startIndex]);
				VarFrequencyTable[startIndex]=0;
				ftmp <- formula(frm1);
				inserted = inserted + 1;
				bestmodel <- modelFitting(ftmp,myTrainSample,type,TRUE);
				startIndex = startIndex + 1;
			}
			
			while (inherits(bestmodel, "try-error"))
			{
				frm1 <- cpyformula;
				ftmp <- formula(frm1);
				cat("Update Formula 2: ",frm1,"\n")
				theTrainSamples <- sample(1:nsize, nsize, replace=TRUE);
				myTrainSample <- data[theTrainSamples,];
				myTestSample <- data[-theTrainSamples,];
				myTestSample <- myTestSample[sample(1:nrow(myTestSample), nsize, replace=TRUE),];
				bestmodel <- modelFitting(ftmp,myTrainSample,type,TRUE)
			}
			cpyformula <- frm1;

			bestResiduals <- residualForNeRIs(bestmodel,myTrainSample,Outcome);
			bestTestResiduals <- residualForNeRIs(bestmodel,myTestSample,Outcome);

			for ( i in startIndex:lastTopVariable)
			{
	#				cat(vnames[topvarID[i]],"-",VarFrequencyTable[i],"\n")
				if ((VarFrequencyTable[i]>0) && (topvarID[i]>0) && (termsinserted <= maxTrainModelSize))
				{
					frma <- paste(frm1," + ");
					frma <-paste(frma,vnames[topvarID[i]]);
					ftmp <- formula(frma);
					newmodel <- modelFitting(ftmp,myTrainSample,type,TRUE)
					if ( !inherits(newmodel, "try-error"))
					{
						iprob <- improvedResiduals(bestResiduals,residualForNeRIs(newmodel,myTrainSample,Outcome),testType);
						iprob_t <- improvedResiduals(bestTestResiduals,residualForNeRIs(newmodel,myTestSample,Outcome),testType);
						piri <- max(iprob$p.value,iprob_t$p.value);
						if (is.numeric(piri) && !is.na(piri) && (piri<cthr))
						{
							bestResiduals <- residualForNeRIs(newmodel,myTrainSample,Outcome);
							bestTestResiduals <- residualForNeRIs(newmodel,myTestSample,Outcome);
							frm1 <- paste(frm1," + ",vnames[topvarID[i]]);
							vnames_model <- append(vnames_model,vnames[topvarID[i]]);
							varlist <- append(varlist,topvarID[i]);
							model_ziri <- append(model_ziri,abs(qnorm(piri)));
							changes = changes + 1;
							inserted = inserted + 1;
							termsinserted = termsinserted + 1;
							kins=1
							VarFrequencyTable[i]=0;
						}	
						if (interaction == 2)
						{
							for (nlist in 1:inserted)
							{
								if (termsinserted<=maxTrainModelSize)
								{
									if (kins==1)
									{
										pthrOl=cthr;
										frma <- paste(frm1," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
									}
									else
									{
										frma <- paste(frm1," + ",vnames[topvarID[i]]," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
										pthrOl=pthrO;
									}
									ftmp <- formula(frma);
									newmodel <- modelFitting(ftmp,myTrainSample,type,TRUE)
									if ( !inherits(newmodel, "try-error"))
									{
										iprob <- improvedResiduals(bestResiduals,residualForNeRIs(newmodel,myTrainSample,Outcome),testType);
										iprob_t <- improvedResiduals(bestTestResiduals,residualForNeRIs(newmodel,myTestSample,Outcome),testType);
										piri <- max(iprob$p.value,iprob_t$p.value);
										if (is.numeric(piri) && !is.na(piri) && (piri<pthrOl))
										{
											bestResiduals <- residualForNeRIs(newmodel,myTrainSample,Outcome);
											bestTestResiduals <- residualForNeRIs(newmodel,myTestSample,Outcome);
											frm1 <- frma;
											vnames_model <- append(vnames_model,vnames[topvarID[i]]);
											model_ziri <- append(model_ziri,abs(qnorm(piri)));
											if (kins == 0)
											{
												varlist <- append(varlist,topvarID[i]);
												inserted = inserted + 1;
											}
											termsinserted = termsinserted + 1;
											kins =1
											VarFrequencyTable[i]=0;
											changes = changes + 1;
										}
									}
								}
							}
						}
					}
					kins=0
				}
			}
			loops = loops+1;
		}
		cat(frm1,"\n");
	}

	ftmp <- formula(frm1);
	bestmodel <- modelFitting(ftmp,data,type)
	cat("Update Formula: ",frm1,"\n");
#	print(summary(bestmodel));
	
  	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	z.NeRI=model_ziri,
	loops=loops);
  
	return (result);
}
