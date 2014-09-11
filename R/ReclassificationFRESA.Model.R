ReclassificationFRESA.Model <-
function(size=100,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,variableList,dataframe,maxTrainModelSize=10,type=c("LM","LOGIT","COX"),timeOutcome="Time",selectionType=c("zIDI","zNRI"),loop.threshold=20,interaction=1)
{

	if (is.na(size))
	{
		stop("Number of variables to be used is not defined\n")
	}

	  type <- match.arg(type)
	  seltype <- match.arg(selectionType)
	  
		zthr = abs(qnorm(pvalue)); 
		zthrO = abs(qnorm(pvalue*pvalue));
		zadjs = zthr/zthrO;

		baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome);
	  baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	}

	baseForm = paste(baseForm,paste(" ~ ",covariates));

	casesample = subset(dataframe,get(Outcome)  == 1);
	controlsample = subset(dataframe,get(Outcome) == 0);

	sizecases = nrow(casesample);
	sizecontrol = nrow(controlsample);
	minsize = min(sizecases,sizecontrol);
	 
	totSamples <- as.integer(fraction*minsize+0.49);
	zthr2 = zthr;
	if (fraction<1) 
	{
		zthr2 = zthr*sqrt(fraction);
	}
	if (zthr2<abs(qnorm(0.1))) zthr2 = abs(qnorm(0.1));
	
	 
	mynames <- vector();

	vnames <- as.vector(variableList[,1]);
	ovnames <- as.vector(variableList[,1]);
	if (size > length(vnames)) size = length(vnames);
	formulas <- vector();

	randomorderNames <- sample(1:size, size, replace=FALSE)
	vnames <- vnames[randomorderNames]; #randomize the order of input
	ovnames <- vnames;
	
	
	
	for ( doOver in 1:loops)
	{ 

		vnames <- ovnames;
		
		if (loops > 1)
		{
			samCases <- sample(1:sizecases, totSamples, replace=TRUE)
			samControl <- sample(1:sizecontrol, totSamples, replace=TRUE)
			mysampleCases = casesample[samCases,];
			mysampleControl = controlsample [samControl,];


			mysample = rbind(mysampleCases,mysampleControl);

			myIndCases = casesample[-samCases,];
			myIndControl = controlsample [-samControl,];

			myIndCases <- myIndCases[sample(1:nrow(myIndCases), totSamples, replace=TRUE),];
			myIndControl <- myIndControl[sample(1:nrow(myIndControl), totSamples, replace=TRUE),];

			
			myTestsample = rbind(myIndCases,myIndControl);
		}
		else
		{
			mysample = dataframe;
			myTestsample = dataframe;
		}

		inname = "Inserted";


		frm1 <- baseForm;
		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,mysample,type)
		if ( !inherits(bestmodel, "try-error"))
		{
			bestpredict_train <- predictForFresa(bestmodel,mysample,type = 'prob');
			bestpredict <- predictForFresa(bestmodel,myTestsample,type = 'prob');

			changes = 1;
			maxrec = 0;
			jmax = -1;
			inserted = 0;
			
			varlist <- vector();

			while (changes>0)
			{
				changes = 0;
				maxrec = 0;
				jmax = 0;
				nlistidx = 0;
				iprob <- improveProb(bestpredict,bestpredict,myTestsample[,Outcome]);
				iprob_t <- iprob;            
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
						gfrm1 <- paste(frm1," + ");
						gfrm1 <- paste(gfrm1,vnames[j]);
						ftmp <- formula(gfrm1);
						newmodel <- modelFitting(ftmp,mysample,type)
						singleTrainPredict <- NULL
						singleTestPredict <- NULL
						if ( !inherits(newmodel, "try-error"))
						{
							izecoef = length(newmodel$coef);
							if (!is.na(newmodel$coef[izecoef]))
							{
								singleTrainPredict <- predictForFresa(newmodel,newdata=mysample,type = 'prob')
								singleTestPredict <- predictForFresa(newmodel,newdata=myTestsample,type = 'prob')
								iprob_t <- improveProb(bestpredict_train,singleTrainPredict,mysample[,Outcome]);            
								iprob <- improveProb(bestpredict,singleTestPredict,myTestsample[,Outcome]);            
								if (seltype=="zIDI") 
								{
									zmin = min(iprob$z.idi,iprob_t$z.idi);
								}
								else
								{
									zmin = min(iprob$z.nri,iprob_t$z.nri);
								}
								if ( !is.na(zmin) && is.numeric(zmin) )
								{
									if (zmin  > maxrec)
									{
										 jmax = j;
										 maxrec = zmin;
										 nlistidx = 0;
									}
								}
							}
						}
						# else
						# {
							# cat("Fit Error\n");
						# }
							 
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
								ftmp <- formula(gfrm1);
								newmodel <- modelFitting(ftmp,mysample,type)
								
								if ( !inherits(newmodel, "try-error"))
								{
									izecoef = length(newmodel$coef);
									if (!is.na(newmodel$coef[izecoef]))
									{
										if (jmax != j)
										{
											iprob_t <- improveProb(bestpredict_train,predictForFresa(newmodel,newdata=mysample,type = 'prob'),mysample[,Outcome]);            
											iprob <- improveProb(bestpredict,predictForFresa(newmodel,newdata=myTestsample,type = 'prob'),myTestsample[,Outcome]);            
										}
										else
										{
											iprob_t <- improveProb(singleTrainPredict,predictForFresa(newmodel,newdata=mysample,type = 'prob'),mysample[,Outcome]);            
											iprob <- improveProb(singleTestPredict,predictForFresa(newmodel,newdata=myTestsample,type = 'prob'),myTestsample[,Outcome]);            
										}
										if ( (is.numeric(iprob$z.idi)) && !is.na(iprob$z.idi) && (is.numeric(iprob_t$z.idi)) && !is.na(iprob_t$z.idi))
										{
											if (seltype=="zIDI") 
											{
												zmin = min(iprob$z.idi,iprob_t$z.idi);
											}
											else
											{
												zmin = min(iprob$z.nri,iprob_t$z.nri);
											}
											if (jmax != j)
											{
												zmin = zadjs*zmin;
											}
											if (zmin  > maxrec)
											{
												 nlistidx = nlist;
												 maxrec = zmin;
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
				if ((jmax > 0) && (maxrec > zthr2) && (vnames[jmax] != inname) && (inserted<maxTrainModelSize))   
				{
					gfrm1 <- paste(frm1," + ");
					if ((interaction == 2) && (inserted>0) && (nlistidx>0))
					{
						gfrm1 <- paste(gfrm1,ovnames[jmax]," + I(",ovnames[varlist[nlistidx]],"*",ovnames[jmax],")")
						mynames <- append(mynames,varlist[nlistidx]);
					}
					else
					{
						gfrm1 <- paste(gfrm1,ovnames[jmax]);
					}
					ftmp <- formula(gfrm1);
					bestmodel <- modelFitting(ftmp,mysample,type)
					if ( !inherits(bestmodel, "try-error"))
					{
						bestpredict_train <- predictForFresa(bestmodel,newdata=mysample,type = 'prob');
						bestpredict <- predictForFresa(bestmodel,newdata=myTestsample,type = 'prob');
						frm1 <- gfrm1;
						changes <- changes + 1;
						mynames <- append(mynames,jmax);
						varlist <- append(varlist,jmax);
						jmax = 0;
						maxrec = 0;
						inserted = inserted + 1;
					}
					vnames[jmax] = inname;
				}
			}
			cat(frm1,"\n")

			formulas <- append(formulas,gfrm1);
			topvar <- table(mynames);
			if (length(topvar)>1)
			{
				topvar <- topvar[order(-topvar)];
				barplot(topvar);
				titname <- paste ( "Var Frequency ",doOver);
				title(main=titname);
			}
		}
	}
	if (length(mynames) == 0) mynames <- append(mynames,1);
	vnames <- as.vector(variableList[,1]);
	topvar <- table(mynames);
	if (length(topvar)>1)
	{
		topvar <- topvar[order(-topvar)];
		topvarID <- as.numeric(rownames(topvar));
		tovarnames <- topvarID;
		for (i in 1:length(topvar)) {tovarnames[i]=randomorderNames[topvarID[i]]; }
		rownames(topvar) <- tovarnames;
	}
	topvarID <- as.numeric(rownames(topvar));


	
	frm1 <- baseForm;
	frm1 <- paste(frm1," + ");
	frm1 <- paste(frm1,vnames[topvarID[1]]);
	ftmp <- formula(frm1);
	bestmodel <- modelFitting(ftmp,dataframe,type)
	if ( !inherits(bestmodel, "try-error"))
	{
		bestpredict <- predictForFresa(bestmodel,newdata=dataframe,type = 'prob');

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

				if (doOver > loop.threshold) 
				{
					frec <- topvar[i];
					if (!is.na(frec))
					{
						if ((frec/doOver) < 1.0/(2*loop.threshold+1.0)) 
						{ 
							topvar[i] = 0;
						}
					}
				}

				if ((topvar[i] > 0) && (inserted < maxTrainModelSize))
				{
					kinserted = 0
					kins = 0 
					frma <- paste(frm1," + ");
					frma <- paste(frma,vnames[topvarID[i]]);
					ftmp <- formula(frma);
					newmodel <- modelFitting(ftmp,dataframe,type)
					if ( !inherits(newmodel, "try-error"))
					{
						iprob <- improveProb(bestpredict,predictForFresa(newmodel,newdata=dataframe,type = 'prob'),dataframe[,Outcome]);
						if (seltype=="zIDI") 
						{
							zmin = iprob$z.idi;
						}
						else
						{
							zmin = iprob$z.nri;
						}
						if (is.numeric(zmin) && !is.na(zmin) && (zmin>zthr))
						{
							bestpredict <- predictForFresa(newmodel,newdata=dataframe,type = 'prob');
							frm1 <- frma;
							vnames_model <- append(vnames_model,vnames[topvarID[i]]);
							model_zmin <- append(model_zmin,zmin);
							varlist <- append(varlist,topvarID[i]);
							inserted = inserted + 1;
							kins = 1;
						}	
						if (interaction == 2)
						{
							for (nlist in 1:inserted)
							{
								if (kins==1)
								{
									zthrOl=zthr;
									frma <- paste(frm1," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
								}
								else
								{
									frma <- paste(frm1," + ",vnames[topvarID[i]]," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
									zthrOl=zthrO;
								}
								ftmp <- formula(frma);
								newmodel <- modelFitting(ftmp,dataframe,type)
								if ( !inherits(newmodel, "try-error"))
								{
									iprob <- improveProb(bestpredict,predictForFresa(newmodel,newdata=dataframe,type = 'prob'),dataframe[,Outcome]);
									if (seltype=="zIDI") 
									{
										zmin = iprob$z.idi;
									}
									else
									{
										zmin = iprob$z.nri;
									}
									if (is.numeric(zmin) && !is.na(zmin) && (zmin>zthrOl))
									{
										bestpredict <- predictForFresa(newmodel,newdata=dataframe,type = 'prob');
										frm1 <- frma;
										vnames_model <- append(vnames_model,vnames[topvarID[i]]);
										model_zmin <- append(model_zmin,zmin);
										kinserted = kinserted + 1;
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
			}
		}
	}
	barplot(topvar);
	titname <- paste ( "Var Frequency Completed");
	title(main=titname);
	ftmp <- formula(frm1);
	bestmodel <- modelFitting(ftmp,dataframe,type)

	print(summary(bestmodel));

	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	ranked.var=topvar,
	z.min=model_zmin,
	formula.list=formulas);
	cat ("Final :",frm1,"\n")
	return (result);
}
