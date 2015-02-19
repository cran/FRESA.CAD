updateModel <-
function(Outcome,covariates="1",pvalue=c(0.05,0.02),VarFrequencyTable,variableList,data,type=c("LM","LOGIT","COX"), lastTopVariable= 0,timeOutcome="Time",selectionType=c("zIDI","zNRI"),numberOfModels=3,interaction=1,maxTrainModelSize=0)
{
	type <- match.arg(type)
  	seltype <- match.arg(selectionType)

	vnames <- as.vector(variableList[,1]);
	topvarID <- as.numeric(rownames(VarFrequencyTable));
	
	casesample = subset(data,get(Outcome)  == 1);
	controlsample = subset(data,get(Outcome) == 0);

	sizecases = nrow(casesample);
	sizecontrol = nrow(controlsample);
	minsize = min(sizecases,sizecontrol);

	if (maxTrainModelSize == 0)
	{
		maxTrainModelSize = minsize/5;
	}

	
	baseForm = Outcome;
					
	#For Cox  models 
	if (type == "COX")
	{
		baseForm = paste("Surv(",timeOutcome);
		baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	}
		
	vnames_model <- vector();
	model_zmin <- vector();
	topvarCpy <- topvarID

	
	if (lastTopVariable < 1) lastTopVariable = length(VarFrequencyTable);
	if (lastTopVariable > length(VarFrequencyTable)) lastTopVariable = length(VarFrequencyTable);
	firstVar = 1
	formulaList <- vector();
	varlist <- vector();
	inserted = 0
	loops = 0
	bestmodel <- NULL;
	ftmp <- NULL;
	frm1 <- NULL;
	if (lastTopVariable>1)
	{


		for (nMod in 1:numberOfModels)
		{
#			cat(nMod," Model \n")
			loops = 0;
			inserted = 0;
			termsinserted = 0;
			if (firstVar>0)
			{
				varlist <- vector();
				for (pidx in 1:length(pvalue))
				{
					changes = 1;
					cthr = abs(qnorm(pvalue[pidx]));	
					zthrO = abs(qnorm(pvalue[pidx]*pvalue[pidx]));
					if (cthr<abs(qnorm(0.1))) cthr = abs(qnorm(0.1));
					cat ("Testing at :",cthr,"\n");
					while ((inserted <= maxTrainModelSize)&&(loops<2) || ((changes>0) && (loops<100)))
					{
						changes = 0;

						caseSam <- sample(1:sizecases, minsize, replace=TRUE);
						contSam <- sample(1:sizecontrol, minsize, replace=TRUE);
						mysampleCases = casesample[caseSam,];
						mysampleControl = controlsample [contSam,];
						testSampleCases = casesample[-caseSam,];
						testSampleControl = controlsample [-contSam,];
						testSampleCases <- testSampleCases[sample(1:nrow(testSampleCases), minsize, replace=TRUE),]
						testSampleControl <- testSampleControl[sample(1:nrow(testSampleControl), minsize, replace=TRUE),]
						
						mysample = rbind(mysampleCases,mysampleControl);
						testsample = rbind(testSampleCases,testSampleControl);

						if ((loops == 0)&&(topvarID[firstVar]>0))
						{
							frm1 = paste(baseForm,paste(" ~ ",covariates));
							frm1 <- paste(frm1," + ");
							frm1 <- paste(frm1,vnames[topvarID[firstVar]]);
							varlist <- append(varlist,topvarID[firstVar]);
							vnames_model <- append(vnames_model,vnames[topvarID[firstVar]]);
							model_zmin <- append(model_zmin,NA);
							topvarID[firstVar] = 0;
							inserted = 1
							ftmp <- formula(frm1);
							bestmodel <- modelFitting(ftmp,mysample,type,TRUE)
						}
						else
						{									
							ftmp <- formula(frm1);
							bestmodel <- modelFitting(ftmp,mysample,type,TRUE)
							while ( inherits(bestmodel, "try-error"))
							{
								caseSam <- sample(1:sizecases, minsize, replace=TRUE);
								contSam <- sample(1:sizecontrol, minsize, replace=TRUE);
								mysampleCases = casesample[caseSam,];
								mysampleControl = controlsample [contSam,];
								testSampleCases = casesample[-caseSam,];
								testSampleControl = controlsample [-contSam,];
								testSampleCases <- testSampleCases[sample(1:nrow(testSampleCases), minsize, replace=TRUE),]
								testSampleControl <- testSampleControl[sample(1:nrow(testSampleControl), minsize, replace=TRUE),]
								
								mysample = rbind(mysampleCases,mysampleControl);
								testsample = rbind(testSampleCases,testSampleControl);
								bestmodel <- modelFitting(ftmp,mysample,type,TRUE)
							}
							
						}
						if ( !inherits(bestmodel, "try-error"))
						{
							bestpredict <- predictForFresa(bestmodel,mysample,'prob');
							bestpredict_test <- predictForFresa(bestmodel,testsample,'prob');


							firstVar = 1;
							for ( i in 2:lastTopVariable)
							{
								if ((VarFrequencyTable[i]>0) && (topvarID[i]>0) && (inserted <= maxTrainModelSize))
								{

									kinserted = 0
									frma <- paste(frm1," + ");
									frma <-paste(frma,vnames[topvarID[i]]);
									ftmp <- formula(frma);
									newmodel <- modelFitting(ftmp,mysample,type,TRUE)
									if ( !inherits(newmodel, "try-error"))
									{
										iprob_t <- improveProb(bestpredict,predictForFresa(newmodel,mysample,'prob'),mysample[,Outcome]);
										iprob <- improveProb(bestpredict_test,predictForFresa(newmodel,testsample,'prob'),testsample[,Outcome]);
										if (seltype=="zIDI") 
										{
											zmin = min(iprob$z.idi,iprob_t$z.idi);
										}
										else
										{
											zmin = min(iprob$z.nri,iprob_t$z.nri);
										}
										if (is.numeric(zmin) && !is.na(zmin) && (zmin>cthr))
										{
											bestpredict <-predictForFresa(newmodel,mysample,'prob');
											bestpredict_test <-predictForFresa(newmodel,testsample,'prob');

											frm1 <- frma;
											vnames_model <- append(vnames_model,vnames[topvarID[i]]);
											model_zmin <- append(model_zmin,zmin);
											varlist <- append(varlist,topvarCpy[i]);
											changes = changes + 1;
											inserted = inserted + 1;
											termsinserted = termsinserted + 1;
											topvarID[i] = 0
										}
										if (is.numeric(zmin) && !is.na(zmin) && (zmin<=cthr))
										{
											if (firstVar == 1) firstVar = i;
										}
										if (interaction == 2)
										{
											chkin <- (topvarID[i] > 0);
											for (nlist in 1:inserted)
											{
												if (termsinserted<=maxTrainModelSize)
												{
													if (topvarID[i] == 0)
													{
														frma <- paste(frm1," + I(",vnames[varlist[nlist]],"*",vnames[topvarCpy[i]],")")
														zthrOl = cthr;
													}
													else
													{
														frma <- paste(frm1," + ",vnames[topvarCpy[i]]," + I(",vnames[varlist[nlist]],"*",vnames[topvarCpy[i]],")")
														zthrOl = zthrO;
													}
													ftmp <- formula(frma);
													newmodel <- modelFitting(ftmp,mysample,type,TRUE)
													if ( !inherits(newmodel, "try-error"))
													{
														iprob_t <- improveProb(bestpredict,predictForFresa(newmodel,mysample,'prob'),mysample[,Outcome]);
														iprob <- improveProb(bestpredict_test,predictForFresa(newmodel,testsample,'prob'),testsample[,Outcome]);
														


														if (seltype=="zIDI") 
														{
															zmin = min(iprob$z.idi,iprob_t$z.idi);
														}
														else
														{
															zmin = min(iprob$z.nri,iprob_t$z.nri);
														}
														if (is.numeric(zmin) && !is.na(zmin) && (zmin>zthrOl))
														{
															bestpredict <- predictForFresa(newmodel,mysample,'prob');
															bestpredict_test <-predictForFresa(newmodel,testsample,'prob');
															
															frm1 <- frma;
															vnames_model <- append(vnames_model,vnames[topvarCpy[i]]);
															model_zmin <- append(model_zmin,zmin);
															changes = changes + 1;
															kinserted = kinserted + 1;
															termsinserted = termsinserted + 1;
															topvarID[i] = 0
														}
													}
												}
											}
											if ((kinserted > 0) && chkin ) 
											{
												varlist <- append(varlist,topvarCpy[i]);
												inserted = inserted + 1;
												if (firstVar == 1) firstVar = i;
											}
										}
									}
								}
							}
#							cat (loops," Form: ",frm1,"\n")
							loops = loops+1;
						}
					}
					cat (frm1,"\n")
				}
				formulaList <- append(formulaList,frm1)
				ftmp <- formula(frm1);
			}
		}
	}
	if (length(formulaList)==0)
	{
		frm1 = paste(baseForm," ~ ",covariates," + ",vnames[topvarID[1]]);
		cat ("Update Formula: ",frm1,"\n")
		formulaList <- append(formulaList,frm1)
		ftmp <- formula(formulaList[1]);
		bestmodel <- modelFitting(ftmp,data,type)
	}
	else
	{
#		cat("Top Formula: \n")
		ftmp <- formula(formulaList[1]);
		cat ("Update Formula: ",formulaList[1],"\n")
		bestmodel <- modelFitting(ftmp,data,type)
	}
#	print(summary(bestmodel));

	
  	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	z.selectionType=model_zmin,
	loops=loops,
	formula.list=formulaList);
  
	return (result);
}
