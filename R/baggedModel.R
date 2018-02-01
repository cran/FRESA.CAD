baggedModel <-
function(modelFormulas,data,type=c("LM","LOGIT","COX"),Outcome=NULL,timeOutcome=NULL,frequencyThreshold=0.025,univariate=NULL,useFreq=TRUE,n_bootstrap=1)
{
	type <- match.arg(type)

	loops <- length(modelFormulas);
	formulaLoops <- 1
	if (is.numeric(useFreq))
	{
		formulaLoops <- useFreq;
		useFreq <- TRUE;
	}
	observations <- nrow(data);
	avgZvalues <- NULL;
	if (length(modelFormulas)>0)
	{
		iscompletef <- (gregexpr(pattern ='~',modelFormulas[1])[1]>0)
		if (iscompletef && is.null(Outcome))
		{
			varsmod <- all.vars(formula(modelFormulas[1]));
			if (substr(modelFormulas[1],1,5)!="Surv(")
			{
				Outcome <- varsmod[1];
			}
			else
			{
				Outcome <- varsmod[2];
				timeOutcome <- varsmod[1];
			}
		}

		theoutcome <- data[,Outcome];
		varOutcome <- var(theoutcome);
		binoutcome <- (length(table(theoutcome))==2) && (min(theoutcome)==0);
		predtype="linear";
		if (binoutcome) predtype="prob";
		if ( (type=="LM") && (binoutcome==TRUE) )
		{
			data[,Outcome] = 2*theoutcome-1.0;
			predtype="linear";
		}
		
		if (type!="COX")
		{
			baseForm <- paste(Outcome,"~ 1");
		}
		else
		{
			baseForm <- paste("Surv(",timeOutcome,",",Outcome,")~ 1");
		}
		EquTrainSet <- data;
		minTrainSamples <- nrow(data);
		maxTrainSamples = minTrainSamples;
		casesample  <- NULL;
		controlsample <- NULL;
		noequalSets <- FALSE;
		nrowcases <- minTrainSamples
		nrowcontrols <- minTrainSamples
		if (type != "LM")
		{
			casesample = subset(data,get(Outcome)  == 1);
			controlsample = subset(data,get(Outcome) == 0);
			nrowcases <- nrow(casesample);
			nrowcontrols <- nrow(controlsample);
			
			minTrainSamples <- min(c(nrowcases,nrowcontrols));
			maxTrainSamples <- max(c(nrowcases,nrowcontrols));
			noequalSets <- (minTrainSamples < (0.6*maxTrainSamples));
	#		cat(nrowcases,":",nrowcontrols,":",noequalSets,"\n")
	#		cat(minTrainSamples,":",maxTrainSamples,":",noequalSets,"\n")
		}

		
		
		
		#	cat("Bagging\n")
		theterms <- list();
		
		features <- vector();
		for (n in 1:loops)
		{
			if (iscompletef)	
			{
				modelFormulas[n] <- unlist(strsplit(as.character(modelFormulas[n]),"[~]"))[2];
			}
			if (nchar(modelFormulas[n])>0)
			{
				modelFormulas[n] <- paste(baseForm,"+",modelFormulas[n]);
			}
			else
			{
				modelFormulas[n] <- baseForm;
			}
			termList <- str_replace_all(attr(terms(formula(modelFormulas[n])),"term.labels"),":","\\*");
			features <- append(features,termList);
			theterms[[n]] <- termList;
		}
		VarFrequencyTable <- table(features);


	#	print(modelFormulas);
			
		if (length(VarFrequencyTable)>1) 
		{
			oVarFrequencyTable <- VarFrequencyTable;
			if (formulaLoops > 1 )
			{
#				barplot(VarFrequencyTable[1:20]);
				cycles <- loops/formulaLoops;
#				cat("Length:",loops,"loops:",formulaLoops,"Cycles:",cycles,"\n")
				VarFrequencyTable <- 0*VarFrequencyTable;
				vnames <- rownames(VarFrequencyTable);
				for (lo in 1:cycles)
				{
#					cat(lo,"\n")
					features <- vector();
					offs <- (lo-1)*formulaLoops;
					for (m in 1:formulaLoops)
					{
						features <- append(features,theterms[[m+offs]]);
					}
					tmpVarFrequencyTable <- table(features);
					thetnames <- vnames %in% rownames(tmpVarFrequencyTable)
					if (!is.null(thetnames))
					{
						otablenames <- rownames(VarFrequencyTable[thetnames]);
						if (!is.null(otablenames))
						{
							freq1 <- as.numeric(VarFrequencyTable[thetnames]);
							freq2 <- as.numeric(tmpVarFrequencyTable[otablenames]);
							VarFrequencyTable[thetnames] <- pmax(freq1,freq2);
						}
					}
				}
				fmax <- max(oVarFrequencyTable);
				VarFrequencyTable <- VarFrequencyTable+(oVarFrequencyTable-1.0)/fmax;
#				barplot(VarFrequencyTable[1:20]);
			}
			if (!is.null(univariate))
			{
				VarFrequencyTable <- VarFrequencyTable[order(-univariate[rownames(VarFrequencyTable),"ZUni"])];			
			}
			if (useFreq)
			{
				VarFrequencyTable <- VarFrequencyTable[order(-VarFrequencyTable)]
			}
		}
		vnames <- rownames(VarFrequencyTable);
		formulaNetwork <- matrix(0,nrow=length(VarFrequencyTable),ncol = length(VarFrequencyTable))
		dimnames(formulaNetwork) <- list(names(VarFrequencyTable),names(VarFrequencyTable))
		m <- formulaNetwork
		Jaccard.SM <- 0;
		tota <- 0;
		for (n in 1:loops)
		{
			feat <- theterms[[n]];
			lft <- length(feat);
			m[,] <- 0;
			if (lft>0)
			{
				m[feat,feat]=1;
				if (n<loops)
				{
					for (i in (n+1):loops)
					{
						feat2 <- theterms[[i]];
						if (length(feat2) > 0)
						{
							Jaccard.SM = Jaccard.SM+sum(duplicated(c(feat2,feat)))/length(unique(c(feat2,feat)));
							tota = tota + 1;
	#						print(feat2)
	#						print(feat)
	#						cat("Dup:",sum(duplicated(c(feat2,feat)))," U:",length(unique(c(feat2,feat)))," JI:",Jaccard.SM,"\n")
						}
					}
				}
			}
			formulaNetwork <- formulaNetwork+m
		}
		if (tota>1) Jaccard.SM = Jaccard.SM/tota;
		formulaNetwork <- round(formulaNetwork/loops,digits = 3);

	#	cat("Size :",nrow(data)," Features :",length(VarFrequencyTable))
		
		nsize <- nrow(data)
		
		lastTopVariable = length(VarFrequencyTable);
#		if (lastTopVariable >= 2*(nsize-2)) lastTopVariable <- 2*(nsize-2); #The largest model size
		
		frma <- baseForm;
		enterlist <- vector();
		toRemove <- vector();
		bmodelsize <- 1;
		removed <- 0;
		avgsize <- 0;
		coefEvolution <- NULL;
		if (length(vnames)>0)
		{
			thrsfreq <- as.integer(frequencyThreshold*VarFrequencyTable[1]+0.5);
			fistfreq <- VarFrequencyTable[1];
			for ( i in 1:length(vnames))
			{
				if ((vnames[i] != " ") && (vnames[i] != ""))
				{
					enterlist <- append(enterlist,vnames[i]);
					if ((i<=lastTopVariable)&&(VarFrequencyTable[i] > thrsfreq))  # Only features with a given frequency
					{
						if ((fistfreq == loops)&&(VarFrequencyTable[i] > (loops/3)))
						{
							fistfreq <- VarFrequencyTable[i];
							thrsfreq <- as.integer(frequencyThreshold*fistfreq+0.5);
						}
						frma <- paste(frma,"+",vnames[i]);	
						bmodelsize = bmodelsize + 1;
					}
					else
					{
						toRemove <- append(toRemove,paste(" ",vnames[i]," ",sep=""));
						removed = removed+1;
					}
				}
			}
			cat("\nNum. Models:",loops," To Test:",length(vnames)," TopFreq:",fistfreq," Thrf:",thrsfreq," Removed:",removed,"\n")
			model <- modelFitting(frma,data,type=type,fitFRESA=TRUE);
			thevars <- all.vars(formula(frma));
			data <- data[,thevars];
			if (noequalSets)
			{
				casesample <- casesample[,thevars]
				controlsample <- controlsample[,thevars]
				trainCaseses <- casesample;
				trainControls <- controlsample;
		#		print(thevars);
			}
			else
			{
				EquTrainSet <- data;
			}
			
			
		#	print(toRemove);
		#	print(model$coefficients);
		#	cat(frma,"\n");
			if (inherits(model, "try-error"))
			{
#				cat("Fitting Error\n");
				warning(frma," Warning Bagging Fitting error\n")
			}
			msize <- length(model$coefficients)
			basecoef <- abs(model$coefficients)+1e-6;
			names(basecoef) <- names(model$coefficients);
			
			if ((type=="COX")&&(class(model)!="fitFRESA"))
			{
				avgZvalues <- numeric(length(model$coefficients));
				names(avgZvalues) <- names(model$coefficients);
			}
			else
			{
				avgZvalues <- numeric(length(model$coefficients)-1);
				names(avgZvalues) <- names(model$coefficients)[-1];
			}
			addedZvalues <- avgZvalues;
			baggingAnalysis <- list();
			baggingAnalysis$uMS_values <- avgZvalues;
			baggingAnalysis$rMS_values <- avgZvalues;
			baggingAnalysis$NeRI_values <- avgZvalues;
			baggingAnalysis$pt_values <- avgZvalues;
			baggingAnalysis$pWilcox_values <- avgZvalues;
			baggingAnalysis$pF_values <- avgZvalues;
			baggingAnalysis$pBin_values <- avgZvalues;
			baggingAnalysis$mMSE_values <- avgZvalues;
			baggingAnalysis$uAcc_values <- avgZvalues;
			baggingAnalysis$rAcc_values <- avgZvalues;
			baggingAnalysis$uAUC_values <- avgZvalues;
			baggingAnalysis$rAUC_values <- avgZvalues;
			baggingAnalysis$idi_values <- avgZvalues;
			baggingAnalysis$nri_values <- avgZvalues;
			baggingAnalysis$zidi_values <- avgZvalues;
			baggingAnalysis$znri_values <- avgZvalues;
			baggingAnalysis$mAUC_values <- avgZvalues;
			baggingAnalysis$mACC_values <- avgZvalues;
			baggingAnalysis$coefstd <- avgZvalues;
			baggingAnalysis$coefficients <- avgZvalues;
			baggingAnalysis$wts <- avgZvalues;
		#	print(basecoef);
			avgsize <- msize-1;
			mado <- NA;
			rnames <- 0;
			nrep <- 1+2*(noequalSets);
			if ((msize > 1)&&(loops>1))
			{
				model$type=type;
				onames <- names(model$coefficients);
				mmult <- 1+1*(type=="COX");
				model$estimations <- rep(0,mmult*msize); 
				wts <- 0;
				model$coefficients <- rep(0,msize);
				names(model$coefficients) <- onames;
				modelmeans <- model$coefficients;
				coefEvolution <- c(0,model$coefficients);
				names(coefEvolution) <- c("Weight",names(model$coefficients));
			#	cat("\n");
				tot_cycles <- 0;
				b_casesample <- casesample;
				b_controlsample <- controlsample;
				for (m in 1:n_bootstrap)
				{
					if (n_bootstrap>1)
					{
						if (type!="LM")
						{
							b_casesample <- casesample[sample(1:nrowcases,nrowcases,replace=TRUE),]
							b_controlsample <- controlsample[sample(1:nrowcontrols,nrowcontrols,replace=TRUE),]						
							EquTrainSet <- rbind(b_casesample,b_controlsample)
						}
						else
						{
							EquTrainSet <- data[sample(1:nrow(data),nrow(data),replace=TRUE),];
						}
						theoutcome <- EquTrainSet[,Outcome];
						varOutcome <- var(theoutcome);
						observations <- nrow(EquTrainSet);
					}
					avgsize = 0;
					for (n in 1:loops)
					{
						if ((n %% 10) == 0) cat(".");
						feat <- theterms[[n]]
						avgsize = avgsize+length(feat);
						if (m==1)
						{
				#			cat(modelFormulas[n],"\n");

							if (length(toRemove)>0)
							{
								modelFormulas[n] <- paste(modelFormulas[n],"  ",sep="");
								for (rml in 1:length(toRemove))
								{
									modelFormulas[n] <- sub(toRemove[rml]," 1 ",modelFormulas[n],fixed=TRUE);
								}
			#					cat("After Rem:",modelFormulas[n],"\n");
								feat <- attr(terms(formula(modelFormulas[n])),"term.labels");
							}
						}
						if (length(feat)>0)
						{
							for (replicates in 1:nrep)
							{
								if (noequalSets)
								{
									if (maxTrainSamples > nrowcases)  trainCaseses <- b_casesample[sample(1:nrowcases,maxTrainSamples,replace=TRUE),]
									if (maxTrainSamples > nrowcontrols)  trainControls <- b_controlsample[sample(1:nrowcontrols,maxTrainSamples,replace=TRUE),]
									EquTrainSet <- rbind(trainCaseses,trainControls)
									theoutcome <- EquTrainSet[,Outcome];
									varOutcome <- var(theoutcome);
									observations <- nrow(EquTrainSet);
								}
								out <- modelFitting(formula(modelFormulas[n]),EquTrainSet,type,fitFRESA=TRUE);
								coef_Zanalysis <- NULL;
								if (!inherits(out, "try-error")) 
								{
									osize <- length(out$coefficients)					
									if (osize > 1)
									{
										tot_cycles = tot_cycles+1;
										curprediction <- predict.fitFRESA(out,EquTrainSet,predtype)
										residual <- as.vector(abs(curprediction-theoutcome));
										onames <- names(out$coefficients);
										znames <- onames;
										if ((type!="COX")||(class(out)=="fitFRESA")) znames <- onames[-1];
										if (predtype=="linear")
										{
	#										Rwts <- (varOutcome-sum(residual^2)/observations)/varOutcome; #Correlation
											gvar <- getVar.Res(out,data=EquTrainSet,Outcome=Outcome,type=type)
											coef_Zanalysis <- -qnorm(gvar$FP.value);
											baggingAnalysis$uMS_values[znames] <- baggingAnalysis$uMS_values[znames] + gvar$unitrainMSE;
											baggingAnalysis$rMS_values[znames] <- baggingAnalysis$rMS_values[znames] + gvar$redtrainMSE;
											baggingAnalysis$NeRI_values[znames] <- baggingAnalysis$NeRI_values[znames] + gvar$NeRIs;
											baggingAnalysis$pF_values[znames] <- baggingAnalysis$pF_values[znames] + gvar$FP.value;
											baggingAnalysis$pt_values[znames] <- baggingAnalysis$pt_values[znames] + gvar$tP.value;
											baggingAnalysis$pBin_values[znames] <- baggingAnalysis$pBin_values[znames] + gvar$BinP.value;
											baggingAnalysis$pWilcox_values[znames] <- baggingAnalysis$pWilcox_values[znames] + gvar$WilcoxP.value;
											baggingAnalysis$mMSE_values[znames] <- baggingAnalysis$mMSE_values[znames] + gvar$FullTrainMSE;
										}
										else
										{
	#										Rwts <- 2.0*(sum(1.0*(residual<0.5))/observations - 0.5); # 2*(ACC-0.5)
											gvar <- getVar.Bin(out,data=EquTrainSet,Outcome=Outcome,type=type)
											coef_Zanalysis <- gvar$z.IDIs;
											baggingAnalysis$uAcc_values[znames] <- baggingAnalysis$uAcc_values[znames] + gvar$uniTrainAccuracy;
											baggingAnalysis$rAcc_values[znames] <- baggingAnalysis$rAcc_values[znames] + gvar$redtrainAccuracy;
											baggingAnalysis$uAUC_values[znames] <- baggingAnalysis$uAUC_values[znames] + gvar$uniTrainAUC;
											baggingAnalysis$rAUC_values[znames] <- baggingAnalysis$rAUC_values[znames] + gvar$redtrainAUC;
											baggingAnalysis$idi_values[znames] <- baggingAnalysis$idi_values[znames] + gvar$IDIs;
											baggingAnalysis$nri_values[znames] <- baggingAnalysis$nri_values[znames] + gvar$NRIs;
											baggingAnalysis$zidi_values[znames] <- baggingAnalysis$zidi_values[znames] + gvar$z.IDIs;
											baggingAnalysis$znri_values[znames] <- baggingAnalysis$znri_values[znames] + gvar$z.NRIs;
											baggingAnalysis$mAUC_values[znames] <- baggingAnalysis$mAUC_values[znames] + gvar$fullTrainAUC;
											baggingAnalysis$mACC_values[znames] <- baggingAnalysis$mACC_values[znames] + gvar$fullTrainAccuracy;
										}
										infnum <- is.infinite(coef_Zanalysis)
										if (sum(infnum)>0)
										{	
											coef_Zanalysis[infnum] <- 20.0;
										}
										avgZvalues[znames] <- avgZvalues[znames] + coef_Zanalysis;
										coef_Zanalysis[coef_Zanalysis<0] <- 0.0;
										coef_Zanalysis[coef_Zanalysis>20] <- 20.0;
										Rwts <- sum(coef_Zanalysis)
										if (Rwts<=0) Rwts <- 1.0e-4;
										Rwts <- Rwts*Rwts;
			#							print(basecoef[onames]);
	#									print(out$coefficients);
	#									print(coef_Zanalysis)
	#									Rwts <- as.vector(Rwts/sum(abs(out$coefficients/basecoef[onames]))) #by sum of weights
										rnames <- append(rnames,tot_cycles)
										outmeans <- out$coefficients;
										wts = wts + Rwts;
										model$coefficients[onames] <- model$coefficients[onames] + Rwts*out$coefficients[onames];
										baggingAnalysis$coefficients[znames] <- baggingAnalysis$coefficients[znames] + out$coefficients[znames];
										baggingAnalysis$coefstd[znames] <- baggingAnalysis$coefstd[znames] + out$coefficients[znames]^2;
										baggingAnalysis$wts[znames] <- baggingAnalysis$wts[znames] + rep(Rwts,length(znames));
										coefEvolution <- rbind(coefEvolution,c(Rwts,model$coefficients/wts));
										addedZvalues[znames] <- addedZvalues[znames] + rep(1,length(znames));
										if (type=="COX")
										{
											fullmodelmeans <- abs(modelmeans[onames]);
											names(fullmodelmeans) <- onames 
											for (ei in 1:osize) 
											{
												outmeans[ei] <- out$estimations[osize+ei];
											}
											for (ei in onames) 
											{
												if (fullmodelmeans[ei]>0)
												{
													modelmeans[ei] <- 0.5*(modelmeans[ei] + outmeans[ei]); 
												}
												else
												{
													modelmeans[ei] <- outmeans[ei]; 
												}
											}
										}
				#						print(model$coefficients)
									}
								}
								else
								{
									cat("+");
				#					print(out$coef);
								}
							}
						}
					}
					avgsize = avgsize/loops;
#					cat("*");
				}
				if( wts>0)
				{
#					print(baggingAnalysis$coefficients^2);
#					print(baggingAnalysis$coefstd);
#					print(addedZvalues);
					model$coefficients <- model$coefficients/wts;
					baggingAnalysis$coefficients <- baggingAnalysis$coefficients/addedZvalues;
					baggingAnalysis$coefstd <- sqrt(baggingAnalysis$coefstd/addedZvalues-(baggingAnalysis$coefficients)^2);
					coefEvolution <- as.data.frame(coefEvolution);
					rownames(coefEvolution) <- rnames
					if (type == "COX")
					{
						model$estimations <- c(model$coefficients,modelmeans);
					}
					else
					{
						model$estimations <- model$coefficients;
					}
					avgZvalues <- avgZvalues/addedZvalues;
					baggingAnalysis$formula.list <- modelFormulas;
					baggingAnalysis$uMS_values <- baggingAnalysis$uMS_values/addedZvalues;
					baggingAnalysis$rMS_values <- baggingAnalysis$rMS_values/addedZvalues;
					baggingAnalysis$NeRI_values <- baggingAnalysis$NeRI_values/addedZvalues;
					baggingAnalysis$pt_values <- baggingAnalysis$pt_values/addedZvalues;
					baggingAnalysis$pWilcox_values <- baggingAnalysis$pWilcox_values/addedZvalues;
					baggingAnalysis$pF_values <- baggingAnalysis$pF_values/addedZvalues;
					baggingAnalysis$pBin_values <- baggingAnalysis$pBin_values/addedZvalues;
					baggingAnalysis$uAcc_values <- baggingAnalysis$uAcc_values/addedZvalues;
					baggingAnalysis$rAcc_values <- baggingAnalysis$rAcc_values/addedZvalues;
					baggingAnalysis$uAUC_values <- baggingAnalysis$uAUC_values/addedZvalues;
					baggingAnalysis$rAUC_values <- baggingAnalysis$rAUC_values/addedZvalues;
					baggingAnalysis$idi_values <- baggingAnalysis$idi_values/addedZvalues;
					baggingAnalysis$nri_values <- baggingAnalysis$nri_values/addedZvalues;
					baggingAnalysis$zidi_values <- baggingAnalysis$zidi_values/addedZvalues;
					baggingAnalysis$znri_values <- baggingAnalysis$znri_values/addedZvalues;
					baggingAnalysis$mAUC_values <- baggingAnalysis$mAUC_values/addedZvalues;
					baggingAnalysis$mACC_values <- baggingAnalysis$mACC_values/addedZvalues;
					baggingAnalysis$mMSE_values <- baggingAnalysis$mMSE_values/addedZvalues;
					baggingAnalysis$wts <- baggingAnalysis$wts/addedZvalues;
					baggingAnalysis$RelativeFrequency <- VarFrequencyTable/loops;
					baggingAnalysis$Jaccard.SM <- Jaccard.SM;
					baggingAnalysis$n_bootstrap <- n_bootstrap;

					model$baggingAnalysis <- baggingAnalysis;
					model$linear.predictors <- predict(model);
				}
				else
				{
					
					model <- modelFitting(formula(frma),data,type=type,fitFRESA=TRUE)
		#			print(model$coefficients)
					model$coefficients[is.nan(model$coefficients)] <- 0.0;
					model$coefficients[is.na(model$coefficients)] <- 0.0;
					model$estimations[is.nan(model$estimations)] <- 0.0;
					model$estimations[is.na(model$estimations)] <- 0.0;
				}
			}
			else
			{
				model <- modelFitting(formula(frma),data,type=type,fitFRESA=TRUE)
				model$coefficients[is.nan(model$coefficients)] <- 0.0;
				model$coefficients[is.na(model$coefficients)] <- 0.0;
				model$estimations[is.nan(model$estimations)] <- 0.0;
				model$estimations[is.na(model$estimations)] <- 0.0;
			}
		}
		else
		{
			model <- modelFitting(formula(frma),data,type=type,fitFRESA=TRUE);
		}
		# print(model$coefficients);
		environment(model$formula) <- globalenv()
		environment(model$terms) <- globalenv()
		
	}
	else
	{
		model <- NULL;
		frma <- NULL;
		VarFrequencyTable <- NULL;
		avgsize <- 0;
		formulaNetwork <- NULL;
		Jaccard.SM <- 0;
		coefEvolution <- NULL;
	}
	

  	result <- list(bagged.model=model,
				   formula=frma,
				   frequencyTable=VarFrequencyTable,
				   averageSize=avgsize,
				   formulaNetwork=formulaNetwork,
				   Jaccard.SM = Jaccard.SM,
				   coefEvolution=coefEvolution,
				   avgZvalues=avgZvalues
				   );
  
	return (result);
}
