univariateRankVariables <-
function(varList,baseModel,Outcome,dataframe,categorizationType=c("Raw","Categorical","ZCategorical","RawZCategorical","RawZTail","RawTail"),FitType=c("LOGIT","LM","COX"),rankingTest=c("zIDI","zNRI","IDI","NRI","NeRI","Ztest","AUC","CStat","Kendall"),cateGroups=c(0.1,0.9),raw.dataFrame=NULL,description=".",uniType=c("Binary","Regression")) 
{

	if (is.null(raw.dataFrame))  raw.dataFrame <- dataframe;
	FitType <- match.arg(FitType);
	uniType <- match.arg(uniType);
	categorizationType <- match.arg(categorizationType);
	rankingTest <- match.arg(rankingTest);
	colnamesList <- as.vector(varList[,1]);


	if (description == ".")
	{
		descripList <- colnamesList;
	}
	else
	{
		descripList <- as.vector(varList[,description]);
	}
	
	if (uniType=="Binary")
	{
		caserawsample <- subset(raw.dataFrame,get(Outcome)  == 1);
		controlrawsample <- subset(raw.dataFrame,get(Outcome) == 0);

		caseZsample <- subset(dataframe,get(Outcome)  == 1);
		controlZsample <- subset(dataframe,get(Outcome) == 0);

		sizecaseZsample <- nrow(caseZsample);
		sizecontrolZsample <- nrow(controlZsample);
	}
	else
	{
		caserawsample <- NULL;
		controlrawsample <- NULL;

		caseZsample <- NULL;
		controlZsample <- NULL;

		sizecaseZsample <- NULL;
		sizecontrolZsample <- NULL;
	}

	size = length(colnamesList);
	Name <- vector();

	parent <- vector();
	descrip <- vector();


	IDI <- vector();
	NRI <- vector();
	zIDI <- vector();
	zNRI <- vector();
	ROCAUC <- vector();
	ZGLM <- vector();

	NeRI <- vector();
	BinRes.p <- vector();
	WilcoxRes.p <- vector();
	TstudentRes.p <- vector();
	FRes.p <- vector();


	cohortMean <- vector();
	cohortStd <- vector();
	cohortKSD <- vector();
	cohortKSP <- vector();
	cohortZKSP <- vector();
	cohortZKSD <- vector();

	
	caseMean <- vector();
	caseStd <- vector();
	caseKSD <- vector();
	caseKSP <- vector();
	caseZKSP <- vector();
	caseZKSD <- vector();
	caseN_Z_Low_Tail <- vector();
	caseN_Z_Hi_Tail <- vector();

	controlMean <- vector();
	controlStd <- vector();
	controlKSD <- vector();
	controlKSP <- vector();
	controlZKSP <- vector();
	controlZKSD <- vector();
	controlN_Z_Low_Tail <- vector();
	controlN_Z_Hi_Tail <- vector();
	kendall.r <- vector();
	spearman.r <- vector();
	pearson.r <- vector();
	kendall.p <- vector();
	cStatCorr <- vector();
	
	t.Rawvalue <- vector();
	t.Zvalue <- vector();
	wilcox.Zvalue <- vector();
	
	frm1 <- baseModel;
	ftmp <- formula(frm1);
	bmodel <- modelFitting(ftmp,dataframe,FitType)
	baseResiduals <- residualForNeRIs(bmodel,newdata=dataframe,Outcome);
	basepredict <- predictForFresa(bmodel,dataframe,type = 'prob');
	if (FitType=="COX")
	{
		termslist <- attr(terms(ftmp),"term.labels");
		if (length(termslist)==0)
		{
			ftmp <- formula(paste(frm1,"+dummy"));
			dummy <-  rnorm(nrow(dataframe));
			coxframe <- cbind(dataframe,dummy);
			bmodel <- modelFitting(ftmp,coxframe,FitType)
			basepredict <- predictForFresa(bmodel,coxframe,type = 'prob');
			baseResiduals <- residualForNeRIs(bmodel,newdata=coxframe,Outcome);
		}
	}
	
	for (j in 1:size)
	{
		frm1 = baseModel;
		categories = 1;
		catlist <- vector();
		caseCount1 = 0;
		caseCount2 = 0;
		controlCount1 = 0;
		controlCount2 = 0;
		stddf=0;
		kendcor <- NA;
		pearcor <- NA;
		speacor <- NA;
		cstat <- NA;
		kstZdf <- NA;
		kstdf <- NA;
		meCa <- NA;
		stdCa <- NA;
		kstCa <- NA;
		kstZCa <- NA;
		meCo <- NA;
		stdCo <- NA;
		kstCo <- NA;
		kstZCo <- NA;
		rtt <- NA;
		ztt <- NA;
		if (uniType=="Binary")
		{
			wtt <- -qnorm(wilcox.test(controlrawsample[,colnamesList[j]],caserawsample[,colnamesList[j]],na.action=na.exclude)$p.value,0,1);
		}
		else
		{
			wtt <- NA;
		}
		stddf <- sd(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
		if (stddf>0)
		{
			kendcor <- try(cor.test(dataframe[,colnamesList[j]],dataframe[,Outcome],method="kendall",na.action=na.omit));
			pearcor <- try(cor.test(dataframe[,colnamesList[j]],dataframe[,Outcome],method="pearson",na.action=na.omit));
			speacor <- try(cor.test(dataframe[,colnamesList[j]],dataframe[,Outcome],method="spearman",na.action=na.omit));
			if (kendcor$estimate > 0)
			{
				cstat <- rcorr.cens(dataframe[,colnamesList[j]],dataframe[,Outcome], outx=FALSE)
			}
			else
			{
				cstat <- rcorr.cens(dataframe[,colnamesList[j]],-dataframe[,Outcome], outx=FALSE)
			}
		}
		if (length(table(dataframe[,colnamesList[j]]))>2)
		{
			if (uniType=="Binary")
			{
				meCa <- mean(caseZsample[,colnamesList[j]],na.rm = TRUE);
				stdCa <- sd(caseZsample[,colnamesList[j]],na.rm = TRUE);
				kstZCa <- ks.test(caseZsample[,colnamesList[j]],"pnorm",meCa,stdCa);
				
				meCo <- mean(controlZsample[,colnamesList[j]],na.rm = TRUE);
				stdCo <- sd(controlZsample[,colnamesList[j]],na.rm = TRUE);
				kstZCo <- ks.test(controlZsample[,colnamesList[j]],"pnorm",meCo,stdCo);
				
				meCa <- mean(caserawsample[,colnamesList[j]],na.rm = TRUE);
				stdCa <- sd(caserawsample[,colnamesList[j]],na.rm = TRUE);
				kstCa <- ks.test(caserawsample[,colnamesList[j]],"pnorm",meCa,stdCa);
				
				meCo <- mean(controlrawsample[,colnamesList[j]],na.rm = TRUE);
				stdCo <- sd(controlrawsample[,colnamesList[j]],na.rm = TRUE);
				kstCo <- ks.test(controlrawsample[,colnamesList[j]],"pnorm",meCo,stdCo);
				
				rtt <- try(t.test(controlrawsample[,colnamesList[j]],caserawsample[,colnamesList[j]],na.action=na.exclude));
				ztt <- try(t.test(controlZsample[,colnamesList[j]],caseZsample[,colnamesList[j]],na.action=na.exclude));

				if (!is.na(cateGroups[1]))
				{
					zthr = sprintf("%5.3f",qnorm(cateGroups[1]));
					subtestt <- paste(Outcome,paste("== 1 & ",paste(colnamesList[j],paste("<",zthr))));
					caseCount1 = eval(parse(text = paste("nrow(subset(dataframe,",paste(subtestt,"))")))); 
					subtestt <- paste(Outcome,paste("== 0 & ",paste(colnamesList[j],paste("<",zthr))));
					controlCount1 = eval(parse(text = paste("nrow(subset(dataframe,",paste(subtestt,"))")))); 

					categories=length(cateGroups);
					if (!is.na(cateGroups[categories]))
					{
						zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
						subtestt <- paste(Outcome,paste("== 1 & ",paste(colnamesList[j],paste(">",zthr))));
						caseCount2 = eval(parse(text = paste("nrow(subset(dataframe,",paste(subtestt,"))")))); 
						subtestt <- paste(Outcome,paste("== 0 & ",paste(colnamesList[j],paste(">",zthr))));
						controlCount2 = eval(parse(text = paste("nrow(subset(dataframe,",paste(subtestt,"))")))); 
					}
					else
					{
						zthr = sprintf("%5.3f",1-qnorm(cateGroups[1]));
						subtestt <- paste(Outcome,paste("== 1 & ",paste(colnamesList[j],paste(">",zthr))));
						caseCount2 = eval(parse(text = paste("nrow(subset(dataframe,",paste(subtestt,"))")))); 
						subtestt <- paste(Outcome,paste("== 0 & ",paste(colnamesList[j],paste(">",zthr))));
						controlCount2 = eval(parse(text = paste("nrow(subset(dataframe,",paste(subtestt,"))")))); 
					}
				}

			}
			medf <- mean(dataframe[,colnamesList[j]],na.rm = TRUE);
			stddf <- sd(dataframe[,colnamesList[j]],na.rm = TRUE);
			kstZdf <- ks.test(dataframe[,colnamesList[j]],"pnorm",medf,stddf);
			medf <- mean(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
			stddf <- sd(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
			kstdf <- ks.test(raw.dataFrame[,colnamesList[j]],"pnorm",medf,stddf);

			switch(categorizationType,
				Raw =
				{
					categories=1;
					catlist <- append(catlist,colnamesList[j]);				
				},
				Categorical =
				{
					categories=length(cateGroups);

					zthr = sprintf("%5.3f",qnorm(cateGroups[1]));


					for (n in 1:categories)
					{
						if (n==1)
						{
							zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
							catvar = paste("I(",colnamesList[j]);
							catvar = paste(catvar," < ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,")");
							catlist <- append(catlist,catvar);
						}
						else
						{
							zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
							zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
							catvar = paste("I((",colnamesList[j]);
							catvar = paste(catvar," >= ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,") & (");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," < ");
							catvar = paste(catvar,zthr2);
							catvar = paste(catvar,"))");
							catlist <- append(catlist,catvar);
						}
					}
					zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
					catvar = paste("I(",colnamesList[j]);
					catvar = paste(catvar," >= ");
					catvar = paste(catvar,zthr);
					catvar = paste(catvar,")");
					catlist <- append(catlist,catvar);
					categories = categories + 1;
				},
				ZCategorical =
				{
					categories=length(cateGroups);
					zthr = sprintf("%5.3f",qnorm(cateGroups[1]));

					for (n in 1:categories)
					{

						if (n==1)
						{
							zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
							catvar = paste("I(",colnamesList[j]);
							catvar = paste(catvar,"* (");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," < ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,"))");
							catlist <- append(catlist,catvar);
						}
						else
						{
							zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
							zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
							catvar = paste("I(",colnamesList[j]);
							catvar = paste(catvar,"* ((");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," >= ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,") & (");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," < ");
							catvar = paste(catvar,zthr2);
							catvar = paste(catvar,")))");
							catlist <- append(catlist,catvar);
						}
					}
					zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
					catvar = paste("I(",colnamesList[j]);
					catvar = paste(catvar,"* (");
					catvar = paste(catvar,colnamesList[j]);
					catvar = paste(catvar," >= ");
					catvar = paste(catvar,zthr);
					catvar = paste(catvar,"))");
					catlist <- append(catlist,catvar);
					categories = categories+1;
				},
				RawZCategorical =
				{
					categories=length(cateGroups);

					zthr = sprintf("%5.3f",qnorm(cateGroups[1]));

					catlist <- append(catlist,colnamesList[j]);				
					for (n in 1:categories)
					{

						if (n==1)
						{
							zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
							catvar = paste("I(",colnamesList[j]);
							catvar = paste(catvar,"* (");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," < ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,"))");
							catlist <- append(catlist,catvar);
						}
						else
						{
							zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
							zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
							catvar = paste("I(",colnamesList[j]);
							catvar = paste(catvar,"* ((");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," >= ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,") & (");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," < ");
							catvar = paste(catvar,zthr2);
							catvar = paste(catvar,")))");
							catlist <- append(catlist,catvar);
						}
					}
					zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
					catvar = paste("I(",colnamesList[j]);
					catvar = paste(catvar,"* (");
					catvar = paste(catvar,colnamesList[j]);
					catvar = paste(catvar," >= ");
					catvar = paste(catvar,zthr);
					catvar = paste(catvar,"))");
					catlist <- append(catlist,catvar);
					categories = categories+2;
				},			
				RawZTail =
				{
					categories = 1;
					catlist <- append(catlist,colnamesList[j]);				
					zthr = sprintf("%5.3f",qnorm(cateGroups[1]));

					f1= caseCount1/sizecaseZsample;
					f2= controlCount1/sizecontrolZsample;
					if ((f1>f2)&&(f1>0.1))
					{
						catvar = paste("I(",colnamesList[j]);
						catvar = paste(catvar,"* (");
						catvar = paste(catvar,colnamesList[j]);
						catvar = paste(catvar," < ");
						catvar = paste(catvar,zthr);
						catvar = paste(catvar,"))");
						catlist <- append(catlist,catvar);
						categories = categories+1;
					}				


					zthr = sprintf("%5.3f",qnorm(1.0-cateGroups[1]));
					f1= caseCount2/sizecaseZsample;
					f2= controlCount2/sizecontrolZsample;
					if ((f1>f2)&&(f1>0.1))
					{
						catvar = paste("I(",colnamesList[j]);
						catvar = paste(catvar,"* (");
						catvar = paste(catvar,colnamesList[j]);
						catvar = paste(catvar," > ");
						catvar = paste(catvar,zthr);
						catvar = paste(catvar,"))");
						catlist <- append(catlist,catvar);
						categories = categories+1;
					}
				},			
				RawTail =
				{
					categories = 1;
					catlist <- append(catlist,colnamesList[j]);				
					zthr = sprintf("%5.3f",qnorm(cateGroups[1]));

					f1 = caseCount1/sizecaseZsample;
					f2 = controlCount1/sizecontrolZsample; 
					if ((f1>f2)&&(f1>0.1))   # will add only if fraction is greater
					{
						catvar = paste("I(",colnamesList[j]);
						catvar = paste(catvar," < ");
						catvar = paste(catvar,zthr);
						catvar = paste(catvar,")");
						catlist <- append(catlist,catvar);
						categories = categories+1;
					}				


					zthr = sprintf("%5.3f",qnorm(1.0-cateGroups[1]));
					f1= caseCount2/sizecaseZsample;
					f2= controlCount2/sizecontrolZsample;
					if ((f1>f2)&&(f1>0.1)) # will add only if fraction is greater
					{
						catvar = paste("I(",colnamesList[j]);
						catvar = paste(catvar," > ");
						catvar = paste(catvar,zthr);
						catvar = paste(catvar,")");
						catlist <- append(catlist,catvar);
						categories = categories+1;
					}
				},			
				{
					categories=1;
					catlist <- append(catlist,colnamesList[j]);
				}
			)
		}
		else
		{
			categories=1;
			catlist <- append(catlist,colnamesList[j]);	
			medf = table(dataframe[,colnamesList[j]])[1];			
			stddf = table(dataframe[,colnamesList[j]])[2];			
			if (uniType=="Binary")
			{
				meCa <- table(caseZsample[,colnamesList[j]])[1];
				stdCa <- table(caseZsample[,colnamesList[j]])[2];
				
				meCo <- table(controlZsample[,colnamesList[j]])[1];
				stdCo <- table(controlZsample[,colnamesList[j]])[2];
			}
		}
		for (n in 1:categories)
		{
			termName <- str_replace_all(catlist[n]," ","");
			termName <- str_replace_all(termName,"<"," < ");
			termName <- str_replace_all(termName,">"," > ");
			termName <- str_replace_all(termName,"&"," & ");
			termName <- str_replace_all(termName,"=","= ");
			termName <- str_replace_all(termName,fixed("> ="),">=");
			termName <- str_replace_all(termName,fixed("*")," * ");
			cat(termName,"\n");
			frmg <- paste( baseModel,paste(" + ",termName));
			ftmg <- formula(frmg);
			if (FitType=="COX") 
			{
				zcol=4;
			}
			else
			{
				zcol=3;
			}
			lmodel <- modelFitting(ftmg,dataframe,FitType)
			if (!inherits(lmodel, "try-error"))
			{
				modcoef <- summary(lmodel)$coefficients;
				sizecoef <- length(lmodel$coef);
			}
			else
			{
				modcoef <- NULL;
				sizecoef <- NULL;
			}
				Name <- append(Name,termName);
				parent <- append(parent,colnamesList[j])
				descrip <- append(descrip,descripList[j])
				if (uniType=="Binary")
				{
					caseMean <- append(caseMean,meCa); 
					caseStd <- append(caseStd,stdCa);
					if (!is.na(kstCa[[1]])) 
					{
						caseKSD <- append(caseKSD,kstCa$statistic);
						caseKSP <- append(caseKSP,kstCa$p.value);
						caseZKSP <- append(caseZKSP,kstZCa$p.value);
						caseZKSD <- append(caseZKSD,kstZCa$statistic);
					}
					else
					{
						caseKSD <- append(caseKSD,NA);
						caseKSP <- append(caseKSP,NA);
						caseZKSP <- append(caseZKSP,NA);
						caseZKSD <- append(caseZKSD,NA);
					}
					controlMean <- append(controlMean,meCo); 
					controlStd <- append(controlStd,stdCo);
					if (!is.na(kstCo[[1]]))
					{
						controlKSD <- append(controlKSD,kstCo$statistic);
						controlKSP <- append(controlKSP,kstCo$p.value);
						controlZKSP <- append(controlZKSP,kstZCo$p.value);
						controlZKSD <- append(controlZKSD,kstZCo$statistic);
					}
					else
					{
						controlKSD <- append(controlKSD,NA);
						controlKSP <- append(controlKSP,NA);
						controlZKSP <- append(controlZKSP,NA);
						controlZKSD <- append(controlZKSD,NA);
					}
					if (!is.na(rtt[[1]]))
					{
						if ( !inherits(rtt, "try-error"))
						{
							t.Rawvalue <- append(t.Rawvalue,rtt$statistic);
						}
						else
						{
							t.Rawvalue <- append(t.Rawvalue,NA);
						}
						if ( !inherits(ztt, "try-error"))
						{
							t.Zvalue <- append(t.Zvalue,ztt$statistic);
						}
						else
						{
							t.Zvalue <- append(t.Zvalue,NA);
						}
					}
					else
					{
						t.Rawvalue <- append(t.Rawvalue,NA);
						t.Zvalue <- append(t.Zvalue,NA);
					}
					if (!is.na(wtt[[1]]))
					{
						wilcox.Zvalue <- append(wilcox.Zvalue,wtt);
					}
					else
					{
						wilcox.Zvalue <- append(wilcox.Zvalue,NA);
					}
				}
			
			cohortMean <- append(cohortMean,medf); 
			cohortStd <- append(cohortStd,stddf);
			if (!is.na(kstdf[[1]]))
			{
				cohortKSD <- append(cohortKSD,kstdf$statistic);
				cohortKSP <- append(cohortKSP,kstdf$p.value);
				cohortZKSP <- append(cohortZKSP,kstZdf$p.value);
				cohortZKSD <- append(cohortZKSD,kstZdf$statistic);
			}
			else
			{
				cohortKSD <- append(cohortKSD,NA);
				cohortKSP <- append(cohortKSP,NA);
				cohortZKSP <- append(cohortZKSP,NA);
				cohortZKSD <- append(cohortZKSD,NA);
			}
			
			if (!is.na(kendcor[[1]])) 
			{
				kendall.r <- append(kendall.r,kendcor$estimate);
				kendall.p <- append(kendall.p,kendcor$p.value);
				pearson.r <- append(pearson.r,pearcor$estimate);
				spearman.r <- append(spearman.r,speacor$estimate);
				cStatCorr <- append(cStatCorr,cstat[1]);
			}
			else
			{
				kendall.r <- append(kendall.r,NA);
				kendall.p <- append(kendall.p,NA);
				pearson.r <- append(pearson.r,NA);
				spearman.r <- append(spearman.r,NA);
				cStatCorr <- append(cStatCorr,NA);
			}
			
			
			if (is.null(sizecoef) || is.na(lmodel$coef[sizecoef])) 
			{
				test=NA;
				if (uniType=="Binary")
				{
					IDI <- append(IDI,test);
					NRI <- append(NRI,test);
					zIDI <- append(zIDI,test);
					zNRI <- append(zNRI,test);
					ROCAUC <- append(ROCAUC,test);
					caseN_Z_Low_Tail <- append(caseN_Z_Low_Tail,test);
					caseN_Z_Hi_Tail <- append(caseN_Z_Hi_Tail,test);
					controlN_Z_Low_Tail <- append(controlN_Z_Low_Tail,test);
					controlN_Z_Hi_Tail <- append(controlN_Z_Hi_Tail,test);
				}
				ZGLM <- append(ZGLM,test);
				NeRI <- append(NeRI,test);
				BinRes.p <- append(BinRes.p,test);
				WilcoxRes.p <- append(WilcoxRes.p,test);
				TstudentRes.p <- append(TstudentRes.p,test);
				FRes.p <- append(FRes.p,test);
			}
			else
			{

				if (uniType=="Binary")
				{
					spredict <- predictForFresa(lmodel,dataframe,type = 'prob');
					iprob <- improveProb(basepredict,spredict,dataframe[,Outcome]);
					IDI <- append(IDI,iprob$idi);
					NRI <- append(NRI,iprob$nri);
					zIDI <- append(zIDI,iprob$z.idi);
					zNRI <- append(zNRI,iprob$z.nri);
					if (length(dataframe[,Outcome])==length(spredict))
					{
						ROCAUC <- append(ROCAUC,roc( dataframe[,Outcome], spredict,plot=FALSE,auc=TRUE)$auc[1]);
					}
					else 
					{
						ROCAUC <- append(ROCAUC,NA);
					}
					caseN_Z_Low_Tail <- append(caseN_Z_Low_Tail,caseCount1);
					caseN_Z_Hi_Tail <- append(caseN_Z_Hi_Tail,caseCount2);
					controlN_Z_Low_Tail <- append(controlN_Z_Low_Tail,controlCount1);
					controlN_Z_Hi_Tail <- append(controlN_Z_Hi_Tail,controlCount2);
				}

				varResiduals <- residualForNeRIs(lmodel,newdata=dataframe,Outcome);
				rprob <- improvedResiduals(baseResiduals,varResiduals);

				ZGLM  <- append(ZGLM,abs(modcoef[sizecoef,zcol]));
				NeRI <- append(NeRI,rprob$NeRI);
				BinRes.p <- append(BinRes.p,rprob$p.value);
				WilcoxRes.p <- append(WilcoxRes.p,rprob$wilcox.pValue);
				TstudentRes.p <- append(TstudentRes.p,rprob$t.test.pValue);
				FRes.p <- append(FRes.p,rprob$F.test.pValue);
			}
		}
	}
	
	if (uniType=="Binary")
	{
		orderframe <- data.frame(Name,parent,descrip,cohortMean,cohortStd,cohortKSD,cohortKSP,caseMean,
		caseStd,caseKSD,caseKSP,caseZKSD,caseZKSP,controlMean,controlStd,controlKSD,controlKSP,controlZKSD,
		controlZKSP,t.Rawvalue,t.Zvalue,wilcox.Zvalue,ZGLM,zNRI,zIDI,ROCAUC,cStatCorr,NRI,IDI,NeRI,kendall.r,
		kendall.p,BinRes.p,TstudentRes.p,WilcoxRes.p,FRes.p,caseN_Z_Low_Tail,caseN_Z_Hi_Tail,controlN_Z_Low_Tail,controlN_Z_Hi_Tail);
		switch(rankingTest,
			zIDI=
			{
				orderframe <- with(orderframe,orderframe[order(-zIDI),]);
			},
			zNRI=
			{
				orderframe <- with(orderframe,orderframe[order(-zNRI),]);
			},
			IDI=
			{
				orderframe <- with(orderframe,orderframe[order(-IDI),]);
			},
			NRI=
			{
				orderframe <- with(orderframe,orderframe[order(-NRI),]);
			},
			NeRI=
			{
				orderframe <- with(orderframe,orderframe[order(-NeRI),]);
			},
			Ztest=
			{
				orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
			},
			AUC=
			{
				orderframe <- with(orderframe,orderframe[order(-ROCAUC),]);
			},
			Kendall=
			{
				orderframe <- with(orderframe,orderframe[order(kendall.p),]);
			},
			{
				orderframe <- with(orderframe,orderframe[order(-ROCAUC),]);
			}
		)

	}
	else
	{
		orderframe <- data.frame(Name,parent,descrip,cohortMean,cohortStd,cohortKSD,cohortKSP,cohortZKSD,cohortZKSP,ZGLM,
		NeRI,cStatCorr,spearman.r,pearson.r,kendall.r,kendall.p,BinRes.p,TstudentRes.p,WilcoxRes.p,FRes.p);
		switch(rankingTest,
			NeRI=
			{
				orderframe <- with(orderframe,orderframe[order(-NeRI),]);
			},
			Ztest=
			{
				orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
			},
			CStat=
			{
				orderframe <- with(orderframe,orderframe[order(-cStatCorr),]);
			},
			Kendall=
			{
				orderframe <- with(orderframe,orderframe[order(kendall.p),]);
			},
			{
				orderframe <- with(orderframe,orderframe[order(-cStatCorr),]);
			}
		)
	}

	row.names(orderframe) <- orderframe$Name;
	return (orderframe);
	
}
