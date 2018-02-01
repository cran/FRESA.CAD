reportEquivalentVariables <-
function (object,pvalue=0.05,data,variableList,Outcome="Class",timeOutcome=NULL, type = c("LOGIT", "LM","COX"),eqPGain=1.0,description=".",method="BH",osize=0,fitFRESA=TRUE) 
{
    type <- match.arg(type);
  
  
	varsList <- unlist(as.list(attr(terms(object),"variables")))
	termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")
	vnames <- as.vector(variableList[,1]);
	
#	print(termList)
#	print(vnames)
#	print(colnames(data))
    
	orgsize <- length(termList);
	if (osize==0) osize=max(ncol(data),nrow(variableList));
	if (osize<length(vnames)) osize=length(vnames);
    

	outCome = paste(varsList[2]," ~ 1");
	auxmodel <- object;
	

	if (type != "LM")
	{
		orgpVAL <- 1.0-pnorm(as.vector(getVar.Bin(object,data,Outcome,type)$testData.z.IDIs));
	}
	else
	{
		orgpVAL <- as.vector(getVar.Res(object,data,Outcome,type)$FP.value);
	}

#	print (orgpVAL,digits=6);
	adj.pVAL <- eqPGain * orgpVAL;
	adj.pVAL[adj.pVAL>pvalue] <- pvalue;
#	print (orgpVAL,digits=6);

	formulaList <- as.character();
	pvalues <- numeric(length(vnames));
	tImprovement <- numeric(length(vnames));
	tunitMetric <- numeric(length(vnames));
	tReducedMetric <- numeric(length(vnames));
	tFullMetric <- numeric(length(vnames));
	formuls <- character(length(vnames));
	theLocus <- numeric(length(vnames));
	pvalueList <- list();
	thename <- character();
	theparent <- character();
	thedescription <- character();
	theImprovement <- numeric();
	theUniPerformance <- numeric();
	theDeltaPerformance <- numeric();
	theFullPerformance <- numeric();
	thepvalue <- numeric();
	smcoff <- NULL;
	varoutcome <- var(as.vector(data[,Outcome]));
	coff=1;
	if (type=="COX") coff=0;
#	print(osize)
#	print(termList)
#	print(vnames)
	if (length( termList)>0)
	{
		for ( i in 1:length(termList))
		{

	#		cat(as.character(termList[i]),"->");	
			cat(":");
			namevector <- character();
			pVALlist <- list();
			the_locus=1;
			while ((the_locus<i)&&(vnames[the_locus]!=termList[i])) the_locus=the_locus+1;
			for (j in 1:length(vnames))
			{
				frm1 = outCome;
				if ( j%%100 == 0) cat(".");
				if ( j%%500 == 0) cat("#");
				for ( n in 1:length(termList))
				{
					if (n == i)
					{
						frm1 <- paste(frm1,paste("+",vnames[j]));
					}
					else
					{
						frm1 <- paste(frm1,paste("+",termList[n]));
					}
				}
	#			print(frm1);
				ftmp <- formula(frm1);
				sdata <- data[,all.vars(ftmp)];
				auxmodel <- modelFitting(ftmp,sdata,type,fitFRESA=TRUE);
				pVAL <- rep(1,length(orgpVAL));
				unitMetric <- rep(NA,length(orgpVAL));
				redMetric <- rep(NA,length(orgpVAL));
				fullMetric <- 0;
				testImprovement <- rep(NA,length(orgpVAL));
				if (orgsize == length(attr(terms(auxmodel),"term.labels")))
				{
					if (type != "LM")
					{
						varana <- getVar.Bin(auxmodel,sdata,Outcome,type);
						pVAL <- 1.0-pnorm(varana$testData.z.IDIs);
						unitMetric <- varana$uniTestAUC;
						redMetric <- varana$redtestAUC;
						fullMetric <- varana$fullTestAUC;
						testImprovement <- 0.5*varana$testData.NRIs;
					}
					else
					{
						varana <- getVar.Res(auxmodel,sdata,Outcome,type);
						pVAL <- varana$testData.FP.value;
						unitMetric <- (varoutcome-as.vector(varana$unitestMSE))/varoutcome;
						redMetric <- (varoutcome-as.vector(varana$redtestMSE))/varoutcome;
						fullMetric <- (varoutcome-as.vector(varana$FullTestMSE))/varoutcome;
						testImprovement <- varana$NeRIs;
					}
				}
				theLocus[j]=the_locus;
				pvalues[j]=pVAL[i];
				tunitMetric[j]=unitMetric[i];
				tReducedMetric[j]=redMetric[i];
				tFullMetric[j]=fullMetric;
	#			cat("idx=",(i)+1,"\n");
				tImprovement[j]=testImprovement[i];
				formuls[j]=frm1;
	#			print(auxmodel$coefficients);
	#			print(pvalues[j]);
	#			print(tImprovement[j]);
			}
			apvalues <- p.adjust(pvalues,method,osize);
			
			for (j in 1:length(vnames))
			{
	#			cat(orgpVAL[i],":",vnames[j],". ad p:",apvalues[j],". un p:",pvalues[j],". a Locust p:",apvalues[theLocus[j]],". u Locus p:",pvalues[theLocus[j]],"\n");
				if (
					( as.character(vnames[j]) == as.character(termList[i]) ) 
					|| 
					(
						   (( apvalues[j] <= pvalue ) && 
							( apvalues[j] <= eqPGain*apvalues[theLocus[j]]))
						|| (( pvalues[j] <= adj.pVAL[i] ) &&
							(apvalues[theLocus[j]] <= 2.0*pvalue )
							)
					)
				)
				{
					listindx <- length(pVALlist)+1;
					namevector <- append(namevector,vnames[j]);
					pVALlist[[listindx]] <- pvalues[j];
					  
					if (as.character(vnames[j])!=as.character(termList[i])) 
					{
						formulaList <- append(formulaList,formuls[j]);
					}
					else
					{
						if (length(formulaList)==0) formulaList <- append(formulaList,formuls[j]);
					}
					if (description != ".") 
					{
						thedescription <- append(thedescription,variableList[j,description]);
					}
					else
					{
						thedescription <- append(thedescription,paste(vnames[j],termList[i],sep=":"));
					}
					if (!fitFRESA)
					{
						sm <- summary(modelFitting(formula(formuls[j]),data,type,fitFRESA=fitFRESA));
						smcoff <- rbind(smcoff,sm$coefficients[(i)+coff,]);
					}
					
					thename <- append(thename,vnames[j]);
					theparent <- append(theparent,as.character(termList[i]));
					theImprovement <- append(theImprovement,tImprovement[j]);
					thepvalue <- append(thepvalue,pvalues[j]);
					theUniPerformance <- append(theUniPerformance,tunitMetric[j]);
					theDeltaPerformance <- append(theDeltaPerformance,tFullMetric[j]-tReducedMetric[j]);
					theFullPerformance <- append(theFullPerformance,tFullMetric[j]);
				}
			}
			names(pVALlist) <- namevector;
			listindx <- length(pvalueList)+1;
			pvalueList[[listindx]] <- pVALlist;
		}
		names(pvalueList) <- termList;
		equmodel <- baggedModel(formulaList,data,type,Outcome,timeOutcome,frequencyThreshold=0.001)$bagged.model;
	}
	else
	{
		equmodel <- NULL;
	}

#	cat("Equivalent\n")
#	print(formulaList);



	Mresult <- NULL
	Mresult$Name <- thename;
	Mresult$Locus <- theparent;
	Mresult$Extendend_Name <- thedescription;
	Mresult$UniPerformance <- as.numeric(theUniPerformance);
	Mresult$FullPerformance <- as.numeric(theFullPerformance);
	Mresult$DeltaPerformance <- as.numeric(theDeltaPerformance);
	Mresult$ImprovementFraction <- as.numeric(theImprovement);
	Mresult <- as.data.frame(Mresult) 
	if (!is.null(smcoff))
	{
		colnames(smcoff) <- colnames(sm$coefficients);
		Mresult <- cbind(Mresult,smcoff);
	}
	Mresult$p.value <- as.numeric(thepvalue);
	result <- list(pvalueList=pvalueList,equivalentMatrix=Mresult,formula.list=formulaList,equivalentModel=equmodel);
    return (result)
}
