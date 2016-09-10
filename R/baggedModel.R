baggedModel <-
function(modelFormulas,data,type=c("LM","LOGIT","COX"),Outcome=NULL,timeOutcome=NULL,pvalue=0.05,backElimination=FALSE,frequencyThreshold=0.05,removeOutliers=4.0)
{
	type <- match.arg(type)

	
	loops <- length(modelFormulas);
	theoutcome <- data[,Outcome];
	binoutcome <- (length(table(theoutcome))==2) && (min(theoutcome)==0);
	predtype="linear";
	if (binoutcome) predtype="prob";
	if ( (type=="LM") && (binoutcome==TRUE) )
	{
		data[,Outcome] = 2*theoutcome-1.0;
	}
	features <- vector();
	observations <- nrow(data);
	
#	cat(removeOutliers,"\n");

	for (n in 1:loops)
	{
		frm <- gsub("+ ","+",modelFormulas[n],fixed = TRUE);
		frm <- gsub(" +","+",frm,fixed = TRUE);
		frm <- gsub("  ","",frm,fixed = TRUE);
#		print(frm);
		if (gregexpr(pattern ='~',frm)[1]>0)
		{
			feat <- unlist(strsplit(frm,"[~]"));
			feat <- unlist(strsplit(feat[2],"[+]"));
		}
		else
		{
			feat <- unlist(strsplit(frm,"[+]"));
		}
		if (length(feat)>0)
		{
			for (i in 1:length(feat))
			{
				if ((feat[i]!="1")&&(feat[i]!=" ")&&(feat[i]!="  ")&&(feat[i]!="")) features <- append(features,feat[i]);
			}
		}
	}
	VarFrequencyTable <- table(features);
	VarFrequencyTable <- VarFrequencyTable[order(-VarFrequencyTable)]
	vnames <- rownames(VarFrequencyTable);

	
	nsize <- nrow(data)
	
	baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome);
	  baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	  mtype="COX";
	}
	else
	{
		mtype=type;
	}
	


	lastTopVariable = length(VarFrequencyTable);
	if (lastTopVariable >= nsize/3) lastTopVariable <- nsize/3;
	frma <- paste(baseForm,"~ ");
	enterlist <- vector();
	thrsfreq <- frequencyThreshold*loops;
	toRemove <- vector();
	fistfreq <- 0;
	if (length(vnames)>0)
	{
		for ( i in 1:length(vnames))
		{
			if ((vnames[i] != " ") && (vnames[i] != ""))
			{
				enterlist <- append(enterlist,vnames[i]);
				if ((i<=lastTopVariable)&&(VarFrequencyTable[i] > thrsfreq))  # Only features with a given frequency
				{
					if (fistfreq == 0) 
					{
						thrsfreq <- frequencyThreshold*VarFrequencyTable[i];
						fistfreq <- VarFrequencyTable[i];
					}
					frma <- paste(frma,"+",vnames[i]);				
				}
				else
				{
					toRemove <- append(toRemove,paste(" ",vnames[i]," ",sep=""));
				}
			}
		}
	}
	else
	{
		frma <- paste(frma,"+",1);	
	}
	model <- modelFitting(formula(frma),data,mtype)
	if (inherits(model, "try-error"))
	{
		cat("Warning Bagging Fitting error\n")
		model <- modelFitting(formula(frma),data,"LM")
		if (type=="COX")
		{
			model$coefficients <- model$coefficients[-1];
		}			
	}

#	print(summary(model));
	msize <- length(model$coefficients)
	wts <- 0;
	zvalues <- model$coefficients;
	zwts <- model$coefficients;

	for (i in 1:length(model$coefficients))
	{
		model$coefficients[i] <- 0;
		zvalues[i] = 0;
		zwts[i] = 0;
	}
	
	if (type=="COX") 
	{
		zcol=4;
		csize=0;
	}
	else
	{
		zcol=3;
		csize=1;
	}
	
#	cat("\n");
	avgsize = 0;
	varoutcome <- var(theoutcome);
	varistds <- 0.00001+sqrt((colMeans(data^2)-colMeans(data)^2)/varoutcome);
	totresidual <- rep(0,length(theoutcome));
	for (n in 1:loops)
	{
	
		if ((n %% 10) == 0) cat(".");
#		cat(modelFormulas[n],"\n");
		if (length(toRemove)>0)
		{
			modelFormulas[n] <- paste(modelFormulas[n]," +1",sep="");
			for (rml in 1:length(toRemove))
			{
				modelFormulas[n] <- sub(toRemove[rml]," 1 ",modelFormulas[n]);
			}
		}
#		cat(modelFormulas[n],"\n");
		if (gregexpr(pattern ='~',modelFormulas[n])[1]>0)
		{
			ftmp <- formula(modelFormulas[n]);
		}
		else
		{
			ftmp <- formula(paste(baseForm,"~",modelFormulas[n]));
		}
		out <- modelFitting(ftmp,data,mtype);
		if (backElimination)
		{
			if (type=="LM")
			{
				out <- backVarElimination_Res(out,pvalue=pvalue,Outcome=Outcome,data=data,type=type,testType="Ftest")$back.model;				
			}
			else
			{
				out <- backVarElimination_Bin(out,pvalue=pvalue,Outcome=Outcome,data=data,type=type,selectionType="zIDI")$back.model;
			}
		}
		if (!inherits(out, "try-error")) 
		{
			class(model) <- class(out)		
			osize <- length(out$coefficients)
			avgsize = avgsize+osize;
			
			if (osize>csize)
			{				
				outcoef  <- summary(out)$coefficients;
				curprediction <- predictForFresa(out,data,predtype)
				residual <- abs(curprediction-theoutcome);
				varresidual <- mean(residual^2);
				w = (varoutcome-varresidual)/varoutcome; # we will prefer models with small residuals
				if (predtype=='prob') 
				{
					w = 2.0*(mean((curprediction>0.5) == theoutcome) - 0.5); # we will use accuracy for binary outcomes
				}
				if (!is.na(w))
				{
					if (w>0) 
					{
						w = w-(1.0-w)*(osize-1.0)/(nsize-osize); # adjusting for model size.
						if (w<0.01) w= 0.01; # we will keep the minumum weight to 0.01
						if (type == "COX") 
						{
							cs <- abs(out$coefficients*varistds[names(out$coefficients)]) # cs has the ratio of coefficient vs the outcome variance
						}
						else
						{
							cs <- abs(out$coefficients[-1]*varistds[names(out$coefficients[-1])])
						}
						if (predtype=='prob')
						{
							cs <- mean(cs*(cs>=2.0)+(cs<2.0),na.rm=TRUE); # The mean of coeffcients that are larger than 2 times the expected logistic range
						}
						else
						{
							cs <- mean(cs*(cs>=1.0/osize)+(cs<1.0/osize),na.rm=TRUE); # The mean of coeffcients that are larger than the expected linear range
						}
	#					print(varistds[names(out$coefficients)]);
						w <- w/cs; # let us penalize models that have large coefficients
						if (!is.na(w))
						{
							totresidual <- totresidual+w*residual;
							wts = wts + w;
	#						cat(w,", ",cs,", \n");

							for (i in 1:msize)
							{
								for (j in 1:osize)
								{
									if (names(model$coefficients)[i] == names(out$coefficients)[j])
									{
										model$coefficients[i] <- model$coefficients[i] + w*out$coefficients[j]; 
										zvalues[i] <- zvalues[i] + outcoef[j,zcol];
										zwts[i] = zwts[i] + 1;
									}
								}
							}
						}
					}
				}
			}
		}
		else
		{
			cat("+");
		}
	}
	avgsize = avgsize/loops;
	if( wts>0)
	{
		model$coefficients <- model$coefficients/wts;
	}
	for (i in 1:msize)
	{
		if( zwts[i]>0)
		{
			zvalues[i] <- zvalues[i]/zwts[i];
		}
	}
#	print(summary(model));
	baggedPredict <- predictForFresa(model,data,"linear");
	ndata <- as.data.frame(cbind(data[,Outcome],baggedPredict));
	if (predtype == 'linear') # for linear predictior we recalibrate the coefficients for minimum RMSE
	{
		out <- modelFitting(formula("V1~baggedPredict"),ndata,"LM");
		model$coefficients <- model$coefficients*out$coefficients[2];
		model$coefficients[1] <- out$coefficients[1]+model$coefficients[1];
	}
	else
	{
		if (type == 'LOGIT') # for logit predict we set the offset to maximum AUC
		{
			ndata <- as.data.frame(cbind(data[,Outcome],baggedPredict));
			ndata <- ndata[order(baggedPredict),];
			ncases <- sum(data[,Outcome]);
			ncontrol <- nsize-ncases;
			mauc <- 0;
			sen <- 1;
			spe <- 0;
			threshold <- 1;
			for (i in 1:nsize)
			{	
				event <- 1.0*(ndata[i,1]==1);
				nevent <- 1.0*(event==0);
				sen <- sen - event/ncases;
				spe <- spe + nevent/ncontrol;
				auc <- (sen+spe)/2;
				if (auc>mauc)
				{
					mauc <- auc;
					threshold <- i;
				}
			}
			cat ("AUC: ",mauc,"THR: ",threshold/ncases);
			if (threshold<nsize)
			{
				indx <- threshold;
				preindx <- threshold-1;
				postindx <- threshold+1;
				if (preindx<1) preindx=1;
				if (postindx>nsize) postindx=nsize;
				model$coefficients[1] <- model$coefficients[1]-(ndata[preindx,2]+2.0*ndata[indx,2]+ndata[postindx,2])/4;
			}
		}
	}

	cat(" Avg Size: ",avgsize,"\n");

#	print(summary(model));
	reducededset <- data;
	baggedPredict <- predictForFresa(model,data,predtype);
	residual <- baggedPredict-theoutcome;
	mado <- mean(abs(residual));
	rnames <- rownames(data);
	if( wts>0)
	{
		totresidual <- totresidual/wts+abs(residual);
	}

	if ((removeOutliers>0) && (type=="LM") && (binoutcome==FALSE ))
	{
		made <- mean(abs(totresidual),trim = 0.10); # trim @ 10% to remove outliers from estimation of mean MAD
#		cat("Bagged MAD: ", made,"\n");
		gooddatapoints <- abs(totresidual) < removeOutliers*made;
		if (sum(gooddatapoints)!=nrow(data))
		{
			reducededset <- data[gooddatapoints,];
			upmodel <- baggedModel(modelFormulas,reducededset,type,Outcome,timeOutcome,pvalue,backElimination,frequencyThreshold,removeOutliers=0);
			if (!any(is.na(upmodel$bagged.model$coefficients)))
			{
				if (length(model$coefficients) == length(upmodel$bagged.model$coefficients))
				{
					model$coefficients <- upmodel$bagged.model$coefficients;
				}
				else
				{
					for (n in 1:length(upmodel$bagged.model$coefficients))
					{
						model$coefficients[n] <- upmodel$bagged.model$coefficients[n];
					}
				}
				frma <- upmodel$formula;
				zvalues <- upmodel$zvalues;
				baggedPredict <- predictForFresa(model,data,predtype);
				residual <- baggedPredict-theoutcome;
				made <- mean(abs(residual));
				if (!any(is.na(residual)))
				{
					residual <- abs(residual/made);
					residual <- residual*(residual<10)+10*(residual>10);
					plot(residual~theoutcome,main="Outliers");
					for (i in 1:length(gooddatapoints))
					{
						if (!gooddatapoints[i])
						{
							text(theoutcome[i],residual[i],rnames[i],pos=1,cex=0.7);
						}
					}
				}
				else
				{
					cat("Warning NA residual in Bagging\n");
					print(residual);
					print(summary(model))
				}
				cat("Start MAD:", mado," Reduced MAD:",upmodel$MAD," Final MAD:",made,", Removed: ",rownames(data[!gooddatapoints,]),"\n")
				mado <- made;
			}
			else
			{
				cat("Warning NA coefficents in Bagging\n");
			}
		}
	}
	
  	result <- list(bagged.model=model,
				   formula=frma,
				   frequencyTable=VarFrequencyTable,
				   averageSize=avgsize,
				   zvalues=zvalues,
				   reducedDataSet= reducededset,
				   MAD=mado
				   );
  
	return (result);
}
