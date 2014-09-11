timeSerieAnalysis <-
function(varList,baseModel,dataframe,timevar="time",contime=".",Outcome=".",...,description=".",Ptoshow = c(1),plegend= c("p"),timesign="-",catgo.names=c("Control", "Case"))
{
# it will do the time analysis for all variables in the var list. It will create a data frame with the results of the analysis
# the analysis will  be done using GLS (nlme package) 

if (!require(nlme)) {
  install.packages("nlme", dependencies = TRUE)
  library(nlme)
}

	frm1 <- NULL

	vnames <- as.vector(varList[,1]);
	if (description != ".")
	{
		plottitles <- as.vector(varList[,description]);
	}
	else
	{
		plottitles <- vnames;
	}
	sigma <- vector();
	size = length(vnames);
	anames <- vector();


	if (contime==".")
	{
	  contime=timevar;
	}

	t <- tapply(dataframe[,contime], dataframe[,timevar], mean,na.rm=TRUE)
	
	if (timesign=="-")
	{
		timeorder <- dataframe[order(-dataframe[,contime]),];
	}
	else
	{
		timeorder <- dataframe[order(dataframe[,contime]),];
	}

	if (Outcome != ".")
	{
		timeorderCases <- subset(timeorder,get(Outcome) == 1);
		timeorderControl <- subset(timeorder,get(Outcome) == 0);

		tCase <- tapply(timeorderCases[,contime], timeorderCases[,timevar], mean,na.rm=TRUE)
		tControl <- tapply(timeorderControl[,contime], timeorderControl[,timevar], mean,na.rm=TRUE)
	}


	for (j in 1:size)
	{
    
		frm1 <- paste(vnames[j]," ~ ",baseModel);
		cat(frm1,"\n");
		
		obj_s <- eval(parse(text=paste("try(gls(formula(",frm1,"),dataframe,na.action=na.exclude,...))")))
		if ( inherits(obj_s, "try-error"))
		{
			cat("Getting the gls without parameters\n")
			obj_s <- eval(parse(text=paste("try(gls(formula(",frm1,"),dataframe,na.action=na.exclude))")))
		}

		predCases <- NULL
		predControl <- NULL
		reg <- NULL
		if ( !inherits(obj_s, "try-error"))
		{
			reg <- summary(obj_s);	
			if (Outcome != ".")
			{
				predCases <- predict(obj_s,timeorderCases,na.action=na.exclude)
				predControl <- predict(obj_s,timeorderControl,na.action=na.exclude)
			}
		}


		mval <- tapply(dataframe[,vnames[j]], dataframe[,timevar], mean,na.rm=TRUE)
		sdval <- tapply(dataframe[,vnames[j]], dataframe[,timevar], sd,na.rm=TRUE)
		size <- tapply(dataframe[,vnames[j]], dataframe[,timevar], length)
		delta <- sdval / sqrt( size)
		miny = min(mval-1.5*sdval);
		maxy = max(mval+sdval);
		
		maxtime = max(t);
		mintime = min(t);
		deltatime = maxtime-mintime;
			
		
		if (Outcome != ".")
		{

		  
      
			if (timesign=="-")
			{
				mval <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], mean,na.rm=TRUE)
				sdval <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], sd,na.rm=TRUE)
				size <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], length)
				delta <- sdval / sqrt( size)
				errbar(-tCase,mval,mval-delta,mval+delta,col="red",ylab=vnames[j],xlab="t",type="b",ylim=c(miny,maxy))
				lines(-timeorderCases[,contime],predCases,col="red",lty=3)
      
				mval <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], mean,na.rm=TRUE)
				sdval <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], sd,na.rm=TRUE)
				size <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], length)
				delta <- sdval / sqrt( size)
				errbar(-tControl,mval,mval-delta,mval+delta,add=TRUE,col="blue",type="b")
      
				lines(-timeorderControl[,contime],predControl,col="blue",lty=2)
				legend(-maxtime+0.1*deltatime,miny + 0.2*(maxy-miny), catgo.names, col=c("blue","red"), lty = 2:3,bty="n")
			
				if (!is.null(reg))
				{
					for (lg in 1:length(Ptoshow))
					{
						legend(-maxtime+0.65*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),sprintf("p(%s) = %5.3f",plegend[lg],reg$tTable[Ptoshow[lg],4]),bty="n");
						if (reg$tTable[Ptoshow[lg],4]<0.0001)
						{
							legend(-maxtime+0.9*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),"****",bty="n");
						}
						else
						{
							if (reg$tTable[Ptoshow[lg],4]<0.001)
							{
								legend(-maxtime+0.9*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),"***",bty="n");
							}
							else
							{
								if (reg$tTable[Ptoshow[lg],4]<0.05)
								{
									legend(-maxtime+0.9*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),"**",bty="n");
								}
								else
								{
									if (reg$tTable[Ptoshow[lg],4]<0.05)
									{
										legend(-maxtime+0.9*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),"*",bty="n");
									}
								}
							}
						}
					}
				}
			}
			else
			{
				mval <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], mean,na.rm=TRUE)
				sdval <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], sd,na.rm=TRUE)
				size <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], length)
				delta <- sdval / sqrt( size)

				errbar(tCase,mval,mval-delta,mval+delta,col="red",ylab=vnames[j],xlab="t",type="b",ylim=c(miny,maxy))
				lines(timeorderCases[,contime],predCases,col="red",lty=3)

				mval <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], mean,na.rm=TRUE)
				sdval <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], sd,na.rm=TRUE)
				size <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], length)
				delta <- sdval / sqrt( size)

				errbar(tControl,mval,mval-delta,mval+delta,add=TRUE,col="blue",type="b")
      
				lines(timeorderControl[,contime],predControl,col="blue",lty=2)
				legend(maxtime-0.9*deltatime,miny + 0.2*(maxy-miny), catgo.names, col=c("blue","red"), lty = 2:3,bty="n")
		
				if (!is.null(reg))
				{
					for (lg in 1:length(Ptoshow))
					{
						legend(maxtime-0.45*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),sprintf("p(%s) = %5.3f",plegend[lg],reg$tTable[Ptoshow[lg],4]),bty="n");
						if (reg$tTable[Ptoshow[lg],4]<0.0001)
						{
							legend(maxtime-0.1*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),"****",bty="n");
						}
						else
						{
							if (reg$tTable[Ptoshow[lg],4]<0.001)
							{
								legend(maxtime-0.1*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),"***",bty="n");
							}
							else
							{
								if (reg$tTable[Ptoshow[lg],4]<0.05)
								{
									legend(maxtime-0.1*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),"**",bty="n");
								}
								else
								{
									if (reg$tTable[Ptoshow[lg],4]<0.05)
									{
										legend(maxtime-0.1*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),"*",bty="n");
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
			if (timesign=="-")
			{			
			  errbar(-t,mval,mval-delta,mval+delta,ylab=vnames[j],xlab="t",type="b",ylim=c(miny,maxy))
			  lines(-timeorder[,contime],predict(obj_s,timeorder,na.action=na.exclude),col="blue",lty=2)      
			  for (lg in 1:length(Ptoshow))
				{
					legend(-maxtime+0.65*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),sprintf("p(%s) = %5.3f",plegend[lg],reg$tTable[Ptoshow[lg],4]),bty="n");
				}
			}
			else
			{
				errbar(t,mval,mval-delta,mval+delta,ylab=vnames[j],xlab="t",type="b",ylim=c(miny,maxy))
				lines(timeorder[,contime],predict(obj_s,timeorder,na.action=na.exclude),col="blue",lty=2)      
				for (lg in 1:length(Ptoshow))
				{
					legend(maxtime-0.65*deltatime,miny + (0.3-0.1*lg)*(maxy-miny),sprintf("p(%s) = %5.3f",plegend[lg],reg$tTable[Ptoshow[lg],4]),bty="n");
				}
			}
		}
		
		title(main=plottitles[j])
		if (!is.null(reg))
		{
			if (j>1) 
			{
				coeff <- rbind (coeff,as.vector(reg$tTable[,1]));
				std.Error <- rbind (std.Error,as.vector(reg$tTable[,2]));
				t.value <- rbind (t.value,as.vector(reg$tTable[,3]));
				p.value <- rbind (p.value,as.vector(reg$tTable[,4]));
			}
			else
			{
				coeff <- rbind (as.vector(reg$tTable[,1]));
				std.Error <- rbind (as.vector(reg$tTable[,2]));
				t.value <- rbind (as.vector(reg$tTable[,3]));
				p.value <- rbind (as.vector(reg$tTable[,4]));
			}
			sigma <- append(sigma,reg$sigma);
			anames <- append(anames,vnames[j]);
		}
		
	}

	rownames(coeff) <- anames;
	rownames(std.Error) <- anames;
	rownames(t.value) <- anames;
	rownames(p.value) <- anames;
	names(sigma) <- rownames(p.value);

	if (!is.null(reg))
	{
		colnames(coeff) <- rownames(reg$tTable);
		colnames(std.Error) <- rownames(reg$tTable);
		colnames(t.value) <- rownames(reg$tTable);
		colnames(p.value) <- rownames(reg$tTable);
	}


	result  <- list(coeff=coeff,
	std.Errors=std.Error,
	t.values=t.value,
	p.values=p.value,
	sigmas=sigma);
	return (result);
}