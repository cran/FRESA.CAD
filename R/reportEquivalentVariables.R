reportEquivalentVariables <-
function (object,pvalue=0.05,dataframe,searchVarList,Outcome="Class", type = c("LOGIT", "LM","COX"),eqFrac=0.9,description=".") 
{
    type <- match.arg(type);
	cthr = abs(qnorm(pvalue));
  
	varsList <- as.list(attr(terms(object),"variables"))
	vnames <- as.vector(searchVarList[,1]);
    
	orgsize <- length(varsList);
    

	outCome = paste(varsList[2]," ~ ");
	startlist = 3;
	auxmodel <- object;
	indexZidi = startlist - 1;
	
	outvarnames <- vector();
	outrelated <- vector();
	outzIDI <- vector();

	orgzIDI <- as.vector(getVarReclassification(object,dataframe,Outcome,type)$z.IDIs);
	print (orgzIDI,digits=3);
	orgzIDI <- eqFrac * orgzIDI;
	print (orgzIDI,digits=3);

	 for ( i in startlist:length(varsList))
	{
		outvarnames <- append(outvarnames,as.character(varsList[i]));
		cat(as.character(varsList[i]),"\n");	
		namelist ="{"; 
		zIDIlist ="{";
		for (j in 1:length(vnames))
		{
			frm1 = outCome;
			for ( n in startlist:length(varsList))
			{
				if (n == i)
				{
					frm1 <- paste(frm1,paste(" + ",vnames[j]));
				}
				else
				{
					frm1 <- paste(frm1,paste(" + ",varsList[n]));
				}
			}
			ftmp <- formula(frm1);
			auxmodel <- modelFitting(ftmp,dataframe,type)
      
			if (orgsize == length(as.list(attr(terms(auxmodel),"variables"))))
			{
				zIDI <- as.vector(getVarReclassification(auxmodel,dataframe,Outcome,type)$z.IDIs);
		##        cat(frm1," : ");
		##        print(zIDI,digits = 2 );
				if ((zIDI[i-indexZidi] > cthr) && ( zIDI[i-indexZidi] > orgzIDI[i-indexZidi] ))
				{
				  namelist <- paste(namelist,vnames[j]);
				  if (description != ".") namelist <- paste(namelist,"[",searchVarList[j,description],"]");
				  zIDIlist <- paste(zIDIlist,sprintf("%5.3f",zIDI[i-indexZidi]));
				}        
			}
		}
		namelist <- paste(namelist,"}");
		zIDIlist <- paste(zIDIlist,"}");
		cat(namelist,"\n\n");
		outrelated <- append(outrelated,namelist);
		outzIDI <- append(outzIDI,zIDIlist);
	}
	result <- cbind(outvarnames);
	result <- cbind(result,outrelated);
	result <- cbind(result,outzIDI);

    return (result)
}
