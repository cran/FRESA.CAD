rankInverseNormalDataFrame <-
function(varList,dataframe,referenceframe) 
{
  
  rankInverseNormal <- function(obs,vname,min_value,max_value,dataframe) 
  {
    size=length(dataframe[,vname]);
    mrank=0;
    
    if (obs<=min_value) {mrank=0.0001;}
    else
    {
      if (obs>=max_value) {mrank=0.9999;}
      else
      {
        i=1;
        while ((obs >= dataframe[i,vname]) & (i <= size))
        {
          i=i+1; 
        }
        i = i-1;
        if (i==0) 
        {
          mrank= (obs-min_value)/(dataframe[1,vname]-min_value)/(1.0+size);
        }
        else
          if (i == size)
          {
            mrank= (size+(obs-dataframe[size,vname])/(max_value-dataframe[size,vname]))/(1.0+size);
          }
        else
        {
          mrank= (i+(obs-dataframe[i,vname])/(dataframe[i+1,vname]-dataframe[i,vname]))/(1.0+size);
        }
      }
    }
    if (mrank==0) mrank=0.0001;
    if (mrank==1) mrank=0.9999;
    
    return (qnorm(mrank,0,1));
  }
  
  
	colnamesList <- as.vector(varList[,1]);
	size = length(colnamesList);
	SortedCtr = referenceframe;
	nrowsCtr =  nrow(referenceframe);
	nrows = nrow(dataframe);
	zRankInverseFrame = dataframe;
	for (i in 1:size)
	{ 
		if (length(table(dataframe[,colnamesList[i]]))>2)
		{

			SortedCtr[,colnamesList[i]]<-SortedCtr[order(referenceframe[,colnamesList[i]]),colnamesList[i]];
			std <- sd (SortedCtr[,colnamesList[i]],na.rm = TRUE);
			minvalue <- SortedCtr[1,colnamesList[i]] - 2*std;
			maxvalue <- SortedCtr[nrowsCtr,colnamesList[i]] + 2*std;
			cat(" Variable: ",colnamesList[i],"Min: ",minvalue," Max: ",maxvalue,"\n")
			for (j in 1:nrows)
			{
				zRankInverseFrame[j,colnamesList[i]]=rankInverseNormal(dataframe[j,colnamesList[i]],colnamesList[i],minvalue,maxvalue,SortedCtr);
			}
		}
	}
	return (zRankInverseFrame);
}
