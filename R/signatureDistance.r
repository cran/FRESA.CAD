signatureDistance <- 
function (template, data=NULL, method = c("pearson","spearman","kendall","RMS","MAN"))
{

#given the template: mean,median,sample, etc....;signatureDistance it will return the distance between the template to each row of the dataframe
#the template is a named numeric vector
#the data is a colnamed data frame
#methods:
# RMS: Root Mean Square
# MAN: (Manhattan distance)/IQR
# pearson: 1-Pearson correlation coefficient
# spearman: 1-spearman correlation coefficient
# kendall: 1-kendall correlation coefficient

	method <- match.arg(method)
	theProbs <- c(0.025,0.16, 0.5,0.84, 0.975);
	wvalues <- 1.0/abs(qnorm(theProbs))

	
	if (class(template)=="matrix")
	{
		vnames <- colnames(template)
	}
	else
	{
		vnames <- names(template)
	}
	datasubset <- as.matrix(data[,vnames]);
	switch(method, 
		RMS = 
		{ 
			if (class(template)=="matrix")
			{
				tem <- template[3,];
				ld <- tem-0.5*(wvalues[1]*template[1,]+wvalues[2]*template[2,]);
				mv <- min(ld[ld>0]);
				if (is.numeric(mv))	{ ld[ld==0] <- mv }
				else { ld[ld==0] <- 1.0; }
				ud <- 0.5*(wvalues[4]*template[4,]+wvalues[5]*template[5,])-tem;
				mv <- min(ud[ud>0]);
				if (is.numeric(mv))	{ ud[ud==0] <- mv }
				else { ud[ud==0] <- 1.0; }
			}
			else
			{
				tem <- template;
				ld <- sd(template);
				ud <- ld;
			}
			RMSDistance <- function (x,template,ld,ud) 
			{
				md <- x-template
				md <- sqrt(mean(pmax(md/ud,-md/ld)^2,na.rm=TRUE));
				return (md)
			}
			metric <- apply(datasubset,1,RMSDistance,tem,ld,ud);
		},
		MAN = 
		{ 
			if (class(template)=="matrix")
			{
				tem <- template[3,];
				ld <- tem-0.5*(wvalues[1]*template[1,]+wvalues[2]*template[2,]);
				mv <- min(ld[ld>0]);
				if (is.numeric(mv))	{ ld[ld==0] <- mv }
				else { ld[ld==0] <- 1.0; }
				ud <- 0.5*(wvalues[4]*template[4,]+wvalues[5]*template[5,])-tem;
				mv <- min(ud[ud>0]);
				if (is.numeric(mv))	{ ud[ud==0] <- mv }
				else { ud[ud==0] <- 1.0; }
			}
			else
			{
				tem <- template;
				ld <- sd(template);
				ud <- ld;
			}
			manDistance <- function (x,template,ld,ud) 
			{
				md <- x-template
				md <- mean(pmax(md/ud,-md/ld),na.rm=TRUE);
				return (md)
			}
			metric <- apply(datasubset,1,manDistance,tem,ld,ud);
	  },
		{
			if (class(template)=="matrix")
			{
				tem <- (0.25*template[2,]+0.5*template[3,]+0.25*template[4,]);
			}
			else
			{
				tem <- template;
			}
			corDistance <- function (x,template,method) {md <- 1.0-cor(x,template,method=method,use="pairwise.complete.obs"); return (md)}
			metric <- apply(datasubset,1,corDistance,template=tem,method=method);
		}
	)
	names(metric) <- rownames(data);
	metric[is.na(metric)] <- 1.0e10;
	
  result <- metric
	return (result);
}
