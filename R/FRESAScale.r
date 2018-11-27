FRESAScale <- function(data,refFrame=NULL,method=c("Norm","Order","RankInv"),refMean=NULL,refDisp=NULL)
{
	usedFeatures <- colnames(data);

	outs <- sapply(data,table);
	outl <- numeric(length(outs));
	for (i in 1:length(outs)) outl[i] <- length(outs[[i]]);
	usedFeatures <- usedFeatures[outl > 10];
	if (is.null(refFrame))
	{
		refFrame <- data;
	}
	method <- match.arg(method);
	switch(method,
			Norm =
			{
				if (is.null(refMean))
				{
					refMean <- apply(refFrame,2,mean, na.rm = TRUE);
					refDisp <- apply(refFrame,2,sd, na.rm = TRUE);
					refDisp[refDisp == 0] <- 1.0;
				}

				meanmat <- matrix(rep(refMean,nrow(data)),nrow=nrow(data),ncol=ncol(data),byrow = TRUE);
				sdmat <- matrix(rep(refDisp,nrow(data)),nrow=nrow(data),ncol=ncol(data),byrow = TRUE);
				data <- (data-meanmat)/sdmat;
			},
			Order =
			{
				if (is.null(refMean))
				{
					refSD <- apply(refFrame,2,sd, na.rm = TRUE);
					refSD[refSD == 0] <- 1.0;
					refMean <- apply(refFrame,2,median, na.rm = TRUE);
					refDisp <- apply(refFrame,2,IQR, na.rm = TRUE);
					refDisp[refDisp == 0] <- refSD[refDisp == 0];
				}

				meanmat <- matrix(rep(refMean,nrow(data)),nrow=nrow(data),ncol=ncol(data),byrow = TRUE);
				sdmat <- matrix(rep(refDisp,nrow(data)),nrow=nrow(data),ncol=ncol(data),byrow = TRUE);
				data <- (data-meanmat)/sdmat;
			},
			RankInv =
			{
				data <- rankInverseNormalDataFrame(usedFeatures,data,refFrame); 			
			}
		)
	result <- list(scaledData=data,refMean=refMean,refDisp=refDisp);
	return (result);
}
