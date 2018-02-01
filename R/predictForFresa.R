#' @method predict fitFRESA
predict.fitFRESA <-
function (object,...) 
{

	impute=FALSE;
	parameters <- list(...);
	if (length(parameters)==3)
	{
		if (!is.null(parameters$impute))
		{
			impute=parameters$impute;
		}
	}
	if (length(parameters)==2)
	{
		testData <- parameters[[1]];
		predictType <- parameters[[2]];
	}
	else
	{
		if (length(parameters)==1)
		{
			testData <- parameters[[1]];
			predictType <- "linear";
		}
		else
		{
			if (is.null(parameters$testData))
			{
				testData <-object$model;
			}
			else
			{
				testData <- parameters$testData
			}
			if (is.null(parameters$predictType))
			{
				predictType <- "linear";
			}
			else
			{
				predictType <- parameters$predictType
			}
		}
	}
	
	frm <- formula(terms(object));
#	print(as.character(frm))
#	cat(as.character(frm),":",nrow(testData),":",predictType,":",class(object),"\n");

	classlen=length(class(object))
	
	cobj <- substr(class(object)[classlen], 1, 2);
	switch(cobj,
		co =
		{
			pobj <- object;
			switch(predictType[1], 
				linear = 
					{		
						out <- predict(pobj,testData,type = 'lp');
					},
				prob = 
					{
						out <- 1.0/(1.0+exp(-predict(pobj,testData,type = 'lp')));
					},
					{
						out <- predict(pobj,testData,type = 'lp');
					}
			)
		},
		fi =
		{
			cf <- object$estimations;
			if (object$type=="COX")
			{			
				mm <- model.matrix(frm, testData)
			}
			else
			{
				mm <- as.matrix(model.frame(frm, testData,na.action=na.pass))
				mm[,1] <- rep(1,nrow(mm));
			}
			if ((any(is.na(mm))) && impute)
			{
				mm[,-1]=nearestneighborimpute(tobeimputed=mm[,-1],referenceSet = object$model[,-1]);
			}
			pred <-.Call("predictForFresa",cf,mm,predictType[1],object$type);
			out <- as.vector(pred$prediction);
			names(out) <- rownames(testData);
		},
		tr =
		{
			cf <- object$estimations;
			if (object$type=="COX")
			{			
				mm <- model.matrix(frm, testData)
			}
			else
			{
				mm <- as.matrix(model.frame(frm, testData,na.action=na.pass))
				mm[,1] <- rep(1,nrow(mm));
			}
			if (!is.null(cf))
			{
				warning("Warning: Fitting error. Object: ",class(object),"NA coff set to zero\n");
#				print(cf);
				pred <-.Call("predictForFresa",cf,mm,predictType[1],object$type);
				out <- as.vector(pred$prediction);
			}
			else
			{
				warning("Warning: Fitting error. Object: ",class(object),"All predictions set to NA \n");
				out <- rep(NA,nrow(testData));
			}
			names(out) <- rownames(testData);
		},
		sv =
		{
			out <- predict(object,testData);
		},
		{
			cf <- coef(object)
			if (is.null(cf))
			{
				warning("Warning: Null object in predict. Object: ",cobj,"\n");
				out <- rep(NA,nrow(testData));
			}
			else
			{
				cf <- object$coefficients;
				s <- is.na(cf);
				if (any(s)) 
				{
					cf[s] <- 0;
				}
				if (object$type=="COX")
				{			
					mm <- model.matrix(frm, testData)
				}
				else
				{
					mm <- as.matrix(model.frame(frm, testData,na.action=na.pass))
					mm[,1] <- rep(1,nrow(mm));
				}
				
				switch(predictType[1], 
					linear = 
						{
						   out <- as.vector(mm  %*% cf);
						}, 
					prob = 
						{
						  out <- as.vector(mm  %*% cf);
						  out <- 1.0/(1.0+exp(-out));
						}, 
						{
						  out <- as.vector(mm  %*% cf);
						}
				)
			}
			names(out) <- rownames(testData);
		}
	)
	s <- is.na(out);
	if (any(s)) 
	{
		warning("Warning NA predict.fitFRESA \n");
	}
	if (length(out)!=nrow(testData))
	{
		warning("Different number of rows:",length(out),"(",nrow(testData),"). Setting to NA missing values\n");
		tout <- out;
		out=rep(NA,nrow(testData));
		names(out) <- rownames(testData);
		out[names(tout)] <- tout;
	}
    return (out)
}
