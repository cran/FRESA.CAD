modelFitting <-
function (model.formula,data,type = c("LOGIT", "LM","COX","SVM"),fitFRESA=TRUE,...) 
{

	type <- match.arg(type);
	if (type=="SVM")
	{
		if (!requireNamespace("e1071", quietly = TRUE)) {
		   install.packages("e1071", dependencies = TRUE)
		} 
	}
	
	if (class(model.formula) == "character")
	{
		model.formula <- formula(model.formula);
	}

	if (!fitFRESA)
	{

		switch(type, 
			LM = 
			{ 
			  fittedModel <- try(lm(model.formula,data=data,na.action=na.exclude));
			},
			LOGIT =
			{
	          fittedModel <- try(glm(model.formula,data=data,na.action=na.exclude,family=binomial(link=logit),...));
			},
			COX =
			{
#			  cat(as.character(model.formula),"\n")
			  fittedModel <- try(survival::coxph(model.formula,data=data,na.action=na.exclude,model=TRUE,...));
			},
			SVM =
			{
			  fittedModel <- try(e1071::svm(model.formula,data=data));
			},
			{
	          fittedModel <- try(glm(model.formula,data=data,na.action=na.exclude,...));
			}
		)
		if (!inherits(fittedModel, "try-error"))
		{
			if (!is.null(fittedModel$coefficients))
			{
				s <- is.na(fittedModel$coefficients) | is.nan(fittedModel$coefficients);
				if (any(s))
				{	
					fittedModel$coefficients[s] <- 0;				
					class(fittedModel) <- c(class(fittedModel),"try-error");
				}
			}
			else
			{
				class(fittedModel) <- c(class(fittedModel),"try-error");
			}
		}
		else
		{
			if (!is.null(fittedModel$coefficients))
			{
				s <- is.na(fittedModel$coefficients) | is.nan(fittedModel$coefficients);
				if (any(s))
				{	
					fittedModel$coefficients[s] <- 0;				
				}
			}
		}
		fittedModel$type <- type;
	}
	else
	{	
		fittedModel <- NULL;
		fittedModel$type <- type;		
		fittedModel$formula <- model.formula;
		varsmod <- all.vars(model.formula);
		if (length(varsmod)>1)
		{
			modelFrame <- data[,varsmod];
			if (type=="SVM")
			{
				fittedModel <- try(e1071::svm(formula=model.formula,data=modelFrame));
			}
			else
			{
				modelMat <- model.matrix(model.formula,modelFrame);
				
				maxterms <- ncol(modelMat)-1;
				if (type=="COX")
				{
					response <- as.matrix(modelFrame[1:2]);
					maxterms <- ncol(modelMat)-2;
				}
				else
				{
					response <- as.matrix(cbind(modelFrame[,1],modelFrame[,1]));
				}
				if (nrow(modelMat)>=maxterms)
				{
					fittedModel <-.Call("modelFittingCpp",response,modelMat,type);
				}
				else
				{
					warning("More features than data rows. Skiping fitting\n");
#					cat("More features than data rows. Skiping fitting\n");
					if (type=="COX")
					{
						fittedModel$coefficients <- numeric(2*ncol(modelMat))
						fittedModel$estimations <- numeric(2*ncol(modelMat))
#						fittedModel$coefficients[(ncol(modelMat)+1):2*ncol(modelMat)] <- colMeans(modelMat);
					}
					else
					{
						fittedModel$coefficients <- numeric(ncol(modelMat))
					}
				}
				class(fittedModel) <- "fitFRESA";
				fittedModel$estimations <- as.vector(fittedModel$coefficients);
				fittedModel$family <- "FRESA.CAD";
				fittedModel$type <- type;
				fittedModel$formula <- model.formula;
				fittedModel$call <- match.call();
				fittedModel$terms <- terms(model.formula);
				fittedModel$model <- modelFrame;
				fittedModel$response <- response;
				if (type=="COX")
				{
					fittedModel$coefficients <- fittedModel$estimations[1:(length(fittedModel$estimations)/2)];
				}
				fittedModel$coefficients <- as.vector(fittedModel$coefficients);
				if (!is.null(fittedModel$coefficients))
				{
					names(fittedModel$coefficients) <- colnames(modelMat);
					s <- is.na(fittedModel$coefficients) | is.nan(fittedModel$coefficients);
					if (any(s))
					{
#						cat("Fitting NA\n");
						fittedModel$coefficients[s] <- 0;
						fittedModel$estimations[s] <- 0;
						class(fittedModel) <- c(class(fittedModel),"try-error");
					}
				}
				else
				{
					class(fittedModel) <- c(class(fittedModel),"try-error");
				}
			}
		}
		else
		{
			if (length(varsmod)>0)
			{
				modelFrame <- data[,varsmod];
				response <- as.matrix(cbind(modelFrame,modelFrame));
				fittedModel$coefficients <- mean(modelFrame,na.rm = TRUE);
				fittedModel$estimations <- c(fittedModel$coefficients,fittedModel$coefficients);
				class(fittedModel) <- "fitFRESA";
				class(fittedModel) <- c(class(fittedModel),"Constant");
				fittedModel$family <- "mean";
			}
			else
			{
			    warning("Warning Zero Length",as.character(model.formula),"\n");
				modelFrame <- data[,1];
				response <- as.matrix(cbind(modelFrame,modelFrame));
				fittedModel$coefficients <- mean(modelFrame,na.rm = TRUE);
				fittedModel$estimations <- c(fittedModel$coefficients,fittedModel$coefficients);
				class(fittedModel) <- "fitFRESA";
				class(fittedModel) <- c(class(fittedModel),"try-error");
				fittedModel$family <- "mean";
			}
			fittedModel$formula <- model.formula;
			fittedModel$call <- match.call();
			fittedModel$terms <- terms(model.formula);
			fittedModel$model <- modelFrame;
			fittedModel$response <- response;
		}
	}
    return (fittedModel)
}
