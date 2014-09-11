modelFitting <-
function (model,dataframe,type = c("LOGIT", "LM","COX"),...) 
{

if (!require(speedglm)) {
	speedglm <- function(...) {}
  install.packages("speedglm", dependencies = TRUE)
  library(speedglm)
}

    type <- match.arg(type);
	switch(type, 
	    LM = 
        { 
          fittedModel <- try(lm(model,data=dataframe,na.action=na.exclude));
         },
        LOGIT =
        {
#          fittedModel <- try(glm(model,data=dataframe,na.action=na.exclude,family=binomial(link=logit),...));
		  fittedModel <- try(speedglm(model,data=dataframe,na.action=na.exclude,family=binomial(link=logit)));
        },
		COX =
		{
		  fittedModel <- try(coxph(model,data=dataframe,,na.action=na.exclude,model=TRUE,...));
		},
        {
#          fittedModel <- try(glm(model,data=dataframe,na.action=na.exclude,...));
		  fittedModel <- try(speedglm(model,data=dataframe,na.action=na.exclude,family=binomial(link=logit)));
        }
	)
	if (!inherits(fittedModel, "try-error"))
	{
		s <- is.na(fittedModel$coefficients);
		if (!any(s))
		{		
			if ((max(fittedModel$coefficients)>1.0e10) || (min(fittedModel$coefficients)< -1.0e10))
			{
				class(fittedModel) <- c(class(fittedModel),"try-error");
			}
		}
		else
		{
			class(fittedModel) <- c(class(fittedModel),"try-error");
		}
	}

    return (fittedModel)
}
