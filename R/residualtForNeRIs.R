residualForNeRIs <-
function (object,newdata,outcome,eta=0.05) 
{
#Set eta to zero to get prediction residuals for cox models. eta to 1 get the Martingale residuals
	classlen=length(class(object))
	
	cobj <- substr(class(object)[classlen], 1, 2);
	switch(cobj,
		co =
		{
			if (classlen==1)
			{
				lpp <- 1.0/(1.0+exp(-predict(object,newdata=newdata,type = 'lp',na.action=na.omit)));
				s <- is.na(lpp);
				if (any(s)) 
				{
					lpp[s] = 10; # set to a large residual
				}
				lppres <- lpp - newdata[,outcome];
				matingale <- newdata[,outcome]-predict(object,newdata=newdata,type = 'expected',na.action=na.omit);
				s <- is.na(matingale);
				if (any(s))
				{
					matingale[s] = 100; # set to a large residual
				}
				out <- (1-eta)*lppres - eta*matingale;			
			}
			else
			{
				out <- (1 - 2*newdata[,outcome]);
			}
		},
		lm =
		{
			out <- predictForFresa(object,newdata=newdata,type = 'linear') - newdata[,outcome];
		},
		{
			if (object$family[1] == "binomial")
			{
				out <- predictForFresa(object,newdata=newdata,type = 'prob') - newdata[,outcome];
			}
			else
			{
				out <- predictForFresa(object,newdata=newdata,type = 'linear') - newdata[,outcome];
			}
#			out <- predict(object,newdata=newdata,type='response',na.action=na.omit) - newdata[,outcome];
		}
	)	
	s <- is.na(out);
	if (any(s)) 
	{
		cat("Warning NA predictFor NeRIs \n");
		out[s] <- 1.0e10; # Set a large residual
	}
    return (out)
}
