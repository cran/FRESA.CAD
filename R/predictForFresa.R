predictForFresa <-
function (object,newdata, type = c("prob", "linear")) 
{

	frm <- formula(terms(object));
	classlen=length(class(object))
	
	cobj <- substr(class(object)[classlen], 1, 2);
	switch(cobj,
		co =
		{
			pobj <- object;
			switch(type, 
				linear = 
					{		
						out <- predict(pobj,newdata=newdata,type = 'lp');
					},
				prob = 
					{
						out <- 1.0/(1.0+exp(-predict(pobj,newdata=newdata,type = 'lp')));
					},
					{
						out <- predict(pobj,newdata=newdata,type = 'lp');
					}
			)
		},
		{
			type <- match.arg(type)
			cf <- coef(object)
			s <- is.na(cf);
			if (any(s)) 
			{
				cf[s] <- 0;
			}
			
			switch(type, 
				linear = 
					{
					   out <- as.vector(model.matrix(frm, newdata)  %*% cf);
					}, 
				prob = 
					{
					  out <- as.vector(model.matrix(frm, newdata)  %*% cf);
					  out <- 1.0/(1.0+exp(-out));
					}, 
					{
					  out <- as.vector(model.matrix(frm, newdata)  %*% cf);
					}
			)
		}
	)
	s <- is.na(out);
	if (any(s)) 
	{
		cat("Warning NA predictForFresa \n");
		out[s] <- 0;
	}
    return (out)
}
