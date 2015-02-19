#include "FRESAcommons.h"

#define EPS 1e-4		

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()


#define MAX_TIES 1000


extern "C" SEXP modelFittingCpp(SEXP _ymat,SEXP _X,SEXP _type)
{
	std::string type = Rcpp::as<std::string>(_type);
	Rcpp::NumericMatrix rymat(_ymat);
	Rcpp::NumericMatrix rX(_X);
	mat ymat(rymat.begin(), rymat.rows(), rymat.cols(),false);
	mat X(rX.begin(), rX.rows(), rX.cols(),false);
	vec betas=modelFittingFunc(ymat,X,type);
	vec linearPredictors = predictForFresaFunc(betas,X,"linear",type);
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("coefficients")=Rcpp::wrap(betas.t()),
										   Rcpp::Named("linear.predictors")=Rcpp::wrap(linearPredictors)); 
	return result;
}

extern "C" SEXP improvedResidualsCpp(SEXP _oldResiduals,SEXP _newResiduals,SEXP _testType)
{
	std::string testType = Rcpp::as<std::string>(_testType);
	Rcpp::NumericVector roldResiduals(_oldResiduals);
	Rcpp::NumericVector rnewResiduals(_newResiduals);
	vec oldResiduals(roldResiduals.begin(), roldResiduals.size(),false);
	vec newResiduals(rnewResiduals.begin(), rnewResiduals.size(),false);
	improvedRes impred = improvedResidualsFunc(oldResiduals,newResiduals,testType);

	Rcpp::List result = Rcpp::List::create(Rcpp::Named("p1")=Rcpp::wrap(impred.p1),
										   Rcpp::Named("p2")=Rcpp::wrap(impred.p2),
										   Rcpp::Named("NeRI")=Rcpp::wrap(impred.NeRI),
										   Rcpp::Named("p.value")=Rcpp::wrap(impred.pvalue),
										   Rcpp::Named("BinP.value")=Rcpp::wrap(impred.binom_pValue),
										   Rcpp::Named("WilcoxP.valu")=Rcpp::wrap(impred.wilcox_pValue),
										   Rcpp::Named("tP.value")=Rcpp::wrap(impred.t_test_pValue),
										   Rcpp::Named("FP.value")=Rcpp::wrap(impred.F_test_pValue)
										   );
	return result;
}

extern "C" SEXP improveProbCpp(SEXP _x1,SEXP _x2,SEXP _y)
{
	Rcpp::NumericVector rx1(_x1);
	Rcpp::NumericVector rx2(_x2);
	Rcpp::NumericVector ry(_y);
	vec x1(rx1.begin(), rx1.size(),false);
	vec x2(rx2.begin(), rx2.size(),false);
	vec y(ry.begin(), ry.size(),false);
	
	vec imp=improveProbFunc(x1,x2,y);
	
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("z_idi")=Rcpp::wrap(imp[0]),
										   Rcpp::Named("z_nri")=Rcpp::wrap(imp[1]),
										   Rcpp::Named("idi")=Rcpp::wrap(imp[2]),
										   Rcpp::Named("nri")=Rcpp::wrap(imp[3])
										   );
	return result;
}

extern "C" SEXP predictForFresaCpp(SEXP _cf,SEXP _newdata,SEXP _typ,SEXP _opc)
{
	std::string typ = Rcpp::as<std::string>(_typ);
	std::string opc = Rcpp::as<std::string>(_opc);
	Rcpp::NumericVector rcf(_cf);
	Rcpp::NumericMatrix rnewdata(_newdata);
	vec cf(rcf.begin(), rcf.size(),false);
	mat newdata(rnewdata.begin(), rnewdata.rows(), rnewdata.cols(), false);
	vec prediction = predictForFresaFunc(cf,newdata,typ,opc);
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("prediction")=Rcpp::wrap(prediction));
	return result;
}

                                                          
double qnorm(double p, double mu, double sigma)
{
	return R::qnorm(p, mu, sigma, 1, 0);
}


/*
** This code is based on the chinv2 function from survival package 
** Reference: Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, http://CRAN.R-project.org/package=survival.
*/
void chinv2(mat &matrix , int n)
{
	register double temp;
	register int i,j,k;
	for (i=0; i<n; i++)
	{
		if (matrix(i,i) >0) 
		{
		  	matrix(i,i) = 1/matrix(i,i);  
		 	for (j= (i+1); j<n; j++) 
		 	{
		   		matrix(i,j) = -matrix(i,j);
		   		for (k=0; k<i; k++)   
					matrix(k,j) += matrix(i,j)*matrix(k,i);
		   	}
		}
	}
	for (i=0; i<n; i++) 
	{
		if (matrix(i,i)==0) 
		{ 
			for (j=0; j<i; j++) matrix(i,j)=0;
			for (j=i; j<n; j++) matrix(j,i)=0;
		}
		else 
		{
			for (j=(i+1); j<n; j++) 
			{
				temp = matrix(i,j)*matrix(j,j);
				if (j!=i) matrix(j,i) = temp;
				for (k=i; k<j; k++)
				matrix(k,i) += temp*matrix(k,j);
			}
		}
	}
}

/*
** This code is based on the cholesky2 function from survival package 
** Reference: Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, http://CRAN.R-project.org/package=survival.
*/

int cholesky2(mat &matrix, int n, double toler)
{
    double temp;
    int  i,j,k;
    double eps, pivot;
    int rank;
    int nonneg;
    nonneg=1;
    eps =0;
    for (i=0; i<n; i++) 
    {
		if (matrix(i,i)> eps)  eps = matrix(i,i);
		for (j=(i+1); j<n; j++)  matrix(i,j) = matrix(j,i);
	}
    eps *= toler;

    rank =0;
    for (i=0; i<n; i++) 
    {
		pivot = matrix(i,i);
		if (pivot < eps) 
		{
		    matrix(i,i) =0;
		    if (pivot < -8*eps) nonneg= -1;
	    }
		else  
		{
		    rank++;
		    for (j=(i+1); j<n; j++) 
		    {
				temp = matrix(i,j)/pivot;
				matrix(i,j) = temp;
				matrix(j,j) -= temp*temp*pivot;
				for (k=(j+1); k<n; k++) matrix(j,k) -= temp*matrix(i,k);
			}
	    }
	}
    return(rank * nonneg);
}

/*
** This code is based on the coxfit6 function from survival package 
** Reference: Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, http://CRAN.R-project.org/package=survival.
** calls functions:  cholesky2, chsolve2, chinv2.
*/
void chsolve2(mat &matrix, int n, vec &y)
	{
	register int i,j;
	register double temp;

	for (i=0; i<n; i++) 
	{
		temp = y(i);
		for (j=0; j<i; j++)
		   temp -= y(j) * matrix(j,i);
		y(i) = temp ;
	}
	for (i=(n-1); i>=0; i--) 
	{
		if (matrix(i,i)==0)  y(i) =0;
		else 
		{
			temp = y(i)/matrix(i,i);
			for (j= i+1; j<n; j++)
				temp -= y(j)*matrix(i,j);
			y(i) = temp;
	  	}
	}
}
/*
** This code is based on the coxfit6 function from survival package 
** Reference: Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, http://CRAN.R-project.org/package=survival.
** calls functions:  cholesky2, chsolve2, chinv2.
*/
vec coxfit(int  maxiter,  vec &time,  vec &status, mat &covar, vec &offset, vec &weights, vec &strata,   int method, double eps, double toler,    vec &beta,    int doscale) 
{
    int i,j,k, obs;  
    double  wtave;
    double  denom=0, zbeta, risk;
    double  temp, temp2;
    int     ndead;  
    double  tdeath=0;  
    double  newlk=0;
    double  dtime, d2;
    double  deadwt;  
    double  efronwt;
    int     halving;   
    int     nrisk;   
    int     nobs=offset.n_elem,  nvar  = covar.n_cols;
    vec a(nvar);
    vec a2(nvar);
    vec maxbeta(nvar);
    vec newbeta(nvar);
    vec scale(nvar);
    mat cmat(nvar,nvar);
    mat cmat2(nvar,nvar);
    mat imat(nvar,nvar);
    vec means(nvar);
    vec u(nvar);
    vec loglik(2);
//    double sctest;
//   int flag;
    int iter;
    tdeath=0; temp2=0;
    for (i=0; i<nobs; i++) 
    {
		temp2 += weights(i);
		tdeath += weights(i) * status(i);
    }	
    for (i=0; i<nvar; i++) 
    {
		temp=0;
		for (obs=0; obs<nobs; obs++) 
	    	temp += weights(obs) * covar(obs,i);
		temp /= temp2;
		means(i) = temp;
		for (obs=0; obs<nobs; obs++) covar(obs,i) -=temp;
		if (doscale==1) 
		{ 
	    	temp =0;
	    	for (obs=0; obs<nobs; obs++) 
	    	{    
				temp += weights(obs) * fabs(covar(obs,i));
	    	}
	    	if (temp > 0) temp = temp2/temp;   
	    	else temp=1.0; 
	    	scale(i) = temp;
	    	for (obs=0; obs<nobs; obs++)  covar(obs,i) *= temp;
	    }
	}
    if (doscale==1) 
    {
		for (i=0; i<nvar; i++) beta(i) /= scale(i); 
	}
    else 
    {
		for (i=0; i<nvar; i++) scale(i) = 1.0;
	}
    strata(nobs-1) =1;
    loglik(1) =0;
    for (i=0; i<nvar; i++)
    {
		u(i) =0;
		a2(i) =0;
		for (j=0; j<nvar; j++) 
		{
		    imat(i,j) =0 ;
		    cmat2(i,j) =0;
	    }
	}
    for (obs=nobs-1; obs>=0; ) 
    {
		if (strata(obs) == 1) 
		{
		    nrisk =0 ;  
		    denom = 0;
		    for (i=0; i<nvar; i++) 
		    {
				a(i) = 0;
				for (j=0; j<nvar; j++) cmat(j,i) = 0;
			}
	    }
		dtime = time(obs);
		ndead =0; 
		deadwt =0;  
		efronwt=0;  
		while(obs >=0 &&time(obs)==dtime) 
		{
		    nrisk++;
		    zbeta = offset(obs);    
		    for (i=0; i<nvar; i++)
				zbeta += beta(i)*covar(obs,i);
		    risk = exp(zbeta) * weights(obs);
		    denom += risk;
		    for (i=0; i<nvar; i++) 
		    {
				a(i) += risk*covar(obs,i);
				for (j=0; j<=i; j++)
			    	cmat(j,i) += risk*covar(obs,i)*covar(obs,j);
		    }

		    if (status(obs)==1.0) 
		    {
				ndead++;
				deadwt += weights(obs);
				efronwt += risk;
				loglik(1)+= weights(obs)*zbeta;

				for (i=0; i<nvar; i++) 
				    u(i) += weights(obs)*covar(obs,i);
				if (method==1) 
				{ 
				    for (i=0; i<nvar; i++) 
				    {
						a2(i) +=  risk*covar(obs,i);
						for (j=0; j<=i; j++)
						    cmat2(j,i) += risk*covar(obs,i)*covar(obs,j);
			        }
				}
		    }
		    obs--;
		    if (obs>=0)
		    if (strata(obs)==1) break;  

	    }

		if (ndead >0) 
		{ 
		    if (method==0) 
		    { 
				loglik(1) -= deadwt* log(denom);
		   
				for (i=0; i<nvar; i++) 
				{
				    temp2= a(i)/ denom; 
				    u(i) -=  deadwt* temp2;
				    for (j=0; j<=i; j++)
						imat(i,j) += deadwt*(cmat(j,i) - temp2*a(j))/denom;
			    }
			}
		    else 
		    { 
				for (k=0; k<ndead; k++) 
				{
				    temp = (double)k/ ndead;
				    wtave = deadwt/ndead;
				    d2 = denom - temp*efronwt;
				    loglik(1) -= wtave* log(d2);
				    for (i=0; i<nvar; i++) 
				    {
						temp2 = (a(i) - temp*a2(i))/ d2;
						u(i) -= wtave *temp2;
						for (j=0; j<=i; j++)
						    imat(i,j) +=  (wtave/d2) * ((cmat(j,i) - temp*cmat2(j,i)) - temp2*(a(j)-temp*a2(j)));
				    }
			    }
		
				for (i=0; i<nvar; i++) 
				{
				    a2(i)=0;
				    for (j=0; j<nvar; j++) cmat2(i,j)=0;
		    	}
			}
	    }
	}  
    loglik(0) = loglik(1);
    for (i=0; i<nvar; i++) 
		maxbeta(i) = 20* sqrt(imat(i,i)/tdeath);
    for (i=0; i<nvar; i++) 
		a(i) = u(i);

//    flag= cholesky2(imat, nvar, toler);
    cholesky2(imat, nvar, toler);
    chsolve2(imat,nvar,a);       
    temp=0;
    for (i=0; i<nvar; i++) temp +=  u(i)*a(i);
//    sctest = temp;  
    for (i=0; i<nvar; i++)
    {
		newbeta(i) = beta(i) + a(i);
	}
    if (maxiter==0) 
    {
		chinv2(imat,nvar);
		for (i=0; i<nvar; i++) 
		{
		    beta(i) *= scale(i); 
		    u(i) /= scale(i);
		    imat(i,i) *= scale(i)*scale(i);
		    for (j=0; j<i; j++) 
		    {
				imat(i,j) *= scale(i)*scale(j);
				imat(j,i) = imat(i,j);
			}
    	}
		goto finish;
	}
    halving =0 ;            
    for (iter=1; iter<= maxiter; iter++) 
    {
		newlk =0;
		for (i=0; i<nvar; i++) 
		{
		    u(i) =0;
		    for (j=0; j<nvar; j++)
			imat(j,i) =0;
	    }
		for (obs=nobs-1; obs>=0; ) 
		{
		    if (strata(obs) == 1) 
		    { 
				denom = 0;
				nrisk =0;
				for (i=0; i<nvar; i++) 
				{
				    a(i) = 0;
				    for (j=0; j<nvar; j++) cmat(j,i) = 0;
			    }
			}
		    dtime = time(obs);
		    deadwt =0;
		    ndead =0;
		    efronwt =0;
		    while(obs>=0 && time(obs)==dtime) 
		    {
				nrisk++;
				zbeta = offset(obs);
				for (i=0; i<nvar; i++)
				    zbeta += newbeta(i)*covar(obs,i);
				risk = exp(zbeta) * weights(obs);
				denom += risk;

				for (i=0; i<nvar; i++) 
				{
				    a(i) += risk*covar(obs,i);
				    for (j=0; j<=i; j++)
				    cmat(j,i) += risk*covar(obs,i)*covar(obs,j);
			    }
				if (status(obs)==1) 
				{
				    ndead++;
				    deadwt += weights(obs);
				    newlk += weights(obs) *zbeta;
				    for (i=0; i<nvar; i++) 
						u(i) += weights(obs) *covar(obs,i);
				    if (method==1) 
				    { 
						efronwt += risk;
						for (i=0; i<nvar; i++) 
						{
						    a2(i) +=  risk*covar(obs,i);
						    for (j=0; j<=i; j++)
								cmat2(j,i) += risk*covar(obs,i)*covar(obs,j);
						}   
			        }
		  	    }
				obs--;
				if (obs>=0)
				if (strata(obs)==1) break; 
		    }
		    if (ndead >0) 
		    {  
				if (method==0) 
				{ 
				    newlk -= deadwt* log(denom);
				    for (i=0; i<nvar; i++) 
				    {
						temp2= a(i)/ denom;  
						u(i) -= deadwt* temp2;
						for (j=0; j<=i; j++)
						  	imat(i,j) +=  (deadwt/denom)* (cmat(j,i) - temp2*a(j));
			        }
			    }
				else  
				{ 
				    for (k=0; k<ndead; k++) 
				    {
						temp = (double)k / ndead;
						wtave= deadwt/ ndead;
						d2= denom - temp* efronwt;
						newlk -= wtave* log(d2);
						for (i=0; i<nvar; i++) 
						{
						    temp2 = (a(i) - temp*a2(i))/ d2;
						    u(i) -= wtave*temp2;
						    for (j=0; j<=i; j++)
								imat(i,j) +=  (wtave/d2)*((cmat(j,i) - temp*cmat2(j,i)) - temp2*(a(j)-temp*a2(j)));
			            }
			        }

				    for (i=0; i<nvar; i++) 
				    { 
						a2(i) =0;
						for (j=0; j<nvar; j++) cmat2(j,i) =0;
			        }
		        }
			}
	    }   
//		flag = cholesky2(imat, nvar, toler);
		cholesky2(imat, nvar, toler);
		if (fabs(1-(loglik(1)/newlk))<= eps && halving==0) 
		{ 
		    loglik(1) = newlk;
		    chinv2(imat, nvar);     
		    for (i=0; i<nvar; i++)
		    {
				beta(i) = newbeta(i)*scale(i);
				u(i) /= scale(i);
				imat(i,i) *= scale(i)*scale(i);
				for (j=0; j<i; j++) 
				{
				    imat(i,j) *= scale(i)*scale(j);
				    imat(j,i) = imat(i,j);
			    }
		    }
		    goto finish;
		}
		if (iter== maxiter) break;  
		if (newlk < loglik(1))   
		{   
			halving =1;
			for (i=0; i<nvar; i++)
			    newbeta(i) = (newbeta(i) + beta(i)) /2;
		}
		else 
		{
		    halving=0;
		    loglik(1) = newlk;
		    chsolve2(imat,nvar,u);
		    j=0;
		    for (i=0; i<nvar; i++) 
		    {
				beta(i) = newbeta(i);
				newbeta(i) = newbeta(i) +  u(i);
				if (newbeta(i) > maxbeta(i)) newbeta(i) = maxbeta(i);
				else if (newbeta(i) < -maxbeta(i)) newbeta(i) = -maxbeta(i);
	    	}
	    }
	} 
    loglik(1) = newlk;
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) 
    {
		beta(i) = newbeta(i)*scale(i);
		u(i) /= scale(i);
		imat(i,i) *= scale(i)*scale(i);
		for (j=0; j<i; j++) 
		{
		    imat(i,j) *= scale(i)*scale(j);
		    imat(j,i) = imat(i,j);
	    }
	}
//    flag = 1000;
	finish:
    return join_cols(beta,means);
}

/* 
** This code is based on the logit_link function from 
** file src/library/stats/src/family.c
** Part of the R package, http://www.R-project.org
*/

vec logit_link(vec mu)
{
    int i, n = mu.n_elem;
    vec rans = mu;
    for (i = 0; i < n; i++)
    {
//		rans[i]=log(mu[i]/(1.0 - mu[i]));
		if (mu[i] < 1.0)
		{
			if(mu[i]> 0.0) 
			{
				rans[i]=log(mu[i]/(1.0 - mu[i]));
			}
			else
			{
				rans[i] = -1.0e100;
			}
		}
		else 
		{
			rans[i] = 1.0e100;
		}
    }
    return rans;
}
/* 
** This code is based on the logit_linkinv function from 
** file src/library/stats/src/family.c
** Part of the R package, http://www.R-project.org
*/

vec logit_linkinv(vec eta)
{
    int i, n = eta.n_elem;
    vec rans = eta;
	double expeta;
    for (i = 0; i < n; i++) 
    {
				expeta = exp(eta[i]);
				rans[i] = expeta/(1.0 + expeta); 
    }
    return rans;
}

/* 
** This code is based on the binomial_dev_resids function from 
** file src/library/stats/src/family.c
** Part of the R package, http://www.R-project.org
*/
vec logit_mu_eta(vec eta)
{
    int i, n = eta.n_elem; 
    vec rans = eta;
	double expi,opexp;
    for (i = 0; i < n; i++)
    {
		expi =  exp(eta[i]);
		opexp = 1.0 + expi;
		rans[i] = expi/(opexp * opexp);
    }
	return rans;
}

/* 
** This code is based on the binomial_dev_resids function from 
** file src/library/stats/src/family.c
** Part of the R package, http://www.R-project.org
*/
vec binomial_dev_resids(vec y, vec mu, vec wt)
{
    int i, n = y.n_elem, lmu = mu.n_elem, lwt = wt.n_elem;
    double mui, yi; 
    vec rans = y;
    if(lmu > 1.0) 
    {
    	for (i = 0; i < n; i++)
    	{
		      mui = mu[i];
		      yi = y[i];
		      rans[i] = 2.0 * wt[lwt > 1.0 ? i : 0] * ( ((yi) ? (yi * log(yi/mui)) : 0.0 ) + ((1-yi) ? ( (1-yi) * log((1.0-yi)/(1.0-mui))) : 0.0 ));
		}
    } 
    else 
    {
	  mui = mu[0];
	  for (i = 0; i < n; i++) 
	  {
	    yi = y[i];
		rans[i] = 2.0 * wt[lwt > 1.0 ? i : 0] * ( ((yi) ? (yi * log(yi/mui)) : 0.0 ) + ((1-yi) ? ( (1-yi) * log((1.0-yi)/(1.0-mui))) : 0.0 ));
	  }
    }
    return rans;
}
       
	  
	   
vec modelFittingFunc(mat ymat, mat X,std::string type)
{
  int  nobs = ymat.n_rows;
  mat WX = X;
  vec start = zeros<vec>(X.n_cols); 
    if (X.n_cols==1)  
    {
        if (X.n_elem==0) 
        {
        	return start;
    	}
    }
  vec offset = zeros<vec>(nobs); 
  vec weights = ones<vec>(nobs); 
  vec newstrat = zeros<vec>(nobs);
  vec minmu(nobs);
  for (int i=0; i<nobs;i++) minmu[i]=DOUBLEEPS;
  vec strata;
if (type=="COX")
{
	vec status = ymat.col(1);
    vec stime = ymat.col(0);
    uvec sorted;
    if (strata.n_elem==0) 
    {
       sorted = sort_index(stime);
    }
    else {
          sorted=sort_index(stime);
          sorted=sort_index(strata(sorted));
          strata=strata(sorted);
         }
         for (unsigned int i=0;i<X.n_cols;i++)
         {
           vec sor=X.col(i);
            X.col(i)=sor(sorted);
         }
         stime=stime(sorted);
         status=status(sorted);
    vec betas=coxfit(20,stime,status,X,offset,weights, newstrat,1,1e-9,std::pow(2.220446e-16,0.75),start,1);
    return betas;
}
  vec y = ymat.col(1);
  vec mustart=ones<vec>(nobs);  
  vec mu=ones<vec>(nobs);  
  try{
  vec eta=ones<vec>(nobs); 
  int iter =0;
  int maxit=25;
  double acc=1e-08;
  double  dev=0.0;
  vec  n=ones<vec>(nobs);
  if (type=="LOGIT")
  {
    mustart = (weights % y + 0.5)/(weights + 1.0);
    eta=logit_link(mustart);
    mu = mustart;
    mu = logit_linkinv(eta + offset);
    dev=sum(binomial_dev_resids(y, mu, weights));
  }
 	vec varmu(mu.size());
	 vec mu_eta_val(eta.size());
  if (type=="LM")
  {
    maxit=1;
    mustart=y;
    eta=mustart;
    mu=mustart;
 	dev=sum(weights%(square(y-mu)));  
	varmu.ones();
	mu_eta_val.ones();  
   }
     double tol = 1.0;
	 double dev0;
	 vec z=ones<vec>(nobs);
	 vec W=ones<vec>(nobs);
	 mat XTX;
	 vec XTz;
    while((tol>acc) && (iter<maxit))
	{
	    iter = iter + 1;
	    dev0 = dev;
	   
	    if (type=="LOGIT")
	    {
			varmu = mu % (1.0 - mu);
			mu_eta_val = logit_mu_eta(eta);
			for (int n=0;n<nobs;n++)
			{
				W[n] = (varmu[n]>0) ?  ((weights[n]*mu_eta_val[n]*mu_eta_val[n])/varmu[n]) : 0.0 ;
			}
			z =  W % ((eta - offset) + (y - mu)/mu_eta_val);
			WX = X;
			WX.each_col() %= W;
	    }
		else
		{
			z=y;
		}
		
		XTz = X.t()*z;  
		XTX = X.t()*WX;    
	
		if (!((is_finite(XTX)) && (is_finite(XTz))) )  
		{
			start[0]=NAN;
			return start;
		}
	    try
		{ 
			start = pinv(XTX)*XTz;
		}
		catch(std::exception const& e)
		{
			Rcpp::Rcout<<"Exception::: "<<e.what()<<"\n";
			start = zeros<vec>(X.n_cols);
			start[0]=NAN; // set to NaN
			return start;
		}
		if (iter<maxit)
		{
		    eta=X*start;
		    if (type=="LOGIT")
		    {
		          mu = logit_linkinv(eta + offset);
		          dev=sum(binomial_dev_resids(y, mu, weights));
		    }
		    if (type=="LM")
		    {
		        mu=eta; 
		        dev=sum(weights%(square(y-mu)));
		    }   
			tol = (std::abs(dev0-dev)/(std::abs(dev)+0.1));
		 }
  	}
  return start;
  }
    catch(std::runtime_error const& e)
{
    Rcpp::Rcout<<"Exception: "<<e.what()<<"\n";
	start[0]=NAN; // set to NaN
    return start;
} 
}

vec predictForFresaFunc(vec cf,mat newdata,std::string typ, std::string opc) 
{
	vec out;
	if(opc=="COX")
	{
	    vec coef = cf.subvec(0,newdata.n_cols-1);
	    vec means = cf.subvec(newdata.n_cols,cf.n_elem-1);
		out =newdata*coef- sum(coef%means);
	     if(typ =="prob")  
	    {
	      for(unsigned int i=0;i<out.n_elem;i++)
	       out(i) = 1.0/(1.0+exp(-out(i)));
	    }
	}
	else
	{	 
		out =newdata*cf;  
	    if(typ =="prob")  
		{
		  for(unsigned int i=0;i<out.n_elem;i++)
		 	 out(i) = 1.0/(1.0+exp(-out(i)));
		}
	}
	 	
    return out;
}

/* 
** This code is based on the improveProb function from Hmisc package 
** Reference: Frank E Harrell Jr. URL: http://biostat.mc.vanderbilt.edu/wiki/Main/Hmisc 
** R package version 2.37-7, http://CRAN.R-project.org/package=Hmisc.
*/

vec improveProbFunc(vec x1, vec x2, vec y)
  {
    unsigned int  n = y.n_elem;
    int n_ev=0;
    int n_ne=0;
    vec d  = x2 - x1;
    vec da=d.elem(find(y==1));
    vec db=d.elem(find(y==0));
    int nup_ev=0;
    int nup_ne =0; 
    int ndown_ev =0;
    int ndown_ne=0;
    double sum_ev=0.0;
    double sum_ne=0.0;
    double sse_ev=0.0;
    double sse_ne=0.0;
    for(unsigned int i=0;i<n;i++)
    {
    	if (y(i)==1)
    	{
    		n_ev++;
    		sum_ev+=d(i);
    		sse_ev+=(d(i)*d(i));
	    	if (d(i)>0) nup_ev++;
	    	if (d(i)<0) ndown_ev++;
		}
		if (y(i)==0)
    	{
	    	n_ne++;
	    	sum_ne+=d(i);
	    	sse_ne+=(d(i)*d(i));
	    	if (d(i)>0) nup_ne++;
	    	if (d(i)<0) ndown_ne++;
		}	
    }
    double pup_ev   = (double)nup_ev/(double)n_ev;
    double pup_ne   = (double)nup_ne/(double)n_ne;
    double pdown_ev = (double)ndown_ev/(double)n_ev;
    double pdown_ne = (double)ndown_ne/(double)n_ne;
    
    double v_nri_ev = (double)(nup_ev + ndown_ev)/(double)std::pow(n_ev,2) - ((double)std::pow(nup_ev - ndown_ev,2)/(double)std::pow(n_ev,3));
    
    double v_nri_ne = (double)(ndown_ne + nup_ne)/(double)((n_ne*n_ne) - ((double)std::pow(ndown_ne - nup_ne,2)/(double)std::pow(n_ne,3)));

    double nri = pup_ev - pdown_ev - (pup_ne - pdown_ne);
    double se_nri = std::sqrt(v_nri_ev + v_nri_ne);
    double z_nri  = nri/se_nri;   

    double improveSens =  sum_ev/(double)n_ev;
    double improveSpec = sum_ne/(double)n_ne;
    double idi = improveSens - (improveSpec);
    double var_ev = ((sse_ev)/(double)(n_ev))-(improveSens*improveSens);//var(a)/n_ev
    double var_ne = ((sse_ne)/(double)(n_ne))-(improveSpec*improveSpec);//var(b)/n_ne
    double se_idi = std::sqrt(var_ev/(double)(n_ev) + var_ne/(double)(n_ne));
    double z_idi = idi/se_idi;
    z_idi=(mean(da)-mean(db))/std::sqrt((var(da)/da.n_elem)+(var(db)/db.n_elem));
    vec out(4);
    	out(0)=z_idi;
    	out(1)=z_nri;
    	out(2)=idi;
    	out(3)=nri;
    return out;
  }
  
getVReclass getVarReclassificationFunc(mat dataframe,std::string type, mat independentFrame) 
{
	const int z_idi=0;
	const int z_nri=1;
	const int idi=2;
	const int nri=3;
	unsigned int i=3;
	int nzero=0;
	int n_var=dataframe.n_cols-3;
	if (type=="COX")
	{
		n_var=dataframe.n_cols-2;
		i=2;
		nzero=1;
	}
	if (independentFrame.n_rows<=1)
	{
		independentFrame = dataframe;
	}
	vec model_zidi(n_var);
	vec model_znri(n_var);
	vec model_idi(n_var);
	vec model_nri(n_var);
	vec t_model_zidi(n_var);
	vec t_model_znri(n_var);
	vec t_model_idi(n_var);
	vec t_model_nri(n_var);
	vec fullModel = modelFittingFunc(dataframe.cols(0,1),dataframe.cols(2,dataframe.n_cols-1),type);
	{
		vec fullPredict_train = predictForFresaFunc(fullModel,dataframe.cols(2,dataframe.n_cols-1),"prob",type);
		vec fullPredict = predictForFresaFunc(fullModel,independentFrame.cols(2,independentFrame.n_cols-1),"prob",type);
		for (int  j=0; i<dataframe.n_cols; i++,j++)
		{
			mat dataframe_sh=dataframe;
			mat independentFrame_sh=independentFrame;
			if (n_var>nzero)
			{
				dataframe_sh.shed_col(i);
				independentFrame_sh.shed_col(i);
			}
			vec redModel= modelFittingFunc(dataframe_sh.cols(0,1),dataframe_sh.cols(2,dataframe_sh.n_cols-1),type);
			vec iprob = improveProbFunc(predictForFresaFunc(redModel,independentFrame_sh.cols(2,independentFrame_sh.n_cols-1),"prob",type),fullPredict,independentFrame.col(1));
			vec iprob_t = improveProbFunc(predictForFresaFunc(redModel,dataframe_sh.cols(2,dataframe_sh.n_cols-1),"prob",type), fullPredict_train ,dataframe.col(1));
			model_zidi(j)=iprob(z_idi);
			model_idi(j)=iprob(idi);
			model_nri(j)=iprob(nri);
			model_znri(j)=iprob(z_nri);
			t_model_zidi(j)=iprob_t(z_idi);
			t_model_idi(j)=iprob_t(idi);
			t_model_nri(j)=iprob_t(nri);
			t_model_znri (j)=iprob_t(z_nri);
		}
	}
	getVReclass result;
	 result.z_IDIs=model_zidi;
	 result.z_NRIs=model_znri;
	 result.IDIs=model_idi;
	 result.NRIs=model_nri;
	 result.tz_IDIs=t_model_zidi;
	 result.tz_NRIs=t_model_znri;
	 result.tIDIs=t_model_idi;
	 result.tNRIs=t_model_nri;
    return result;
}

double rocAUC( vec response, vec predictor, std::string direction,std::string r) 
{
	int dir = 1;
	if (direction == "auto")
	{
		vec controls=predictor(find(response==0));
		vec cases=predictor(find(response==1));
		if (median(controls) > median(cases)) 	dir = 0;
	}
	else dir= 1*(direction==">");
    	
	double ncases = sum(response);
	double ncontrol = predictor.n_elem-ncases;
	
	uvec sindex = arma::sort_index(predictor,dir);
	double auc=0.0;
	double sen=0.0;
	double spe=0.0;
	double sena=0.0;
	double spea=0.0;
	for (unsigned int  i=0;i<predictor.n_elem;i++)
	{
		sen = sen + response[sindex[i]]/ncases;
		spe = spe + (response[sindex[i]]==0)/ncontrol;
		auc = auc + (sen+(sen-sena)/2)*(spe-spea);
		sena = sen;
		spea = spe;
	}

	return(auc);
}


vec Fresarank(vec _x)
{
	vec so=sort(_x);
	int n=so.n_elem;
    uvec si_x=sort_index(_x);
    vec sortRankx(n);
    int s;
	for(unsigned i=0;i<so.n_elem;)
	{
		s=i+1;
		int h=1;
		if (i+h<so.n_elem)
		while(so(i)==so(i+h))
		{
			s=s+(i+1+h);
			h++;
			if (i+h==so.n_elem) break;
		}
		double in=s/(double)h;
		unsigned int t=i+h;
		for (;i<t;i++)
			sortRankx(i)=in;
	}
    vec x(n);
	for(int i=0;i<n;i++)
		x(si_x(i))=sortRankx(i);
	return(x);
}

vec residualForNeRIsFunc(vec cf,mat newdata,std::string typ, std::string type,vec outcome) 
{
	vec out;
	if(type=="COX")
	{
		vec lpp =predictForFresaFunc(cf,newdata,"prob",type);
		lpp.elem(find_nonfinite(lpp)).fill(1000.0);
		out = lpp - outcome;
	}
    if(type=="LM")
	{
		out=predictForFresaFunc(cf,newdata,"linear",type)- outcome;
	}
    if(type=="LOGIT")
    {
		if (typ == "prob")
		{
			out=predictForFresaFunc(cf,newdata,"prob",type)- outcome;
		}
		else
		{
			out=predictForFresaFunc(cf,newdata,"linear",type)- outcome;
		}
	}
	if (!(out.is_finite())) 
	{
		Rcpp::Rcout<<"Warning NA predictFor NeRIs \n";
		out.elem(find_nonfinite(out)).fill(1.0e10);
	}
    return (out);
}

improvedRes improvedResidualsFunc(vec oldResiduals,vec newResiduals, std::string testType)
{
	double pwil = 1.0;
	double pbin = 1.0;
	double ptstu = 1.0;
	double f_test = 1.0;
	double p1=0.0;
	double p2=0.0;
	double pvalue = 0.0;
	double size = oldResiduals.n_elem;
	double size1 = size - 1;
	vec oldres = abs(oldResiduals);
	vec newres = abs(newResiduals);
	vec delta = newres - oldres;
	double reduction = 0.0;
	double increase  = 0.0;
	for (int i=0; i< size;i++)
	{
		reduction += (delta[i]<0.0);
		increase += (delta[i]>0.0);
	}
	double improved = (reduction-increase)/size; 									//			# the net improvement in residuals
	p1 = reduction/size;															//	#proportion of subjects with improved residuals
	p2 = increase/size;																//#proportion of subjects with worst residuals
	double rss1 = sum(square(oldResiduals));
	double rss2 = sum(square(newResiduals))/size1;
	f_test = 1-R::pf(rss1/rss2-size1,1.0,size1,1,0);
	bool paired = TRUE;
	bool var_equal = TRUE;
	std::string tail="greater";
	double mu=0.0;
	double t_pvalue = ttest(oldres, newres,mu, paired,var_equal, tail);
	if (is_finite(t_pvalue))
	{
		ptstu = t_pvalue;   	
	}
	else 
	{
		ptstu = 1.0;
	}
	if (improved>=0)
	{
		pbin = binomtest(reduction,size,0.5,tail);					// Lets do a sign test to test a significant improvement in residual variance
		pwil = wilcoxtest(oldres, newres,mu, paired,tail ,TRUE);  		// let compute that the probability that the new residuals are better than the old residuals via wilcoxon
	}
	if(testType=="Binomial") pvalue = pbin;
	if(testType=="Wilcox")   pvalue = pwil;
	if(testType=="tStudent") pvalue = ptstu;
	if(testType=="Ftest")    pvalue = f_test;
    improvedRes result;
	result.p1=p1;
	result.p2=p2;
	result.NeRI=improved;
	result.pvalue = pvalue;
	result.binom_pValue = pbin;
	result.wilcox_pValue = pwil;
	result.t_test_pValue = ptstu;
	result.F_test_pValue = f_test;
	return (result);
}

/* 
** This code is based on the t.test function from 
** file src/library/stats/R/t.test.R
** Part of the R package, http://www.R-project.org
*/

double ttest(vec x, vec y , double mu, bool paired, bool var_equal, std::string tail)
{
    if (paired) 
    {
		x = x-y;
		y.clear();
    }
    double nx =(double)x.n_elem;
    double mx = mean(x);
    double vx = var(x);
    double df;
    double tstat;
    double stderrr; 
    if(y.empty()) 
    {
        if(nx < 2.0) stop("not enough 'x' observations");
		df = nx-1.0;
		stderrr = std::sqrt(vx/nx);
	        if(stderrr < (10.0 * DOUBLEEPS* std::abs(mx)))
			{
				
//	            stop("data are essentially constant");
				tstat = 0.0;
			}
		else
		{
			tstat = (mx-mu)/stderrr;
		}
    } 
    else 
    {
		double ny = (double)y.n_elem;
		double my = mean(y);
		double vy = var(y);
		
		if(var_equal) 
		{
		    df = nx+ny-2;
	            double v = 0.0;
	            if(nx > 1) v = v + (nx-1.0)*vx;
	            if(ny > 1) v = v + (ny-1.0)*vy;
		    v = v/df;
		    stderrr = std::sqrt(v*(1/nx+1/ny));
		} 
		else 
		{
		    double stderrx = std::sqrt(vx/nx);
		    double stderry = std::sqrt(vy/ny);
		    stderrr = std::sqrt(stderrx*stderrx + stderry*stderry);
		    df = std::pow(stderrr,4)/(std::pow(stderrx,4)/(nx-1.0) + std::pow(stderry,4)/(ny-1.0));
		}
	        if(stderrr < 10.0 *DOUBLEEPS * std::max(std::abs(mx), std::abs(my)))
			{
//	            stop("data are essentially constant");
				tstat=0.0;
			}
			else
			{
				tstat = (mx - my - mu)/stderrr;
			}
    }
    double pval;
    if (tail == "less") 
    {
		pval = R::pt(tstat, df,1,0);
    }
    else 
    	if (tail == "greater") 
    	{
			pval = R::pt(tstat, df, 0,0);
    	}
   		else 
   		{
			pval = 2 * R::pt(-std::abs(tstat), df,1,0);
		}
    return(pval);
}
/* 
** This code is based on the wilcox.test function from 
** file src/library/stats/R/wilcox.test.R
** Part of the R package, http://www.R-project.org
*/

double wilcoxtest(vec x, vec y , double mu, bool paired, std::string tail,bool correct)
{ 	
    if(paired) 
    {
        x = x - y;
        y.clear(); 
    }
    double pvalue;
    double  correc = 0.0;
    if(y.empty() ) 
    {
        x = x - mu;
        bool zeros = any(x == 0.0);
        if(zeros)
            x = x.elem(find(x != 0));
        double n = (double)x.n_elem;
        bool exact=FALSE;
        exact = (n < 50);
        vec r = Fresarank(abs(x));
        double stat = sum(r.elem(find(x > 0)));
        bool ties=FALSE;
        vec ur=unique(r);
        if (r.n_elem != ur.n_elem) ties = TRUE;

        if(exact && !ties && !zeros) 
        {
           if(tail=="greater") pvalue=R::psignrank(stat - 1.0, n,0,0);
           else if(tail=="less") pvalue=R::psignrank(stat, n,1,0);
           else
            {
            	double p;
                if(stat > (n * (n + 1)/4.0)) p =R::psignrank(stat - 1, n, 0,0);
               	else p =R::psignrank(stat, n,1,0);
               	pvalue= std::min(2.0 * p, 1.0);
            }
        } 
        else 
        { 
             std::map <double, double> rtiesm; 	
                for(unsigned int i=0;i<r.n_elem;i++)
                {
                	rtiesm[r(i)]++;
                }
                vec rtiesv(rtiesm.size());
                std::map<double,double>::iterator ii=rtiesm.begin();
                 for(unsigned int i=0 ; ii!=rtiesm.end(); ++ii)
				   {
				   	   rtiesv(i)=(*ii).second;
				   	   i++;
				   }
            double z = stat - n * (n + 1)/4.0;
            double sigma = std::sqrt(n * (n + 1) * (2 * n + 1)/24.0 - sum(pow(rtiesv,3) - rtiesv)/48.0);
            if(correct) 
            {
               if(tail=="greater") correc=0.5;
	           else if(tail=="less") correc=-0.5;
	           else 
	           {
		           	if (z!=0)
		            	correc=(z/std::abs(z)) * 0.5;
		            else
		            	correc=0.5;
		        }
	           
            }
	     	z = (z - correc)/sigma;
	  
		   if(tail=="less")pvalue= R::pnorm5(z,0.0,1.0,1,0);
		   else if(tail=="greater")pvalue = R::pnorm5(z,0.0,1.0,0,0);
		   else pvalue = 2 * std::min(R::pnorm5(z,0.0,1.0,1,0),R::pnorm5(z,0.0,1.0,0,0));
           
			
		}
    }
    else 
    { 
    //-------------------------- 2-sample case ---------------------------
        vec r = Fresarank(x - y);
        double nx = (double)x.n_elem;
        double ny = (double)y.n_elem;
        bool exact=FALSE;
             if((nx < 50) && (ny < 50))
             		exact=TRUE;

        double stat =  sum(r) - nx * (nx + 1) / 2.0;
        bool ties=FALSE;
         vec ur=unique(r);
        if (r.n_elem != ur.n_elem) ties = TRUE;
     
        if(exact && !ties) 
        {
        	 if(tail=="greater")
                pvalue =R::pwilcox(stat - 1, nx, ny, 0,0);
             else if(tail=="less") pvalue =R::pwilcox(stat, nx, ny,1,0);
             else
                   {
                   	double p;
                       if(stat > (nx * ny / 2.0))
                           p =R::pwilcox(stat - 1, nx, ny, 0,0);
                       else
                           p =R::pwilcox(stat, nx, ny,1,0);
                        pvalue =std::min(2.0 * p, 1.0);
                   }
        }
        else 
        {
            std::map <double, double> rtiesm; 	
                for(unsigned int i=0;i<r.n_elem;i++)
                {
                	rtiesm[r(i)]++;
                }
                vec rtiesv(rtiesm.size());
                std::map<double,double>::iterator ii=rtiesm.begin();
                 for(unsigned int i=0 ; ii!=rtiesm.end(); ++ii)
				   {
				   	   rtiesv(i)=(*ii).second;
				   	   i++;
				   }
            double z = stat - nx * ny/2.0;
            double sigma = std::sqrt((nx * ny/12.0) * ((nx + ny + 1)- sum(pow(rtiesv,3) - rtiesv)/((nx + ny) * (nx + ny - 1.0))));
            if(correct) 
            {
               if(tail=="greater") correc=0.5;
	           else if(tail=="less") correc=-0.5;
	           else 
	           {
		           	if (z!=0)
		            	correc=(z/std::abs(z)) * 0.5;
		            else
		            	correc=0.5;
		        }
            }
		    z = (z - correc) / sigma;
		    if(tail=="less")pvalue= R::pnorm5(z,0.0,1.0,1,0);
		    else if(tail=="greater")pvalue = R::pnorm5(z,0.0,1.0,0,0);
		    else pvalue = 2 * std::min(R::pnorm5(z,0.0,1.0,1,0),R::pnorm5(z,0.0,1.0,0,0));
        }
    }
    return(pvalue);
}

double binomtest(double x, double n, double p , std::string tail)
{
    x = round(x);
    n = round(n);
    double pvalue=1.0;
    if (tail=="less") pvalue =R::pbinom(x, n, p,1,0);
    else if (tail=="greater") pvalue =R::pbinom(x - 1, n, p, 0,0);
return(pvalue);
}

gvarNeRI getVarNeRIFunc(mat dataframe, std::string type,mat testdata) 
{
	if (testdata.n_elem<=1)
	{
		testdata = dataframe;
	}
	unsigned int i=3;
	int nzero=0;
	int nvar=dataframe.n_cols-3;
	if (type=="COX")
	{
		nvar=dataframe.n_cols-2;
		nzero=1;
		i=2;
	}
	vec fullModel = modelFittingFunc(dataframe.cols(0,1),dataframe.cols(2,dataframe.n_cols-1),type);
	vec fullResiduals = residualForNeRIsFunc(fullModel,dataframe.cols(2,dataframe.n_cols-1)," ",type,dataframe.cols(1,1));
	vec testResiduals = residualForNeRIsFunc(fullModel,testdata.cols(2,dataframe.n_cols-1)," ",type,testdata.cols(1,1));
	vec model_tpvalue(nvar);
	vec model_bpvalue(nvar);
	vec model_neri(nvar);
	vec model_wpvalue(nvar);
	vec model_fpvalue(nvar);
	vec testmodel_tpvalue(nvar);
	vec testmodel_bpvalue(nvar);
	vec testmodel_neri(nvar);
	vec testmodel_wpvalue(nvar);
	vec testmodel_fpvalue(nvar);
	if (nvar>1)
	{
		for (int  j=0; i<dataframe.n_cols; i++,j++)
		{
			mat dataframe_sh=dataframe;
			mat testdata_sh=testdata;
		
		if (nvar>nzero)
		{
			dataframe_sh.shed_col(i);
			testdata_sh.shed_col(i);
		}
			vec redModel = modelFittingFunc(dataframe_sh.cols(0,1),dataframe_sh.cols(2,dataframe_sh.n_cols-1),type);

			if ( !is_finite(redModel))
			{
				redModel = fullModel;
			}
			vec redResiduals = residualForNeRIsFunc(redModel,dataframe_sh.cols(2,dataframe_sh.n_cols-1)," ",type,dataframe_sh.cols(1,1));
			vec redTestResiduals = residualForNeRIsFunc(redModel,testdata_sh.cols(2,testdata_sh.n_cols-1)," ",type,testdata_sh.cols(1,1));
			improvedRes iprob = improvedResidualsFunc(redResiduals,fullResiduals," ");
			improvedRes testiprob = improvedResidualsFunc(redTestResiduals,testResiduals," ");
			model_tpvalue(j) = iprob.t_test_pValue;
			model_bpvalue(j) = iprob.binom_pValue;
			model_wpvalue(j) = iprob.wilcox_pValue;
			model_fpvalue(j) = iprob.F_test_pValue;
			model_neri(j) = iprob.NeRI;
			testmodel_tpvalue(j) = testiprob.t_test_pValue;
			testmodel_bpvalue(j) = testiprob.binom_pValue;
			testmodel_wpvalue(j) = testiprob.wilcox_pValue;
			testmodel_fpvalue(j) = testiprob.F_test_pValue;
			testmodel_neri(j) = testiprob.NeRI;
		}
	}
	 gvarNeRI result;
	 result.tP_value=model_tpvalue;
	 result.BinP_value=model_bpvalue;
	 result.WilcoxP_value=model_wpvalue;
	 result.FP_value=model_fpvalue;
	 result.NeRIs=model_neri;
	 result.testData_tP_value=testmodel_tpvalue;
	 result.testData_BinP_value=testmodel_bpvalue;
	 result.testData_WilcoxP_value=testmodel_wpvalue;
	 result.testData_FP_value=testmodel_fpvalue;
	 result.testData_NeRIs=testmodel_neri;
     return (result);
}


