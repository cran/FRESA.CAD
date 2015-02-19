#ifndef Fresa_commons_h
#define Fresa_commons_h


#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R_ext/Utils.h>
#include <exception>
#include <iostream>
#include <iomanip>
#include <map>
#ifdef _WIN32
#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#elif MACOS
#include <errno.h>
#include <unistd.h>
#include <dlfcn.h>
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <dlfcn.h>
#include <unistd.h>
#endif
#define HAVE_UINTPTR_T
#define CSTACK_DEFNS 7
#include "Rinterface.h"
using namespace arma;
using namespace Rcpp;
static const double THRESH = 36.04365;
static const double MTHRESH = -36.04365;
static const double INVEPS = 1/2.220446e-16;
static const double DOUBLEEPS = 2.220446e-16;

struct getVReclass{
  	 vec z_IDIs;
	 vec z_NRIs;
	 vec IDIs;
	 vec NRIs;
	 vec tz_IDIs;
	 vec tz_NRIs;
	 vec tIDIs;
	 vec tNRIs;

  };
 struct gvarNeRI{
	 vec tP_value;
	 vec BinP_value;
	 vec WilcoxP_value;
	 vec FP_value;
	 vec NeRIs;
	 vec testData_tP_value;
	 vec testData_BinP_value;
	 vec testData_WilcoxP_value;
	 vec testData_FP_value;
	 vec testData_NeRIs;
	};
 struct improvedRes {
	double p1;
	double p2;
	double NeRI;
	double pvalue;
	double binom_pValue;
	double wilcox_pValue;
	double t_test_pValue;
	double F_test_pValue;
  };
#include <stdio.h>
#include <math.h>
double qnorm(double p, double mu, double sigma);
void chinv2(mat &matrix , int n);
int cholesky2(mat &matrix, int n, double toler);
void chsolve2(mat &matrix, int n, vec &y);
vec coxfit(int  maxiter,  vec &time,  vec &status, mat &covar, vec &offset, vec &weights, vec &strata,   int method, double eps, double toler,    vec &beta,    int doscale);
vec logit_link(vec mu);
vec logit_linkinv(vec eta);
vec logit_mu_eta(vec eta);
vec binomial_dev_resids(vec y, vec mu, vec wt);
vec modelFittingFunc(mat ymat, mat X,std::string type);
vec predictForFresaFunc(vec cf,mat newdata,std::string typ, std::string opc);		
vec improveProbFunc(vec x1, vec x2, vec y); 
getVReclass getVarReclassificationFunc(mat dataframe,std::string type, mat independentFrame);
double rocAUC( vec controls, vec cases, std::string direction,std::string r);
vec Fresarank(vec _x);
vec residualForNeRIsFunc(vec cf,mat newdata,std::string typ, std::string type,vec outcome);
improvedRes improvedResidualsFunc(vec oldResiduals,vec newResiduals, std::string testType);
double ttest(vec x, vec y , double mu, bool paired, bool var_equal, std::string tail);
double wilcoxtest(vec x, vec y , double mu, bool paired, std::string tail,bool correct);
double binomtest(double x, double n, double p , std::string tail);
gvarNeRI getVarNeRIFunc(mat dataframe, std::string type,mat testdata);
#endif
// END

