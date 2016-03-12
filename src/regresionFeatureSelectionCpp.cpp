/* FRESA.CAD: utilities for building and testing formula-based 
	models (linear, logistic or COX) for Computer Aided Diagnosis/Prognosis 
	applications.  Utilities include data adjustment, univariate analysis, 
	model building, model-validation, longitudinal analysis, reporting and visualization.. 

   This program is free software under the terms of the 
   GPL Lesser General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
   
   Jose Tamez and Israel Alanis
  
*/

#include "FRESAcommons.h"


struct NeRIbootVal{
	mat bootmodel;
	mat bootmodelmeans;
	mat NeRi;
	mat tpvalue;
	mat wpvalue;
	mat spvalue;
	mat Fpvalue;
	mat test_tpvalue;
	mat test_wpvalue;
	mat test_spvalue;
	mat test_Fpvalue;
 }NeRIredBoot;

int redCurmodel_S_lastRemovedNeRI;


//********************************************************************************************
//**==========================================================================================
//********************************************************************************************
extern "C" SEXP bootstrapValidationResCpp(SEXP _fraction,SEXP _loops,SEXP _dataframe,SEXP _type,SEXP _response)
{
	try 
	{  //R_CStackLimit=(uintptr_t)-1;
		double fraction = Rcpp::as<double>(_fraction);
		int loops = Rcpp::as<int>(_loops);
		std::string type = Rcpp::as<std::string>(_type);
		Rcpp::NumericMatrix dataf(_dataframe);
		Rcpp::NumericMatrix resp(_response);
	    mat dataframe;
	    mat Response(resp.rows(), 2);
		mat dat(dataf.begin(), dataf.rows(), dataf.cols(), false);
		vec dtime=ones<vec>(Response.n_rows);
		if (resp.cols()==1) 
			{
				Response.col(0)=dtime;
				Response.col(1)=Rcpp::as<vec>(Rcpp::NumericVector(_response));
			}
			else 
			{
			 mat Res(resp.begin(), resp.rows(), resp.cols(), false);
			 Response=Res;
			}
		vec Outcome=Response.col(1);
		dtime=Response.col(0);
        int n_var2=dat.n_cols-1;
        if (type=="COX")
		{
			dataframe=dat.cols(1,dat.n_cols-1);
		}else dataframe=dat;
	    mat casesample=join_rows(Response,dataframe);
		int sizecases   =dataframe.n_rows;
		int totSamples = (int)(fraction*sizecases);
		std::vector<double>trainOutcome;
		std::vector<double>trainPrediction;
		std::vector<double>trainResiduals;
		vec trainSampledRMSE(loops);
		std::vector<double>testOutcome;
		std::vector<double>testPrediction;
		std::vector<double>testResiduals;
		vec testSampledRMSE(loops);
		mat bcoef(loops,dataframe.n_cols); 	  
		mat bmeans(loops,dataframe.n_cols); 
		mat	NeRi(loops,n_var2); 
		mat	tpvalue=ones<mat>(loops,n_var2);  
		mat	wpvalue=ones<mat>(loops,n_var2);  
		mat	spvalue=ones<mat>(loops,n_var2);  
		mat	Fpvalue=ones<mat>(loops,n_var2);  
		mat	test_tpvalue=ones<mat>(loops,n_var2); 
		mat	test_wpvalue=ones<mat>(loops,n_var2); 
		mat	test_spvalue=ones<mat>(loops,n_var2); 
		mat	test_Fpvalue=ones<mat>(loops,n_var2);
		int lastInserted=0;		
#pragma omp parallel for schedule(dynamic) ordered shared(lastInserted,trainOutcome,trainPrediction,trainResiduals,testOutcome,testPrediction,testResiduals,trainSampledRMSE,testSampledRMSE,bcoef,bmeans,NeRi,tpvalue,wpvalue,spvalue,Fpvalue,test_tpvalue,test_wpvalue,test_spvalue,test_Fpvalue)
		for (int doOver=0;doOver<loops;doOver++)
		{ 
			mat trainingSample;
			uvec ntestSample;
			unsigned int mintestSample = 0.2*totSamples;
			do
			{
				uvec samCases= randi<uvec>(totSamples, distr_param(0,sizecases-1));
				vec auxcasesample=zeros<vec>(sizecases);
				trainingSample = casesample.row(samCases(0));
				auxcasesample[samCases(0)]=1;
				for (int i =1;i<totSamples;i++)
				{
					trainingSample=join_cols(trainingSample,casesample.row(samCases(i)));
					auxcasesample[samCases(i)]=1;
				}
				ntestSample=find(auxcasesample==0);
			} while (ntestSample.n_elem<mintestSample);			
			vec trainmodel = modelFittingFunc(trainingSample.cols(0,1),trainingSample.cols(2,trainingSample.n_cols-1),type);
			if (is_finite(trainmodel))
			{
				mat testSample = casesample.rows(ntestSample);
				vec coef = trainmodel;
				vec coxmean;
				if (type == "COX") 
				{
					 coef = trainmodel.subvec(0,(trainmodel.n_elem/2)-1);
					 coxmean = trainmodel.subvec((trainmodel.n_elem/2),trainmodel.n_elem-1);
				}

				vec residualTrain = residualForFRESAFunc(trainmodel,trainingSample.cols(2,trainingSample.n_cols-1),"",type,trainingSample.cols(0,1));
				vec residualTest = residualForFRESAFunc(trainmodel,testSample.cols(2,testSample.n_cols-1),"",type,testSample.cols(0,1));
				gvarNeRI modelReclas = getVarResFunc(trainingSample,type,testSample);

#pragma omp critical
{
				for (unsigned int i=0;i<residualTrain.n_elem;i++)
				{
					trainResiduals.push_back(residualTrain(i));
				}

				for (unsigned int i=0;i<residualTest.n_elem;i++)
				{
					testOutcome.push_back(testSample(i,1));
					testPrediction.push_back(testSample(i,1)+residualTest(i));
					testResiduals.push_back(residualTest(i));
				}
				trainSampledRMSE(lastInserted) = std::sqrt(mean(residualTrain%residualTrain));
				testSampledRMSE(lastInserted) = std::sqrt(mean(residualTest%residualTest));
				bcoef.row(lastInserted) = coef.t();	
				if (type == "COX") bmeans.row(lastInserted) = coxmean.t();
				NeRi.row(lastInserted) = modelReclas.NeRIs.t();
				tpvalue.row(lastInserted) = modelReclas.tP_value.t();
				wpvalue.row(lastInserted) = modelReclas.WilcoxP_value.t();
				spvalue.row(lastInserted) = modelReclas.BinP_value.t();
				Fpvalue.row(lastInserted) = modelReclas.FP_value.t();
				test_tpvalue.row(lastInserted) = modelReclas.testData_tP_value.t();
				test_wpvalue.row(lastInserted) = modelReclas.testData_WilcoxP_value.t();
				test_spvalue.row(lastInserted) = modelReclas.testData_BinP_value.t();
				test_Fpvalue.row(lastInserted) = modelReclas.testData_FP_value.t();
				lastInserted = lastInserted + 1;
}
			}
		}
		if ((lastInserted>0)&&(lastInserted<loops))
		{
			trainSampledRMSE.resize(lastInserted);
			testSampledRMSE.resize(lastInserted);
			bcoef.resize(lastInserted,dataframe.n_cols);	
			if (type == "COX") bmeans.resize(lastInserted,dataframe.n_cols);
			NeRi.resize(lastInserted,n_var2);
			if ((3*lastInserted)<loops) Rcout<<"Warning: only "<<lastInserted<<" samples of "<<loops<<" were inserted\n";
		}

	 	rowvec bootmodel = mean (bcoef);
		if (type == "COX") 
		{ 
			rowvec CoxMmodel = mean (bmeans);
			rowvec aux(bootmodel.n_elem*2);
			for (unsigned int i=0;i<bootmodel.n_elem;i++)
			{
				aux[i] = bootmodel[i];
				aux[bootmodel.n_elem+i] = CoxMmodel[i];
			}
			bootmodel.resize(aux.n_elem);
			bootmodel = aux;
		}
	 	vec basemodel = modelFittingFunc(casesample.cols(0,1),casesample.cols(2,casesample.n_cols-1),type);
		double startRMSE = std::sqrt(mean(square(residualForFRESAFunc(basemodel,casesample.cols(2,casesample.n_cols-1),"",type,casesample.cols(0,1)))));
		double bootRMSE = std::sqrt(mean(square(residualForFRESAFunc(bootmodel.t(),casesample.cols(2,casesample.n_cols-1),"",type,casesample.cols(0,1)))));
		Rcpp::List result = Rcpp::List::create(Rcpp::Named("NeRi")=Rcpp::wrap(NeRi), \
						Rcpp::Named("tpvalue")=Rcpp::wrap(tpvalue), \
						Rcpp::Named("wpvalue")=Rcpp::wrap(wpvalue), \
						Rcpp::Named("bcoef")=Rcpp::wrap(bcoef), \
						Rcpp::Named("spvalue")=Rcpp::wrap(spvalue), \
						Rcpp::Named("Fpvalue")=Rcpp::wrap(Fpvalue), \
						Rcpp::Named("test_tpvalue")=Rcpp::wrap(test_tpvalue), \
						Rcpp::Named("test_wpvalue")=Rcpp::wrap(test_wpvalue), \
						Rcpp::Named("test_spvalue")=Rcpp::wrap(test_spvalue), \
						Rcpp::Named("test_Fpvalue")=Rcpp::wrap(test_Fpvalue), \
						Rcpp::Named("trainSampledRMSE")=Rcpp::wrap(trainSampledRMSE), \
						Rcpp::Named("testOutcome")=Rcpp::wrap(testOutcome), \
						Rcpp::Named("testPrediction")=Rcpp::wrap(testPrediction), \
						Rcpp::Named("testResiduals")=Rcpp::wrap(testResiduals), \
						Rcpp::Named("testSampledRMSE")=Rcpp::wrap(testSampledRMSE),\
						Rcpp::Named("startRMSE")=Rcpp::wrap(startRMSE), \
						Rcpp::Named("bootRMSE")=Rcpp::wrap(bootRMSE), \
						Rcpp::Named("trainResiduals")=Rcpp::wrap(trainResiduals) 
					);
						

		return (result);
	} 
	catch( std::exception &ex ) 
	{ 
		forward_exception_to_r( ex );
		return 0;
	}	
	catch(...) 
	{ 
		::Rf_error( "c++ exception (unknown reason)" );
		return 0;
	}
}

std::string residualFowardSelection(const unsigned int size,const int totSamples,const double pthr2,const std::string testType,
									const int loops,const std::string Outcome,const std::string timeOutcome,const std::string type,const unsigned int maxTrainModelSize,
									const mat &dataframe, const std::vector < std::string > &ovnames,
									std::map<std::string, int> &lockup,const std::string baseForm,std::vector <int> &mynamesLoc, 
									const std::vector<std::string> &covari)
{
	std::string frm1 = baseForm;
	unsigned int inserted = 0;
	double iprob=0;
	unsigned int sizecases = dataframe.n_rows;

	

	arma::mat mysample = dataframe;
	arma::mat myTestsample = dataframe;

	arma::mat mysampleMatbase;
	arma::mat myTestsampleMatbase;
	vec bestmodel;
	vec newmodel;
	vec bestTrainResiduals;
	vec bestResiduals;
	std::vector < std::string > vnames = ovnames;
	std::string inname = "Inserted";

	int outinx = 0;
	vec randOutput=zeros<vec>(sizecases);
	if (Outcome=="RANDOM")
	{
		if (type=="LM")
		{
			randOutput = randu<vec>(sizecases);
		}
		else
		{
			randOutput = randi<vec>(sizecases, distr_param(0,1));
		}
	}
	else
	{
		outinx=lockup[Outcome];
	}
	
	int timeinx	= 0;
	if (type=="COX")
	{
		timeinx	= lockup[timeOutcome];
	}

	for(unsigned  int j=0;j<covari.size();j++)
		for(unsigned int i=0;i<vnames.size();i++)
			if (covari[j]==vnames[i]) vnames[i]=inname;
	if (loops > 1)
	{
		int testelements=0;
		int mitestsize=0.2*sizecases;
		do
		{
			uvec samCases= randi<uvec>(totSamples, distr_param(0,sizecases-1));
			vec auxcasesample=zeros<vec>(sizecases);
			mysample = dataframe.row(samCases(0));
			auxcasesample[samCases(0)]=1;
			if (Outcome=="RANDOM") mysample.col(outinx)(0) = randOutput(samCases(0));
			for (int i = 1;i<totSamples;i++)
			{
				mysample=join_cols(mysample,dataframe.row(samCases(i)));
				if (Outcome=="RANDOM") mysample.col(outinx)(i) = randOutput(samCases(i));
				auxcasesample[samCases(i)]=1;
			}
			uvec notinserted = find(auxcasesample==0);
			testelements = notinserted.n_elem;
			myTestsample = dataframe.rows(notinserted);
			if (Outcome=="RANDOM") myTestsample.col(outinx) = randOutput(notinserted);
		}	
		while (testelements<mitestsize);

	}

	mat train_outcome = join_rows(mysample.col(timeinx),mysample.col(outinx));
	mat test_outcome = join_rows(myTestsample.col(timeinx),myTestsample.col(outinx));

    mat myoutcome(mysample.n_rows,2);
	myoutcome.col(1) = train_outcome.col(1);
	mysampleMatbase.set_size(mysample.n_rows);
	myTestsampleMatbase.set_size(myTestsample.n_rows);
	mysampleMatbase.fill(1.0);
	myTestsampleMatbase.fill(1.0);
	if (type=="COX")
	{
		myoutcome.col(0)=mysample.col(lockup[timeOutcome]);
	}
	if (covari[0]!="1")
		for (unsigned int i=0;i<covari.size();i++)
	  	{
			mysampleMatbase=join_rows(mysampleMatbase,mysample.col(lockup[covari[i]]));
	  		myTestsampleMatbase=join_rows(myTestsampleMatbase,myTestsample.col(lockup[covari[i]]));		  
  	    }
	if((mysampleMatbase.n_cols>1)or(type=="COX"))
		bestmodel =modelFittingFunc(myoutcome,mysampleMatbase,type);
	else
    {
		bestmodel=mean(myoutcome.col(1));
       	if (type=="LOGIT")
		{
			if ((bestmodel(0)>0)&&(bestmodel(0)<1))
			{
				bestmodel=log(bestmodel/(1.0-bestmodel));
			}
			else
			{
				if (bestmodel(0)==1) bestmodel(0) = THRESH;
				else bestmodel(0) = -THRESH;
			}
		}
    } 
  	bestTrainResiduals=residualForFRESAFunc(bestmodel,mysampleMatbase,"",type,train_outcome);
  	if (loops > 1) bestResiduals=residualForFRESAFunc(bestmodel,myTestsampleMatbase,"",type,test_outcome);
	int changes = 1;
	double minpiri = 1;
	int jmax = -1;
  	std::string gfrm1;
	mat myTestsampleMat;
	mat mysampleMat;
	vec testResiduals;
	vec trainResiduals;
	while (changes>0)
	{
		changes = 0;
	  	minpiri = 1;
	  	jmax = -1;
	  	for (unsigned int j=0; j<size;j++)
	  	{
			if (vnames[j] != inname)
	  		{
				mysampleMat=join_rows(mysampleMatbase,mysample.col(lockup[vnames[j]]));
				if (loops > 1) myTestsampleMat=join_rows(myTestsampleMatbase,myTestsample.col(lockup[vnames[j]]));
	  			gfrm1 = frm1+" + "+vnames[j];
             	if((mysampleMat.n_cols>1)or(type=="COX"))
					newmodel =modelFittingFunc(myoutcome,mysampleMat,type);
			    else
			    {
					newmodel=mean(myoutcome.col(1));
        			if (type=="LOGIT")
					{
						if ((newmodel(0)>0)&&(newmodel(0)<1))
						{
							newmodel=log(newmodel/(1.0-newmodel));
						}
						else
						{
							if (newmodel(0)==1) newmodel(0) = THRESH;
							else newmodel(0) = -THRESH;
						}
					}
        		}
	  			if (is_finite(newmodel))
	  			{
					trainResiduals = residualForFRESAFunc(newmodel,mysampleMat,"",type,train_outcome);
			   		double piri = improvedResidualsFunc(bestTrainResiduals,trainResiduals,testType).pvalue;            
					if (loops > 1) 
					{
						testResiduals = residualForFRESAFunc(newmodel,myTestsampleMat,"",type,test_outcome);
						iprob = improvedResidualsFunc(bestResiduals,testResiduals,testType,totSamples).pvalue;
						if ((is_finite(iprob)) and (is_finite(piri)))
						{
							piri =std::max(iprob,piri);
						}
					}						
					if (is_finite(piri)&&(piri  < minpiri))
					{
						jmax = j;
						minpiri = piri;
					}
	  			}
	  		}
	  	}
		if ((jmax >= 0) and (minpiri < pthr2) and (inserted<maxTrainModelSize))   
	  	{
			if (vnames[jmax] != inname)
			{
				gfrm1 = frm1+" + ";
				gfrm1 = gfrm1+ovnames[jmax];
				mysampleMatbase = join_rows(mysampleMatbase,mysample.col(lockup[vnames[jmax]]));
				bestmodel = modelFittingFunc(myoutcome,mysampleMatbase,type);
				if (loops > 1) 
				{
					myTestsampleMatbase = join_rows(myTestsampleMatbase,myTestsample.col(lockup[vnames[jmax]]));
					bestResiduals = residualForFRESAFunc(bestmodel,myTestsampleMatbase,"",type,test_outcome);
				}
				bestTrainResiduals = residualForFRESAFunc(bestmodel,mysampleMatbase,"",type,train_outcome);
				frm1 = gfrm1;	
				changes = changes + 1;
				mynamesLoc.push_back(jmax);
				vnames[jmax] = inname;
				inserted = inserted + 1;
			}
	  	} 
	}
	return frm1;
}

extern "C" SEXP ForwardResidualModelCpp(SEXP _size, SEXP _fraction,SEXP _pvalue, SEXP _loops,  SEXP _covariates, SEXP _Outcome, SEXP _variableList, SEXP _maxTrainModelSize, SEXP _type, SEXP _timeOutcome, SEXP _testType,SEXP _loop_threshold, SEXP _interaction, SEXP _dataframe,SEXP _colnames,SEXP _cores)
{
try {   //R_CStackLimit=(uintptr_t)-1;
		std::string type = as<std::string>(_type);
		std::string testType = as<std::string>(_testType);
#ifdef _OPENMP
//		Rcout<<"::============ Available CPU Cores(Threads) : "<<omp_get_num_procs()<<" ===Requested Threads : " << as<unsigned int>(_cores) <<endl;
		omp_set_num_threads(as<unsigned int>(_cores));
#endif	
		unsigned int size = as<unsigned int>(_size); 
		double fraction = as<double>(_fraction);
		double pvalue = as<double>(_pvalue);
		int loops = as<int>(_loops);
		std::string Outcome = as<std::string>(_Outcome); 
		CharacterVector colnames(_colnames);
	    int maxTrainModelSize = Rcpp::as<int>(_maxTrainModelSize);
		std::string timeOutcome = as<std::string>(_timeOutcome);
		Rcpp::NumericMatrix DF(_dataframe);
    	arma::mat dataframe(DF.begin(), DF.rows(), DF.cols(), false);
		std::vector<std::string>covari=as<std::vector<std::string> >(_covariates);
		std::string baseForm=Outcome;
		std::map<std::string, int> lockup;
		std::vector<std::string> formulas;
		std::vector<int> mynames;
		std::vector<std::string>ovnames=as<std::vector<std::string> >(CharacterVector(_variableList)); 
		double pthr = pvalue;
		for(int i=0;i<CharacterVector(colnames).size();i++)
		{
			lockup[std::string(CharacterVector(colnames)[i])]=i;
		}
		if (type=="COX")
		{
		  baseForm = "Surv("+timeOutcome+","+Outcome+")";
		}
		baseForm = baseForm + " ~";
		for(unsigned int i=0;i<covari.size();i++)
		{
			baseForm = baseForm+" + "+covari[i];
		}
		int sizecases   =dataframe.n_rows;
		int totSamples = (int)(fraction*sizecases);
		double pthr2 = 1.0-R::pnorm(std::sqrt(fraction)*std::abs(qnorm(pthr,0.0,1.0)),0.0,1.0,1,0);
	    if (pthr2>0.1) pthr2 = 0.1;
//		Rcout<<pvalue<<"pv: "<<pthr2<<" : "<<baseForm<< "\n";
		if (size > ovnames.size()) size = ovnames.size();

		
#pragma omp parallel for schedule(dynamic) ordered shared(mynames,formulas)
    	for (int doOver=0; doOver<loops;doOver++)
	 	{   
			std::vector <int> mynamesLoc;
			std::string  frm;
			mynamesLoc.clear();
			frm = residualFowardSelection(size,totSamples,pthr2,testType,loops,
			Outcome,timeOutcome,type,maxTrainModelSize,dataframe,
			ovnames,lockup,baseForm,mynamesLoc,covari);
			if ((doOver%10)==1)
			{				
#pragma omp critical
				Rcout << ".";	
//				Rcout<<doOver<<" : "<<frm<<std::endl;	
			}
#pragma omp critical
{
				mynames.insert(mynames.end(), mynamesLoc.begin(), mynamesLoc.end());
				formulas.push_back(frm);
}
  		}
  		if (mynames.size() == 0) mynames.push_back(0);
	    List result = List::create(Named("mynames")=wrap(mynames),Named("formula.list")=wrap(formulas));
	    return (result);
	} 
	catch( std::exception &ex ) 
	{ 
		forward_exception_to_r( ex );
		return 0;
	}	
	catch(...) 
	{ 
		::Rf_error( "c++ exception (unknown reason)" );
		return 0;
	}

}
