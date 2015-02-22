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


struct bootVal{
  	 mat bootmodel;
  	 mat bootmodelmeans;
	 double BlindAccuracy;
	 double BlindSensitivity;
	 double BlindSpecificity;
	 std::vector<double> testoutcome;
	 std::vector<double> testprediction;
	 mat zNRI;
	 mat zIDI;
	 mat test_zIDI;
	 mat test_zNRI;
  }redBoot;

int redCurmodel_S_lastRemoved;


//********************************************************************************************
//**==========================================================================================
//********************************************************************************************
extern "C" SEXP bootstrapValidationCpp(SEXP _fraction,SEXP _loops,SEXP _dataframe,SEXP _type,SEXP _response)
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
		int n_var=dataframe.n_cols;
	    mat casesample=join_rows(Response.rows(find(Outcome==1)),dataframe.rows(find(Outcome==1)));
	    mat controlsample=join_rows(Response.rows(find(Outcome==0)),dataframe.rows(find(Outcome==0)));

		int sizecases   =casesample.n_rows;
		int sizecontrol =controlsample.n_rows;
		int minsize = std::min(sizecases,sizecontrol);
		int totSamples = (int)(fraction*minsize);

		std::vector<double> testoutcome;
		std::vector<double> testprediction;

		vec trainRoc=ones<vec>(loops);
		vec accuracy(loops);
		vec sensitivity(loops);
		vec specificity(loops);
		vec taccuracy(loops);
		vec tsensitivity(loops);
		vec tspecificity(loops);
		mat bcoef(loops,n_var); 	  
		mat zNRI(loops,n_var2); 
		mat zIDI(loops,n_var2); 
		mat test_zNRI(loops,n_var2); 
		mat test_zIDI(loops,n_var2); 
		mat NRI(loops,n_var2); 
		mat IDI(loops,n_var2); 
		vec sumwtdcf=zeros<vec>(n_var);
		vec sumwts = sumwtdcf;
		double gacc = 0.0; 
		double gsen = 0.0; 
		double gspe = 0.0; 
		int smaptot = 0;
		int smapsen = 0;
		int smapspe = 0;
		int lastInserted=0;

		vec rawbetas = modelFittingFunc(Response,dataframe,type);

		
		#pragma omp parallel for schedule(dynamic) ordered shared(lastInserted,sumwtdcf,sumwts,testoutcome,testprediction,gacc,gsen,gspe,smaptot,smapsen,smapspe,trainRoc,bcoef,accuracy,sensitivity,specificity,taccuracy,tsensitivity,tspecificity,zNRI,zIDI,test_zNRI,test_zIDI,NRI,IDI)  
		for (int doOver=0;doOver<loops;doOver++)
		{ 
			mat trainingSample;
			mat testSample;
			mat myTestCases;
			mat myTestControl;
			uvec samCases= randi<uvec>(totSamples, distr_param(0,sizecases-1));
			uvec samControl=randi<uvec>(totSamples, distr_param(0,sizecontrol-1));
			vec auxcasesample=zeros<vec>(sizecases);
			vec auxcontrolsample=zeros<vec>(sizecontrol);
			for (int i =0;i<totSamples;i++)
			{
				trainingSample=join_cols(trainingSample,casesample.row((samCases(i))));
				trainingSample=join_cols(trainingSample,controlsample.row((samControl(i))));
				auxcasesample[samCases(i)]=1;
				auxcontrolsample[samControl(i)]=1;
			}	 

		    vec trainmodel = modelFittingFunc(trainingSample.cols(0,1),trainingSample.cols(2,trainingSample.n_cols-1),type);
			if (is_finite(trainmodel))
			{
				// the model beta coefficients
				vec coef=trainmodel;
				if (type == "COX") 
				{
					  vec coefaux = trainmodel.subvec(0,(trainmodel.n_elem/2)-1);
   					  coef=coefaux;
				}

				// Getting the test set
				myTestCases=casesample.rows(find(auxcasesample==0));
				myTestControl=controlsample.rows(find(auxcontrolsample==0));
				testSample=join_cols(myTestCases,myTestControl);
				mat indpendentSample= join_cols(myTestCases.rows(randi<uvec>(totSamples, distr_param(0,myTestCases.n_rows-1)))  ,  myTestControl.rows(randi<uvec>(totSamples, distr_param(0,myTestControl.n_rows-1))));

				// predicting the test
				vec p = predictForFresaFunc(trainmodel,testSample.cols(2,testSample.n_cols-1),"prob",type);
				// predicting the training
				vec trainmodel_linear_predictors = predictForFresaFunc(trainmodel,trainingSample.cols(2,trainingSample.n_cols-1),"linear",type);

				// performance AUC and NRI/IDI
				double auul = rocAUC(trainingSample.col(1), trainmodel_linear_predictors,"auto","auc");
				getVReclass modelReclas = getVarReclassificationFunc(trainingSample,type,indpendentSample);
				
				double acc = 0.0;
				double sen = 0.0;
				double spe = 0.0;
				for(unsigned int i=0;i<testSample.n_rows;i++)
				{
					if (testSample(i,1)>0)
					{
				    	if (p(i)>=0.5) sen++;
					}
					if (testSample(i,1)==0)
					{
				    	if (p(i)<0.5) spe++;
					}	
				}
				acc = sen+spe;

				double trsen=0;
				double trspe=0;
				double tracc=0;
				for(unsigned int i=0;i<trainmodel_linear_predictors.n_rows;i++)
				{
					if (trainingSample(i,1)>0)
					{
					 	if (trainmodel_linear_predictors(i)>=0) trsen++;
					}
					if (trainingSample(i,1)==0)
					{
				    	if (trainmodel_linear_predictors(i)<0) trspe++;
					}	
				}
				tracc = trsen+trspe;
				tracc = tracc/static_cast<double>(trainmodel_linear_predictors.n_rows);
				trsen = trsen/static_cast<double>(totSamples);
				trspe = trspe/static_cast<double>(totSamples);

				vec wtsIdi = modelReclas.z_IDIs;
				for (unsigned int n=0;n<wtsIdi.n_elem;n++)
				{
					wtsIdi(n) = 0;
					if (is_finite(modelReclas.tz_IDIs(n)) && is_finite(modelReclas.IDIs(n)))
				 	{
						if (modelReclas.tz_IDIs(n) > 0) 
						{
							wtsIdi(n) = (modelReclas.IDIs(n)-std::abs(modelReclas.IDIs(n)-modelReclas.tIDIs(n)))/modelReclas.tIDIs(n);
							if (wtsIdi(n) < -0.1) 
					 		{
								wtsIdi(n) = -0.1;	// just down-weight bad prediction
							}
						}
					}
				}
				if (type != "COX") 
				{
					wtsIdi.insert_rows(0,ones<vec>(1));
				}
				
				
#pragma omp critical
{
				for (unsigned int i=0;i<p.n_elem;i++)
				{
					testoutcome.push_back(testSample(i,1));
					testprediction.push_back(p(i));
				}
				gacc = gacc + acc;
				smaptot = smaptot+testSample.n_rows;
				gsen = gsen + sen; 
				smapsen = smapsen + myTestCases.n_rows;
				gspe = gspe + spe; 
				smapspe = smapspe + myTestControl.n_rows;
				acc = acc/static_cast<double>(testSample.n_rows);
				sen = sen/static_cast<double>(myTestCases.n_rows);
				spe = spe/static_cast<double>(myTestControl.n_rows);
				sumwtdcf = sumwtdcf+(wtsIdi % coef);
				sumwts = sumwts+wtsIdi;

				trainRoc(lastInserted) = auul;
				bcoef.row(lastInserted) = coef.t();	
				accuracy(lastInserted) = acc;
				sensitivity(lastInserted) = sen;
				specificity(lastInserted) = spe;
				taccuracy(lastInserted)=tracc;
				tsensitivity(lastInserted)=trsen; 
				tspecificity(lastInserted)=trspe;
				zNRI.row(lastInserted) = modelReclas.tz_NRIs.t();
				zIDI.row(lastInserted) = modelReclas.tz_IDIs.t();
				test_zNRI.row(lastInserted) = modelReclas.z_NRIs.t();
				test_zIDI.row(lastInserted) = modelReclas.z_IDIs.t();
				NRI.row(lastInserted) = modelReclas.NRIs.t();
				IDI.row(lastInserted) = modelReclas.IDIs.t();
				lastInserted=lastInserted+1;
}
			}
		}
		if (lastInserted<loops)
		{
			trainRoc.resize(lastInserted);
			accuracy.resize(lastInserted);
			sensitivity.resize(lastInserted);
			specificity.resize(lastInserted);
			taccuracy.resize(lastInserted);
			tsensitivity.resize(lastInserted);
			tspecificity.resize(lastInserted);
			bcoef.resize(lastInserted,n_var);
			zNRI.resize(lastInserted,n_var2);
			zIDI.resize(lastInserted,n_var2);
			test_zNRI.resize(lastInserted,n_var2);
			test_zIDI.resize(lastInserted,n_var2);
			NRI.resize(lastInserted,n_var2);
			IDI.resize(lastInserted,n_var2);
			Rcout<<"Warning: only "<<lastInserted<<" samples of "<<loops<<" where inserted\n";
		}
		
		double blindAUC = rocAUC(testoutcome,testprediction,"auto","auc");
		double BlindAccuracy = gacc/static_cast<double>(smaptot);
		double BlindSensitivity = gsen/static_cast<double>(smapsen);
		double BlindSpecificity = gspe/static_cast<double>(smapspe);
//		Rcout<<"Blind AUC: "<<blindAUC<<"\n";
		Rcpp::List resul =Rcpp::List::create(
			Rcpp::Named("test_zNRI")=Rcpp::wrap(test_zNRI), 
			Rcpp::Named("test_zIDI")=Rcpp::wrap(test_zIDI),
			Rcpp::Named("rawbetas")=Rcpp::wrap(rawbetas), 
			Rcpp::Named("blindAUC")=Rcpp::wrap(blindAUC) 
			);
		Rcpp::List result = Rcpp::List::create(Rcpp::Named("BlindAccuracy")=Rcpp::wrap(BlindAccuracy), 
						Rcpp::Named("BlindSensitivity")=Rcpp::wrap(BlindSensitivity), 
						Rcpp::Named("BlindSpecificity")=Rcpp::wrap(BlindSpecificity), 
						Rcpp::Named("bcoef")=Rcpp::wrap(bcoef), 
						Rcpp::Named("zNRI")=Rcpp::wrap(zNRI), 
						Rcpp::Named("zIDI")=Rcpp::wrap(zIDI), 
						Rcpp::Named("trainRoc")=Rcpp::wrap(trainRoc), 
						Rcpp::Named("NRI")=Rcpp::wrap(NRI), 
						Rcpp::Named("IDI")=Rcpp::wrap(IDI), 
						Rcpp::Named("sumwtdcf")=Rcpp::wrap(sumwtdcf), 
						Rcpp::Named("sumwts")=Rcpp::wrap(sumwts), 
						Rcpp::Named("taccuracy")=Rcpp::wrap(taccuracy), 
						Rcpp::Named("tsensitivity")=Rcpp::wrap(tsensitivity), 
						Rcpp::Named("tspecificity")=Rcpp::wrap(tspecificity), 
						Rcpp::Named("accuracy")=Rcpp::wrap(accuracy), 
						Rcpp::Named("sensitivity")=Rcpp::wrap(sensitivity), 
						Rcpp::Named("specificity")=Rcpp::wrap(specificity), 
						Rcpp::Named("testprediction")=Rcpp::wrap(testprediction), 
						Rcpp::Named("testoutcome")=Rcpp::wrap(testoutcome),
						Rcpp::Named("resul")=resul
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

std::string binaryFowardSelection(const unsigned int size,const int totSamples,const double zthr2,const std::string selType,
								const int loops,std::string Outcome,const std::string timeOutcome,const std::string type,const unsigned int maxTrainModelSize,
								const mat &casesample,const mat &controlsample,const  std::vector < std::string > &ovnames,
								std::map<std::string, int> &lockup,const std::string baseForm,std::vector <int> &mynamesLoc, 
								const std::vector<std::string> &covari)
{
	unsigned int inserted = 0;
	double zmin = 0;
	double maxrec = 0;
	arma::mat mysample;
	arma::mat myTestsample;
	arma::mat mysampleMat;
	arma::mat myTestsampleMat;
	arma::mat mysampleMatbase;
	arma::mat myTestsampleMatbase;
	std::string frm1 = baseForm;
	vec bestmodel;
	vec bestmodelCoef;
	vec newmodel;
	vec bestpredict_train;
	vec bestpredict;
	vec singleTrainPredict;
	vec singleTestPredict;
	std::vector < std::string > vnames = ovnames;
	std::string gfrm1;
  	const int z_idi=0;
  	const int z_nri=1;
	int sizecases = casesample.n_rows;
	int sizecontrol = controlsample.n_rows;
	if (loops > 1)
	{
		uvec samCases= randi<uvec>(totSamples, distr_param(0,sizecases-1));
		uvec samControl=randi<uvec>(totSamples, distr_param(0,sizecontrol-1));
		vec auxcasesample=zeros<vec>(sizecases);
   		vec auxcontrolsample=zeros<vec>(sizecontrol);
  		mat myTestCases;
  		mat myTestControl;
		for (int i =0;i<totSamples;i++)
		{
			mysample=join_cols(mysample,casesample.row((samCases(i))));
			mysample=join_cols(mysample,controlsample.row((samControl(i))));
			auxcasesample[(samCases(i))]=1;
			auxcontrolsample[(samControl(i))]=1;
		}				
		myTestCases=casesample.rows(find(auxcasesample==0));
		myTestControl=controlsample.rows(find(auxcontrolsample==0));
		myTestsample=join_cols(myTestCases,myTestControl);
	}
	else
	{
		mysample=join_cols(casesample,controlsample);
  	  	myTestsample=join_cols(casesample,controlsample);
	}
			
    mat myoutcome(mysample.n_rows,2);
	myoutcome.col(1)=mysample.col(lockup[Outcome]);
       
	mysampleMatbase.set_size(mysample.n_rows);
	myTestsampleMatbase.set_size(myTestsample.n_rows);
	mysampleMatbase.fill(1.0);
	myTestsampleMatbase.fill(1.0);
	if (type == "COX")
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
	{
        bestmodel =modelFittingFunc(myoutcome,mysampleMatbase,type);
	}
    else
    {
       bestmodel=mean(myoutcome.col(1));
       if (type=="LOGIT") bestmodel=log(bestmodel/(1-bestmodel));
	}
  	bestpredict_train=predictForFresaFunc(bestmodel,mysampleMatbase,"prob",type);
  	bestpredict=predictForFresaFunc(bestmodel,myTestsampleMatbase,"prob",type);
	int changes = 1;
	unsigned int jmax = -1;
  	unsigned int j;
	std::string inname = "Inserted";
	while (changes>0)
	{
		changes = 0;
	  	maxrec = 0;
	  	jmax = -1;
	  	for ( j=0; j<size;j++)
	  	{
			vec iprob;
			vec iprob_t;
	  		if (vnames[j] != inname)
	  		{
				mysampleMat=join_rows(mysampleMatbase,mysample.col(lockup[vnames[j]]));
				myTestsampleMat=join_rows(myTestsampleMatbase,myTestsample.col(lockup[vnames[j]]));
	  			gfrm1 = frm1+" + "+vnames[j];
	  			singleTrainPredict.clear();
	  			singleTestPredict.clear();
             	if((mysampleMat.n_cols>1)or(type=="COX"))
             	newmodel =modelFittingFunc(myoutcome,mysampleMat,type);
			    else
			    {
        			newmodel=mean(myoutcome.col(1));
        			if (type=="LOGIT") 	newmodel=log(newmodel/(1-newmodel));
        		}
				if (is_finite(newmodel)) 
	  			{
  					singleTrainPredict=predictForFresaFunc(newmodel,mysampleMat,"prob",type);
					singleTestPredict=predictForFresaFunc(newmodel,myTestsampleMat,"prob",type);
					iprob_t = improveProbFunc(bestpredict_train,singleTrainPredict,mysample.col(lockup[Outcome]));            
  					iprob = improveProbFunc(bestpredict,singleTestPredict,myTestsample.col(lockup[Outcome]));            
  					if ((is_finite(iprob)) and (is_finite(iprob_t)))
					{
						if (selType=="zIDI") 
		  				{
		  					zmin = std::min(iprob(z_idi),iprob_t(z_idi));
		  				}
		  				else
		  				{
		  					zmin = std::min(iprob(z_nri),iprob_t(z_nri));
		  				}
	  					if (zmin  > maxrec)
	  					{
	  						jmax = j;
							maxrec = zmin;
	  					}
					}
	  			}
			}
	  	}
		if ((jmax >= 0) and (maxrec > zthr2) and (vnames[jmax] != inname) and (inserted<maxTrainModelSize))   
	  	{
			gfrm1 = frm1+" + ";
	  		gfrm1 = gfrm1+ovnames[jmax];
	  		mysampleMatbase=join_rows(mysampleMatbase,mysample.col(lockup[vnames[jmax]]));
			myTestsampleMatbase=join_rows(myTestsampleMatbase,myTestsample.col(lockup[vnames[jmax]]));
           	bestmodel =modelFittingFunc(myoutcome,mysampleMatbase,type);
	  		bestpredict_train=predictForFresaFunc(bestmodel,mysampleMatbase,"prob",type);
			bestpredict=predictForFresaFunc(bestmodel,myTestsampleMatbase,"prob",type);
           	frm1 = gfrm1;	
	  		changes = changes + 1;
	  		mynamesLoc.push_back(jmax);
	  		vnames[jmax] = inname;
	  		jmax = 0;
	  		maxrec = 0;
	  		inserted = inserted + 1;
	  	} 
	}
	return frm1;
}

extern "C" SEXP ReclassificationFRESAModelCpp(SEXP _size, SEXP _fraction,SEXP _pvalue, SEXP _loops,  SEXP _covariates, SEXP _Outcome, SEXP _variableList, SEXP _maxTrainModelSize, SEXP _type, SEXP _timeOutcome, SEXP _selType,SEXP _loop_threshold, SEXP _interaction, SEXP _dataframe,SEXP _colnames,SEXP _cores)
{
try {   // R_CStackLimit=(uintptr_t)-1;
 		#ifdef _OPENMP
		Rcout<<"::============ Available CPU Cores(Threads) : "<<omp_get_num_procs()<<" ===Requested Threads : "<< as<unsigned int>(_cores) << endl;
	    omp_set_num_threads(as<unsigned int>(_cores));
		#endif	
		unsigned int size = as<unsigned int>(_size); 
		double fraction = as<double>(_fraction);
		double pvalue = as<double>(_pvalue);
		int loops = as<int>(_loops);
		std::string Outcome = as<std::string>(_Outcome); 
		CharacterVector colnames(_colnames);
	    unsigned int maxTrainModelSize = Rcpp::as<unsigned int>(_maxTrainModelSize);
		std::string type = as<std::string>(_type);
		std::string timeOutcome = as<std::string>(_timeOutcome);
		std::string selType = as<std::string>(_selType);
		Rcpp::NumericMatrix DF(_dataframe);
    	arma::mat dataframe(DF.begin(), DF.rows(), DF.cols(), false);
		std::vector<std::string>covari=as<std::vector<std::string> >(_covariates);
		std::string baseForm=Outcome;
		std::map<std::string, int> lockup;
		std::vector<std::string> formulas;
		std::vector<int> mynames;
		std::vector<std::string> ovnames=as<std::vector<std::string> >(CharacterVector(_variableList)); 
		double zthr = abs(qnorm(pvalue,0.0,1.0)); //double zthr = abs(Rf_qnorm5(pvalue,0.0, 1.0, 1, 0));
		for(int i=0;i<CharacterVector(colnames).size();i++)
		{
			lockup[std::string(CharacterVector(colnames)[i])]=i;
		}
		if (type=="COX")
		{
		  baseForm = "Surv("+timeOutcome+","+Outcome+")";
		}
		baseForm = baseForm+" ~ 1";
		if (covari[0]!="1")
			for(unsigned int i=0;i<covari.size();i++)
			    baseForm = baseForm+"+"+covari[i];
		mat casesample=dataframe.rows(find(dataframe.col(lockup[Outcome])==1));
	    mat controlsample=dataframe.rows(find(dataframe.col(lockup[Outcome])==0));
		int sizecases   =casesample.n_rows;
		int sizecontrol =controlsample.n_rows;
		int minsize = std::min(sizecases,sizecontrol);
		int totSamples = (int)(fraction*minsize+0.49);
		double zthr2 = zthr;
		if (fraction<1) 
		{
			zthr2 = zthr*sqrt(fraction);
		}
		if (zthr2<abs(qnorm(0.1,0,1))) zthr2 = abs(qnorm(0.1,0,1)); 

		Rcout<<pvalue<<" <-pv : z-> "<<zthr2<<" Form: "<<  baseForm << "\n";

		if (size > ovnames.size()) size = ovnames.size();
		std::string inname = "Inserted";
		for(unsigned int j=0;j<covari.size();j++)
		    for(unsigned int i=0;i<ovnames.size();i++)
					 if (covari[j]==ovnames[i]) ovnames[i]=inname;  	
#pragma omp parallel for schedule(dynamic) ordered shared (mynames,formulas)  
    	for (int doOver=0; doOver<loops;doOver++)
	 	{   
			std::vector <int> mynamesLoc;
			std::string  frm = binaryFowardSelection(size,totSamples,zthr2,selType,loops,
			Outcome,timeOutcome,type,maxTrainModelSize,casesample,controlsample,
			ovnames,lockup,baseForm,mynamesLoc,covari);

			if ((doOver%10)==0) 
				{
#pragma omp critical
					Rcout << doOver <<" : "<< frm << std::endl;	
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

