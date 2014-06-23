/*
 * TimmeParameterInferenceEngine.cpp
 *
 *  Created on: 12 Jan 2014
 *      Author: rusty
 */

#include "TimmeParameterInferenceEngine.h"


TimmeParameterInferenceEngine::TimmeParameterInferenceEngine(int ID, int unitCount, int nTimesteps, int unitDimension, int parametersPerUnit, IDynamicalUnit** units, double timestepSize){
	this->ID = ID;
	this->unitCount = unitCount;
	this->nTimesteps = nTimesteps;
	this->unitDimension = unitDimension;
	this->parametersPerUnit = parametersPerUnit;
	this->timeseriesFileName = "null";

	X = zeros<mat>(1,1);
	G = zeros<mat>(1,1);

	this->jHat = new double*[unitCount];
	for (int i=0; i < this->unitCount; i++){jHat[i] = new double[unitCount+(parametersPerUnit*unitCount)];}

	this->fullTimeseries = new double*[nTimesteps];
	int columnCount = unitCount * unitDimension + 1;
	for (int i = 0; i < nTimesteps; i++) {
		fullTimeseries[i] = new double[columnCount];
		for (int j=0; j<columnCount; j++){
			fullTimeseries[i][j] = 0;
		}
	}

	this->timestepSize = timestepSize;
	this->units =units;


}

TimmeParameterInferenceEngine::TimmeParameterInferenceEngine(int unitCount, int nTimesteps, int unitDimension, int parametersPerUnit, const char* timeseriesFileName, IDynamicalUnit** dummyUnits){

		this->unitCount = unitCount;
		this->nTimesteps = nTimesteps;
		this->unitDimension = unitDimension;
		this->parametersPerUnit = parametersPerUnit;
		this->timeseriesFileName = timeseriesFileName;
//		this->parameters = parameters;

		X = zeros<mat>(1,1);
		G = zeros<mat>(1,1);

		this->jHat = new double*[unitCount];
		this->fullTimeseries = new double*[nTimesteps];

		for (int i=0; i < this->unitCount; i++){jHat[i] = new double[unitCount+(parametersPerUnit*unitCount)];}

		readTimseriesToArray();

		this->timestepSize = 0;
		this->units = dummyUnits;

	}

TimmeParameterInferenceEngine::TimmeParameterInferenceEngine(int unitCount, int nTimesteps, int unitDimension, int parametersPerUnit, vector<vector <double> >& populationTimeseries, IDynamicalUnit** dummyUnits){

	this->ID = 0;
	this->unitCount = unitCount;
//	this->nTimesteps = nTimesteps;
	this->nTimesteps = (int)populationTimeseries.size();
	this->unitDimension = unitDimension;
	this->parametersPerUnit = parametersPerUnit;
	this->timeseriesFileName = "not needed";

	X = zeros<mat>(1,1);
	G = zeros<mat>(1,1);

	this->jHat = new double*[unitCount];
	for (int i=0; i < this->unitCount; i++){jHat[i] = new double[unitCount+(parametersPerUnit*unitCount)];}

	this->fullTimeseries = new double*[nTimesteps];
	int columnCount = unitCount + 1;
	for (int i = 0; i <= nTimesteps; i++) {
		fullTimeseries[i] = new double[columnCount];
		for (int j=0; j<columnCount; j++){
			fullTimeseries[i][j] = populationTimeseries.at(i).at(j);
		}
	}


	this->timestepSize = 0;
	this->units = dummyUnits;

}

bool TimmeParameterInferenceEngine::runInference(int sampleStepLength){
		// fill timeseries:
//		readUnitHistoriesToArray();

	// REPLACE THIS!! DATA BODGE..
//	double replaceZero = 100;
//	int columnCount = unitCount * unitDimension + 1;
//	for (int i = 0; i<nTimesteps; i++){
//		for (int j=0; j<columnCount; j++){
//			if (fullTimeseries[i][j]==0){
//				fullTimeseries[i][j] = replaceZero;
//			}
//		}
//	}
//

		if(!fillJhat(sampleStepLength)){
			return false;
		}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!switch!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		saveResult();

		return true;
	}

bool TimmeParameterInferenceEngine::saveResult(){

	string fname = "_GM_result.log";
	stringstream str_ID;
	str_ID << ID;
	fname.insert(0, str_ID.str());

	ofstream ofile;
	ofile.open(fname.c_str());
	if(!ofile){
		cerr << "GM result log file not open" << endl;
		return false;
	}
	ofile << "recontructed adjacency (Jhat) = " << endl << endl;
	for (int i=0;i<unitCount;i++){

		for (int j=0;j<unitCount;j++){
			ofile << jHat[i][j] << "        ";
		}
		ofile << endl;
	}
	ofile << endl;

	// this should be a separate function:
	ofile << "recontructed parameters = " << endl << endl;
	for (int i=0; i<unitCount; i++){ofile << jHat[i][unitCount] << endl;}

	ofile << "*********************************************************" << endl;
	ofile.close();


	return true;
}

bool TimmeParameterInferenceEngine::fillJhat(int sampleStepLength){

			for (int i=0;i<unitCount;i++){
				mat temp = zeros(1,unitCount+(parametersPerUnit));

				if(!timmeJhatI(i, sampleStepLength,temp)){
					return false;
				}

				for (int j=0;j<unitCount+(parametersPerUnit);j++){
					jHat[i][j] = temp(0,j);
				}
			}

			return true;

		}

bool TimmeParameterInferenceEngine::timmeJhatI(int unitID, int sampleStepLength, mat &JhatI){

			int t = 0;
			int t2;
			int M = (nTimesteps/sampleStepLength) - 1;
			mat X = zeros<mat>(1,M);
			mat G = zeros<mat>(unitCount+(parametersPerUnit),M);
			double xtau[unitCount];

			double xdotTau;

			for (int m=0;m<M;m++){

				t2 = t + sampleStepLength;

				xtau[unitID] = (fullTimeseries[t2][1+unitID] + fullTimeseries[t][1+unitID])/2;
				xdotTau = (fullTimeseries[t2][1+unitID] - fullTimeseries[t][1+unitID])/((fullTimeseries[t2][0] - fullTimeseries[t][0]));
				//xtau[unitID] = (histories[unitID]->getHistoryDiTi(0,t2)+histories[unitID]->getHistoryDiTi(0,t))/2;
				//xdotTau = (histories[unitID]->getHistoryDiTi(0,t2)-histories[unitID]->getHistoryDiTi(0,t))/(thistories[t2]-thistories[t]);

				X(0,m) = xdotTau;// - (xtau[unitID]*parameters[unitID]); // !! naughty workaround now
				//X(0,m) = xdotTau - (xtau[unitID]*Units[unitID]->getParamI(0)); // !! need unit parameters...


				for(int i=0;i<unitCount;i++){
					xtau[i] = (fullTimeseries[t2][1+i] + fullTimeseries[t][1+i])/2;

					//G(i,m) = sin(xtau[i]);
					G(i,m) = xtau[i] * xtau[unitID];
				}

				G(unitCount,m) = xtau[unitID];
				t = t2;

			}

			//mat Jhat = X*trans(G)*inv(G*trans(G));
			// INSERT CHECK FOR IF THIS OPERATION IS SUCCESSFUL..
			mat Iv_temp;
			if(!inv(Iv_temp, G*trans(G))){
				cerr << "inference failure: sampled matrix singular" << endl;
				return false;
			}

			JhatI = X*trans(G)*Iv_temp;
			return true;
		}

double** TimmeParameterInferenceEngine::getJhat(){return jHat;}

void TimmeParameterInferenceEngine::printJhat(){
			cout << "recontructed adjacency (Jhat) = " << endl << endl;
			for (int i=0;i<unitCount;i++){

				for (int j=0;j<unitCount;j++){
					cout << jHat[i][j] << "        ";
				}
				cout << endl;
			}
			cout << endl;

			// this should be a separate function:
			cout << "recontructed parameters = " << endl << endl;
			for (int i=0; i<unitCount; i++){cout << jHat[i][unitCount] << endl;}

			cout << "*********************************************************" << endl;
		}


//		void printQ(double accuracy, INetwork *net){
//			cout << "quality of reconstruction = " << qualityOfReconstruction(accuracy, net);
//			cout << "  at accuracy: " << accuracy << endl << endl;
//		}


bool TimmeParameterInferenceEngine::originalQualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements){
// this only considers inter-specific interactions, and implements the oringal TImme metric.
			quality = 0;
			double deltaJij = 0;
			for (int i = 0;i<unitCount;i++){
				for (int j = 0;j<unitCount;j++){
					if (i!=j){
						deltaJij = abs(jHat[i][j] - IMelements.at((i*2)+j))/(2*max(abs(jHat[i][j]), abs(IMelements.at((2*i)+j))));
						quality = quality + heaviside((1-accuracy)-deltaJij);
					}
				}
			}
			quality = quality / (unitCount*unitCount - unitCount);


			return true;
		}


bool TimmeParameterInferenceEngine::qualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements){

			quality = 0;
			double deltaJij = 0;
			for (int i = 0;i<unitCount;i++){
				for (int j = 0;j<unitCount;j++){

					deltaJij = abs(jHat[i][j] - IMelements.at((i*2)+j))/(2*max(abs(jHat[i][j]), abs(IMelements.at((2*i)+j))));
					//quality = quality + heaviside((1-accuracy)-deltaJij);
					// just using normalised elemnt wise difference:
					quality = quality + deltaJij;
				}
			}
			quality = quality / (unitCount*unitCount);

			// metric that only considers the inter-specific interactions...
//			quality = 0;
//			double deltaJij = 0;
//			for (int i = 0;i<unitCount;i++){
//				for (int j = 0;j<unitCount;j++){
//					if (i!=j){
//						deltaJij = abs(jHat[i][j] - IMelements.at((i*2)+j))/(2*max(abs(jHat[i][j]), abs(IMelements.at((2*i)+j))));
//						//quality = quality + heaviside((1-accuracy)-deltaJij);
//						// just using normalised elemnt wise difference:
//						quality = quality + deltaJij;
//					}
//				}
//			}
//			quality = quality / (unitCount);


			return true;
		}

bool TimmeParameterInferenceEngine::percentageQualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements){
// for two species, only considering inter-specific interactions:
			quality = 0;
			double deltaJij = 0;
			for (int i = 0;i<unitCount;i++){
				for (int j = 0;j<unitCount;j++){
					if (i!=j){
						deltaJij = abs(jHat[i][j] - IMelements.at((i*2)+j))/(abs(IMelements.at((2*i)+j)));
						//quality = quality + heaviside((1-accuracy)-deltaJij);
						// just using normalised elemnt wise difference:
						quality = quality + deltaJij;
					}
				}
			}
			quality = quality / (unitCount);
			return true;
		}


bool TimmeParameterInferenceEngine::unweightedQualityOfReconstruction(double threshold, double& quality, vector< vector<double > > adjacency){


	if (unitCount!=(signed int)adjacency.size()){
		cerr << "cannot calculate quality: adjacency given is of wrong dimension" << endl;
		return false;
	}
	// do thresholding here or in inference method??
	quality = 0;
	double estimate = 0;
	for (unsigned int i =0; i< adjacency.size(); i++){
		for (unsigned int j=0; j<adjacency.size(); j++){
			if (i!=j){
				if (abs(jHat[i][j])<threshold){estimate=0;}
				else {estimate=jHat[i][j];}
				if (estimate*adjacency.at(i).at(j)>0){quality++;}
//				else if (adjacency.at(i).at(j)==0 && estimate==0){quality++;}
				else if (abs(adjacency.at(i).at(j))<=threshold && estimate==0){quality++;}  //!!!
			}
		}
	}
	quality = quality / (unitCount*unitCount-unitCount);
	return true;
}

double TimmeParameterInferenceEngine::heaviside(double x){
			if (x >= 0){return 1.0;}
			else {return 0.0;}
		}

double TimmeParameterInferenceEngine::getSampleRate(int sampleStepLength){return 1/(fullTimeseries[sampleStepLength][0] - fullTimeseries[0][0]);}



