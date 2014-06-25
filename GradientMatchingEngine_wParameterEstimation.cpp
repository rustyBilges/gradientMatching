/*
 * GradientMatchingEngine_wParameterEstimation.cpp
 *
 *  Created on: 23 Jun 2014
 *      Author: rusty
 */

#include "IInferenceEngine.h"

GradientMatchingEngine_wParameterEstimation::GradientMatchingEngine_wParameterEstimation(int unitCount, int timestepCount, int parametersPerUnit, vector<vector <double> >* interpolatedTimeseries, int simID){

	this->unitCount=unitCount;
	this->timestepCount = timestepCount;
	this->parametersPerUnit = parametersPerUnit;
	this->interpolatedTimeseries = interpolatedTimeseries;

	this->simID=simID;


	X = zeros<mat>(1,1);
	G = zeros<mat>(1,1);

	jHat.resize(unitCount, vector<double>(unitCount + parametersPerUnit*unitCount, 0));
//
//	this->jHat = new double*[unitCount];
//	for (int i=0; i < this->unitCount; i++){jHat[i] = new double[unitCount+(parametersPerUnit*unitCount)];}
//
//	this->fullTimeseries = new double*[nTimesteps];
//	int columnCount = unitCount + 1;
//	for (int i = 0; i <= nTimesteps; i++) {
//		fullTimeseries[i] = new double[columnCount];
//		for (int j=0; j<columnCount; j++){
//			fullTimeseries[i][j] = populationTimeseries.at(i).at(j);
//		}
//	}

}


bool GradientMatchingEngine_wParameterEstimation::runInference(int sampleStepLength){

		cout << interpolatedTimeseries->size() << endl;

		if ((int)interpolatedTimeseries->size()<=unitCount+1){
			cerr << "cannot run gradient matching. Too few observations -> system is not over-constrained. Please re-do interpolation." << endl;
			return false;
		}

		if(!fillJhat(sampleStepLength)){
			return false;
		}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!switch!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		saveResult();

		return true;
	}

// WRITE THIS:****************************************************
bool GradientMatchingEngine_wParameterEstimation::saveResult(){

//	string fname = "_GM_result.log";
//	stringstream str_ID;
//	str_ID << ID;
//	fname.insert(0, str_ID.str());
//
//	ofstream ofile;
//	ofile.open(fname.c_str());
//	if(!ofile){
//		cerr << "GM result log file not open" << endl;
//		return false;
//	}
//	ofile << "recontructed adjacency (Jhat) = " << endl << endl;
//	for (int i=0;i<unitCount;i++){
//
//		for (int j=0;j<unitCount;j++){
//			ofile << jHat[i][j] << "        ";
//		}
//		ofile << endl;
//	}
//	ofile << endl;
//
//	// this should be a separate function:
//	ofile << "recontructed parameters = " << endl << endl;
//	for (int i=0; i<unitCount; i++){ofile << jHat[i][unitCount] << endl;}
//
//	ofile << "*********************************************************" << endl;
//	ofile.close();


	return true;
}

bool GradientMatchingEngine_wParameterEstimation::fillJhat(int sampleStepLength){

			for (int i=0;i<unitCount;i++){
				mat temp = zeros(1,unitCount+(parametersPerUnit));

				if(!timmeJhatI(i, sampleStepLength,temp)){
					return false;
				}

				for (int j=0;j<unitCount+(parametersPerUnit);j++){
					jHat.at(i).at(j) = temp(0,j);
				}
			}

			return true;

		}

bool GradientMatchingEngine_wParameterEstimation::timmeJhatI(int unitID, int sampleStepLength, mat &JhatI){

			int t = 0;
			int t2;
//			int M = (timestepCount/sampleStepLength) - 1;
			int M = (int)interpolatedTimeseries->size();
			mat X = zeros<mat>(1,M);
			mat G = zeros<mat>(unitCount+(parametersPerUnit),M);

			double xtau[unitCount];

			double xdotTau;

			for (int m=0;m<M;m++){

				t2 = t + sampleStepLength;

				xtau[unitID] = interpolatedTimeseries->at(m).at(unitID+1);
				xdotTau      = interpolatedTimeseries->at(m).at(unitID + 1 + unitCount);


//				xtau[unitID] = (fullTimeseries[t2][1+unitID] + fullTimeseries[t][1+unitID])/2;
//				xdotTau = (fullTimeseries[t2][1+unitID] - fullTimeseries[t][1+unitID])/((fullTimeseries[t2][0] - fullTimeseries[t][0]));
				//xtau[unitID] = (histories[unitID]->getHistoryDiTi(0,t2)+histories[unitID]->getHistoryDiTi(0,t))/2;
				//xdotTau = (histories[unitID]->getHistoryDiTi(0,t2)-histories[unitID]->getHistoryDiTi(0,t))/(thistories[t2]-thistories[t]);

				X(0,m) = xdotTau;// - (xtau[unitID]*parameters[unitID]); // !! naughty workaround now
				//X(0,m) = xdotTau - (xtau[unitID]*Units[unitID]->getParamI(0)); // !! need unit parameters...


				for(int i=0;i<unitCount;i++){
					xtau[i] = interpolatedTimeseries->at(m).at(i+1);
							//							(fullTimeseries[t2][1+i] + fullTimeseries[t][1+i])/2;

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


void GradientMatchingEngine_wParameterEstimation::printJhat(){
	cout << "recontructed adjacency (Jhat) = " << endl << endl;
	for (int i=0;i<unitCount;i++){

		for (int j=0;j<unitCount;j++){
			cout << jHat.at(i).at(j) << "        ";
		}
		cout << endl;
	}
	cout << endl;

	// this should be a separate function:
	cout << "recontructed parameters = " << endl << endl;
	for (int i=0; i<unitCount; i++){cout << jHat.at(i).at(unitCount) << endl;}

	cout << "*********************************************************" << endl;
}

bool GradientMatchingEngine_wParameterEstimation::getEstimatedIM(vector<vector<double> >& estIM){

	if ((int)estIM.size()!=unitCount || (int)estIM.at(0).size()!=unitCount){
		cerr << "estimated interaction matrix could not be obtained: size mismatch." << endl;
		return false;
	}

	for (int i=0; i<unitCount; i++){
		for (int j=0; j<unitCount; j++){
			estIM.at(i).at(j) = jHat.at(i).at(j);
		}
	}

	return true;
}
bool GradientMatchingEngine_wParameterEstimation::getEstimatedParams(vector<double>& estParams){

	if ((int)estParams.size()!=unitCount){
		cerr << "estimated parameters could not be obtained: size mismatch." << endl;
		return false;
	}

	for (int i=0; i<unitCount; i++){estParams.at(i) = jHat.at(i).at(unitCount);}
	return true;
}


bool GradientMatchingEngine_wParameterEstimation::qualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements){
	return true;
}

bool GradientMatchingEngine_wParameterEstimation::percentageQualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements){
	return true;
}

bool GradientMatchingEngine_wParameterEstimation::originalQualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements){
	return true;
}

bool GradientMatchingEngine_wParameterEstimation::unweightedQualityOfReconstruction(double threshold, double& quality, vector< vector<double > > adjacency){
	return true;
}

double GradientMatchingEngine_wParameterEstimation::getSampleRate(int sampleStepLength){
	return 0;
}


