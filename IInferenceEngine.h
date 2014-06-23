/*
 * IInferenceEngine.h
 *
 *  Created on: 20 Dec 2013
 *      Author: rusty
 */

#ifndef IINFERENCEENGINE_H_
#define IINFERENCEENGINE_H_

#include "INetwork.h"

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <armadillo>
// depricated use of time.h, sstream, new

using namespace std;
using namespace arma;


class IInferenceEngine {

public:
		virtual ~IInferenceEngine(){}

		virtual bool runInference(int sampleStepLength)=0;
		//virtual double ** runInference(const char* timeseriesFileName, int nTimesteps, int unitDimension, int sampleStepLength) = 0;
		// sampleStepLength gives frequency to sample timeseries at (=1  => every measurement. Maximum.)

//		virtual void buildJhat() = 0;
//
		virtual bool fillJhat(int sampleStepLength) = 0;
//
//		virtual mat timmeJhatI(double *thistories, UnitHistory *histories[], int Tsteps, int dim, int SpeciesID, int sampleStep, IDynamicalUnit **Units) = 0;
//
//		virtual double **getJhat() = 0;
//
		virtual void printJhat() = 0;
//
//		virtual void printQ(double accuracy, INetwork *net) = 0;

//		updated version:
		virtual bool qualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements) = 0;

		virtual bool percentageQualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements)=0;

		virtual bool originalQualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements)=0;

		virtual bool unweightedQualityOfReconstruction(double threshold, double& quality, vector< vector <double > > adjacency)=0;

		virtual double getSampleRate(int sampleStepLength) = 0;
};

class GradientMatchingEngine_wParameterEstimation : public IInferenceEngine
{

private:
	int simID;
	int unitCount;
	int timestepCount;
	int parametersPerUnit;  // assumed to be 1?

	vector <vector <double> >* interpolatedTimeseries;
//	double **fullTimeseries;

	double **jHat;   // following Timme notation: jHat is the inferred adjacency matrix.
	mat X;
	mat G;


	bool saveResult();

public:
	~GradientMatchingEngine_wParameterEstimation(){
//		for (int i = 0; i < timestepCount; i++) {
//			delete[] fullTimeseries[i];
//		}
//		delete[] fullTimeseries;
		for (int i=0; i < this->unitCount; i++){
			delete[] jHat[i];
		}
		delete[] jHat;
	}

	// this constructor takes a timeseires vector to run inference on. THIS SHOULD BE USED IN ALL CASES!
	GradientMatchingEngine_wParameterEstimation(int unitCount, int timestepCount, int parametersPerUnit, vector<vector <double> >* interpolatedTimeseries, int simID=0);

	bool runInference(int sampleStepLength);

	bool fillJhat(int sampleStepLength);

	void printJhat();

	// updated from:
	//double qualityOfReconstruction(double accuracy, INetwork *net);
	// to:
	bool qualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements);

	bool percentageQualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements);

	bool originalQualityOfReconstruction(double accuracy, double& quality, vector<double > IMelements);

	bool unweightedQualityOfReconstruction(double threshold, double& quality, vector< vector<double > > adjacency);

	double getSampleRate(int sampleStepLength);

private:
	// helper functions, not accessible by base class pointer:
//	void runInferenceNoLoops(int sampleStepLength);

	bool timmeJhatI(int unitID, int sampleStepLength, mat &JHatI);

//	double **getJhat();

//	double heaviside(double x);
//
//	bool movingWindow(int windowLength, int unitID, double **result){
//		if (windowLength<(this->unitCount+1)){return 1;}
//
//		int M = (timestepCount) - windowLength;
//		mat X = zeros<mat>(1,windowLength);
//		mat G = zeros<mat>(unitCount+1,windowLength);
//		double xtau[unitCount];
//
//		double xdotTau;
//
//		mat jHat;
//
//		// iterate over all windows:
//		for (int i=0;i<M;i++){
//			cout << "window " << i << endl;
//			// iterate of entries in window:
//			for (int w=0; w<windowLength; w++){
//				xdotTau = (fullTimeseries[i+1+w][1+unitID] - fullTimeseries[i+w][1+unitID])/((fullTimeseries[i+1+w][0] - fullTimeseries[i+w][0]));
//				X(0,w) = xdotTau;
//
//				for(int j=0;j<unitCount;j++){
//					xtau[j] = (fullTimeseries[i+1+w][1+j] + fullTimeseries[i+w][1+j])/2;
//
//					G(j,w) = xtau[j] * xtau[unitID];
//				}
//				G(unitCount,w) = xtau[unitID];
//			}
//			jHat = X*trans(G)*inv(G*trans(G));
//
//			for (int j=0;j<unitCount+1;j++){
//				result[i][j] = jHat[j];
//			}
//		}
//
//
//		return true;
//	}
};

#endif /* IINFERENCEENGINE_H_ */
