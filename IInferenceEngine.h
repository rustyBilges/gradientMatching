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



#endif /* IINFERENCEENGINE_H_ */
