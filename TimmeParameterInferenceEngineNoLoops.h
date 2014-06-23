/*
 * TimmeParameterInferenceEngineNoLoops.h
 *
 *  Created on: 2 May 2014
 *      Author: rusty
 */

#ifndef TIMMEPARAMETERINFERENCEENGINENOLOOPS_H_
#define TIMMEPARAMETERINFERENCEENGINENOLOOPS_H_

//CURRENTLY ONLY WORKS WITH READING FROM UNITS..

#include "IInferenceEngine.h"
#include "ReaderFileToArray.h"
// from simulationResources
#include "IDynamicalUnit.h"
#include "UnitHistory.h"

class TimmeParameterInferenceEngineNoLoops : public IInferenceEngine
{

private:
	int ID;

	int unitCount;
	int unitDimension;
	int nTimesteps;
	int parametersPerUnit;

	double **fullTimeseries;

	double **jHat;   // following Timme notation: jHat is the inferred adjacency matrix.
	mat X;
	mat G;

	const char* timeseriesFileName;
	// for reading directly from units:
	double timestepSize;
	IDynamicalUnit** units;

	void readTimseriesToArray() {
		int columnCount = unitCount * unitDimension + 1;

		ReaderFileToArray reader(nTimesteps, columnCount, timeseriesFileName);

		reader.readFile();
		cout << "Timseries File Read." << endl;
		double** tempStore = reader.getArray();

		for (int i = 0; i < nTimesteps; i++) {
			fullTimeseries[i] = new double[columnCount];
			for (int j = 0; j < columnCount; j++) {
				fullTimeseries[i][j] = tempStore[i][j];
				//cout << i << endl;
			}
		}
		//reader.~ReaderFileToArray();  // why does it not like this??
		cout << "Timeseries stored in local array" << endl;
	}

	void readUnitHistoriesToArray(){
		int columnCount = unitCount * unitDimension + 1;

		for (int i = 0; i < nTimesteps; i++) {
			//fullTimeseries[i] = new double[columnCount];
			fullTimeseries[i][0] = timestepSize*i;
			for (int j = 1; j < columnCount; j++) {
				fullTimeseries[i][j] = units[j-1]->getHistory()->getHistoryTi(i);
				//cout << i << endl;
			}
		}

//		cout << "Timeseries stored in local array" << endl;

	}

	bool saveResult();

public:
	~TimmeParameterInferenceEngineNoLoops(){
		for (int i = 0; i < nTimesteps; i++) {
			delete[] fullTimeseries[i];
		}
		delete[] fullTimeseries;
		for (int i=0; i < this->unitCount; i++){
			delete[] jHat[i];
		}
		delete[] jHat;
	}

	// this constructor reads time series from units that are still active (instead of from file):
	TimmeParameterInferenceEngineNoLoops(int ID, int unitCount, int nTimesteps, int unitDimension, int parametersPerUnit, IDynamicalUnit** units , double timestepSize);

	//TimmeParameterInferenceEngine(int unitCount, int nTimesteps, int unitDimension, int parametersPerUnit, const char* timeseriesFileName, double *parameters, IDynamicalUnit** dummyUnits);

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
	void runInferenceNoLoops(int sampleStepLength);

	bool timmeJhatI(int unitID, int sampleStepLength, mat &JHatI);

	double **getJhat();

	double heaviside(double x);

	bool movingWindow(int windowLength, int unitID, double **result){
		if (windowLength<(this->unitCount+1)){return 1;}

		int M = (nTimesteps) - windowLength;
		mat X = zeros<mat>(1,windowLength);
		mat G = zeros<mat>(unitCount+1,windowLength);
		double xtau[unitCount];

		double xdotTau;

		mat jHat;

		// iterate over all windows:
		for (int i=0;i<M;i++){
			cout << "window " << i << endl;
			// iterate of entries in window:
			for (int w=0; w<windowLength; w++){
				xdotTau = (fullTimeseries[i+1+w][1+unitID] - fullTimeseries[i+w][1+unitID])/((fullTimeseries[i+1+w][0] - fullTimeseries[i+w][0]));
				X(0,w) = xdotTau;

				for(int j=0;j<unitCount;j++){
					xtau[j] = (fullTimeseries[i+1+w][1+j] + fullTimeseries[i+w][1+j])/2;

					G(j,w) = xtau[j] * xtau[unitID];
				}
				G(unitCount,w) = xtau[unitID];
			}
			jHat = X*trans(G)*inv(G*trans(G));

			for (int j=0;j<unitCount+1;j++){
				result[i][j] = jHat[j];
			}
		}


		return true;
	}
};



#endif /* TIMMEPARAMETERINFERENCEENGINENOLOOPS_H_ */
