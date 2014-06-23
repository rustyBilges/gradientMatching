/*
 * RangeSample.h
 *
 *  Created on: 8 Apr 2014
 *      Author: rusty
 *
 *      Class to handle a range sample of a full dynamical timeseries, based on the range of one of the species populations invloved
 *
 *      Each range is contains subsetCount subsets of the original data, that may be discontinuous in time.
 *      i.e. The population of the rangeSpecies may fall within the specified range several times over the course of the original timeseries.
 *
 *      originalTimepoints: the timepoints at which the original timeseries was observed. This is 2D. There is one 1D array for each subset in the range.
 *      observations: the populations of each species at orignialTimepoints. This is 3D. For each species there is one 1D array for each subset in the range.
 *
 *		subsetCount: the number of subsets in the range. (used for accessing the observations)
 *      sampleCount: the number of samples resulting from interpolation
 *      rangeSpeies: the ID of the species on which the range is generated (by default =0 for prey)
 */

#ifndef RANGESAMPLE_H_
#define RANGESAMPLE_H_

#include <vector>
using namespace std;

class RangeSample{

private:

	//int subsetCount,
	int sampleCount;
	int rangeID;
	int speciesCount;
	int rangeSpecies;
	double rangeMin, rangeMax;

	vector<double> tau;                       // interpolated timepoints. lenght=sampleCount;
	vector<vector<double> > interpPopulation; // interpolated estimates of x fore each species at each tau
	vector<vector<double> > interpRates;      // interpolated estimates of xdot fore each species at each tau

public:
	RangeSample(int rangeID, vector< vector<double> > originalTimepoints, vector<vector<vector<double> > > observations, double rangeMin, double rangeMax, int rangeSpecies=0){

		this->rangeID =rangeID;
		this->speciesCount = observations.size();
		//this->subsetCount = subsetCount;

		// check that RangeSample is of suitable dimension for inference.
		sampleCount = 0;

		if(!countSamples(originalTimepoints, sampleCount)){
			cerr << "Failure in RangeSample contstructor. Error in: " << rangeID << endl;
			throw 60;
		}
		// resize the vector data memebers according to the number of samples to be generated: (initialise all elements to zero)
		tau.resize(sampleCount, 0);
		interpPopulation.resize(speciesCount , vector<double>( sampleCount , 0 ) );
		interpRates.resize(speciesCount , vector<double>( sampleCount , 0 ) );
		// THEN FILL THESE BY CALL TO INTERPLATEOBSERVATIONS
		// DO THIS WITHIN THE CONSTRUCTOR?? - depends on ifthe exception handling is working...
		try {
			if(!interpolateObservations(originalTimepoints, observations)){
				throw 70;
			}
		}
		catch (int x){
			cerr << "Interpolation Failed." << endl;
		}


		this->rangeSpecies = rangeSpecies;
		this->rangeMin = rangeMin;
		this->rangeMax = rangeMax;
	}

	bool countSamples(vector< vector<double> > originalTimepoints, int& sampleCount){
		// before interpolation, counts how many samples will be generated.

		// iterate over subsets in the range:
		for (unsigned int i=0; i<originalTimepoints.size(); i++){
			// increment sampleCount by number of contributions resulting from interpolation in this subset:
			sampleCount += originalTimepoints.at(i).size()-1;
		}
		if (sampleCount<= speciesCount+1){
			cerr << "System is underconstrained. Too few samples resulting from interpolation of range: " << rangeID << endl;
			return false;
		}

		return true;
	}

	bool interpolateObservations(vector< vector<double> > originalTimepoints, vector<vector<vector<double> > > observations){

		int sampleIndex = 0;

		double tau0 = 0;
		double tau1 = 0;
		double xtau0 = 0;
		double xtau1 = 0;

		// iterate over subsets:
		for (unsigned int i=0; i<originalTimepoints.size(); i++){
			// iteratate over points within subset
			for (unsigned int j=0; j < originalTimepoints.at(i).size()-1; j++){

				tau0 = originalTimepoints.at(i).at(j);
				tau1 = originalTimepoints.at(i).at(j+1);
				tau.at(sampleIndex) = (tau1 + tau0)/2;

				for (int k = 0; k<speciesCount; k++){

					xtau0 = observations.at(k).at(i).at(j);
					xtau1 = observations.at(k).at(i).at(j+1);

					interpPopulation.at(k).at(sampleIndex) = (xtau0 + xtau1)/2;
					interpRates.at(k).at(sampleIndex) = (xtau1 - xtau0) / (tau1 - tau0);
				}

				sampleIndex++;
				if (sampleIndex > sampleCount){
					cerr << "Error in interpolation: sampleCount exceeded. Range ID: " << rangeID << endl;
					return false;
				}
			}
		}

		return true;
	}

	friend ostream& operator<< (ostream &out, RangeSample &sample);


};

ostream& operator<< (ostream &out, RangeSample &sample)
{
    // Since operator<< is a friend of the RangeSample class, we can access
    // sample's members directly.
    out << "range min, " << sample.rangeMin << endl;
    out << "range max, " << sample.rangeMax << endl;

    out << "sample count, " << sample.sampleCount << endl;
    out << "species count, " << sample.speciesCount << endl;

    // output interpolated times:
    out << "tau, ";
    for (int i=0; i< sample.sampleCount; i++){
    	out << sample.tau.at(i) << ", ";
    }
    out << endl;
    // output interpolated populations:
    for (int n=0; n<sample.speciesCount; n++){

    	out << "x" << n <<", ";
		for (int i=0; i< sample.sampleCount; i++){
			out << sample.interpPopulation.at(n).at(i) << ", ";
		}
		out << endl;
    }
    // output interpolated rates:
    for (int n=0; n<sample.speciesCount; n++){

    	out << "xdot" << n <<", ";
		for (int i=0; i< sample.sampleCount; i++){
			out << sample.interpRates.at(n).at(i) << ", ";
		}
		out << endl;
    }


    return out;
}


#endif /* RANGESAMPLE_H_ */
