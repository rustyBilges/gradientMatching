/*
 * dataBin.h
 *
 *  Created on: 25 Jun 2014
 *      Author: rusty
 *
 *      Class to use in bining of population timeseries data according to prey densities
 *
 *      NOTES:
 *      > if bining according to a single species density then the dimension=1
 *      > initially we will do naive bining: uniform grid. However should implement kd-tree algorithm
 *
 */

#ifndef DATABIN_H_
#define DATABIN_H_

#include <vector>
using namespace std;

class dataBin{

private:
	int binID, dimension, memberCount;
	vector <vector<double> > boundaries;   //size: dimension x 2 (i.e. min&max for each dimension)
	vector<double> binCentre;
//	vector<double> binMean;                // mean calculated from data members

	vector<int> memeberIDs;                // corresponding to the points in the timeseries which are members of this bin

public:

	dataBin(int binID, int dimension, vector <vector<double> >& boundaries):binID(binID), dimension(dimension), memberCount(0), boundaries(boundaries){

//		boundaries.resize(dimension, vector<double>(dimension, 2));

		binCentre.resize(dimension,0);
		for (int d=0;d<dimension; d++){
			binCentre.at(d) = (boundaries.at(d).at(0)+boundaries.at(d).at(1))/2;
		}
	}

	void store(int ID){
		memeberIDs.push_back(ID);
		memberCount++;
	}

	bool runInference(vector<vector <double> >& fullPopulationTimeseries, vector<vector <double> >& estimatedIM){

		int unitCount = (int)fullPopulationTimeseries.at(0).size() - 1;

		vector<vector <double> > interpolatedTimeseries;
		if(!interpolate(fullPopulationTimeseries, interpolatedTimeseries)){
			return false;
		}

		IInferenceEngine* IE = new GradientMatchingEngine_wParameterEstimation(unitCount, 0, 1, &interpolatedTimeseries);
		IE->runInference(1);
//		IE->printJhat();
		IE->getEstimatedIM(estimatedIM);
		delete IE;


		return true;
	}

	void printCentreDi(int dimension, ofstream& os){
		os << binCentre.at(dimension);
	}
private:

	bool isMember(int t){

		// determine if t is in memeberIDs. could use a faster search (binary tree?)
		for (int m=0; m<memberCount; m++){
			if(t==memeberIDs.at(m)){return true;}
		}

		return false;
	}

	bool interpolate(vector<vector <double> >& populationTimeseries, vector<vector <double> >& interpolatedTimeseries){


		unsigned int unitCount = populationTimeseries.at(0).size() - 1;
//		interpolatedTimeseries.resize(populationTimeseries.size()-1, vector<double>(unitCount*2 + 1, 0));
//		interpolatedTimeseries.resize(1, vector<double>(unitCount*2 + 1, 0));
		interpolatedTimeseries.clear(); // ensure it is empty
		vector<double> tempStore(unitCount*2 + 1, 0);

		int subsetMemeberCount = 0;

		for (unsigned int t=0; t<populationTimeseries.size(); t++){

			if(isMember(t)){
				subsetMemeberCount++;

				if (subsetMemeberCount>1){
					tempStore.at(0) = (populationTimeseries.at(t-1).at(0)+populationTimeseries.at(t).at(0))/2;

					for (int i=0; i<(int)unitCount; i++){
						tempStore.at(i+1) = (populationTimeseries.at(t-1).at(i+1)+populationTimeseries.at(t).at(i+1))/2;
						tempStore.at(i+1+unitCount) = (populationTimeseries.at(t).at(i+1)-populationTimeseries.at(t-1).at(i+1))/(populationTimeseries.at(t).at(0)-populationTimeseries.at(t-1).at(0));
					}
					interpolatedTimeseries.push_back(tempStore);
				}
			}
			else{
				subsetMemeberCount=0;
			}

		}

		return true;
	}

};


#endif /* DATABIN_H_ */
