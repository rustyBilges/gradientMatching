/*
 * ReaderFileToArray.h
 *
 *  Created on: 18 Dec 2013
 *      Author: rusty
 */

#include "ReaderFileToArray.h"

ReaderFileToArray::~ReaderFileToArray(){
		delete [] array;
	}
ReaderFileToArray::ReaderFileToArray(int rowCount, int columnCount, const char* inputFileName){
		 this -> rowCount = rowCount;
		 this -> columCount = columnCount;
		 this -> inputFileName = inputFileName;


		 array = new double*[rowCount];
		 for (int i=0; i<rowCount; i++){
			 array[i]   = new double[columnCount];
		 }
	 }

void ReaderFileToArray::readFile(){
		int countRows = 0;

		ifstream inputFileStream;
		 inputFileStream.open(this->inputFileName);
		 if(!inputFileStream){cout<< "file not open" << endl; exit(1);}
		 //else                     {cout<< "file is open!!" << endl;}

		 double i = 1000;
		 int indexI = 0;
		 int indexJ = 0;
		 //statesFileStream >> i;
		 while (inputFileStream.good()){
			 inputFileStream >> i;
			 //cout << i << endl;
			 array[indexI][indexJ] = i;
			 indexJ++;
			 //cout << i << endl;
			 if (indexJ==columCount){indexJ=0; indexI++; countRows++;}
		 }
		 if (countRows!=rowCount){cout<< "Number of Rows read not equal to number specified." << endl;}

		 inputFileStream.close();

	 }

double** ReaderFileToArray::getArray(){return array;}
double*  ReaderFileToArray::getRow(int ID){return array[ID];}




