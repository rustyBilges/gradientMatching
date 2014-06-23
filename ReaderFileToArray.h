/*
 * ReaderFileToArray.h
 *
 *  Created on: 18 Dec 2013
 *      Author: rusty
 */

#ifndef READERFILETOARRAY_H_
#define READERFILETOARRAY_H_


#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>
using namespace std;

class ReaderFileToArray
{


private:

	int rowCount , columCount;
	double** array;
	const char* inputFileName;

public:
	~ReaderFileToArray();
	 ReaderFileToArray(int rowCount, int columnCount, const char* inputFileName);

	 void     readFile();
	 double** getArray();
	 double*  getRow(int ID);



};


#endif /* READERFILETOARRAY_H_ */
