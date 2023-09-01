#pragma once 

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
//#include <algorithm>

#ifdef HAVE_OPENMP
	#include <omp.h> 
#endif

// Define typical error and exit log message.
#define ERROR( msg ) Terminate( __FUNCTION__, __FILE__, __LINE__, msg )



// Type definition of double.
typedef double as3double;

// Vector of data type: 1D.
template<typename T>
using as3vector1d = std::vector<T>;

// Vector of data type: 2D.
template<typename T>
using as3vector2d = std::vector<std::vector<T>>;

// Vector of data type: 3D.
template<typename T>
using as3vector3d = std::vector<std::vector<std::vector<T>>>;

// Vector of data type: 4D.
template<typename T>
using as3vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>;


// Magic number in the restart file when using the AS3 format.
const int AS3_MagicNumber  = 535532;
// Length of the strings to store the variables names in the restart files.
const int CGNS_STRING_SIZE = 33;

// Total number of variables per DOF read. 
// Note, these are fixed to the conservative variables in 2D:
// i.e. the density, momentum and total energy.
const unsigned short nVar = 4;

// Specific heat ratio.
const as3double GAMMA = 1.4;
// Abbreviation involving gamma.
const as3double GAMMA_MINUS_ONE = GAMMA - 1.0;

// The expected number of elements, according to the initial header read.
extern unsigned long  nElemExpected;
// The expected number of nodes, according to the initial header read.
extern unsigned short nNodeExpected;

// Function that pre-processes the data information and returns the 
// number of chunks of data distribution with their starting/ending indices.
as3vector2d<unsigned long> PreprocessDataDistribution(const char         *dirA,
		                                                  const char         *dirB,
																		                  const char         *filename,
																		                  const unsigned long nData);

// Function, which checks if all the specified files exist.
void CheckFilesExist(const char   *directory,
		                 const char   *filename,
										 unsigned long nFile);

// Function, which reads from a reference input file the 
// necessary header information. 
// This also returns nElem*nNode in this current directory.
unsigned long DeduceHeaderInformation(const char    *directory,
		                                  const char    *filename,
														          unsigned long  iFile);


// Function that terminates program immediately.
void Terminate(const char        *functionName,
               const char        *fileName,
               const int          lineNumber,
               const std::string &errorMessageI);

// Function that swaps bytes.
void SwapBytes(void   *buffer,
		           size_t  nBytes,
							 int     nItems);

