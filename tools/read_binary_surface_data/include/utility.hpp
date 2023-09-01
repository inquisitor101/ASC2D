#pragma once

#include <sstream>
#include <iostream>
#include <vector>


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


// Function that terminates program immediately.
void Terminate(const char        *functionName,
               const char        *fileName,
               const int          lineNumber,
               const std::string &errorMessageI);

// Function that swaps bytes.
void SwapBytes(void   *buffer,
		           size_t  nBytes,
							 int     nItems);
