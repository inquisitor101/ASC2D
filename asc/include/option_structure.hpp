#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>

// Include file for runtime NaN catching.
#ifdef ENABLE_NAN_CHECK
#include <fenv.h>
#endif


// Forward declaration of CData to avoid compiler problems
// when aliasing as3element.
class CData;

// Type definition of double.
typedef double as3double;

// Element data structure.
using as3element = std::vector<CData*>;

// Vector of data type: 1D.
template<typename T>
using as3vector1d = std::vector<T>;

// Vector of data type: 2D.
template<typename T>
using as3vector2d = std::vector<std::vector<T>>;

// Vector of data type: 3D.
template<typename T>
using as3vector3d = std::vector<std::vector<std::vector<T>>>;


// Vector of pointer data type: 1D.
template<typename T>
using as3data1d = std::vector<T*>;

// Vector of pointer data type: 2D.
template<typename T>
using as3data2d = std::vector<std::vector<T*>>;

// Vector of pointer data type: 3D.
template<typename T>
using as3data3d = std::vector<std::vector<std::vector<T*>>>;


// Number of variables.
const unsigned short nVar = 4;
// Number of dimensions.
const unsigned short nDim = 2;
// Number of faces per element (for now all elements are assumed to be quadrilateral).
const unsigned short nFace = 4;
// Value of pi.
const as3double PI_CONSTANT = 4.0*atan(1.0);
// Epsilon value used as a threshold.
const as3double EPS_VALUE = 1e-15;
// Ratio of specific heats.
const as3double GAMMA = 1.4;
// Abbreviation for gamma-1.
const as3double GAMMA_MINUS_ONE = GAMMA-1.0;
// Gas constant.
const as3double GAS_CONSTANT = 287.058;
// Specific heat at constant volume.
const as3double CV_CONSTANT = GAS_CONSTANT/GAMMA_MINUS_ONE;
// Specific heat at constant pressure.
const as3double CP_CONSTANT = GAMMA*CV_CONSTANT;
// For vectorization, the vector length used must be a multiple of vecLen1D.
const size_t vecLen1D = 4;


// Enumerated type for face indices on an element.
enum FACE_ELEM_ID {
  IDX_UNKNOWN = 99, // unknown face index.
	IDX_SOUTH   = 0,  // south-most face.
	IDX_NORTH   = 1,  // north-most face.
	IDX_WEST    = 2,  // west-most face.
	IDX_EAST    = 3   // east-most face.
};

// Enumerated type for dimension index.
enum DIM_ID {
	XDIM = 0, // x dimension.
	YDIM = 1  // y-dimension.
};

// Enumerated type of nodal DOFs.
enum NODAL_DOFS {
	TYPE_DOF_UNKOWN = 99, // unknown type of nodal points.
	TYPE_DOF_EQD    = 0,  // equidistant nodal points.
	TYPE_DOF_LGL    = 1   // Legendre-Gauss-Lobato nodal points.
};

// Enumerated type of problem-dependent variables.
enum PROBLEM_EQUATIONS_VARIABLE {
  CONT_VAR   = 0,  // continuity variable.
	XMOM_VAR   = 1,  // x-momentum variable.
	YMOM_VAR   = 2,  // y-momentum variable.
	ENER_VAR   = 3,  // energy variable.
  CONTQ1_VAR = 4,  // q1[density] variable.
  XMOMQ1_VAR = 5,  // q1[x-momentum] variable.
  YMOMQ1_VAR = 6,  // q1[x-momentum] variable.
  ENERQ1_VAR = 7,  // q1[energy] variable.
};

// Enumerated type of PDE problem to solve.
enum TYPE_PDE_SOLVER {
	NO_SOLVER      = 99, // solver option not valid.
	SOLVER_LEE     = 0,  // linearized Euler equations.
	SOLVER_EE      = 1,  // nonlinear Euler equations.
	SOLVER_NS      = 2,  // Navier-Stokes.
};

// Enumerated type of riemann solver.
enum RIEMANN_SOLVER_TYPE {
  RIEMANN_UNKNOWN   = 99, // unknown riemann solver.
  RIEMANN_ROE       = 0,  // Roe's riemann solver.
  RIEMANN_RUSANOV   = 1,  // Rusanov's riemann solver.
  RIEMANN_ROEISMAIL = 2   // Roe-Ismail riemann solver.
};

// Enumerated type of initial conditions.
enum INITIAL_CONDITION_TYPE {
	IC_UNKNOWN     			   = 99, // unknown type of IC.
	IC_IGNORE      			   = 0,  // no IC.
	IC_FREE_STREAM 			   = 1,  // free-stream IC.
	IC_GAUSSIAN_PRESSURE   = 2,  // gaussian pressure IC.
  IC_ISENTROPIC_VORTEX   = 3,  // isentropic vortex IC.
  IC_ENTROPY_WAVE        = 4,  // entropy wave IC.
  IC_VORTEX_ROLLUP       = 5,  // vortex roll-up IC.
  IC_ACOUSTIC_PLANE_WAVE = 6 // acoustic plane-wave IC.
};

// Enumerated type of boundary condition.
enum BOUNDARY_CONDITION_TYPE {
	BC_UNKNOWN           = 99, // unknown type of BC.
	BC_INTERFACE         = 0,  // interface/periodic boundary.
	BC_SYMMETRY          = 1,  // inviscid/slip wall.
  BC_CBC_OUTLET        = 2,  // outlet CBC.
  BC_CBC_INLET         = 3,  // inlet CBC.
  BC_STATIC_INLET      = 4,  // static-condition inlet.
  BC_TOTAL_INLET       = 5,  // total-conditions inlet.
  BC_STATIC_OUTLET     = 6,  // static-condition subsonic outlet.
  BC_SUPERSONIC_OUTLET = 7,  // supersonic outlet.
  BC_SUPERSONIC_INLET  = 8   // supersonic inlet.
};

// Enumerated type of temporal scheme used.
enum TEMPORAL_SCHEME {
	TEMPORAL_SCHEME_UNKNOWN = 99, // temporal scheme not valid.
	TEMPORAL_SCHEME_LSRK4   = 0,  // low-storage 4th-order RK.
	TEMPORAL_SCHEME_CRK4    = 1,  // classic 4th-order RK.
  TEMPORAL_SCHEME_SSPRK3  = 2   // strong-stability-preserving RK3.
};

// Enumerated type for LSRK4 options.
enum LSRK4_OPTIONS {
  LSRK4_N_STORAGE = 1, // number of tentative storage.
  LSRK4_N_STAGES  = 5  // number of grid sweeps.
};

// Enumerated type for SSPRK3 options.
enum SSPRK3_OPTIONS {
  SSPRK3_N_STORAGE = 1, // Number of tentative storage.
  SSPRK3_N_STAGES  = 3  // Number of grid sweeps.
};

// Enumerated type of multizone strategies.
enum MULTIZONE_STRATEGY_TYPE {
  MULTIZONE_STRATEGY_MAIN      = 0, // only zone 0 is used.
  MULTIZONE_STRATEGY_WEST      = 1, // use only zone: 0 and 1.
  MULTIZONE_STRATEGY_EAST      = 2, // use only zone: 0 and 2.
  MULTIZONE_STRATEGY_SOUTH     = 3, // use only zone: 0 and 3.
  MULTIZONE_STRATEGY_NORTH     = 4, // use only zone: 0 and 4.
  MULTIZONE_STRATEGY_ALL       = 5, // uses all zones: 0 to 8.
  MULTIZONE_STRATEGY_HORIZONAL = 6, // uses zones: 0, 1 and 2.
  MULTIZONE_STRATEGY_VERTICAL  = 7  // uses zones: 0, 3 and 4.
};

// Enumerated type of zones.
enum ZONE_TYPE {
  ZONE_UNKNOWN  = 99, // unknown zone.
	ZONE_MAIN     = 0,  // zone corresponding to physical region
	ZONE_WEST     = 1,  // zone PML corresponding to west
	ZONE_EAST     = 2,  // zone PML corresponding to east
	ZONE_SOUTH    = 3,  // zone PML corresponding to south
	ZONE_NORTH    = 4,  // zone PML corresponding to north
	ZONE_CORNER_0 = 5,  // zone PML corresponding to SW-corner C0
	ZONE_CORNER_1 = 6,  // zone PML corresponding to SE-corner C1
	ZONE_CORNER_2 = 7,  // zone PML corresponding to NW-corner C2
	ZONE_CORNER_3 = 8   // zone PML corresponding to NE-corner C3
};

// Enumerated type of 3-zone horizontal strategy.
enum HORIZONTAL_ZONE_TYPE {
  HORIZONTAL_ZONE_MAIN = 0, // zone corresponding to physical region.
  HORIZONTAL_ZONE_WEST = 1, // zone corresponding to west.
  HORIZONTAL_ZONE_EAST = 2  // zone corresponding to east.
};

// Enumerated type of 3-zone vertical strategy.
enum VERTICAL_ZONE_TYPE {
  VERTICAL_ZONE_MAIN  = 0, // zone corresponding to physical region.
  VERTICAL_ZONE_SOUTH = 1, // zone corresponding to south.
  VERTICAL_ZONE_NORTH = 2  // zone corresponding to north.
};

// Enumerated type of processing data.
enum PROCSES_DATA {
  PROCESS_NOTHING = 0, // no processing is done.
  PROCESS_IC      = 1  // process data of IC.
};

// Enumerated type of number of nodes per nPoly=1 element.
enum ELEM_POINTS {
	N_POINTS_LINE          = 2, // line
	N_POINTS_TRIANGLE      = 3, // triangle
	N_POINTS_QUADRILATERAL = 4, // quadrilateral
	N_POINTS_TETRAHEDRON   = 4, // tetrahedron
	N_POINTS_HEXAHEDRON    = 8, // hexaheron
	N_POINTS_PYRAMID       = 5, // pyramid
	N_POINTS_PRISM         = 6  // prism
};

// Enumerated type of geometric entities based on VTK nomenclature.
enum GEO_TYPE {
  NONE          = 0,  // Point.
  VERTEX        = 1,  // VTK nomenclature for defining a vertex element.
  LINE          = 3,  // VTK nomenclature for defining a line element.
  TRIANGLE      = 5,  // VTK nomenclature for defining a triangle element.
  QUADRILATERAL = 9,  // VTK nomenclature for defining a quadrilateral element.
  TETRAHEDRON   = 10, // VTK nomenclature for defining a tetrahedron element.
  HEXAHEDRON    = 12, // VTK nomenclature for defining a hexahedron element.
  PRISM         = 13, // VTK nomenclature for defining a prism element.
  PYRAMID       = 14  // VTK nomenclature for defining a pyramid element.
};

// Enumerated type of data filtering.
enum FILTERING_TYPE {
  NO_FILTER          = 0, // No filtering used.
  EXPONENTIAL_FILTER = 1  // Exponential filter.
};

// Enumerated type of buffer layer.
enum BUFFER_LAYER_TYPE {
  NO_LAYER                = 0, // No layer.
  PML_LAYER               = 1, // PML layer.
  SPONGE_LAYER            = 2  // sponge layer.
};

// Enumerated type of modified boundary condition.
enum MODIFIED_BOUNDARY_CONDITION {
  NO_BC_MODIFICATION                  = 0, // No modification.
  GAUSSIAN_PRESSURE_BC_MODIFICATION   = 1, // Gaussian-pressure profile.
  SINUSOIDAL_VELOCITY_BC_MODIFICATION = 2  // Sinusoidal velocity profile.
};

// Function used for printing the type of the input zone.
std::string DisplayTypeZone(unsigned short iTypeZone);
// Function used for printing the side of the input boundary.
std::string DisplayBoundarySide(unsigned short iBoundary);
// Function used for printing the type of solver.
std::string DisplaySolverType(unsigned short iSolver);
// Function used for printing the type of BC.
std::string DisplayBoundaryType(unsigned short iBoundary);


// Function that removes blanks from a string.
void RemoveBlanks(std::string &lineBuf);

// Function that replaces tabs and return characters
// in a string with a blank.
void ReplaceTabsAndReturns(std::string &lineBuf);

// Function that returns a lowered-case version of input string.
void CreateLowerCase(std::string &lineBuf);

// Function that searches for and finds a string vector from input file.
bool FindStringVectorFromKeyword(std::ifstream            &su2file,
                          			 const char               *keyword,
                          			 std::vector<std::string> &valStringVector);

// Function that searches for and finds a string from input file.
bool FindStringFromKeyword(std::ifstream &su2file,
                          const char     *keyword,
                          std::string    &valString);

// Function that terminates program immediately.
void Terminate(const char        *functionName,
               const char        *fileName,
               const int          lineNumber,
               const std::string &errorMessageI);


// Function used for printing, usually debugging only. Works on 1D arrays.
template<class TValueType>
void PrintData
(
 const TValueType  *val,
 const unsigned int nRow,
 const unsigned int nCol,
 const char        *message = nullptr,
 const unsigned int nDigits = 6
)
 /*
	* Function used for printing data.
	*/
{
	if( message ) std::cout << message << std::endl;
	std::cout.precision(nDigits);
	for(unsigned int iRow=0; iRow<nRow; iRow++){
		for(unsigned int iCol=0; iCol<nCol; iCol++)
			std::cout << std::scientific << std::showpos << val[iRow*nCol+iCol] << " ";
		std::cout << std::endl;
	}
}

// Function used for printing, usually debugging only. Works on 2D arrays.
template<class TValueType>
void PrintData
(
 TValueType        **val,
 const unsigned int  nRow,
 const unsigned int  nCol,
 const char         *message = nullptr,
 const unsigned int  nDigits = 6
)
 /*
	* Function used for printing data.
	*/
{
	if( message ) std::cout << message << std::endl;
	std::cout.precision(nDigits);
	for(unsigned int iRow=0; iRow<nRow; iRow++){
		for(unsigned int iCol=0; iCol<nCol; iCol++)
			std::cout << std::scientific << std::showpos << val[iRow][iCol] << " ";
		std::cout << std::endl;
	}
}

// Function used for printing, usually debugging only. Works on 2D vectors.
template<class TValueType>
void PrintData
(
 as3vector2d<TValueType> &val,
 const char        			 *message = nullptr,
 const unsigned int			  nDigits = 6
)
 /*
	* Function used for printing data.
	*/
{
	if( message ) std::cout << message << std::endl;
	std::cout.precision(nDigits);
	for(unsigned int iRow=0; iRow<val.size(); iRow++){
		for(unsigned int iCol=0; iCol<val[iRow].size(); iCol++)
			std::cout << std::scientific << std::showpos << val[iRow][iCol] << " ";
		std::cout << std::endl;
	}
}


// Function that determines a boolean input.
void AddBoolOption(std::ifstream &inputFile,
                   const char    *keyword,
                   bool          &value,
                   std::string   &defaultValue,
                   bool           ResetLine = false);


// Templated function used in reading "scalar" data, takes default parameter.
template<class TValueType>
void AddScalarOption
(
  std::ifstream  &inputFile,
  const char     *keyword,
  TValueType     &value,
	TValueType      defaultValue,
  bool            ResetLine = false
 )
  /*
   * Function that assigns a value from input file, otherwise uses default.
   */
 {
  // Initialize string.
  std::string valString;
  // Initialize line number.
  int CurrentLine;
  // Record current line.
  if( ResetLine ) CurrentLine = inputFile.tellg();

  // Use default value, if string could not be found.
  if( !FindStringFromKeyword(inputFile, keyword, valString) ){
		// Assign default value.
    value = defaultValue;
    // Go back to previous line location.
    if( ResetLine ) { inputFile.clear(); inputFile.seekg(CurrentLine); }
    return;
  }

  // Assign string value.
  std::istringstream istr(valString);
  istr >> value;

	// Check if any bad character is input.
	if( istr.fail() )
		Terminate("AddScalarOption", __FILE__, __LINE__,
							"Wrong character input inside expected parameter type!");

  // Go back to previous line location.
  if( ResetLine ) inputFile.seekg(CurrentLine);
}


// Templated function used in reading "scalar" data.
template<class TValueType>
void AddScalarOption
(
  std::ifstream  &inputFile,
  const char     *keyword,
  TValueType     &value,
  bool            ResetLine = false
 )
  /*
   * Function that assigns a value from input file, otherwise exits.
   */
 {
  // Initialize string.
  std::string valString;
  // Initialize line number.
  int CurrentLine;
  // Record current line.
  if( ResetLine ) CurrentLine = inputFile.tellg();

  // Use default value, if string could not be found.
  if( !FindStringFromKeyword(inputFile, keyword, valString) ){
    // Output message.
		std::string message = "Keyword: ";
		message += keyword;
		message += " not found!";
		// Exit program.
		Terminate("AddScalarOption", __FILE__, __LINE__, message);
  }

  // Assign string value.
  std::istringstream istr(valString);
  istr >> value;

	// Check if any bad character is input.
	if( istr.fail() )
		Terminate("AddScalarOption", __FILE__, __LINE__,
							"Wrong character input inside expected parameter type!");

  // Go back to previous line location.
  if( ResetLine ) inputFile.seekg(CurrentLine);
}


// Function that works as a template to read vector data, otherwise assigns default value.
template<class TValueType>
void AddVectorOption
(
 std::ifstream           &inputFile,
 const char              *keyword,
 std::vector<TValueType> &value,
 std::vector<TValueType>  defaultValue,
 bool                     ResetLine = false
)
 /*
	* Function that searched for keyword from input file and assigns it to
	* the templated value, otherwise uses default value.
	*/
{
  // Initialize string.
	std::vector<std::string> valStringVector;
  // Initialize line number.
  int CurrentLine;
  // Record current line.
  if( ResetLine ) CurrentLine = inputFile.tellg();

  // Use default value, if string could not be found.
  if( !FindStringVectorFromKeyword(inputFile, keyword, valStringVector) ){
		// Assign default value.
		value = defaultValue;
    // Go back to previous line location.
    if( ResetLine ) { inputFile.clear(); inputFile.seekg(CurrentLine); }
    return;
  }

	// Assign string value.
	for(unsigned short iData=0; iData<valStringVector.size(); iData++){
		TValueType tmp;
		std::istringstream istr(valStringVector[iData]);
		istr >> tmp;

		// Check if any bad character is input.
		if( istr.fail() )
			Terminate("AddVectorOption", __FILE__, __LINE__,
								"Wrong character input inside expected parameter type!");

		// Add to input value vector.
		value.push_back(tmp);
	}

  // Go back to previous line location.
  if( ResetLine ) inputFile.seekg(CurrentLine);
}


// Function that works as a template to read vector data, otherwise exits.
template<class TValueType>
void AddVectorOption
(
 std::ifstream           &inputFile,
 const char              *keyword,
 std::vector<TValueType> &value,
 bool                     ResetLine = false
)
 /*
	* Function that searched for keyword from input file and assigns it to
	* the templated value, otherwise exits.
	*/
{
  // Initialize string.
	std::vector<std::string> valStringVector;
  // Initialize line number.
  int CurrentLine;
  // Record current line.
  if( ResetLine ) CurrentLine = inputFile.tellg();

  // Use default value, if string could not be found.
  if( !FindStringVectorFromKeyword(inputFile, keyword, valStringVector) ){
    // Output message.
		std::string message = "Keyword: ";
		message += keyword;
		message += " not found!";
		// Exit program.
		Terminate("AddStringOption", __FILE__, __LINE__, message);
  }

	// Assign string value.
	for(unsigned short iData=0; iData<valStringVector.size(); iData++){
		TValueType tmp;
		std::istringstream istr(valStringVector[iData]);
		istr >> tmp;

		// Check if any bad character is input.
		if( istr.fail() )
			Terminate("AddVectorOption", __FILE__, __LINE__,
								"Wrong character input inside expected parameter type!");

		// Add to input value vector.
		value.push_back(tmp);
	}

  // Go back to previous line location.
  if( ResetLine ) inputFile.seekg(CurrentLine);
}


// Function that pads remaining entries in a vector.
template<class TValueType>
void PadEntriesVectorData
(
  std::vector<TValueType> &data,
  std::string        			 keyword,
  unsigned short           nExpected,
  unsigned short           nCondition1 = 0,
  unsigned short           nCondition2 = 0,
  unsigned short           nCondition3 = 0
)
 /*
  * Function that pds the remaining entries in a given vector of templated value.
  */
{
  // Size of current data vector.
  unsigned short nData = data.size();

  // Check first condition.
  if( nData != nExpected ){

    // Data flag.
    bool PadData = false;
    // Check if there need be any modification, if so flag the data.
    if( nData == nCondition1 ) PadData = true;
    if( nData == nCondition2 ) PadData = true;
    if( nData == nCondition3 ) PadData = true;

    // Pad data accordingly.
    if( PadData ) for(unsigned short i=nData; i<nExpected; i++) data.push_back( data[i-1] );
  }

  // Consistency check.
  if( data.size() != nExpected ){
    std::string message = "Keyword: ";
    message += keyword;
    message += " must be of size: ";
    message += std::to_string(nExpected);
    Terminate("PadEntriesVectorData", __FILE__, __LINE__, message);
  }
}


// Function that pads remaining entries in a vector based on default reference data.
template<class TValueType>
void PadEntriesVectorDefaultData
(
  std::vector<TValueType> &data,
  std::vector<TValueType> &reference,
  std::string        			 keyword
)
 /*
  * Function that pds the remaining entries in a given vector of templated value
  * based on a default reference value.
  */
{
  // Size of current data vector.
  unsigned short nData     = data.size();
  // Size of expected default reference data vector.
  unsigned short nExpected = reference.size();

  // Check if padding is needed.
  if( nData != nExpected ) for(unsigned short i=nData; i<nExpected; i++) data.push_back( reference[i] );

  // Consistency check.
  if( data.size() != nExpected ){
    std::string message = "Keyword: ";
    message += keyword;
    message += " must be of size: ";
    message += std::to_string(nExpected);
    Terminate("PadEntriesVectorDefaultData", __FILE__, __LINE__, message);
  }
}



// Function that determines if there are any floating point errors.
#ifdef ENABLE_NAN_CHECK
void CheckFloatingError(void);
#endif

// Include remaining classes.
#include "TensorProduct.hpp"
