#pragma once

/*!
 * @file option_structure.hpp
 * @brief The file containing all the global option functionalities.
 */

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

typedef double as3double;                                       ///< Type definition of double.
using as3element = std::vector<CData*>;                         ///< Element data structure.

template<typename T>
using as3vector1d = std::vector<T>;                             ///< Vector of data type: 1D.

template<typename T>
using as3vector2d = std::vector<std::vector<T>>;                ///< Vector of data type: 2D.

template<typename T>
using as3vector3d = std::vector<std::vector<std::vector<T>>>;   ///< Vector of data type: 3D.


template<typename T>
using as3data1d = std::vector<T*>;                              ///< Vector of pointer data type: 1D.

template<typename T>
using as3data2d = std::vector<std::vector<T*>>;                 ///< Vector of pointer data type: 2D.

template<typename T>
using as3data3d = std::vector<std::vector<std::vector<T*>>>;    ///< Vector of pointer data type: 3D.

const int AS3_MagicNumber       = 535532;                       ///< Magic number in the restart file when using the AS3 format.
const int CGNS_STRING_SIZE      = 33;                           ///< Length of the strings to store the variables names in the restart files.
const unsigned short nVar       = 4;                            ///< Number of working variables.
const unsigned short nDim       = 2;                            ///< Number of spatial dimensions.
const unsigned short nFace      = 4;                            ///< Number of faces per element (for now all elements are assumed to be quadrilateral).
const as3double PI_CONSTANT     = 4.0*atan(1.0);                ///< Value of the constant \pi.
const as3double EPS_VALUE       = 1e-15;                        ///< Epsilon value used as a threshold.
const as3double GAMMA           = 1.4;                          ///< Ratio of specific heats.
const as3double GAMMA_MINUS_ONE = GAMMA-1.0;                    ///< Abbreviation for: gamma - 1.
const as3double GAS_CONSTANT    = 287.058;                      ///< Gas constant.
const as3double CV_CONSTANT     = GAS_CONSTANT/GAMMA_MINUS_ONE; ///< Specific heat at constant volume.
const as3double CP_CONSTANT     = GAMMA*CV_CONSTANT;            ///< Specific heat at constant pressure.
const size_t vecLen1D           = 4;                            ///< For vectorization, the vector length used must be a multiple of vecLen1D.


/*!
 * @brief Enumerated type for face indices on an element.
 */
enum FACE_ELEM_ID {
  IDX_UNKNOWN = 99, ///< Unknown face index.
	IDX_SOUTH   = 0,  ///< South-most face.
	IDX_NORTH   = 1,  ///< North-most face.
	IDX_WEST    = 2,  ///< West-most face.
	IDX_EAST    = 3   ///< East-most face.
};

/*!
 * @brief Enumerated type for dimension index.
 */
enum DIM_ID {
	XDIM = 0, ///< Index of x-dimension.
	YDIM = 1  ///< Index of y-dimension.
};

/*!
 * @brief Enumerated type of nodal DOFs.
 */
enum NODAL_DOFS {
	TYPE_DOF_UNKOWN = 99, ///< Unknown type of nodal points.
	TYPE_DOF_EQD    = 0,  ///< Equidistant nodal points.
	TYPE_DOF_LGL    = 1   ///< Legendre-Gauss-Lobato nodal points.
};

/*!
 * @brief Enumerated type of problem-dependent variables.
 */
enum PROBLEM_EQUATIONS_VARIABLE {
  CONT_VAR   = 0,  ///< Continuity variable.
	XMOM_VAR   = 1,  ///< x-momentum variable.
	YMOM_VAR   = 2,  ///< y-momentum variable.
	ENER_VAR   = 3,  ///< energy variable.
  CONTQ1_VAR = 4,  ///< q1[density] variable.
  XMOMQ1_VAR = 5,  ///< q1[x-momentum] variable.
  YMOMQ1_VAR = 6,  ///< q1[x-momentum] variable.
  ENERQ1_VAR = 7,  ///< q1[energy] variable.
};

/*!
 * @brief Enumerated type of PDE problem to solve.
 */
enum TYPE_PDE_SOLVER {
	NO_SOLVER      = 99, ///< Solver option not valid.
	SOLVER_LEE     = 0,  ///< Linearized Euler equations.
	SOLVER_EE      = 1,  ///< Non-linear Euler equations.
	SOLVER_NS      = 2,  ///< Navier-Stokes.
};

/*!
 * @brief Enumerated type of riemann solver.
 */
enum RIEMANN_SOLVER_TYPE {
  RIEMANN_UNKNOWN   = 99, ///< Unknown Riemann solver.
  RIEMANN_ROE       = 0,  ///< Roe's Riemann solver.
  RIEMANN_RUSANOV   = 1,  ///< Rusanov's Riemann solver.
  RIEMANN_ROEISMAIL = 2   ///< Ismail-Roe Riemann solver.
};

/*!
 * @brief Enumerated type of initial conditions.
 */
enum INITIAL_CONDITION_TYPE {
	IC_UNKNOWN     			      = 99, ///< Unknown type of IC.
	IC_IGNORE      			      = 0,  ///< No IC.
	IC_FREE_STREAM 			      = 1,  ///< Free-stream IC.
	IC_GAUSSIAN_PRESSURE      = 2,  ///< Gaussian pressure IC.
  IC_ISENTROPIC_VORTEX      = 3,  ///< Isentropic vortex IC.
  IC_ENTROPY_WAVE           = 4,  ///< Entropy wave IC.
  IC_VORTEX_ROLLUP          = 5,  ///< Vortex roll-up IC.
  IC_ACOUSTIC_PLANE_WAVE    = 6,  ///< Acoustic plane-wave IC.
	IC_GAUSSIAN_PRESSURE_1D_X = 7,  ///< 1D (x-)pressure pulse IC.
  IC_GAUSSIAN_PRESSURE_1D_Y = 8   ///< 1D (y-)pressure pulse IC.
};

/*!
 * @brief Enumerated type of boundary condition.
 */
enum BOUNDARY_CONDITION_TYPE {
	BC_UNKNOWN           = 99, ///< Unknown type of BC.
	BC_INTERFACE         = 0,  ///< Interface/periodic boundary.
	BC_SYMMETRY          = 1,  ///< Inviscid/slip wall.
  BC_CBC_OUTLET        = 2,  ///< Outlet CBC.
  BC_CBC_INLET         = 3,  ///< Inlet CBC.
  BC_CBC_TOTAL_INLET   = 4,  ///< Total-condition CBC inlet.
	BC_CBC_STATIC_INLET  = 5,  ///< Static-condition CBC inlet.
	BC_STATIC_INLET      = 6,  ///< Static-condition inlet.
  BC_TOTAL_INLET       = 7,  ///< Total-conditions inlet.
  BC_STATIC_OUTLET     = 8,  ///< Static-condition subsonic outlet.
  BC_SUPERSONIC_OUTLET = 9,  ///< Supersonic outlet.
  BC_SUPERSONIC_INLET  = 10  ///< Supersonic inlet.
};

/*!
 * @brief Enumerated type of temporal scheme used.
 */
enum TEMPORAL_SCHEME {
	TEMPORAL_SCHEME_UNKNOWN = 99, ///< Temporal scheme not valid.
	TEMPORAL_SCHEME_LSRK4   = 0,  ///< Low-storage 4th-order RK.
	TEMPORAL_SCHEME_CRK4    = 1,  ///< Classic 4th-order RK.
  TEMPORAL_SCHEME_SSPRK3  = 2   ///< Strong-stability-preserving RK3.
};

/*!
 * @brief Enumerated type for LSRK4 options.
 */
enum LSRK4_OPTIONS {
  LSRK4_N_STORAGE = 1, ///< Number of tentative storage.
  LSRK4_N_STAGES  = 5  ///< Number of grid sweeps.
};

/*!
 * @brief Enumerated type for SSPRK3 options.
 */
enum SSPRK3_OPTIONS {
  SSPRK3_N_STORAGE = 1, ///< Number of tentative storage.
  SSPRK3_N_STAGES  = 3  ///< Number of grid sweeps.
};

/*!
 * @brief Enumerated type of multizone strategies.
 */
enum MULTIZONE_STRATEGY_TYPE {
  MULTIZONE_STRATEGY_MAIN      = 0, ///< Only zone 0 is used.
  MULTIZONE_STRATEGY_WEST      = 1, ///< Use only zone: 0 and 1.
  MULTIZONE_STRATEGY_EAST      = 2, ///< Use only zone: 0 and 2.
  MULTIZONE_STRATEGY_SOUTH     = 3, ///< Use only zone: 0 and 3.
  MULTIZONE_STRATEGY_NORTH     = 4, ///< Use only zone: 0 and 4.
  MULTIZONE_STRATEGY_ALL       = 5, ///< Uses all zones: 0 to 8.
  MULTIZONE_STRATEGY_HORIZONAL = 6, ///< Uses zones: 0, 1 and 2.
  MULTIZONE_STRATEGY_VERTICAL  = 7  ///< Uses zones: 0, 3 and 4.
};

/*!
 * @brief Enumerated type of zones.
 */
enum ZONE_TYPE {
  ZONE_UNKNOWN  = 99, ///< Unknown zone.
	ZONE_MAIN     = 0,  ///< Zone corresponding to physical region.
	ZONE_WEST     = 1,  ///< Zone PML corresponding to west.
	ZONE_EAST     = 2,  ///< Zone PML corresponding to east.
	ZONE_SOUTH    = 3,  ///< Zone PML corresponding to south.
	ZONE_NORTH    = 4,  ///< Zone PML corresponding to north.
	ZONE_CORNER_0 = 5,  ///< Zone PML corresponding to SW-corner C0.
	ZONE_CORNER_1 = 6,  ///< Zone PML corresponding to SE-corner C1.
	ZONE_CORNER_2 = 7,  ///< Zone PML corresponding to NW-corner C2.
	ZONE_CORNER_3 = 8   ///< Zone PML corresponding to NE-corner C3.
};

/*!
 * @brief Enumerated type of 3-zone horizontal strategy.
 */
enum HORIZONTAL_ZONE_TYPE {
  HORIZONTAL_ZONE_MAIN = 0, ///< Zone corresponding to physical region.
  HORIZONTAL_ZONE_WEST = 1, ///< Zone corresponding to west.
  HORIZONTAL_ZONE_EAST = 2  ///< Zone corresponding to east.
};

/*!
 * @brief Enumerated type of 3-zone vertical strategy.
 */
enum VERTICAL_ZONE_TYPE {
  VERTICAL_ZONE_MAIN  = 0, ///< Zone corresponding to physical region.
  VERTICAL_ZONE_SOUTH = 1, ///< Zone corresponding to south.
  VERTICAL_ZONE_NORTH = 2  ///< Zone corresponding to north.
};

/*!
 * @brief Enumerated type of processing data.
 */
enum PROCSES_DATA {
  PROCESS_NOTHING      = 0, ///< No processing is done.
  PROCESS_RATIO_P_U    = 1, ///< Process: ratio( p / u ).
	PROCESS_RATIO_P_M    = 2, ///< Process: ratio( p / M ).
	PROCESS_RATIO_WM_WP  = 3, ///< Process: ratio( w(-) / w(+) ). 
	PROCESS_RATIO_LM_LP  = 4, ///< Process: ratio( L(-) / L(+) ).
	PROCESS_WAVE_ENTROPY = 5  ///< Process: entropy wave.
};

/*!
 * @brief Enumerated type of processing location.
 */
enum PROCESS_DATA_LOCATION {
	PROCESS_LOCATION_UNKNOWN = 0, ///< Unknown location.
	PROCESS_LOCATION_XMIN    = 1, ///< X-min boundary.
	PROCESS_LOCATION_XMAX    = 2, ///< Ximax boundary.
	PROCESS_LOCATION_YMIN    = 3, ///< Yimin boundary.
	PROCESS_LOCATION_YMAX    = 4, ///< Yimax boundary.
	PROCESS_LOCATION_DOMAIN  = 5  ///< Total (main-zone) domain.
};

/*!
 * @brief Enumerated type of variables to probe.
 */
enum PROBE_VARIABLES_DATA {
  PROBE_NOTHING     = 0, ///< No variables specified.
  PROBE_DENSITY     = 1, ///< Density.
	PROBE_XMOMENTUM   = 2, ///< X-momentum.
	PROBE_YMOMENTUM   = 3, ///< U-momentum.
	PROBE_TOTALENERGY = 4, ///< Total energy.
	PROBE_XVELOCITY   = 5, ///< U-velocity.
	PROBE_YVELOCITY   = 6, ///< V-velocity.
	PROBE_PRESSURE    = 7  ///< Pressure.
};

/*!
 * @brief Enumerated type of VTK variables written.
 */
enum VTK_VARIABLE_WRITTEN {
	UNKNOWN_VTK_VARIABLE     = 0, ///< Unknown variable.
	VTK_VARIABLE_DENSITY     = 1, ///< Scalar: density.
	VTK_VARIABLE_MOMENTUM    = 2, ///< Vector: momentum.
	VTK_VARIABLE_TOTALENERGY = 3, ///< Scalar: total energy.
	VTK_VARIABLE_PRESSURE    = 4, ///< Scalar: pressure.
	VTK_VARIABLE_VELOCITY    = 5, ///< Vector: velocity.
	VTK_VARIABLE_VORTICITY   = 6, ///< Scalar: vorticity.
	VTK_VARIABLE_MACH        = 7, ///< Scalar: Mach number.
	VTK_VARIABLE_TEMPERATURE = 8, ///< Scalar: temperature.
	VTK_VARIABLE_ENTROPY     = 9  ///< Scalar: specific entropy.
};

/*!
 * @brief Enumerated type of number of nodes per nPoly=1 element.
 */
enum ELEM_POINTS {
	N_POINTS_LINE          = 2, ///< Line.
	N_POINTS_TRIANGLE      = 3, ///< Triangle.
	N_POINTS_QUADRILATERAL = 4, ///< Quadrilateral.
	N_POINTS_TETRAHEDRON   = 4, ///< Tetrahedron.
	N_POINTS_HEXAHEDRON    = 8, ///< Hexaheron.
	N_POINTS_PYRAMID       = 5, ///< Pyramid.
	N_POINTS_PRISM         = 6  ///< Prism.
};

/*!
 * @brief Enumerated type of geometric entities based on VTK nomenclature.
 */
enum GEO_TYPE {
  NONE          = 0,  ///< Point.
  VERTEX        = 1,  ///< VTK nomenclature for defining a vertex element.
  LINE          = 3,  ///< VTK nomenclature for defining a line element.
  TRIANGLE      = 5,  ///< VTK nomenclature for defining a triangle element.
  QUADRILATERAL = 9,  ///< VTK nomenclature for defining a quadrilateral element.
  TETRAHEDRON   = 10, ///< VTK nomenclature for defining a tetrahedron element.
  HEXAHEDRON    = 12, ///< VTK nomenclature for defining a hexahedron element.
  PRISM         = 13, ///< VTK nomenclature for defining a prism element.
  PYRAMID       = 14  ///< VTK nomenclature for defining a pyramid element.
};

/*!
 * @brief Enumerated type of data filtering.
 */
enum FILTERING_TYPE {
  NO_FILTER          = 0, ///< No filtering used.
  EXPONENTIAL_FILTER = 1  ///< Exponential filter.
};

/*!
 * @brief Enumerated type of buffer layer.
 */
enum BUFFER_LAYER_TYPE {
  NO_LAYER                = 0, ///< No layer.
  PML_LAYER               = 1, ///< PML layer.
  SPONGE_LAYER            = 2  ///< sponge layer.
};

/*!
 * @brief Enumerated type of modified boundary condition.
 */
enum MODIFIED_BOUNDARY_CONDITION {
  NO_BC_MODIFICATION                  = 0, ///< No modification.
  GAUSSIAN_PRESSURE_BC_MODIFICATION   = 1, ///< Gaussian-pressure profile.
  SINUSOIDAL_VELOCITY_BC_MODIFICATION = 2  ///< Sinusoidal velocity profile.
};

/*!
 * @brief Enumerated type of file format.
 */
enum FORMAT_FILE {
	UNKNOWN_FORMAT = 0, ///< Unknown file format. 
	ASCII_FORMAT   = 1, ///< ASCII  format.
	BINARY_FORMAT  = 2  ///< Binary format.
};


/*!
 * @brief Function used for printing the type of the input zone.
 *
 * @param[in] iTypeZone type of input zone.
 *
 * @return string name for type of zone
 */
std::string DisplayTypeZone(unsigned short iTypeZone);
/*!
 * @brief Function used for printing the side of the input boundary.
 *
 * @param[in] iBoundary input boundary ID.
 *
 * @return string name for boundary ID
 */
std::string DisplayBoundarySide(unsigned short iBoundary);
/*!
 * @brief Function used for printing the type of solver.
 *
 * @param[in] iSolver input solver type.
 *
 * @return string name for type of solver
 */
std::string DisplaySolverType(unsigned short iSolver);
/*!
 * @brief Function used for printing the type of BC.
 *
 * @param[in] iBoundary input boundary ID.
 *
 * @return string name for type of boundary
 */
std::string DisplayBoundaryType(unsigned short iBoundary);


/*!
 * @brief Function that removes blanks from a string.
 *
 * @param[in] lineBuf reference to input string line.
 */
void RemoveBlanks(std::string &lineBuf);

/*!
 * @brief Function that replaces tabs and return characters in a string with a blank.
 *
 * @param[in] lineBuf reference to input string line.
 */
void ReplaceTabsAndReturns(std::string &lineBuf);

/*!
 * @brief Function that returns a lowered-case version of input string.
 *
 * @param[in] lineBuf reference to input string line.
 */
void CreateLowerCase(std::string &lineBuf);

/*!
 * @brief Function that searches for and finds a string vector from input file.
 *
 * @param[in] inputFile reference to input file.
 * @param[in] keyword pointer to keyword looked for.
 * @param[out] valStringVector reference to extracted string vector. 
 *
 * @return indicator on whether the extraction worked or not
 */
bool FindStringVectorFromKeyword(std::ifstream            &inputFile,
                          			 const char               *keyword,
                          			 std::vector<std::string> &valStringVector);

/*!
 * @brief Function that searches for and finds a string from input file.
 *
 * @param[in] inputFile reference to input file.
 * @param[in] keyword pointer to keyword looked for.
 * @param[out] valString reference to extracted string. 
 *
 * @return indicator on whether the extraction worked or not
 */
bool FindStringFromKeyword(std::ifstream &inputFile,
                          const char     *keyword,
                          std::string    &valString);

/*!
 * @brief Function that terminates program immediately.
 *
 * @param[in] functionName pointer to the function name.
 * @param[in] fileName pointer to the file name.
 * @param[in] lineNumber input line number.
 * @param[in] errorMessageI reference to the output error message displayed.
 */
void Terminate(const char        *functionName,
               const char        *fileName,
               const int          lineNumber,
               const std::string &errorMessageI);

/*!
 * @brief Function that swaps bytes.
 *
 * @param[in,out] buffer pointer to the data buffer being swapped.
 * @param[in] nBytes input number of bytes.
 * @param[in] nItems input number of items in the buffer.
 */
void SwapBytes(void   *buffer,
		           size_t  nBytes,
							 int     nItems);

/*!
 * @brief Function used for printing, usually debugging only. Works on 1D arrays.
 *
 * @param[in] val pointer to data (matrix/array) written.
 * @param[in] nRow number of rows of input data matrix.
 * @param[in] nCol number of columns of input data matrix.
 * @param[in] message pointer to the header message displayed.
 * @param[in] nDigits input precision of values displayed.
 */
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

/*!
 * @brief Function used for printing, usually debugging only. Works on 2D arrays.
 *
 * @param[in] val pointer to data (matrix/array) written.
 * @param[in] nRow number of rows of input data matrix.
 * @param[in] nCol number of columns of input data matrix.
 * @param[in] message pointer to the header message displayed.
 * @param[in] nDigits input precision of values displayed.
 */
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

/*!
 * @brief Function used for printing, usually debugging only. Works on 2D vectors.
 *
 * @param[in] val reference to vector of data (matrix/array) written.
 * @param[in] message pointer to the header message displayed.
 * @param[in] nDigits input precision of values displayed.
 */
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


/*!
 * @brief Function that determines a boolean input.
 *
 * @param[in] inputFile reference to input file.
 * @param[in] keyword pointer to keyword looked for.
 * @param[out] value reference to value of boolean read.
 * @param[in] defaultValue reference to default value specified.
 * @param[in] ResetLine option to reset to first line of file.
 */
void AddBoolOption(std::ifstream &inputFile,
                   const char    *keyword,
                   bool          &value,
                   std::string   &defaultValue,
                   bool           ResetLine = false);


/*!
 * @brief Function used in reading "scalar" data, takes default parameter.
 *
 * @param[in] inputFile reference to input file.
 * @param[in] keyword pointer to keyword looked for.
 * @param[out] value reference to value read.
 * @param[in] defaultValue reference to default value specified.
 * @param[in] ResetLine option to reset to first line of file.
 */
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


/*!
 * @brief Function used in reading "scalar" data.
 *
 * @param[in] inputFile reference to input file.
 * @param[in] keyword pointer to keyword looked for.
 * @param[out] value reference to value read.
 * @param[in] ResetLine option to reset to first line of file.
 */
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


/*!
 * @brief Function that works as a template to read vector data, otherwise assigns default value.
 *
 * @param[in] inputFile reference to input file.
 * @param[in] keyword pointer to keyword looked for.
 * @param[out] value reference to value read.
 * @param[in] defaultValue reference to default value specified.
 * @param[in] ResetLine option to reset to first line of file.
 */
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


/*!
 * @brief Function that works as a template to read vector data, otherwise exits.
 *
 * @param[in] inputFile reference to input file.
 * @param[in] keyword pointer to keyword looked for.
 * @param[out] value reference to value read.
 * @param[in] ResetLine option to reset to first line of file.
 */
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


/*!
 * @brief Function that pads remaining entries in a vector.
 *
 * @param[in,out] data reference to data that is padded.
 * @param[in] keyword input option string name.
 * @param[in] nExpected expected size of input/output vector.
 * @param[in] nCondition1 input condition on whether to pad. 
 * @param[in] nCondition2 input condition on whether to pad. 
 * @param[in] nCondition3 input condition on whether to pad. 
 */
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


/*!
 * @brief Function that pads remaining entries in a vector based on default reference data.
 *
 * @param[in,out] data reference to data that is padded.
 * @param[in] reference reference to default value in padding.
 * @param[in] keyword input option string name.
 */
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
