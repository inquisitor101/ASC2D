#pragma once

/*!
 * @file variable_structure.hpp
 * @brief The file containing all the variable extraction expressions.
 */

#include "option_structure.hpp"


/* * * * * * * *
 * Scalar data.
 * * * * * * * * */


/*!
 * @brief An interface structure used for extracting a generic scalar variable.
 */
struct CScalarVariable {
 	/*!
	 * @brief Default constructor of CScalarVariable, which does nothing.
	 */
	CScalarVariable(void){}
  
	/*!
	 * @brief Destructor, which does nothing.
	 */
  virtual ~CScalarVariable(void){}

	/*!
	 * @brief Pure virtual function that returns a scalar value. 
	 * Must be overriden by a derived class.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the relevant scalar value
	 */
	virtual as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) = 0;
};


/*!
 * @brief A structure used for extracting a density scalar variable from conservative variables.
 */
struct CDensityConservative : public CScalarVariable {
 	/*!
	 * @brief Default constructor of CDensityConservative, which does nothing.
	 */
	CDensityConservative(void){}
  
  /*!
	 * @brief Destructor, which does nothing.
	 */
	~CDensityConservative(void){}

	/*!
	 * @brief Function that returns the density, based on conservative working variables.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the density scalar value
	 */
	inline as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
		// Return density.
		return WorkingVariable[CONT_VAR][index];
	}
};


/*!
 * @brief A structure used for extracting a total energy scalar variable from conservative variables.
 */
struct CEnergyConservative : public CScalarVariable {
 	/*!
	 * @brief Default constructor of CEnergyConservative, which does nothing.
	 */
	CEnergyConservative(void){}
  
  /*!
	 * @brief Destructor, which does nothing.
	 */
	~CEnergyConservative(void){}

	/*!
	 * @brief Function that returns the total energy, based on conservative working variables.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the total energy scalar value
	 */
	inline as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
		// Return energy.
		return WorkingVariable[ENER_VAR][index];
	}
};


/*!
 * @brief A structure used for extracting a pressure scalar variable from conservative variables.
 */
struct CPressureConservative : public CScalarVariable {
 	/*!
	 * @brief Default constructor of CPressureConservative, which does nothing.
	 */
	CPressureConservative(void){}
  
  /*!
	 * @brief Destructor, which does nothing.
	 */
	~CPressureConservative(void){}

	/*!
	 * @brief Function that returns the pressure, based on conservative working variables.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the pressure scalar value
	 */
	inline as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
		// Extract conservative variables.
		const as3double rho  = WorkingVariable[CONT_VAR][index];
		const as3double rhou = WorkingVariable[XMOM_VAR][index];
		const as3double rhov = WorkingVariable[YMOM_VAR][index];
		const as3double rhoE = WorkingVariable[ENER_VAR][index];
		// Compute pressure.
		const as3double p    = (GAMMA_MINUS_ONE)*( rhoE - 0.5*( rhou*rhou + rhov*rhov )/rho );
		// Return pressure.
		return p;
	}
};


/*!
 * @brief A structure used for extracting a temperature scalar variable from conservative variables.
 */
struct CTemperatureConservative : public CScalarVariable {
 	/*!
	 * @brief Default constructor of CTemperatureConservative, which does nothing.
	 */
	CTemperatureConservative(void){}
  
  /*!
	 * @brief Destructor, which does nothing.
	 */
	~CTemperatureConservative(void){}

  /*!
	 * @brief Function that returns the temperature, based on conservative working variables.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the temperature scalar value
	 */
  inline as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
    // Extract conservative variables.
    const as3double rho  = WorkingVariable[CONT_VAR][index];
    const as3double rhou = WorkingVariable[XMOM_VAR][index];
    const as3double rhov = WorkingVariable[YMOM_VAR][index];
    const as3double rhoE = WorkingVariable[ENER_VAR][index];
    // Compute pressure.
    const as3double p    = (GAMMA_MINUS_ONE)*( rhoE - 0.5*( rhou*rhou + rhov*rhov )/rho );
    // Compute temperature.
    const as3double T    = p/(rho*GAS_CONSTANT);
    // Return temperature.
    return T;
  }
};


/*!
 * @brief A structure used for extracting a Mach number scalar variable from conservative variables.
 */
struct CMachNumberConservative : public CScalarVariable {
 	/*!
	 * @brief Default constructor of CMachNumberConservative, which does nothing.
	 */
	CMachNumberConservative(void){}
  
  /*!
	 * @brief Destructor, which does nothing.
	 */
	~CMachNumberConservative(void){}

  /*!
	 * @brief Function that returns the Mach number, based on conservative working variables.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the Mach number scalar value
	 */
  inline as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
    // Extract conservative variables.
    const as3double rho  = WorkingVariable[CONT_VAR][index];
    const as3double rhou = WorkingVariable[XMOM_VAR][index];
    const as3double rhov = WorkingVariable[YMOM_VAR][index];
    const as3double rhoE = WorkingVariable[ENER_VAR][index];
    // Abbreviation.
    const as3double ovrho   = 1.0/rho;
    const as3double rhou2v2 = (rhou*rhou + rhov*rhov)*ovrho;
    // Compute velocity magnitude.
    const as3double umag    = sqrt( rhou2v2*ovrho );
    // Compute pressure.
    const as3double p       = GAMMA_MINUS_ONE*( rhoE - 0.5*rhou2v2 );
    // Compute speed of sound.
    const as3double a       = sqrt( GAMMA*p*ovrho );
    // Compute magnitude of Mach number.
    const as3double Mach    = umag/a;
    // Return magnitude of Mach number.
    return Mach;
  }
};


/*!
 * @brief A structure used for extracting a specific entropy scalar variable from conservative variables.
 */
struct CEntropyConservative : public CScalarVariable {
 	/*!
	 * @brief Default constructor of CEntropyConservative, which does nothing.
	 */
	CEntropyConservative(void){}
  
  /*!
	 * @brief Destructor, which does nothing.
	 */
	~CEntropyConservative(void){}

  /*!
	 * @brief Function that returns the specific entropy, based on conservative working variables.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the specific entropy scalar value
	 */
  inline as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
    // Extract conservative variables.
    const as3double rho  = WorkingVariable[CONT_VAR][index];
    const as3double rhou = WorkingVariable[XMOM_VAR][index];
    const as3double rhov = WorkingVariable[YMOM_VAR][index];
    const as3double rhoE = WorkingVariable[ENER_VAR][index];
    // Abbreviation.
    const as3double ovrho   = 1.0/rho;
    // Compute pressure.
    const as3double p       = GAMMA_MINUS_ONE*( rhoE - 0.5*(rhou*rhou + rhov*rhov)*ovrho );
		// Compute the specitific entropy.
		const as3double s       = logf( p/pow( rho, GAMMA ) );
    // Return magnitude of Mach number.
    return s;
  }
};



/* * * * * * * *
 * Vector data.
 * * * * * * * * */


/*!
 * @brief An interface structure used for extracting a generic vector variable.
 */
struct CVectorVariable {
 	/*!
	 * @brief Default constructor of CVectorVariable, which does nothing.
	 */
	CVectorVariable(void){}
  
  /*!
	 * @brief Destructor, which does nothing.
	 */
	virtual ~CVectorVariable(void){}

	/*!
	 * @brief Pure virtual function that returns a vector value. 
	 * Must be overriden by a derived class.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the relevant vector value
	 */
	virtual as3vector1d<as3double> GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) = 0;
};


/*!
 * @brief A structure used for extracting a momentum vector variable from conservative variables.
 */
struct CMomentumConservative : public CVectorVariable {
 	/*!
	 * @brief Default constructor of CMomentumConservative, which does nothing.
	 */
	CMomentumConservative(void){}
  
  /*!
	 * @brief Destructor, which does nothing.
	 */
	~CMomentumConservative(void){}

	/*!
	 * @brief Function that returns the momentum, based on conservative working variables.
	 *
	 * @param[in] WorkingVariable reference to the input conservatie variables.
	 * @param[in] index input nodal index of the variable.
	 *
	 * @return the value of the momentum vector value
	 */
	inline as3vector1d<as3double> GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
		// Assemble momentum.
		as3vector1d<as3double> momentum = { WorkingVariable[XMOM_VAR][index],
                                        WorkingVariable[YMOM_VAR][index] };
    // Return momentum.
		return momentum;
	}
};


