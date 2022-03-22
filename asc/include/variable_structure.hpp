#pragma once

#include "option_structure.hpp"


/* * * * * * * *
 * Scalar data.
 * * * * * * * * */

struct CScalarVariable {
  // Constructor.
  CScalarVariable(void){}
  // Destructor.
  virtual ~CScalarVariable(void){}

	// Pure virtual function that returns a value. Must be
	// overriden by a derived class.
	virtual as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) = 0;
};


struct CDensityConservative : public CScalarVariable {
  // Constructor.
  CDensityConservative(void){}
  // Destructor.
  ~CDensityConservative(void){}

	// Function that returns the density, based on conservative working variables.
	inline as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
		// Return density.
		return WorkingVariable[CONT_VAR][index];
	}
};


struct CEnergyConservative : public CScalarVariable {
  // Constructor.
  CEnergyConservative(void){}
  // Destructor.
  ~CEnergyConservative(void) {}

	// Function that returns the density, based on conservative working variables.
	inline as3double GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
		// Return energy.
		return WorkingVariable[ENER_VAR][index];
	}
};


struct CPressureConservative : public CScalarVariable {
  // Constructor.
  CPressureConservative(void){}
  // Destructor.
  ~CPressureConservative(void) {}

	// Function that returns the density, based on conservative working variables.
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


struct CTemperatureConservative : public CScalarVariable {
  // Constructor.
  CTemperatureConservative(void){}
  // Destructor.
  ~CTemperatureConservative(void) {}

  // Function that returns the temperature, based on conservative working variables.
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


struct CMachNumberConservative : public CScalarVariable {
  // Constructor.
  CMachNumberConservative(void){}
  // Destructor.
  ~CMachNumberConservative(void) {}

  // Function that returns the temperature, based on conservative working variables.
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


/* * * * * * * *
 * Vector data.
 * * * * * * * * */


struct CVectorVariable {
  // Constructor.
  CVectorVariable(void){}
  // Destructor.
  virtual ~CVectorVariable(void){}

	// Pure virtual function that returns a vector value. Must be
	// overriden by a derived class.
	virtual as3vector1d<as3double> GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) = 0;
};


struct CMomentumConservative : public CVectorVariable {
  // Constructor.
  CMomentumConservative(void){}
  // Destructor.
  ~CMomentumConservative(void) {}

	// Function that returns the momentum, based on conservative working variables.
	inline as3vector1d<as3double> GetValue(const as3data1d<as3double> &WorkingVariable, unsigned long index) {
		// Assemble momentum.
		as3vector1d<as3double> momentum = { WorkingVariable[XMOM_VAR][index],
                                        WorkingVariable[YMOM_VAR][index] };
    // Return momentum.
		return momentum;
	}
};


