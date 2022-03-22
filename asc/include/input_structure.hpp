#pragma once

#include "option_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "element_structure.hpp"



class CInput {

	public:
		// Constructor.
		CInput(CConfig   *config_container,
					 CGeometry *geometry_container);

		// Destructor.
		~CInput(void);

    // Function that reads the solution from a restart file.
    void ReadSolutionRestartFile(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CElement  **element_container,
                                 CSolver   **solver_container,
                                 as3double  &SimTime);

	protected:

	private:
		// Number of zones.
		unsigned short nZone;

    // Interpolating polynomial that converts from the restart solution to the
    // current specified/expected solution polynomial.
    as3vector1d<as3double> lagrange1DTranspose;

};



