#pragma once

#include "option_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"

// Forward declaration to avoid compiler errors.
#include "variable_structure.hpp"


class COutput {

	public:

		// Constructor.
		COutput(CConfig   *config_container,
						CGeometry *geometry_container);

		// Destructor.
		~COutput(void);


		// Function that writes VTK data file.
		void WriteFileVTK(CConfig    *config_container,
											CGeometry  *geometry_container,
											CSolver   **solver_container);

    // Function that writes a data to a file in ASCII format.
    void WriteDataToFile(CConfig                *config_container,
                         CGeometry              *geometry_container,
                         const char             *fileinfo,
                         as3vector1d<as3double> &time,
                         as3vector2d<as3double> &data);

    // Function that writes the data to a file.
    void WriteSolutionToFile(CConfig    *config_container,
                             CGeometry  *geometry_container,
                             CElement  **element_container,
                             CSolver   **solver_container,
                             as3double   SimTime);

    // Function that writes zone data to a file.
    void WriteZoneDataToFile(CConfig    *config_container,
                             CGeometry  *geometry_container,
                             CSolver   **solver_container,
                             as3double   SimTime);

	protected:

	private:
		// Number of zones.
		unsigned short nZone;

		// Output VTK Filename.
		const char *OutputVTKFilename;

		// VTK file number used in writing sequence.
		unsigned long FileNumberVTK;
    // Processed data file number.
    unsigned long FileNumberProcessed;
    // Data solution file number.
    unsigned long FileNumberDataSolution;
    // Zone data solution file number.
    unsigned long FileNumberZoneData;

		// Local indicial connectivity matrix.
		unsigned short **ConnLocal;

		// Type of density scalar variable used.
		CScalarVariable *VariableDensity     = nullptr;
		// Type of energy scalar variable used.
		CScalarVariable *VariableEnergy      = nullptr;
		// Type of pressure scalar variable used.
		CScalarVariable *VariablePressure    = nullptr;
    // Type of temperature scalar variable used.
    CScalarVariable *VariableTemperature = nullptr;
    // Type of Mach number scalar variable used.
    CScalarVariable *VariableMachNumber  = nullptr;
		// Type of momentum vector variable used.
		CVectorVariable *VariableMomentum    = nullptr;


		// Function that writes a scalar parameter.
		void WriteScalar(CConfig         *config_container,
										 CGeometry       *geometry_container,
										 CSolver     	  **solver_container,
										 std::ofstream   &Paraview_File,
										 CScalarVariable *Variable);

		// Function that writes a vector parameter.
		void WriteVector(CConfig         *config_container,
										 CGeometry       *geometry_container,
										 CSolver     	  **solver_container,
										 std::ofstream   &Paraview_File,
										 CVectorVariable *Variable);
    
    // Function that writes a PML auxiliary scalar parameter.
    void WriteScalarAuxPML(CConfig         *config_container,
                           CGeometry       *geometry_container,
                           CSolver     	  **solver_container,
                           std::ofstream   &Paraview_File,
                           unsigned short   iVar);

    // Function that writes a PML auxiliary vector parameter.
    void WriteVectorAuxPML(CConfig         *config_container,
                           CGeometry       *geometry_container,
                           CSolver     	  **solver_container,
                           std::ofstream   &Paraview_File,
                           unsigned short   iVar);
};


