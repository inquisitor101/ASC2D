#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "temporal_structure.hpp"
#include "iteration_structure.hpp"
#include "solver_structure.hpp"
#include "spatial_structure.hpp"
#include "output_structure.hpp"
#include "initial_structure.hpp"
#include "input_structure.hpp"
#include "element_structure.hpp"
#include "process_structure.hpp"



class CDriver {

	public:
		// Constructor.
		CDriver(CConfig *config);

		// Destructor.
		~CDriver(void);

		// Function that starts the solver.
		void StartSolver(void);


		// Simulation start time.
		as3double SimTimeStart;
		// Simulation final time.
		as3double SimTimeFinal;
		// Maximum temporal iterations.
		unsigned long MaxTimeIter;

	protected:
		// Containers needed to drive the solver.
		CConfig 			*config_container;
		CGeometry     *geometry_container;
		CInput        *input_container;
		COutput       *output_container;
		CTemporal     *temporal_container;
		CIteration   **iteration_container;
		CElement     **element_container;
		CSolver      **solver_container;
		CInitial     **initial_container;
		CSpatial     **spatial_container;
    CProcess     **process_container;

	private:
		// Number of zones.
		unsigned short nZone;

    // Total number of elements in all zones combined.
    unsigned long nElemTotal;

    // Map data from [iZone][iElemZone] to [iElem]. This is use for
    // parallelization efficiency. Dimension: [iElem][iData], where:
    // [iData] is: [0]: iZone, [1]: iElemZone.
    as3vector2d<unsigned long> MapGlobalToLocal;

		// Function that initializes all containers to nullptr.
		void SetContainers_Null(void);

		// Function that preprocesses the element container.
		void Element_Preprocessing(CConfig *config_container);

		// Function that preprocesses the geometry container.
		void Geometry_Preprocessing(CConfig   *config_container,
																CElement **element_container);

		// Function that preprocesses the input container.
		void Input_Preprocessing(CConfig   *config_container,
														 CGeometry *geometry_container);

		// Function that preprocesses the output container.
		void Output_Preprocessing(CConfig   *config_container,
															CGeometry *geometry_container);

		// Function that preprocesses the initial container.
		void Initial_Preprocessing(CConfig   	  *config_container,
															 CGeometry 	  *geometry_container,
															 CElement    **element_container);

		// Function that preprocesses the spatial container.
		void Spatial_Preprocessing(CConfig      *config_container,
															 CGeometry    *geometry_container,
															 CElement    **element_container,
                               CInitial    **initial_container);

		// Function that preprocesses the solver container.
		void Solver_Preprocessing(CConfig 		 *config_container,
															CGeometry 	 *geometry_container,
                              CInitial    **initial_container,
															CElement  	**element_container,
															CSpatial    **spatial_container);

		// Function that preprocesses the iteration container.
		void Iteration_Preprocessing(CConfig  	  *config_container,
																 CGeometry	  *geometry_container,
																 CSolver  	 **solver_container,
																 CElement 	 **element_container,
																 CSpatial    **spatial_container);

		// Function that preprocesses the temporal container.
		void Temporal_Preprocessing(CConfig      *config_container,
																CGeometry    *geometry_container,
																CIteration  **iteration_container,
																CSolver     **solver_container,
																CElement    **element_container,
																CSpatial    **spatial_container);

    // Function that preprocesses the process container.
    void Process_Preprocessing(CConfig       *config_container,
															 CGeometry     *geometry_container,
                               COutput       *output_container,
															 CElement     **element_container,
                               CInitial     **initial_container,
                               CSolver      **solver_container,
															 CSpatial     **spatial_container);

    // Function that preprocesses the parallelization set-up.
    void Parallelization_Preprocessing(void);

		// Function that preprocesses the driver once at the start.
		void Preprocess(void);

		// Function that runs the solver.
		void Run(void);

		// Function that outputs the data being monitored and header.
		void MonitorOutput(unsigned long           iIter,
											 as3double               time,
											 as3double               dt,
                       as3vector1d<as3double> &MonitoringData,
											 bool                    MonitorData);

		// Function that estimates the time step.
		as3double ComputeTimeStep(void);
};




