#pragma once

/*!
 * @file driver_structure.hpp
 * @brief The file responsible for setting and and running the solver.
 */

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


/*!
 * @brief A class used for initializing all necessary classes and running the solver.
 */
class CDriver {

	public:
		/*!
		 * @brief Default constructor of CDriver, which initializes the solver.
		 *
		 * @param[in] config input configuration/dictionary container.
		 */
		CDriver(CConfig *config);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CDriver(void);

		/*!
		 * @brief Function that starts the solver.
		 */
		void StartSolver(void);


		as3double     SimTimeStart;         ///< Simulation (physical) starting time.
		as3double     SimTimeFinal;         ///< Simulation (physical) ending time.
		unsigned long MaxTimeIter;          ///< Maximum number of temporal iterations.

	protected:
		CConfig 			*config_container;    ///< Container for the dictionary information. 
		CGeometry     *geometry_container;  ///< Container for the geometry data.
		CInput        *input_container;     ///< Container for the input imported data.
		COutput       *output_container;    ///< Container for the output data written.
		CTemporal     *temporal_container;  ///< Container for the temporal discretization.
		CIteration   **iteration_container; ///< Container for the iteration structure at each time step, per each zone.
		CElement     **element_container;   ///< Container for the standard/reference element, per each zone.
		CSolver      **solver_container;    ///< Container for the solver step, per each zone. 
		CInitial     **initial_container;   ///< Container for the initial condition specified, per each zone.
		CSpatial     **spatial_container;   ///< Container for the spatial discretization, per each zone.
    CProcess     **process_container;   ///< Container for the (post-)processing commands, per each zone.

	private:
		unsigned short             nZone;            ///< Total number of zones.
    unsigned long              nElemTotal;       ///< Total number of elements in all the zones combined.
    as3vector2d<unsigned long> MapGlobalToLocal; ///< Mapping vector used for parallelization efficiency.
																								 ///< This maps data from the format: [iZone][iElemZone] to [iElem].
																								 ///< Dimension: [iElem][iData], where [iData] is [0]: iZone, [1]: iElemZone.

		/*!
		 * @brief Function that initializes all containers to nullptr.
		 */
		void SetContainers_Null(void);

		/*!
		 * @brief Function that preprocesses the element container.
		 *
		 * @param[in] config_container configuration container.
		 */
		void Element_Preprocessing(CConfig *config_container);

		/*!
		 * @brief Function that preprocesses the geometry container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] element_container pointer to the element container.
		 */
		void Geometry_Preprocessing(CConfig   *config_container,
																CElement **element_container);

		/*!
		 * @brief Function that preprocesses the input container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 */
		void Input_Preprocessing(CConfig   *config_container,
														 CGeometry *geometry_container);

		/*!
		 * @brief Function that preprocesses the output container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 */
		void Output_Preprocessing(CConfig   *config_container,
															CGeometry *geometry_container);

		/*!
		 * @brief Function that preprocesses the initial container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] element_container pointer to the element container.
		 */
		void Initial_Preprocessing(CConfig   	  *config_container,
															 CGeometry 	  *geometry_container,
															 CElement    **element_container);

		/*!
		 * @brief Function that preprocesses the spatial container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] element_container pointer to the element container.
		 * @param[in] initial_container pointer to the initial container.
		 */
		void Spatial_Preprocessing(CConfig      *config_container,
															 CGeometry    *geometry_container,
															 CElement    **element_container,
                               CInitial    **initial_container);

		/*!
		 * @brief Function that preprocesses the solver container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the initial container.
		 * @param[in] element_container pointer to the element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 */
		void Solver_Preprocessing(CConfig 		 *config_container,
															CGeometry 	 *geometry_container,
                              CInitial    **initial_container,
															CElement  	**element_container,
															CSpatial    **spatial_container);

		/*!
		 * @brief Function that preprocesses the iteration container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 */
		void Iteration_Preprocessing(CConfig  	  *config_container,
																 CGeometry	  *geometry_container,
																 CSolver  	 **solver_container,
																 CElement 	 **element_container,
																 CSpatial    **spatial_container);

		/*!
		 * @brief Function that preprocesses the temporal container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] iteration_container pointer to the iteration container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 */
		void Temporal_Preprocessing(CConfig      *config_container,
																CGeometry    *geometry_container,
																CIteration  **iteration_container,
																CSolver     **solver_container,
																CElement    **element_container,
																CSpatial    **spatial_container);

    /*!
		 * @brief Function that preprocesses the process container.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] output_container pointer to the output container.
		 * @param[in] element_container pointer to the element container.
		 * @param[in] initial_container pointer to the initial container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] spatial_container pointer to the spatial container.
		 */
    void Process_Preprocessing(CConfig       *config_container,
															 CGeometry     *geometry_container,
                               COutput       *output_container,
															 CElement     **element_container,
                               CInitial     **initial_container,
                               CSolver      **solver_container,
															 CSpatial     **spatial_container);

    /*!
		 * @brief Function that preprocesses the parallelization set-up.
		 */
    void Parallelization_Preprocessing(void);

		/*!
		 * @brief Function that preprocesses the driver once at the start.
		 */
		void Preprocess(void);

		/*!
		 * @brief Function that runs the solver.
		 */
		void Run(void);

		/*!
		 * @brief Function that outputs the data being monitored and header.
		 *
		 * @param[in] iIter current iteration number.
		 * @param[in] time current physical time.
		 * @param[in] dt current time step.
		 * @param[in] MonitoringData vector of data to monitor (passed by reference). 
		 * @param[in] MonitorData option whether to monitor data or not.
		 */
		void MonitorOutput(unsigned long           iIter,
											 as3double               time,
											 as3double               dt,
                       as3vector1d<as3double> &MonitoringData,
											 bool                    MonitorData);

		/*!
		 * @brief Function that estimates the time step.
		 *
		 * @return estimates time step.
		 */
		as3double ComputeTimeStep(void);
};




