#pragma once

/*!
 * @file temporal_structure.hpp
 * @brief The file containing all the temporal discretization functionalities.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "iteration_structure.hpp"
#include "spatial_structure.hpp"
#include "initial_structure.hpp"


/*!
 * @brief A class used for initializing a generic temporal class.
 */
class CTemporal {

	public:
		/*!
		 * @brief Default constructor of CTemporal, which initializes a generic temporal class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] iteration_container pointer to input iteration container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] input_MapGlobalToLocal reference to the indices of the elements.
		 */
		CTemporal(CConfig                    *config_container,
							CGeometry                  *geometry_container,
							CIteration                **iteration_container,
							CSolver                   **solver_container,
							CElement                  **element_container,
							CSpatial                  **spatial_container,
              as3vector2d<unsigned long> &input_MapGlobalToLocal);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */ 
		virtual ~CTemporal(void);

		/*!
		 * @brief Pure virtual function that performs an update over the entire simulation in time.
		 * Must be overriden by one of the derived classes.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] iteration_container pointer to input iteration container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] initial_container pointer to input initial conditions container.
		 * @param[in] physicalTime current physical (simulation) time.
		 * @param[in] dtTime current time step.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
		virtual void TimeMarch(CConfig                *config_container,
									 				 CGeometry              *geometry_container,
									 				 CIteration            **iteration_container,
									 				 CSolver               **solver_container,
									 				 CElement              **element_container,
									 				 CSpatial              **spatial_container,
                           CInitial              **initial_container,
									 				 as3double               physicalTime,
													 as3double               dtTime,
                           as3vector1d<as3double> &MonitoringData) = 0;

	protected:
		unsigned short              nZone;              ///< Total number of zones.
		as3vector1d<unsigned long>  nElemZone;          ///< Number of elements per each zone.
		as3vector1d<unsigned short> nNodeZone;          ///< Number of solution DOFs in each element per each zone.
    unsigned long               nElemTotal;         ///< Total number of elements in all zones combined.
    unsigned short              nWorkingArrayDOFs;  ///< Number of DOFs in the working array.
    unsigned short              nWorkingArrayVar;   ///< Number of variables in the working array.
    unsigned short              nWorkingArrayEntry; ///< Number of data entries in the working array.
    as3vector2d<unsigned long>  MapGlobalToLocal;   ///< Element index data that maps from [iZone][iElemZone] to [iElem].
																										///< This is used for parallelization efficiency.
																										///< Dimension: [iElem][iData], where: [iData] is: [0]: iZone, [1]: iElemZone.

    /*!
		 * @brief Function that initializes and defined the dimension parameters needed in the working array.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] element_container pointer to input standard element container.
		 */
    void InitializeWorkArrayDimension(CConfig   *config_container,
                                      CElement **element_container);

    /*!
		 * @brief Function that initializes the working array by reserving its memory.
		 *
		 * @param[in] work_array reference to the working array.
		 */
    void InitializeWorkArray(as3data1d<as3double> &work_array);

	private:

};


/*!
 * @brief A class used for initializing a low-storage 4th-order Runge-Kutta (LSRK4) temporal class.
 */
class CLSRK4Temporal : public CTemporal {

	public:
		/*!
		 * @brief Default constructor of CLSRK4Temporal, which initializes a LSRK4 temporal class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] iteration_container pointer to input iteration container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] input_MapGlobalToLocal reference to the indices of the elements.
		 */
		CLSRK4Temporal(CConfig                    *config_container,
									 CGeometry                  *geometry_container,
									 CIteration                **iteration_container,
									 CSolver                   **solver_container,
									 CElement                  **element_container,
									 CSpatial                  **spatial_container,
                   as3vector2d<unsigned long> &input_MapGlobalToLocal);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */ 
		~CLSRK4Temporal(void) final;

		/*!
		 * @brief Function that performs an update over the entire simulation in time via LSRK4.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] iteration_container pointer to input iteration container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] initial_container pointer to input initial conditions container.
		 * @param[in] physicalTime current physical (simulation) time.
		 * @param[in] dtTime current time step.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
		void TimeMarch(CConfig                *config_container,
					 				 CGeometry              *geometry_container,
					 				 CIteration            **iteration_container,
					 				 CSolver               **solver_container,
					 				 CElement              **element_container,
					 				 CSpatial              **spatial_container,
                   CInitial              **initial_container,
					 				 as3double               physicalTime,
									 as3double               dtTime,
                   as3vector1d<as3double> &MonitoringData) final;

	protected:

	private:
		unsigned short nStageRK = LSRK4_N_STAGES;   ///< Number of RK stages: this is a 5-stage scheme.
		unsigned short nStorage = LSRK4_N_STORAGE;  ///< Number of storage for tentative data.

		as3vector1d<as3double> rk4a;                ///< LSRK4: a-coefficients.
		as3vector1d<as3double> rk4b;                ///< LSRK4: b-coefficients.
		as3vector1d<as3double> rk4c;                ///< LSRK4: c-coefficients.

		as3data3d<as3double> DataDOFsSolTentative;  ///< Tentative/intermediate solution stored.
																								///< Dimension: [iZone][iElem][iVar][iNode].

    /*!
		 * @brief Function that performs a single stage sweep of a LSRK4.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] iteration_container pointer to input iteration container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] initial_container pointer to input initial conditions container.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[in] dtTime current time step.
		 * @param[in] alpha input LSRK4: a-coefficient.
		 * @param[in] beta input LSRK4: b-coefficient.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
		void UpdateTime(CConfig                *config_container,
									  CGeometry              *geometry_container,
									  CIteration            **iteration_container,
									  CSolver               **solver_container,
									  CElement              **element_container,
									  CSpatial              **spatial_container,
                    CInitial              **initial_container,
									  as3double               localTime,
                    as3double               dtTime,
                    as3double               alpha,
                    as3double               beta,
                    as3vector1d<as3double> &MonitoringData);
};


/*!
 * @brief A class used for initializing a strong stability preserving 3rd-order Runge Kutta (SSPRK3) temporal class.
 */
class CSSPRK3Temporal : public CTemporal {

	public:
		/*!
		 * @brief Default constructor of CSSPRK3Temporal, which initializes a SSPRK3 temporal class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] iteration_container pointer to input iteration container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] input_MapGlobalToLocal reference to the indices of the elements.
		 */
		CSSPRK3Temporal(CConfig                    *config_container,
									  CGeometry                  *geometry_container,
									  CIteration                **iteration_container,
									  CSolver                   **solver_container,
									  CElement                  **element_container,
									  CSpatial                  **spatial_container,
                    as3vector2d<unsigned long> &input_MapGlobalToLocal);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */ 
		~CSSPRK3Temporal(void) final;

		/*!
		 * @brief Function that performs an update over the entire simulation in time via SSPRK3.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] iteration_container pointer to input iteration container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] initial_container pointer to input initial conditions container.
		 * @param[in] physicalTime current physical (simulation) time.
		 * @param[in] dtTime current time step.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
		void TimeMarch(CConfig                *config_container,
					 				 CGeometry              *geometry_container,
					 				 CIteration            **iteration_container,
					 				 CSolver               **solver_container,
					 				 CElement              **element_container,
					 				 CSpatial              **spatial_container,
                   CInitial              **initial_container,
					 				 as3double               physicalTime,
									 as3double               dtTime,
                   as3vector1d<as3double> &MonitoringData) final;

	protected:

	private:
		unsigned short nStageRK = SSPRK3_N_STAGES;   ///< Number of RK stages: this is a 3-stage scheme.
		unsigned short nStorage = SSPRK3_N_STORAGE;  ///< Number of storage for tentative data.

		as3vector1d<as3double> rk4a;                 ///< SSPRK3: a-coefficients.
		as3vector1d<as3double> rk4b;                 ///< SSPRK3: b-coefficients.
		as3vector1d<as3double> rk4c;                 ///< SSPRK3: c-coefficients.

		as3data3d<as3double> DataDOFsSolTentative;   ///< Tentative/intermediate solution stored.
																								 ///< Dimension: [iZone][iElem][iVar][iNode].

    /*!
		 * @brief Function that performs a single stage sweep of a SSPRK3.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] iteration_container pointer to input iteration container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] initial_container pointer to input initial conditions container.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[in] dtTime current time step.
		 * @param[in] alpha input LSRK4: a-coefficient.
		 * @param[in] beta input LSRK4: b-coefficient.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
		void UpdateTime(CConfig                *config_container,
									  CGeometry              *geometry_container,
									  CIteration            **iteration_container,
									  CSolver               **solver_container,
									  CElement              **element_container,
									  CSpatial              **spatial_container,
                    CInitial              **initial_container,
									  as3double               localTime,
                    as3double               dtTime,
                    as3double               alpha,
                    as3double               beta,
                    as3vector1d<as3double> &MonitoringData);
};

