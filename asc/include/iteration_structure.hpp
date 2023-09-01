#pragma once

/*!
 * @file iteration_structure.hpp
 * @brief The file containing all the iteration steps.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "spatial_structure.hpp"
#include "initial_structure.hpp"
#include "blas.hpp"


/*!
 * @brief An interface class used for iterating over the element in a single time step.
 */
class CIteration {

	public:
		/*!
		 * @brief Default constructor of CIteration, which initializes a generic iteration class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */
		CIteration(CConfig  		 *config_container,
							 CGeometry		 *geometry_container,
							 CSolver  		**solver_container,
							 CElement 		**element_container,
							 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		virtual ~CIteration(void);

		/*!
		 * @brief Pure virtual function that preprocesses every iteration.
		 * Must be overriden by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in[ spatial_container pointer to input spatial container.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical time.
		 * @param[in] iElem current element ID.
		 */
		virtual void Preprocess(CConfig              *config_container,
														CGeometry            *geometry_container,
														CSolver             **solver_container,
														CElement            **element_container,
														CSpatial            **spatial_container,
                            as3data1d<as3double> &work_array,
														as3double             localTime,
                            unsigned long         iElem) = 0;

    /*!
		 * @brief Pure virtual function that iterates over an element to compute its residual.
		 * This sweeps through a single element in this zone and updates the residual accordingly.
		 * Must be overriden by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in[ spatial_container pointer to the current zone spatial container.
		 * @param[in] initial_container pointer to the current zone initial container.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical time.
		 * @param[in] iElem current element ID.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
    virtual void ComputeResidual(CConfig                *config_container,
                                 CGeometry              *geometry_container,
                                 CSolver                *solver_container,
                                 CElement               *element_container,
                                 CSpatial               *spatial_container,
                                 CInitial               *initial_container,
                                 as3data1d<as3double>   &work_array,
                                 as3double               localTime,
                                 unsigned long           iElem,
                                 as3vector1d<as3double> &MonitoringData) = 0;

	protected:
		unsigned short zoneID; ///< Current zone ID.
    unsigned long  nElem;  ///< Number of elements in this zone.

	private:

};


/*!
 * @brief A class used for iterating over the Euler-equation (EE) element in a single time step.
 */
class CEEIteration : public CIteration {

	public:
		/*!
		 * @brief Default constructor of CEEIteration, which initializes an Euler-equation(EE) iteration class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */	
		CEEIteration(CConfig  		 *config_container,
							 	 CGeometry		 *geometry_container,
							 	 CSolver  		**solver_container,
							 	 CElement 		**element_container,
							 	 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */		
		~CEEIteration(void) override;

	protected:
		/*!
		 * @brief Function that preprocesses every iteration.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in[ spatial_container pointer to input spatial container.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical time.
		 * @param[in] iElem current element ID.
		 */
		void Preprocess(CConfig              *config_container,
										CGeometry            *geometry_container,
										CSolver             **solver_container,
										CElement            **element_container,
										CSpatial            **spatial_container,
                    as3data1d<as3double> &work_array,
										as3double             localTime,
                    unsigned long         iElem) override;

 		/*!
		 * @brief Function that iterates over an Euler-equation(EE) element to compute its residual.
		 * This sweeps through a single element in this zone and updates the residual accordingly.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in[ spatial_container pointer to the current zone spatial container.
		 * @param[in] initial_container pointer to the current zone initial container.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical time.
		 * @param[in] iElem current element ID.
		 * @param[out] MonitoringData reference to the data monitored.
		 */   
		void ComputeResidual(CConfig                *config_container,
                         CGeometry              *geometry_container,
                         CSolver                *solver_container,
                         CElement               *element_container,
                         CSpatial               *spatial_container,
                         CInitial               *initial_container,
                         as3data1d<as3double>   &work_array,
                         as3double               localTime,
                         unsigned long           iElem,
                         as3vector1d<as3double> &MonitoringData) override;

	private:

};


