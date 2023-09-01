#pragma once

/*!
 * @file process_structure.hpp
 * @brief The file containing all the (post-)processing functionalities.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "spatial_structure.hpp"
#include "output_structure.hpp"
#include "initial_structure.hpp"
#include "probe_structure.hpp"
#include "reflection_structure.hpp"


/*!
 * @brief An interface class used for initializing a generic process class.
 */
class CProcess {

  public:
 		/*!
		 * @brief Default constructor of CProcess, which initializes a generic process class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] output_container pointer to the output container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] initial_container pointer to input initial conditions container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] iZone input zone ID.
		 */
		CProcess(CConfig       *config_container,
             CGeometry     *geometry_container,
             COutput       *output_container,
             CElement     **element_container,
             CInitial     **initial_container,
             CSolver      **solver_container,
             CSpatial     **spatial_container,
             unsigned short iZone);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */   
		virtual ~CProcess(void);

    /*!
		 * @brief Function that filters the solution in this zone.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] spatial_container pointer to current zone spatial container.
		 */
    void FilterSolution(CConfig   *config_container,
                        CGeometry *geometry_container,
                        CElement  *element_container,
                        CSolver   *solver_container,
                        CSpatial  *spatial_container);

    /*!
		 * @brief Pure virtual function that specifies the processing method.
		 * Must be overriden by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] spatial_container pointer to current zone spatial container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] output_container pointer to output container.
		 * @param[in] localTime current physical (simulation) time.
		 */
    virtual void ProcessData(CConfig   *config_container,
                             CGeometry *geometry_container,
                             CElement  *element_container,
                             CSpatial  *spatial_container,
                             CSolver   *solver_container,
                             CInitial  *initial_container,
														 COutput   *output_container,
                             as3double  localTime) = 0;

  protected:
    unsigned short zoneID;  ///< Current zone ID.

  private:

};


/*!
 * @brief A class used for initializing an Euler-equation(EE) process class.
 */
class CEEProcess : public CProcess {

  public:
 		/*!
		 * @brief Default constructor of CEEProcess, which initializes an Euler-equation(EE) process class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] output_container pointer to the output container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] initial_container pointer to input initial conditions container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] iZone input zone ID.
		 */ 
		CEEProcess(CConfig       *config_container,
               CGeometry     *geometry_container,
               COutput       *output_container,
               CElement     **element_container,
               CInitial     **initial_container,
               CSolver      **solver_container,
               CSpatial     **spatial_container,
               unsigned short iZone);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */      
		~CEEProcess(void) override;

    /*!
		 * @brief Function that specifies the processing method.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] spatial_container pointer to current zone spatial container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] output_container pointer to output container.
		 * @param[in] localTime current physical (simulation) time.
		 */  
		void ProcessData(CConfig   *config_container,
                     CGeometry *geometry_container,
                     CElement  *element_container,
                     CSpatial  *spatial_container,
                     CSolver   *solver_container,
                     CInitial  *initial_container,
										 COutput   *output_container,
                     as3double  localTime) override;

  protected:

	private:
    as3data1d<CProbe>      probe_container;      ///< Probe data object.
		as3data1d<CReflection> reflection_container; ///< Reflection coefficient object.
};



