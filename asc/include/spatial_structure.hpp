#pragma once

/*!
 * @file spatial_structure.hpp
 * @brief The file containing all the spatial discretization functionalities.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "riemann_structure.hpp"
#include "data_structure.hpp"
#include "initial_structure.hpp"


/*!
 * @brief A class used for initializing a generic spatial class.
 */
class CSpatial {

	public:
		/*!
		 * @brief Default constructor of CSpatial, which initializes a generic spatial class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] iZone input zone ID.
		 */	
		CSpatial(CConfig   	   *config_container,
						 CGeometry 	   *geometry_container,
						 CElement  	  **element_container,
             CInitial      *initial_container,
						 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */ 
		virtual ~CSpatial(void);

    /*!
		 * @brief Function that applies a filter to the solution.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] data_container pointer to input element data container.
		 */
    void ComputeFilteredSolution(CConfig   *config_container,
                                 CGeometry *geometry_container,
                                 CElement  *element_container,
                                 CData     *data_container);

    /*!
		 * @brief Pure virtual function that initializes the solution needed.
		 * Note, this must be overriden by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] grid_element pointer to input element grid.
		 * @param[in] data_element pointer to input element data container.
		 * @param[in] time current physical (simulation) time.
		 */
    virtual void InitializeSolution(CConfig                *config_container,
                                    CInitial               *initial_container,
                                    CElement               *element_container,
                                    const CGeometryElement *grid_element,
                                    CData                  *data_element,
                                    as3double               time) = 0;

    /*!
		 * @brief Pure virtual function that computes the volume term contribution to the residual. 
		 * Note, this must be overriden by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] data_container pointer to input element data container.
		 * @param[in] geometry_container pointer to input element grid.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
    virtual void ComputeVolumeResidual(CConfig                *config_container,
                                       CGeometry              *geometry_container,
                                       CElement               *element_container,
                                       CInitial               *initial_container,
                                       CData                  *data_container,
                                       const CGeometryElement *geometry_element,
                                       as3data1d<as3double>   &work_array,
                                       as3double               localTime,
                                       as3vector1d<as3double> &MonitoringData) = 0;

    /*!
		 * @brief Pure virtual function that computes the surface term contribution to the residual. 
		 * Note, this must be overriden by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] data_container reference to current zone data container.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[in] iElem input element ID.
		 */
    virtual void ComputeSurfaceResidual(CConfig              *config_container,
                                        CGeometry            *geometry_container,
                                        CElement             *element_container,
                                        as3element           &data_container,
                                        as3data1d<as3double> &work_array,
                                        as3double             localTime,
                                        unsigned long         iElem) = 0;

	protected:
    unsigned short zoneID;                   ///< Current zone ID.
		unsigned short nDOFsSol1D;               ///< Number of solution DOFs per element in this zone (1D).
    unsigned short nDOFsSol2D;               ///< Number of solution DOFs per element in this zone (2D).
		unsigned short nDOFsInt1D;               ///< Number of integration nodes per element in this zone (1D).
    unsigned short nDOFsInt2D;               ///< Number of integration nodes per element in this zone (2D).

    as3double *FilterMatrix      = nullptr;  ///< Filter matrix.
    CRiemann  *riemann_container = nullptr;  ///< Object for the Riemann container.

    /*!
		 * @brief Pure virtual function that preprocesses the Riemann container.
		 * Not, this must be overriden by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] iZone input zone ID.
		 */
    virtual void Riemann_Preprocessing(CConfig       *config_container,
                                       unsigned short iZone) = 0;

    /*!
		 * @brief Function that computes the filter matrix: F = V*C*Vinv.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 */
    void InitializeFilterMatrix(CConfig   *config_container,
                                CGeometry *geometry_container,
                                CElement  *element_container);

	private:

};


/*!
 * @brief A class used for initializing an Euler-equation(EE) spatial class.
 */
class CEESpatial : public CSpatial {

	public:
		/*!
		 * @brief Default constructor of CEESpatial, which initializes an Euler-equation(EE) spatial class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] iZone input zone ID.
		 */	
		CEESpatial(CConfig   	   *config_container,
						 	 CGeometry 	   *geometry_container,
						 	 CElement  	  **element_container,
               CInitial      *initial_container,
						 	 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */ 		
		~CEESpatial(void);

    /*!
		 * @brief Function that initializes the Euler-equation(EE) solution needed.   
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] grid_element pointer to input element grid.
		 * @param[in] data_element pointer to input element data container.
		 * @param[in] time current physical (simulation) time.
		 */
		void InitializeSolution(CConfig                *config_container,
                            CInitial               *initial_container,
                            CElement               *element_container,
                            const CGeometryElement *grid_element,
                            CData                  *data_element,
                            as3double               time) override;

    /*!
		 * @brief Function that computes the EE volume term contribution to the residual.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] data_container pointer to input element data container.
		 * @param[in] geometry_container pointer to input element grid.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
    void ComputeVolumeResidual(CConfig                *config_container,
                               CGeometry              *geometry_container,
                               CElement               *element_container,
                               CInitial               *initial_container,
                               CData                  *data_container,
                               const CGeometryElement *geometry_element,
                               as3data1d<as3double>   &work_array,
                               as3double               localTime,
                               as3vector1d<as3double> &MonitoringData) override;

    /*!
		 * @brief Function that computes the EE surface term contribution to the residual.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] data_container reference to current zone data container.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[in] iElem input element ID.
		 */
    void ComputeSurfaceResidual(CConfig              *config_container,
                                CGeometry            *geometry_container,
                                CElement             *element_container,
                                as3element           &data_container,
                                as3data1d<as3double> &work_array,
                                as3double             localTime,
                                unsigned long         iElem) override;

	protected:
    /*!
		 * @brief Function that preprocesses the Riemann container.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] iZone input zone ID.
		 */
    void Riemann_Preprocessing(CConfig       *config_container,
                               unsigned short iZone) final;

	private:

};


/*!
 * @brief A class used for initializing an Euler-equation(EE) sponge-layer spatial class.
 */
class CEESpongeSpatial : public CEESpatial {

	public:
		/*!
		 * @brief Default constructor of CEESpongeSpatial, which initializes an EE sponge-layer spatial class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] iZone input zone ID.
		 */
		CEESpongeSpatial(CConfig   	    *config_container,
      						 	 CGeometry 	    *geometry_container,
      						 	 CElement   	 **element_container,
                     CInitial       *initial_container,
      						 	 unsigned short  iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CEESpongeSpatial(void);

    /*!
		 * @brief Function that initializes the EE sponge-layer solution needed.   
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] grid_element pointer to input element grid.
		 * @param[in] data_element pointer to input element data container.
		 * @param[in] time current physical (simulation) time.
		 */
		void InitializeSolution(CConfig                *config_container,
                            CInitial               *initial_container,
                            CElement               *element_container,
                            const CGeometryElement *grid_element,
                            CData                  *data_element,
                            as3double               time) final;

    /*!
		 * @brief Function that computes the EE sponge-layer volume term contribution to the residual.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] data_container pointer to input element data container.
		 * @param[in] geometry_container pointer to input element grid.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
		void ComputeVolumeResidual(CConfig                *config_container,
                               CGeometry              *geometry_container,
                               CElement               *element_container,
                               CInitial               *initial_container,
                               CData                  *data_container,
                               const CGeometryElement *geometry_element,
                               as3data1d<as3double>   &work_array,
                               as3double               localTime,
                               as3vector1d<as3double> &MonitoringData) final;

    /*!
		 * @brief Function that computes the EE sponge-layer surface term contribution to the residual.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] data_container reference to current zone data container.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[in] iElem input element ID.
		 */
		void ComputeSurfaceResidual(CConfig              *config_container,
                                CGeometry            *geometry_container,
                                CElement             *element_container,
                                as3element           &data_container,
                                as3data1d<as3double> &work_array,
                                as3double             localTime,
                                unsigned long         iElem) final;

	protected:

	private:
    as3vector1d<bool> ArtificialConvectionFace;  ///< Vector of options for using artificial-convection on each face.
};


/*!
 * @brief A class used for initializing an Euler-equation(EE) perfectly matched layer spatial class.
 */
class CEEPMLSpatial : public CEESpatial {

	public:
		/*!
		 * @brief Default constructor of CEEPMLSpatial, which initializes an EE perfectly matched layer spatial class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] iZone input zone ID.
		 */	
		CEEPMLSpatial(CConfig   	   *config_container,
    						 	CGeometry 	   *geometry_container,
    						 	CElement  	  **element_container,
                  CInitial       *initial_container,
    						 	unsigned short  iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CEEPMLSpatial(void);

    /*!
		 * @brief Function that initializes the EE perfectly matched layer solution needed.   
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] grid_element pointer to input element grid.
		 * @param[in] data_element pointer to input element data container.
		 * @param[in] time current physical (simulation) time.
		 */ 
		void InitializeSolution(CConfig                *config_container,
                            CInitial               *initial_container,
                            CElement               *element_container,
                            const CGeometryElement *grid_element,
                            CData                  *data_element,
                            as3double               time) final;

    /*!
		 * @brief Function that computes the EE perfectly matched layer volume term contribution to the residual.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] data_container pointer to input element data container.
		 * @param[in] geometry_container pointer to input element grid.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[out] MonitoringData reference to the data monitored.
		 */
		void ComputeVolumeResidual(CConfig                *config_container,
                               CGeometry              *geometry_container,
                               CElement               *element_container,
                               CInitial               *initial_container,
                               CData                  *data_container,
                               const CGeometryElement *geometry_element,
                               as3data1d<as3double>   &work_array,
                               as3double               localTime,
                               as3vector1d<as3double> &MonitoringData) final;

    /*!
		 * @brief Function that computes the EE perfectly matched layer surface term contribution to the residual.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] data_container reference to current zone data container.
		 * @param[in] work_array reference to the work array, which temporarily stores data.
		 * @param[in] localTime current physical (simulation) time.
		 * @param[in] iElem input element ID.
		 */
		void ComputeSurfaceResidual(CConfig              *config_container,
                                CGeometry            *geometry_container,
                                CElement             *element_container,
                                as3element           &data_container,
                                as3data1d<as3double> &work_array,
                                as3double             localTime,
                                unsigned long         iElem) final;

	protected:

	private:
    as3double DispersionCorrection;  ///< Dispersion-relation correction coefficient.
    as3double VelocityTransverse;    ///< Background velocity in the transverse direction.
};
