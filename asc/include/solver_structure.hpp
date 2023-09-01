#pragma once

/*!
 * @file solver_structure.hpp
 * @brief The file containing all the solver functionalities.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "data_structure.hpp"
#include "element_structure.hpp"
#include "boundary_structure.hpp"
#include "spatial_structure.hpp"
#include "initial_structure.hpp"
#include "blas.hpp"

// Forward declaration to avoid compiler problems.
class CBoundary;


/*!
 * @brief A class used for initializing a generic solver class.
 */
class CSolver {

	public:
		/*!
		 * @brief Default constructor of CSolver, which initializes a generic solver.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] iZone input zone ID.
		 */	
		CSolver(CConfig   		 *config_container,
						CGeometry 		 *geometry_container,
            CInitial       *initial_container,
						CElement      **element_container,
						CSpatial      **spatial_container,
						unsigned short  iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */ 	
		virtual ~CSolver(void);


    /*!
		 * @brief Pure virtual function that computes the maximum stable time step per element.
		 * Must be overriden by a derived class.
		 *
		 * @param[in] iElem input element ID.
		 * @param[in] hElem reference to input element dimension.
		 *
		 * @return estimated time step
		 */
    virtual as3double ComputeTimeStep(unsigned long                 iElem,
                                      const as3vector1d<as3double> &hElem) = 0;

		/*!
		 * @brief Getter function which returns the value of nDOFsSol2D.
		 *
		 * @return nDOFsSol2D
		 */
		unsigned short GetnDOFsSol2D(void)           const {return nDOFsSol2D;}
		/*!
		 * @brief Getter function which returns the value of nElem.
		 *
		 * @return nElem
		 */
		unsigned long GetnElem(void)             		 const {return nElem;}
    /*!
		 * @brief Getter function which returns the value of InvMassMatrix.
		 *
		 * @return InvMassMatrix
		 */
    const as3double *GetInvMassMatrix(void)      const {return InvMassMatrix;}
		/*!
		 * @brief Getter function which returns the value of data_container.
		 *
		 * @return data_container
		 */
		as3element &GetDataContainer(void) 		             {return data_container;}
		/*!
		 * @brief Getter function which returns the value of data_container[iElem].
		 *
		 * @param[in] iElem input element ID.
		 *
		 * @return data_container[iElem]
		 */
		CData* GetDataContainer(unsigned long iElem) const {return data_container[iElem];}
		/*!
		 * @brief Getter function which returns the value of nBoundary.
		 *
		 * @return boundary_container.size()
		 */
		unsigned short GetnBoundary(void)            const {return boundary_container.size();}
    /*!
		 * @brief Getter function which returns the value of boundary_container[iBoundary].
		 *
		 * @param[in] iBoundary input boundary ID.
		 *
		 * @return boundary_container[iBoundary]
		 */
    CBoundary* GetBoundaryContainer(unsigned short iBoundary) const {
      return boundary_container[iBoundary];
    }

	protected:
		unsigned short zoneID;                         ///< Current zone ID.
		unsigned short nPoly;                          ///< Current solution polynomial order.
		unsigned short nDOFsSol2D;                     ///< Number of solution DOFs in 2D over an element.
		unsigned long  nElem;                          ///< Total number of elements in this zone.
    unsigned short nBoundary;                      ///< Total number of boundary conditions in this zone.

		as3element           data_container;           ///< Object for each element data container.
    as3double           *InvMassMatrix = nullptr;  ///< Inverse mass matrix in this zone.
		as3data1d<CBoundary> boundary_container;       ///< Object for each boundary condition container.

    /*!
		 * @brief Function that preprocesses and initializes the data container.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] iZone input zone ID.
		 */
    void Data_Preprocessing(CConfig       *config_container,
    												CGeometry     *geometry_container,
                            CInitial      *initial_container,
    												CElement      *element_container,
    												unsigned short iZone);

    /*!
		 * @brief Function used to initialize the boundary container.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */
    void Boundary_Preprocessing(CConfig  	    *config_container,
    														CGeometry	    *geometry_container,
                                CInitial      *initial_container,
    														CElement 	   **element_container,
    														unsigned short iZone);
	private:

};


// IDEA
// when extending solver to handle NS, create another EE class names EENS
// and use that one to handle the EE implementation as an intermediate state
// before inheriting it in the NS class.
// Use this one (CEESolver) as a pure Euler solver (i.e. no NS is used) which
// also inherits the CEENSSolver class, but the end node in that class tree is
// CEESolver instead of CEENSSolver. Note, don't forget to move any common code
// from CEESolver to CEENSSolver.
// This should mitigate the tedious condition statements when preprocessing
// certain objects such as boundary_container and data_container, etc...

/*!
 * @brief A class used for initializing an Euler-equation(EE) solver class.
 */
class CEESolver : public CSolver {

	public:
		/*!
		 * @brief Default constructor of CEESolver, which initializes an Euler-equation(EE) solver.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] spatial_container pointer to input spatial container.
		 * @param[in] iZone input zone ID.
		 */		
		CEESolver(CConfig   		 *config_container,
							CGeometry 		 *geometry_container,
              CInitial       *initial_container,
							CElement      **element_container,
							CSpatial      **spatial_container,
							unsigned short  iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */ 		
		~CEESolver(void) override;

    /*!
		 * @brief Function that computes the maximum stable time step per element.
		 *
		 * @param[in] iElem input element ID.
		 * @param[in] hElem reference to input element dimension.
		 *
		 * @return estimated time step
		 */
    as3double ComputeTimeStep(unsigned long                 iElem,
                              const as3vector1d<as3double> &hElem);

	protected:
    /*!
		 * @brief Function that computes the inverse mass matrix.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 */
    void ComputeInvMassMatrix(CConfig   *config_container,
                              CGeometry *geometry_container,
                              CElement  *element_container);

	private:

};


