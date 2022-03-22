#pragma once

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


class CSolver {

	public:
		// Constructor.
		CSolver(CConfig   		 *config_container,
						CGeometry 		 *geometry_container,
            CInitial       *initial_container,
						CElement      **element_container,
						CSpatial      **spatial_container,
						unsigned short  iZone);

		// Destructor.
		virtual ~CSolver(void);


    // Pure virtual function that computes the maximum stable time step per element.
    // Must be overriden by a derived class.
    virtual as3double ComputeTimeStep(unsigned long                 iElem,
                                      const as3vector1d<as3double> &hElem) = 0;

		// Getter: returns nDOFsSol2D.
		unsigned short GetnDOFsSol2D(void)           const {return nDOFsSol2D;}
		// Getter: returns nElem.
		unsigned long GetnElem(void)             		 const {return nElem;}
    // Getter: returns InvMassMatrix.
    const as3double *GetInvMassMatrix(void)      const {return InvMassMatrix;}
		// Getter: returns data_container.
		as3element &GetDataContainer(void) 		             {return data_container;}
		// Getter: returns data_container[iElem].
		CData* GetDataContainer(unsigned long iElem) const {return data_container[iElem];}
		// Getter: returns nBoundary.
		unsigned short GetnBoundary(void)            const {return boundary_container.size();}
    // Getter: returns boundary_container[iBoundary].
    CBoundary* GetBoundaryContainer(unsigned short iBoundary) const {
      return boundary_container[iBoundary];
    }

	protected:
		// Zone ID.
		unsigned short zoneID;
		// Polynomial order.
		unsigned short nPoly;
		// Number of solution DOFs in 2D.
		unsigned short nDOFsSol2D;
		// Number of elements.
		unsigned long  nElem;
    // Number of boundary conditions in this zone.
    unsigned short nBoundary;

		// Data element container.
		as3element data_container;
    // Inverse mass matrix.
    as3double *InvMassMatrix = nullptr;
		// Boundary condition container.
		as3data1d<CBoundary> boundary_container;

    // Function that preprocesses and initializes the data container.
    void Data_Preprocessing(CConfig       *config_container,
    												CGeometry     *geometry_container,
                            CInitial      *initial_container,
    												CElement      *element_container,
    												unsigned short iZone);

    // Function used to initialize the boundary container.
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
class CEESolver : public CSolver {

	public:
		// Constructor.
		CEESolver(CConfig   		 *config_container,
							CGeometry 		 *geometry_container,
              CInitial       *initial_container,
							CElement      **element_container,
							CSpatial      **spatial_container,
							unsigned short  iZone);

		// Destructor.
		~CEESolver(void) override;

    // Function that computes the maximum stable time step per element.
    // Must be overriden by a derived class.
    as3double ComputeTimeStep(unsigned long                 iElem,
                              const as3vector1d<as3double> &hElem);

	protected:
    // Function that computes the inverse mass matrix.
    void ComputeInvMassMatrix(CConfig   *config_container,
                              CGeometry *geometry_container,
                              CElement  *element_container);

	private:

};


