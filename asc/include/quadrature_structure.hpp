#pragma once

#include "option_structure.hpp"



class CQuadrature {

  public:
    // Constructor.
    CQuadrature(void){}

    // Destructor.
    virtual ~CQuadrature(void){}


    // Pure virtual function that computes the quadrature points and weights in 1D.
    // Note, must be overriden by a derived class.
    virtual void InitializeGLQuadrature(as3double  nPoints,
                                        as3double *rQuad,
                                        as3double *wQuad) = 0;
};


class CGaussJacobiQuadrature : public CQuadrature {

	public:
		// Constructor.
		CGaussJacobiQuadrature(void);

		// Destructor.
		~CGaussJacobiQuadrature(void){}


		// Function that computes the quadrature points and weights in 1D.
		void InitializeGLQuadrature(as3double  nPoints,
																as3double *rQuad,
																as3double *wQuad) final;

	protected:


	private:
		// Newton solver configuration.
		unsigned int maxIter;
		as3double    Err;
		as3double    TOL;

		// Resolution of root-finder.
		as3double 	 dx;

		// Total samples for polynomial curve.
		unsigned int nSamples;

		// Function that computes the legendre polynomial.
		as3double PolyLegendre(as3double 		  r,
								  				 unsigned short iDegree);

		// Function that computes the derivative of the legendre polynomial.
		as3double dPolyLegendre(as3double      r,
								 					  unsigned short iDegree);

};


