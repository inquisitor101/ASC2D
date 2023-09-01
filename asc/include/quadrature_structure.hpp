#pragma once

/*!
 * @file quadrature_structure.hpp
 * @brief The file containing all the numerical quadrature implementations.
 */

#include "option_structure.hpp"


/*!
 * @brief An interface class used for initializing a generic numerical quadrature class.
 */
class CQuadrature {

  public:
 		/*!
		 * @brief Default constructor of CQuadrature, which initializes a generic quadrature class.
		 */
		CQuadrature(void){}

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */   
		virtual ~CQuadrature(void){}


    /*!
		 * @brief Pure virtual function that computes the quadrature points and weights in 1D.
		 * Note, must be overriden by a derived class.
		 *
		 * @param[in] nPoints number of quadrature points in 1D.
		 * @param[out] rQuad quadrature points in 1D.
		 * @param[out] wQuad quadrature weights in 1D.
		 */
    virtual void InitializeGLQuadrature(as3double  nPoints,
                                        as3double *rQuad,
                                        as3double *wQuad) = 0;
};


/*!
 * @brief A class used for initializing a Gauss-Jacobi numerical quadrature class.
 */
class CGaussJacobiQuadrature : public CQuadrature {

	public:
		/*!
		 * @brief Default constructor of CGaussJacobiQuadrature, which initializes a Gauss-Jacobi quadrature class.
		 */	
		CGaussJacobiQuadrature(void);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CGaussJacobiQuadrature(void){}


	  /*!
		 * @brief Function that computes the quadrature points and weights in 1D.
		 *
		 * @param[in] nPoints number of quadrature points in 1D.
		 * @param[out] rQuad quadrature points in 1D.
		 * @param[out] wQuad quadrature weights in 1D.
		 */	
		void InitializeGLQuadrature(as3double  nPoints,
																as3double *rQuad,
																as3double *wQuad) final;

	protected:


	private:
		unsigned int maxIter;  ///< Maximum iterations for the Newton-Rhapson solver.
		as3double    Err;      ///< Error of the solution in the Newton solver.
		as3double    TOL;      ///< Tolerance of the error in the Newton solver.
		as3double 	 dx;       ///< Resolution of root-finder.
		unsigned int nSamples; ///< Total samples for the polynomial curve.

		/*!
		 * @brief Function that computes the Legendre polynomial.
		 *
		 * @param[in] r evaluation point.
		 * @param[in] iDegree degree of the polynomial.
		 *
		 * @return Legendre polynomial of degree (iDegree) evaluated at point (r) 
		 */
		as3double PolyLegendre(as3double 		  r,
								  				 unsigned short iDegree);

		/*!
		 * @brief Function that computes the derivative of the Legendre polynomial.
		 *
		 * @param[in] r evaluation point.
		 * @param[in] iDegree degree of the polynomial.
		 *
		 * @return derivative of the Legendre polynomial of degree (iDegree) evaluated at point (r) 
		 */
		as3double dPolyLegendre(as3double      r,
								 					  unsigned short iDegree);

};


