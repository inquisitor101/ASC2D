#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "quadrature_structure.hpp"
#include "blas.hpp"



class CElement {

	public:
		// Constructor.
		CElement(CConfig       *config_container,
						 unsigned short iZone);

		// Destructor.
		~CElement(void);


    // Function that computes the location of the DOFs in 1D.
    void LocationDOFs1D(unsigned short          TypeDOFs,
                        as3vector1d<as3double> &rDOFs1D);

    // Function that computes the 1D lagrange basis functions and first-order derivatives.
    void LagrangeBasisFunctions(const as3vector1d<as3double> &rBasis,
                                const as3vector1d<as3double> &rPoints,
                                as3double                    *lagrangeDOFs1D,
                                as3double                    *derLagrangeDOFs1D,
                                as3double                    *lagrangeDOFs1DTranspose,
                                as3double                    *derLagrangeDOFs1DTranspose,
                                as3double                    *derLagrangeDOFsMinFace1D,
                                as3double                    *derLagrangeDOFsMaxFace1D);

    // Function that computes the mass matrix.
    void ComputeMassMatrix(unsigned short nDOFs,
                           as3double     *MassMatrix);

    // Getter: returns TypeDOFs.
    unsigned short GetTypeDOFs(void)                              const {return TypeDOFs;}
    // Getter: returns nPolySol.
    unsigned short GetnPolySol(void)                              const {return nPolySol;}
    // Getter: returns nDOFsSol1D,
    unsigned short GetnDOFsSol1D(void)              			 				const {return nDOFsSol1D;}
    // Getter: returns nDOFsSol2D.
    unsigned short GetnDOFsSol2D(void)              			 				const {return nDOFsSol2D;}
    // Getter: returns nDOFsInt1D.
    unsigned short GetnDOFsInt1D(void)              			 				const {return nDOFsInt1D;}
    // Getter: returns nDOFsInt2D.
    unsigned short GetnDOFsInt2D(void)              			 				const {return nDOFsInt2D;}
    // Getter: returns rDOFsSol1D.
    const as3vector1d<as3double> &GetrDOFsSol1D(void)             const {return rDOFsSol1D;}
    // Getter: returns rDOFsInt1D.
    const as3vector1d<as3double> &GetrDOFsInt1D(void)             const {return rDOFsInt1D;}
    // Getter: returns wDOFsInt1D.
    const as3vector1d<as3double> &GetwDOFsInt1D(void)             const {return wDOFsInt1D;}
    // Getter: returns wDOFsInt2D.
    const as3vector1d<as3double> &GetwDOFsInt2D(void)             const {return wDOFsInt2D;}
    // Getter: returns Vandermonde1D.
    const as3double *GetVandermonde1D(void)                       const {return Vandermonde1D;}
    // Getter: returns InvVandermonde1D.
    const as3double *GetInvVandermonde1D(void)                    const {return InvVandermonde1D;}
    // Getter: returns lagrangeInt1D.
    const as3double *GetLagrangeInt1D(void)      			          	const {return lagrangeInt1D;}
    // Getter: returns derLagrangeInt1D.
    const as3double *GetDerLagrangeInt1D(void)   			            const {return derLagrangeInt1D;}
    // Getter: returns lagrangeInt1DTranspose.
		const as3double *GetLagrangeInt1DTranspose(void)      			 	const {return lagrangeInt1DTranspose;}
		// Getter: returns derLagrangeInt1DTranspose.
		const as3double *GetDerLagrangeInt1DTranspose(void)   			  const {return derLagrangeInt1DTranspose;}
    // Getter: returns derLagrangeDOFsSol1DMinFace.
		const as3double *GetDerLagrangeDOFsSol1DMinFace(void)       	const {return derLagrangeDOFsSol1DMinFace;}
		// Getter: returns derLagrangeDOFsSol1DMaxFace.
		const as3double *GetDerLagrangeDOFsSol1DMaxFace(void)       	const {return derLagrangeDOFsSol1DMaxFace;}
    // Getter: returns derLagrangeDOFsSol1DMinFace or derLagrangeDOFsSol1DMaxFace.
		const as3double *GetDerLagrangeDOFsSol1DFace(unsigned short iFace) const {
			if( iFace == IDX_SOUTH || iFace == IDX_WEST )
				return derLagrangeDOFsSol1DMinFace;
			else
				return derLagrangeDOFsSol1DMaxFace;
		}
    // Getter: returns IndexDOFsSol.
    const as3vector2d<unsigned short> &GetIndexDOFsSol(void)      const {return IndexDOFsSol;}
    // Getter: returns IndexDOFsSol, per input face.
    const as3vector1d<unsigned short> &GetIndexDOFsSol(unsigned short iFace) const {
      return IndexDOFsSol[iFace];
    }

  protected:

  private:
    // Zone ID.
		unsigned short zoneID;
		// Solution polynomial order.
		unsigned short nPolySol;
		// Number of nodes in 1D, based on solution DOFs.
		unsigned short nDOFsSol1D;
    // Number of nodes in 2D, based on solution DOFs.
		unsigned short nDOFsSol2D;
		// Integration rule in 1D, based on solution.
		unsigned short nDOFsInt1D;
    // Integration rule in 2D, based on solution.
		unsigned short nDOFsInt2D;
    // Type of DOFs used.
    unsigned short TypeDOFs;

    // Indices of SolDOFs on all four sides of the element.
		as3vector2d<unsigned short> IndexDOFsSol;

    // Solution nodes on a reference line, based on the solution polynomial.
    as3vector1d<as3double> rDOFsSol1D;
    // Integration nodes on a reference line, based on the solution polynomial integration rule.
    as3vector1d<as3double> rDOFsInt1D;
    // Integration weights on a reference line, based on the solution polynomial integration rule.
    as3vector1d<as3double> wDOFsInt1D;
    // Integration weights on entire element, based on the solution polynomial integration rule.
    as3vector1d<as3double> wDOFsInt2D;

    // Vandermonde matrix in 1D that converts from the modal-to-nodal variables.
    as3double *Vandermonde1D    = nullptr;
    // Inverse of the Vandermonde matrix in 1D that converts from the nodal-to-modal variables.
    as3double *InvVandermonde1D = nullptr;

    // Lagrange polynomials in 1D: interpolates from SolDOFs to IntDOFs.
    as3double *lagrangeInt1D    = nullptr;
    // Derivative of lagrange polynomials in 1D: interpolates from SolDOFs to IntDOFs.
    as3double *derLagrangeInt1D = nullptr;

    // Lagrange polynomials (transposed) in 1D: interpolates from SolDOFs to IntDOFs.
    as3double *lagrangeInt1DTranspose    = nullptr;
    // Derivative of lagrange polynomials (transposed) in 1D: interpolates from SolDOFs to IntDOFs.
    as3double *derLagrangeInt1DTranspose = nullptr;

    // Lagrange differentiation vector in 1D on min face: SolDOFs.
    as3double *derLagrangeDOFsSol1DMinFace = nullptr;
    // Lagrange differentiation vector in 1D on max face: SolDOFs.
    as3double *derLagrangeDOFsSol1DMaxFace = nullptr;


    // Function that initializes the quadrature integration data.
    void InitializeQuadratureIntegration(as3vector1d<as3double> &r1D,
                                         as3vector1d<as3double> &w1D,
                                         as3vector1d<as3double> &w2D);

    // Function that determines the indices of each of the element faces
    // based on the polynomial degree of this container and zone.
    void IdentifyFaceIndices();

    // Function that computes the Lagrange interpolation entry in 1D for a given basis.
    as3double EvaluateEntryLagrangePolynomial(const as3double               xPoint,
                                              const as3vector1d<as3double> &rBasis,
                                              const unsigned short          iDegree);

    // Function that computes the derivative of lagrange interpolation entry in 1D
    // for a given basis.
    as3double EvaluateEntryDerLagrangePolynomial(const as3double               xPoint,
                                                 const as3vector1d<as3double> &rBasis,
                                                 const unsigned short          iDegree);

    // Function that computes the 1D Vandermonde matrix.
    void ComputeVandermonde1D(const unsigned short          nCol,
                              const as3vector1d<as3double> &r,
                              as3double                    *V);

    // Function that computes the generalized normalized-Jacobi polynomial for an
    // order n and at a given entry x.
    as3double JacobiNormalized(unsigned short n,
                               unsigned short alpha,
                               unsigned short beta,
                               as3double      x);
  };