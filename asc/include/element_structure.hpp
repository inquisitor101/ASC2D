#pragma once

/*!
 * @file element_structure.hpp
 * @brief The file responsible for the definitions of a standard element in reference space.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "quadrature_structure.hpp"
#include "blas.hpp"


/*!
 * @brief A class used for initializing and definining all operations on a standard element in reference space.
 */
class CElement {

	public:
		/*!
		 * @brief Default constructor of CElement, which initializes CElement class.
		 *
		 * @param[in] config input configuration/dictionary container.
		 * @param[in] iZone current zone ID.
		 */
		CElement(CConfig       *config_container,
						 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CElement(void);


    /*!
		 * @brief Function that computes the location of the DOFs in 1D.
		 *
		 * @param[in] TypeDOFs type of DOFs used.
		 * @param[in] rDOFs1D DOFs in reference space in 1D (passed by reference).
		 */
    void LocationDOFs1D(unsigned short          TypeDOFs,
                        as3vector1d<as3double> &rDOFs1D);

    /*!
		 * @brief Function that computes the 1D lagrange basis functions and first-order derivatives.
		 *
		 * param[in] rBasis basis points in 1D (passed by reference).
		 * param[in] rPoints evaluation points in 1D (passed by reference).
		 * param[out] lagrangeDOFs1D pointer to the Lagrange interpolating function in 1D.
		 * param[out] derLagrangeDOFs1D pointer to the Lagrange differential function in 1D.
		 * param[out] lagrangeDOFs1DTranspose pointer to the transpose of the Lagrange interpolating function in 1D.
		 * param[out] derLagrangeDOFs1DTranspose pointer to the transpose of the Lagrange differential function in 1D.
		 * param[out] derLagrangeDOFsMinFace1D pointer to the Lagrange differential function on the min face in 1D.
		 * param[out] derLagrangeDOFsMaxFace1D pointer to the Lagrange differential function on the max face in 1D.
		 */
    void LagrangeBasisFunctions(const as3vector1d<as3double> &rBasis,
                                const as3vector1d<as3double> &rPoints,
                                as3double                    *lagrangeDOFs1D,
                                as3double                    *derLagrangeDOFs1D,
                                as3double                    *lagrangeDOFs1DTranspose,
                                as3double                    *derLagrangeDOFs1DTranspose,
                                as3double                    *derLagrangeDOFsMinFace1D,
                                as3double                    *derLagrangeDOFsMaxFace1D);

    /*!
		 * @brief Function that computes the mass matrix.
		 *
		 * param[in] nDOFs number of DOFs input.
		 * param[in,out] MassMatrix pointer to the matrix that is inverted.
		 */
    void ComputeMassMatrix(unsigned short nDOFs,
                           as3double     *MassMatrix);

    /*! 
		 * @brief Getter function which returns the value of TypeDOFs.
		 *
		 * @return TypeDOFs
		 */
    unsigned short GetTypeDOFs(void)                              const {return TypeDOFs;}
    /*!
		 * @brief Getter function which returns the value of nPolySol.
		 *
		 * @return nPolySol
		 */
    unsigned short GetnPolySol(void)                              const {return nPolySol;}
    /*!
		 * @brief Getter function which returns the value of nDOFsSol1D.
		 *
		 * @return nDOFsSol1D
		 */
    unsigned short GetnDOFsSol1D(void)              			 				const {return nDOFsSol1D;}
    /*!
		 * @brief Getter function which returns the value of nDOFsSol2D.
		 *
		 * @return nDOFsSol2D
		 */
    unsigned short GetnDOFsSol2D(void)              			 				const {return nDOFsSol2D;}
    /*!
		 * @brief Getter function which returns the value of nDOFsInt1D.
		 *
		 * @return nDOFsInt1D
		 */
    unsigned short GetnDOFsInt1D(void)              			 				const {return nDOFsInt1D;}
    /*!
		 * @brief Getter function which returns the value of nDOFsInt2D.
		 *
		 * @return nDOFsInt2D
		 */
    unsigned short GetnDOFsInt2D(void)              			 				const {return nDOFsInt2D;}
    /*!
		 * @brief Getter function which returns the value of rDOFsSol1D.
		 *
		 * @return rDOFsSol1D
		 */
    const as3vector1d<as3double> &GetrDOFsSol1D(void)             const {return rDOFsSol1D;}
    /*!
		 * @brief Getter function which returns the value of rDOFsInt1D.
		 *
		 * @return rDOFsInt1D
		 */
    const as3vector1d<as3double> &GetrDOFsInt1D(void)             const {return rDOFsInt1D;}
    /*!
		 * @brief Getter function which returns the value of wDOFsInt1D.
		 *
		 * @return wDOFsInt1D
		 */
    const as3vector1d<as3double> &GetwDOFsInt1D(void)             const {return wDOFsInt1D;}
    /*! Getter function which returns the value of wDOFsInt2D.
		 *
		 * @return wDOFsInt2D
		 */
    const as3vector1d<as3double> &GetwDOFsInt2D(void)             const {return wDOFsInt2D;}
    /*!
		 * @brief Getter function which returns the value of Vandermonde1D.
		 *
		 * @return Vandermonde1D
		 */
    const as3double *GetVandermonde1D(void)                       const {return Vandermonde1D;}
    /*!
		 * @brief Getter function which returns the value of InvVandermonde1D.
		 *
		 * @return InvVandermonde1D
		 */
    const as3double *GetInvVandermonde1D(void)                    const {return InvVandermonde1D;}
    /*!
		 * @brief Getter function which returns the value of lagrangeInt1D.
		 *
		 * @return lagrangeInt1D
		 */
    const as3double *GetLagrangeInt1D(void)      			          	const {return lagrangeInt1D;}
		/*!
		 * @brief Getter function which returns the value of lagrangeSol1D.
		 *
		 * @return lagrangeSol1D
		 */
		const as3double *GetLagrangeSol1D(void)                       const {return lagrangeSol1D;}
    /*!
		 * @brief Getter function which returns the value of derLagrangeInt1D.
		 *
		 * @return derLagrangeInt1D
		 */
    const as3double *GetDerLagrangeInt1D(void)   			            const {return derLagrangeInt1D;}
    /*!
		 * @brief Getter function which returns the value of lagrangeInt1DTranspose.
		 *
		 * @return lagrangeInt1DTranspose
		 */
		const as3double *GetLagrangeInt1DTranspose(void)      			 	const {return lagrangeInt1DTranspose;}
		/*!
		 * @brief Getter function which returns the value of derLagrangeInt1DTranspose.
		 *
		 * @return derLagrangeInt1DTranspose
		 */
		const as3double *GetDerLagrangeInt1DTranspose(void)   			  const {return derLagrangeInt1DTranspose;}
		/*!
		 * @brief Getter function which returns the value of derLagrangeDOFsSol1D.
		 *
		 * @return derLagrangeDOFsSol1D
		 */
		const as3double *GetDerLagrangeDOFsSol1D(void)                const {return derLagrangeDOFsSol1D;}
		/*!
		 * @brief Getter function which returns the value of derLagrangeDOFsSol1DTranspose.
		 *
		 * @return derLagrangeDOFsSol1DTranspose
		 */
		const as3double *GetDerLagrangeDOFsSol1DTranspose(void)       const {return derLagrangeDOFsSol1DTranspose;}
    /*!
		 * @brief Getter function which returns the value of derLagrangeDOFsSol1DMinFace.
		 *
		 * @return derLagrangeDOFsSol1DMinFace
		 */
		const as3double *GetDerLagrangeDOFsSol1DMinFace(void)       	const {return derLagrangeDOFsSol1DMinFace;}
		/*! 
		 * @brief Getter function which returns the value of derLagrangeDOFsSol1DMaxFace.
		 *
		 * @return derLagrangeDOFsSol1DMaxFace
		 */
		const as3double *GetDerLagrangeDOFsSol1DMaxFace(void)       	const {return derLagrangeDOFsSol1DMaxFace;}
    /*!
		 * @brief Getter function which returns the value of derLagrangeDOFsSol1DMinFace or derLagrangeDOFsSol1DMaxFace.
		 *
		 * @param[in] iFace input face ID.
		 *
		 * @return derLagrangeDOFsSol1DMinFace or derLagrangeDOFsSol1DMaxFace
		 */
		const as3double *GetDerLagrangeDOFsSol1DFace(unsigned short iFace) const {
			if( iFace == IDX_SOUTH || iFace == IDX_WEST )
				return derLagrangeDOFsSol1DMinFace;
			else
				return derLagrangeDOFsSol1DMaxFace;
		}
    /*!
		 * @brief Getter function which returns the value of IndexDOFsSol.
		 *
		 * @return IndexDOFsSol
		 */
    const as3vector2d<unsigned short> &GetIndexDOFsSol(void)      const {return IndexDOFsSol;}
    /*!
		 * @brief Getter function which returns the value of IndexDOFsSol, per input face.
		 *
		 * @param[in] iFace input face ID.
		 *
		 * @return IndexDOFsSol[iFace]
		 */
    const as3vector1d<unsigned short> &GetIndexDOFsSol(unsigned short iFace) const {
      return IndexDOFsSol[iFace];
    }

  protected:

  private:
		unsigned short zoneID;     ///< Current zone ID.
		unsigned short nPolySol;   ///< Current solution polynomial order.
		unsigned short nDOFsSol1D; ///< Number of solution nodes in 1D.
		unsigned short nDOFsSol2D; ///< Number of solution nodes in 2D.
		unsigned short nDOFsInt1D; ///< Number of integration nodes in 1D.
		unsigned short nDOFsInt2D; ///< Number of integration nodes in 2D.
    unsigned short TypeDOFs;   ///< Type of DOFs used.

		as3vector2d<unsigned short> IndexDOFsSol; ///< Indices of the solution nodes (SolDOFs) on all four sides of the element.
    as3vector1d<as3double>      rDOFsSol1D;   ///< Solution nodes on a reference line (1D).
    as3vector1d<as3double>      rDOFsInt1D;   ///< Integration nodes on a reference line (1D).
    as3vector1d<as3double>      wDOFsInt1D;   ///< Integration weights on a reference line (1D).
    as3vector1d<as3double>      wDOFsInt2D;   ///< Integration weights on a reference element (2D).

    as3double *Vandermonde1D                 = nullptr; ///< Vandermonde matrix in 1D that converts from the modal-to-nodal variables.
    as3double *InvVandermonde1D              = nullptr; ///< Inverse of the Vandermonde matrix in 1D that converts from the nodal-to-modal variables.
    as3double *lagrangeInt1D                 = nullptr; ///< Lagrange polynomials in 1D: interpolates from SolDOFs to IntDOFs.
		as3double *lagrangeSol1D                 = nullptr; ///< Lagrange polynomials in 1D: no interpolation, from SolDOFs to SolDOFs.
    as3double *derLagrangeInt1D              = nullptr; ///< Derivative of lagrange polynomials in 1D: interpolates from SolDOFs to IntDOFs.
    as3double *lagrangeInt1DTranspose        = nullptr; ///< Lagrange polynomials (transposed) in 1D: interpolates from SolDOFs to IntDOFs.
    as3double *derLagrangeInt1DTranspose     = nullptr; ///< Derivative of lagrange polynomials (transposed) in 1D: interpolates from SolDOFs to IntDOFs.
		as3double *derLagrangeDOFsSol1D          = nullptr; ///< Derivative of lagrange polynomials in 1D: no interplation: SolDOFs to SolDOFs.
		as3double *derLagrangeDOFsSol1DTranspose = nullptr; ///< Derivative of lagrange polynomials (transposed) in 1D: no interplation: SolDOFs to SolDOFs.
    as3double *derLagrangeDOFsSol1DMinFace   = nullptr; ///< Lagrange differentiation vector in 1D on min face: SolDOFs.
    as3double *derLagrangeDOFsSol1DMaxFace   = nullptr; ///< Lagrange differentiation vector in 1D on max face: SolDOFs.


    /*!
		 * @brief Function that initializes the quadrature integration data.
		 *
		 * @param[in,out] r1D reference to the integration points in 1D.
		 * @param[in,out] w1D reference to the integration weights in 1D.
		 * @param[in,out] w2D reference to the integration weights in 2D.
		 */
    void InitializeQuadratureIntegration(as3vector1d<as3double> &r1D,
                                         as3vector1d<as3double> &w1D,
                                         as3vector1d<as3double> &w2D);

    /*!
		 * @brief Function that determines the indices of each of the element faces based on the polynomial degree of this container and zone.
		 */
    void IdentifyFaceIndices(void);

    /*!
		 * @brief Function that computes the Lagrange interpolation entry in 1D for a given basis.
		 *
		 * @param[in] xPoint point of evaluation.
		 * @param[in] rBasis reference to input basis nodes.
		 * @param[in] iDegree degree of Lagrange polynomial.
		 *
		 * @return Lagrange polynomial of degree (iDegree) using basis (rBasis) at point (xPoint).
		 */
    as3double EvaluateEntryLagrangePolynomial(const as3double               xPoint,
                                              const as3vector1d<as3double> &rBasis,
                                              const unsigned short          iDegree);

    /*! Function that computes the derivative of lagrange interpolation entry in 1D for a given basis.
		 *
		 * @param[in] xPoint point of evaluation.
		 * @param[in] rBasis reference to input basis nodes.
		 * @param[in] iDegree degree of Lagrange polynomial.
		 * 
		 * @return Lagrange differential polynomial of degree (iDegree) using basis (rBasis) at point (xPoint).
		 */
    as3double EvaluateEntryDerLagrangePolynomial(const as3double               xPoint,
                                                 const as3vector1d<as3double> &rBasis,
                                                 const unsigned short          iDegree);

    /*!
		 * @brief Function that computes the 1D Vandermonde matrix.
		 *
		 * @param[in] nCol Number of columns in a 1D Vandermonde matrix.
		 * @param[in] r reference to the solution DOFs in 1D.
		 * @param[out] V pointer to the Vandermonde matrix in 1D.
		 */
    void ComputeVandermonde1D(const unsigned short          nCol,
                              const as3vector1d<as3double> &r,
                              as3double                    *V);

    /*!
		 * @brief Function that computes the generalized normalized-Jacobi polynomial for an order n and at a given entry x.
		 *
		 * @param[in] n order of the Jacobi polynomial.
		 * @param[in] alpha coefficient alpha in the Jacobi polynomial.
		 * @param[in] beta coefficient beta in the Jacobi polynomial.
		 * @param[in] x evaluation point.
		 *
		 * @return normalized Jacobi polynomial for an order (n) at the evaluation point (x) with coefficients alpha and beta.
		 */
    as3double JacobiNormalized(unsigned short n,
                               unsigned short alpha,
                               unsigned short beta,
                               as3double      x);
  };


