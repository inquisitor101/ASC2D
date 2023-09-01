#include "element_structure.hpp"




CElement::CElement
(
 CConfig       *config_container,
 unsigned short iZone
)
 /*
	* Constructor, initializes CElement class.
	*/
{
  // Integration factor in 1D, based on solution.
	unsigned short cAlias = config_container->GetIntegrationFactor();

	// Set zone ID.
	zoneID     = iZone;
	// Set solution polynomial order.
	nPolySol   = config_container->GetnPolySolZone(iZone);
	// Number of nodes in 1D, based on solution DOFs.
	nDOFsSol1D = nPolySol+1;
	// Number of integration points in 1D, used in Gauss-Legendre quadrature.
	nDOFsInt1D = cAlias*nPolySol/2 + 1;
  // Use multiple of vecLen1D.
  // nDOFsInt1D = ( (nDOFsInt1D+vecLen1D-1)/vecLen1D )*vecLen1D;

  // Number of solution points in 2D.
  nDOFsSol2D = nDOFsSol1D*nDOFsSol1D;
  // Number of integration points in 2D, used in Gauss-Legendre quadrature.
	nDOFsInt2D = nDOFsInt1D*nDOFsInt1D;

  // Get type of DOFs used in this zone.
  TypeDOFs   = config_container->GetTypeDOFs(iZone);

	// Initialize nodal solution DOFs on a reference line.
	rDOFsSol1D.resize(nDOFsSol1D);

	// Initialize integration nodes on a reference line.
	rDOFsInt1D.resize(nDOFsInt1D);
	// Initialize integration weights on a reference line.
	wDOFsInt1D.resize(nDOFsInt1D);
  // Initialize integration weights on entire element.
  wDOFsInt2D.resize(nDOFsInt2D);

  // Initialize lagrange interpolation matrix (transposed) in 1D: SolDOFs to IntDOFs.
	lagrangeInt1D    = new as3double[nDOFsInt1D*nDOFsSol1D]();
	// Initialize lagrange differentiation matrix (transposed) in 1D: SolDOFs to IntDOFs.
	derLagrangeInt1D = new as3double[nDOFsInt1D*nDOFsSol1D]();

	// Initialize lagrange interpolation matrix (transposed) in 1D: SolDOFs to IntDOFs.
	lagrangeInt1DTranspose    = new as3double[nDOFsSol1D*nDOFsInt1D]();
	// Initialize lagrange differentiation matrix (transposed) in 1D: SolDOFs to IntDOFs.
	derLagrangeInt1DTranspose = new as3double[nDOFsSol1D*nDOFsInt1D]();

	// Initialize lagrange differentiation vector in 1D on min face: SolDOFS.
	derLagrangeDOFsSol1DMinFace = new as3double[nDOFsSol1D]();
	// Initialize lagrange differentiation vector in 1D on max face: SolDOFS.
	derLagrangeDOFsSol1DMaxFace = new as3double[nDOFsSol1D]();


  // Compute location of nodal solution DOFs.
  LocationDOFs1D(TypeDOFs, rDOFsSol1D);

  // Compute integration DOFs and weights.
  InitializeQuadratureIntegration(rDOFsInt1D, wDOFsInt1D, wDOFsInt2D);

  // Compute lagrange polynomials in 1D that interpolate from SolDOFs to IntDOFs.
  LagrangeBasisFunctions(rDOFsSol1D, rDOFsInt1D,
                         lagrangeInt1D,
                         derLagrangeInt1D,
  											 lagrangeInt1DTranspose,
  											 derLagrangeInt1DTranspose,
  											 derLagrangeDOFsSol1DMinFace,
  											 derLagrangeDOFsSol1DMaxFace);

  // Determine the nodal indices per each element boundary/trace.
  IdentifyFaceIndices();

  // Reserve memory for the Vandermonde matrix in 1D: converts modal-to-nodal.
  Vandermonde1D = new as3double[nDOFsSol1D*nDOFsSol1D]();
  // Compute the Vandermonde1D matrix in 1D: converts modal-to-nodal.
  ComputeVandermonde1D(nDOFsSol1D, rDOFsSol1D, Vandermonde1D);

  // Reserve memory for the inverse of the Vandermonde matrix in 1D: converts nodal-to-modal.
  InvVandermonde1D = new as3double[nDOFsSol1D*nDOFsSol1D]();
  // Compute inverse of Vandermonde1D matrix.
  ComputeInverseMatrix(nDOFsSol1D, Vandermonde1D, InvVandermonde1D);

	// Reserve memory for the lagrange polynomial matrix in 1D: on SolDOFs only (it is an identity matrix).
	lagrangeSol1D = new as3double[nDOFsSol1D*nDOFsSol1D]();

	// Reserve memory for the lagrange derivative matrix in 1D: on SolDOFs only.
	derLagrangeDOFsSol1D = new as3double[nDOFsSol1D*nDOFsSol1D]();

	// Reserve memory for the lagrange derivative matrix (transposed) in 1D: on SolDOFs only.
	derLagrangeDOFsSol1DTranspose = new as3double[nDOFsSol1D*nDOFsSol1D]();

	// Compute lagrange polynomials in 1D that do not interpolate, just act on the SolDOFs.
	LagrangeBasisFunctions(rDOFsSol1D, rDOFsSol1D,
			                   nullptr, 
												 derLagrangeDOFsSol1D, lagrangeSol1D,
												 derLagrangeDOFsSol1DTranspose,
												 nullptr, nullptr);
}


CElement::~CElement
(
 void
)
 /*
	* Destructor for CElement class, frees allocated memory.
	*/
{
  if( Vandermonde1D                 != nullptr ) delete [] Vandermonde1D;
  if( InvVandermonde1D              != nullptr ) delete [] InvVandermonde1D;
  if( lagrangeInt1D                 != nullptr ) delete [] lagrangeInt1D;
	if( lagrangeSol1D                 != nullptr ) delete [] lagrangeSol1D;
  if( derLagrangeInt1D              != nullptr ) delete [] derLagrangeInt1D;
  if( derLagrangeDOFsSol1D          != nullptr ) delete [] derLagrangeDOFsSol1D;
	if( lagrangeInt1DTranspose        != nullptr ) delete [] lagrangeInt1DTranspose;
	if( derLagrangeInt1DTranspose     != nullptr ) delete [] derLagrangeInt1DTranspose;
	if( derLagrangeDOFsSol1DTranspose != nullptr ) delete [] derLagrangeDOFsSol1DTranspose;
	if( derLagrangeDOFsSol1DMinFace   != nullptr ) delete [] derLagrangeDOFsSol1DMinFace;
	if( derLagrangeDOFsSol1DMaxFace   != nullptr ) delete [] derLagrangeDOFsSol1DMaxFace;
}


void CElement::ComputeMassMatrix
(
  unsigned short nDOFs,
  as3double     *MassMatrix
)
 /*
  * Function that computes the mass matrix in this zone.
  */
{
  // Deduce number of rows, and make sure they match the zone data.
  if( nDOFs != nDOFsSol2D )
    Terminate("CElement::ComputeMassMatrix", __FILE__, __LINE__,
              "Dimension of mass matrix does not match zone data.");

  // Cast the transposed lagrange from 1D to 2D.
	const as3double (*ellT)[nDOFsInt1D] = (const as3double (*)[nDOFsInt1D]) lagrangeInt1DTranspose;
  // Cast the mass matrix from 1D to 2D.
  as3double (*M)[nDOFs] = (as3double (*)[nDOFs]) MassMatrix;

  // Initialize mass matrix to zero.
  for(unsigned short i=0; i<nDOFsSol2D; i++)
    for(unsigned short j=0; j<nDOFsSol2D; j++)
      M[i][j] = 0.0;

  // Compute mass matrix.
  unsigned short iRow = 0;
  // The below two loops are the equivalent to the rows of the mass matrix.
  for(unsigned short iRowJ=0; iRowJ<nDOFsSol1D; iRowJ++){
    for(unsigned short iRowI=0; iRowI<nDOFsSol1D; iRowI++){
      unsigned short iCol = 0;
      // The below two loops are the equivalent of the columns of the mass matrix.
      for(unsigned short iColJ=0; iColJ<nDOFsSol1D; iColJ++){
        for(unsigned short iColI=0; iColI<nDOFsSol1D; iColI++){
          unsigned short iQuad = 0;
          // The below two loops are the equivalent to integration in 2D.
          for(unsigned short iQuadJ=0; iQuadJ<nDOFsInt1D; iQuadJ++){
            for(unsigned short iQuadI=0; iQuadI<nDOFsInt1D; iQuadI++){

              // Test function.
              const as3double f1 = ellT[iRowI][iQuadI]*ellT[iRowJ][iQuadJ];
              // Trial function.
              const as3double f2 = ellT[iColI][iQuadI]*ellT[iColJ][iQuadJ];

              // Mass matrix entry.
              M[iRow][iCol] += wDOFsInt2D[iQuad]*f1*f2;

              // Update quadrature index.
              iQuad++;
            }
          }
          // Update mass matrix column index.
          iCol++;
        }
      }
      // Update mass matrix row index.
      iRow++;
    }
  }
}


void CElement::InitializeQuadratureIntegration
(
  as3vector1d<as3double> &r1D,
  as3vector1d<as3double> &w1D,
  as3vector1d<as3double> &w2D
)
 /*
  * Function that computes the integration DOFs and weights.
  */
{
  // Create a quadrature integration object.
  CGaussJacobiQuadrature quadrature_container;

  // Initialize and compute the 1D integration nodes and weights.
  quadrature_container.InitializeGLQuadrature(r1D.size(),
                                              r1D.data(),
                                              w1D.data());

  // Compute the two-dimensional tensor-product of the weights.
  unsigned short idx = 0;
  for(unsigned short j=0; j<w1D.size(); j++)
    for(unsigned short i=0; i<w1D.size(); i++)
      w2D[idx++] = w1D[i]*w1D[j];
}


void CElement::LagrangeBasisFunctions
(
 const as3vector1d<as3double> &rBasis,
 const as3vector1d<as3double> &rPoints,
 as3double                    *lagrangeDOFs1D,
 as3double                    *derLagrangeDOFs1D,
 as3double                    *lagrangeDOFs1DTranspose,
 as3double                    *derLagrangeDOFs1DTranspose,
 as3double                    *derLagrangeDOFsMinFace1D,
 as3double                    *derLagrangeDOFsMaxFace1D
)
 /*
	* Function that computes the lagrange interpolation and its derivatives
	* in 1D based on the rBasis basis, evaluated (or interpolated) on rPoints.
	* Note, the transpose is computed due to computational purposes.
	*/
{
	// Tolerance.
	as3double TOL = 1.0e-10;

  // Number of points for the basis polynomial.
  const unsigned short nBasis = rBasis.size();
  // Number of points for the interpolated points.
  const unsigned short nPoint = rPoints.size();

	// For ease of readability.
	const unsigned short nRow = nBasis;
	const unsigned short nCol = nPoint;


  // Check if transpose lagrange basis is needed.
  if( lagrangeDOFs1DTranspose ){

  	for(unsigned short iRow=0; iRow<nRow; iRow++){

  		// Compute lagrange entries.
  		for(unsigned short iCol=0; iCol<nCol; iCol++){

  			// Point of evaluation.
  			as3double x = rPoints[iCol];

  			// Assign index.
  			const unsigned short idx = iRow*nCol + iCol;

  			// Compute lagrange entry.
  			lagrangeDOFs1DTranspose[idx] = EvaluateEntryLagrangePolynomial(x, rBasis, iRow);
  		}
  	}

  	// Check if Lagrange computations is correct.
  	for(unsigned short iCol=0; iCol<nCol; iCol++){
  		as3double tmp = 0.0;
  		for(unsigned short iRow=0; iRow<nRow; iRow++){

  			// Assign index.
  			const unsigned short idx = iRow*nCol + iCol;
        // Accumulate values.
  			tmp += lagrangeDOFs1DTranspose[idx];
  		}

  		if( fabs(tmp - 1.0) > TOL )
  			Terminate("CElement::LagrangeBasisFunctions", __FILE__, __LINE__,
  								"Error in lagrange polynomial computation.");
  	}
  }


  // Check if lagrange basis is needed.
  if( lagrangeDOFs1D ){

    // Make sure the lagrange computations are actually computed.
    if( !lagrangeDOFs1DTranspose )
      Terminate("CElement::LagrangeBasisFunctions", __FILE__, __LINE__,
                "Transpose of Lagrange polynomial must be computed first");

    for(unsigned short iRow=0; iRow<nRow; iRow++){
      for(unsigned short iCol=0; iCol<nCol; iCol++){

        // Assign indices for transpose and non-transpose.
        const unsigned short idxT = iRow*nCol + iCol;
        const unsigned short idx  = iCol*nRow + iRow;

        // Compute the transpose.
        lagrangeDOFs1D[idx] = lagrangeDOFs1DTranspose[idxT];
      }
    }
  }


  // Check if transpose lagrange differentiation is needed.
  if( derLagrangeDOFs1DTranspose ){

    for(unsigned short iRow=0; iRow<nRow; iRow++){

  		// Compute lagrange entries.
  		for(unsigned short iCol=0; iCol<nCol; iCol++){

  			// Point of evaluation.
  			as3double x = rPoints[iCol];

  			// Assign index.
  			const unsigned short idx = iRow*nCol + iCol;

  			// Compute lagrange differentiation entry.
  			derLagrangeDOFs1DTranspose[idx] = EvaluateEntryDerLagrangePolynomial(x, rBasis, iRow);
  		}
  	}

  	// Check if Lagrange computations is correct.
  	for(unsigned short iCol=0; iCol<nCol; iCol++){
  		as3double tmp = 0.0;
  		for(unsigned short iRow=0; iRow<nRow; iRow++){

  			// Assign index.
  			const unsigned short idx = iRow*nCol + iCol;
        // Accumulate values.
  			tmp += derLagrangeDOFs1DTranspose[idx];
  		}

  		if( fabs(tmp) > TOL )
  			Terminate("CElement::LagrangeBasisFunctions", __FILE__, __LINE__,
  								"Error in derivative of lagrange polynomial computation.");
    }
  }


  // Check if lagrange differentiation is needed.
  if( derLagrangeDOFs1D ){

    // Make sure the lagrange computations are actually computed.
    if( !derLagrangeDOFs1DTranspose )
      Terminate("CElement::LagrangeBasisFunctions", __FILE__, __LINE__,
                "Transpose of Lagrange differentiation must be computed first");

    for(unsigned short iRow=0; iRow<nRow; iRow++){
  		for(unsigned short iCol=0; iCol<nCol; iCol++){

        // Assign indices for transpose and non-transpose.
        const unsigned short idxT = iRow*nCol + iCol;
        const unsigned short idx  = iCol*nRow + iRow;

  			// Compute lagrange differentiation entry.
				derLagrangeDOFs1D[idx] = derLagrangeDOFs1DTranspose[idxT];
  		}
  	}
  }


	// Check if min/max face lagrange differentiation is needed.
	if( derLagrangeDOFsMinFace1D && derLagrangeDOFsMaxFace1D ){

		// Evaluation points at min and max face.
		as3double rMin = -1.0;
		as3double rMax =  1.0;

		// Evaluate lagrange differentiation on extremities of element (min and max face).
		for(unsigned short iRow=0; iRow<nRow; iRow++){
			derLagrangeDOFsMinFace1D[iRow] = EvaluateEntryDerLagrangePolynomial(rMin, rBasis, iRow);
			derLagrangeDOFsMaxFace1D[iRow] = EvaluateEntryDerLagrangePolynomial(rMax, rBasis, iRow);
		}

		// Check if Lagrange face computations is correct.
		as3double tmp_min = 0.0, tmp_max = 0.0;
		for(unsigned short iRow=0; iRow<nRow; iRow++){
			tmp_min += derLagrangeDOFsMinFace1D[iRow];
			tmp_max += derLagrangeDOFsMaxFace1D[iRow];
		}
		if( fabs(tmp_min) > TOL || fabs(tmp_max) > TOL )
			Terminate("CElement::LagrangeBasisFunctions", __FILE__, __LINE__,
								"Error in derivative of lagrange polynomial computation on faces");

		// Consistency check in between the anti-symmetric nature of the face lagrange differentiation.
		for(unsigned short iRow=0; iRow<nRow; iRow++){
			if( fabs(derLagrangeDOFsMinFace1D[iRow] + derLagrangeDOFsMaxFace1D[nRow-1-iRow]) > TOL )
				Terminate("CElement::LagrangeBasisFunctions", __FILE__, __LINE__,
									"Error in min/max face lagrange computation, anti-symmetry not preserved.");
		}
	}
}


as3double CElement::EvaluateEntryLagrangePolynomial
(
 const as3double               xPoint,
 const as3vector1d<as3double> &rBasis,
 const unsigned short          iDegree
)
 /*
	* Function that evaluates a single entry (xPoint) on a
	* given basis (rBasis) using the lagrange polynomial at
	* degree (iDegree).
	*/
{
  // Number of basis points.
  unsigned short nBasis = rBasis.size();
	// Temporary value to store results in.
	as3double ell = 1.0;
	for(unsigned short j=0; j<nBasis; j++)
		if( j != iDegree )
			ell *= (xPoint - rBasis[j])/(rBasis[iDegree]-rBasis[j]);

	// Return lagrange entry.
	return ell;
}


as3double CElement::EvaluateEntryDerLagrangePolynomial
(
 const as3double               xPoint,
 const as3vector1d<as3double> &rBasis,
 const unsigned short          iDegree
)
 /*
	* Function that evaluates a single entry (xPoint) on a
	* given basis (rBasis) using the derivative of the lagrange
	* polynomial at degree (iDegree).
	*/
{
  // Number of basis points.
  unsigned short nBasis = rBasis.size();
  // Reset value.
  as3double Lprime = 0.0;
  // Compute Lagrange derivative.
  for(unsigned short l=0; l<nBasis; l++){
    if(l != iDegree){
      as3double tmp = 1.0/( rBasis[iDegree] - rBasis[l] );
      for(unsigned short m=0; m<nBasis; m++){
        if(m != iDegree && m != l)
          tmp *= ( xPoint - rBasis[m] )/( rBasis[iDegree] - rBasis[m] );
      }
      Lprime += tmp;
    }
  }

	// Return derivative of lagrange entry.
  return Lprime;
}


void CElement::IdentifyFaceIndices
(
 void
)
 /*
	* Function that determines the indices of the element faces,
	* based on the degree of the polynomial used in this container.
	*/
{
  // Reserve memory for indices.
  IndexDOFsSol.resize(nFace);
  for(unsigned short iFace=0; iFace<nFace; iFace++)
    IndexDOFsSol[iFace].resize(nDOFsSol1D);

	// Offset variable.
	unsigned short offset;

	// South face.
	for(unsigned short i=0; i<nDOFsSol1D; i++)
		IndexDOFsSol[IDX_SOUTH][i] = i;

	// North face.
	offset = nDOFsSol1D*(nDOFsSol1D-1);
	for(unsigned short i=0; i<nDOFsSol1D; i++)
		IndexDOFsSol[IDX_NORTH][i] = IndexDOFsSol[IDX_SOUTH][i] + offset;

	// West face.
	for(unsigned short i=0; i<nDOFsSol1D; i++)
		IndexDOFsSol[IDX_WEST][i] = nDOFsSol1D*i;

	// East face.
	offset = nDOFsSol1D-1;
	for(unsigned short i=0; i<nDOFsSol1D; i++)
		IndexDOFsSol[IDX_EAST][i] = IndexDOFsSol[IDX_WEST][i] + offset;
}


void CElement::ComputeVandermonde1D
(
  const unsigned short          nCol,
  const as3vector1d<as3double> &r,
  as3double                    *V
)
 /*
  * Function that computes the Vandermonde matrix in 1D.
  */
{
  // Determine the number or rows of the Vandermonde matrix and check
  // if the dimension of V is correct.
  const unsigned short nRow = r.size();

  // Cast the 1D array into 2D for convenience.
  as3double (*v)[nCol] = (as3double (*)[nCol]) V;

  for(unsigned short i=0; i<nRow; i++)
    for(unsigned short j=0; j<nCol; j++)
      v[i][j] = JacobiNormalized(j, 0, 0, r[i]);
}


as3double CElement::JacobiNormalized
(
 unsigned short n,
 unsigned short alpha,
 unsigned short beta,
 as3double      x
)
 /*
  * Function that computes the normalized Jacobi polynomials for an order n
  * at the evaluation point x and with coefficients alpha and beta.
  */
{
  // Some abbreviations.
  const as3double ap1   = alpha + 1;
  const as3double bp1   = beta  + 1;
  const as3double apb   = alpha + beta;
  const as3double apbp1 = apb + 1;
  const as3double apbp2 = apb + 2;
  const as3double apbp3 = apb + 3;
  const as3double b2ma2 = beta*beta - alpha*alpha;

  // Determine the terms, which involves the gamma function. As the
  // arguments are integers, this term can be computed easily, because
  // Gamma(n+1) = n!.
  as3double Gamap1 = 1.0, Gambp1 = 1.0, Gamapbp2  = 1.0;
  for(int i=2; i<=alpha;          ++i)  Gamap1   *= i;
  for(int i=2; i<=beta;           ++i)  Gambp1   *= i;
  for(int i=2; i<=(alpha+beta+1); ++i)  Gamapbp2 *= i;

  // Initialize the normalized polynomials.
  as3double Pnm1 = sqrt(pow(0.5,apbp1)*Gamapbp2/(Gamap1*Gambp1));
  as3double Pn   = 0.5*Pnm1*(apbp2*x + alpha - beta)*sqrt(apbp3/(ap1*bp1));

  // Take care of the special situation of n == 0.
  if(n == 0) Pn = Pnm1;
  else {
    // The value of the normalized Jacobi polynomial must be obtained
    // via recursion.
    for(int i=2; i<=n; ++i){
      // Compute the coefficients a for i and i-1 and the coefficient bi.
      const int j          = i-1;
      as3double tmp        = 2*j + apb;
      const as3double aim1 = 2.0*sqrt(j*(j+apb)*(j+alpha)*(j+beta)/((tmp-1.0)*(tmp+1.0)))/tmp;

      const as3double bi   = b2ma2/(tmp*(tmp+2.0));

      tmp = 2*i + apb;
      const as3double ai   = 2.0*sqrt( i*(i+apb)*(i+alpha)*(i+beta)/((tmp-1.0)*(tmp+1.0)) )/tmp;

      // Compute the new value of Pn and make sure to store Pnm1 correctly.
      tmp  = Pnm1;
      Pnm1 = Pn;

      Pn = ((x-bi)*Pn - aim1*tmp)/ai;
    }
  }

  // Return Pn.
  return Pn;
}


void CElement::LocationDOFs1D
(
 unsigned short          TypeDOFs,
 as3vector1d<as3double> &rDOFs1D
)
 /*
	* Function that computes the locations of the DOFs in 1D.
	* This is based in computational space: [-1, 1].
	*/
{
	// Relative tolerance.
	as3double TOL = 1.0e-10;

  // Size of basis data.
  unsigned short nDOFs1D = rDOFs1D.size();

	// Determine what type of distribution the DOFs obey.
	switch( TypeDOFs ){

		// Equidistant distribution.
		case(TYPE_DOF_EQD):
		{
			// Step size.
			const as3double dh = 2.0/(nDOFs1D-1);

			// Starting position.
			const as3double r0 = -1.0;

			for(unsigned short i=0; i<nDOFs1D; i++)
				rDOFs1D[i] = r0 + i*dh;

			// Hard-enforce last node.
			rDOFs1D[nDOFs1D-1] = 1.0;

			break;
		}

		// Legendre-Gauss-Lobatto distribution.
		case(TYPE_DOF_LGL):
		{
			switch( nDOFs1D ){
				// Use already precomputed values.
				case 1:
				{
					rDOFs1D[0] = 0.0;
					break;
				}

				case 2:
				{
					rDOFs1D[0] = -1.0; rDOFs1D[1] = 1.0;
					break;
				}

				case 3:
				{
					rDOFs1D[0] = -1.0; rDOFs1D[1] = 0.0; rDOFs1D[2] = 1.0;
					break;
				}

				case 4:
				{
					const as3double t0 = sqrt(1.0/5.0);
					rDOFs1D[0] = -1.0;
					rDOFs1D[1] = -t0; rDOFs1D[2] =  t0;
					rDOFs1D[3] =  1.0;
					break;
				}

				case 5:
				{
					const as3double t0 = sqrt(3.0/7.0);
					rDOFs1D[0] = -1.0;
					rDOFs1D[1] = -t0; rDOFs1D[2] = 0.0; rDOFs1D[3] = t0;
					rDOFs1D[4] =  1.0;
					break;
				}

				case 6:
				{
					const as3double t1 = 2.0*sqrt(7.0)/21.0;
					const as3double t2 = sqrt(1.0/3.0 + t1);
					const as3double t3 = sqrt(1.0/3.0 - t1);
					rDOFs1D[0] = -1.0;
					rDOFs1D[1] = -t2; rDOFs1D[2] = -t3;
					rDOFs1D[3] =  t3; rDOFs1D[4] =  t2;
					rDOFs1D[5] =  1.0;
					break;
				}

				case 7:
				{
					const as3double t1 = 2.0*sqrt(5.0/3.0)/11.0;
					const as3double t2 = sqrt(5.0/11.0 + t1);
					const as3double t3 = sqrt(5.0/11.0 - t1);
					rDOFs1D[0] = -1.0;
					rDOFs1D[1] = -t2; rDOFs1D[2] = -t3;
					rDOFs1D[3] =  0.0;
					rDOFs1D[4] =  t3; rDOFs1D[5] =  t2;
					rDOFs1D[6] =  1.0;

					break;
				}

				case 8:
				{
					rDOFs1D[0] = -1.0;
          rDOFs1D[1] = -0.8717401485096066153375;
          rDOFs1D[2] = -0.5917001814331423021445;
          rDOFs1D[3] = -0.2092992179024788687687;
          rDOFs1D[4] =  0.2092992179024788687687;
          rDOFs1D[5] =  0.5917001814331423021445;
          rDOFs1D[6] =  0.8717401485096066153375;
          rDOFs1D[7] =  1.0;

					break;
				}

				case 9:
				{
          rDOFs1D[0] = -1.0;
          rDOFs1D[1] = -0.8997579954114601573124;
          rDOFs1D[2] = -0.6771862795107377534459;
          rDOFs1D[3] = -0.3631174638261781587108;
          rDOFs1D[4] =  0.0;
          rDOFs1D[5] =  0.3631174638261781587108;
          rDOFs1D[6] =  0.6771862795107377534459;
          rDOFs1D[7] =  0.8997579954114601573124;
          rDOFs1D[8] =  1.0;

					break;
				}

				case 10:
				{
          rDOFs1D[0] = -1.0;
          rDOFs1D[1] = -0.9195339081664588138289;
          rDOFs1D[2] = -0.7387738651055050750031;
          rDOFs1D[3] = -0.4779249498104444956612;
          rDOFs1D[4] = -0.1652789576663870246262;
          rDOFs1D[5] =  0.1652789576663870246262;
          rDOFs1D[6] =  0.4779249498104444956612;
          rDOFs1D[7] =  0.7387738651055050750031;
          rDOFs1D[8] =  0.9195339081664588138289;
          rDOFs1D[9] =  1.0;

					break;
				}

				case 11:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.93400143040805913433227413609938;
					rDOFs1D[2]  = -0.78448347366314441862241781610846;
					rDOFs1D[3]  = -0.56523532699620500647096396947775;
 					rDOFs1D[4]  = -0.29575813558693939143191151555906;
          rDOFs1D[5]  =  0.0;
					rDOFs1D[6]  =  0.29575813558693939143191151555906;
					rDOFs1D[7]  =  0.56523532699620500647096396947775;
					rDOFs1D[8]  =  0.78448347366314441862241781610846;
					rDOFs1D[9]  =  0.93400143040805913433227413609938;
					rDOFs1D[10] =  1.0;

					break;
				}

				case 12:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.94489927222288222340758013830322;
					rDOFs1D[2]  = -0.81927932164400667834864158171690;
					rDOFs1D[3]  = -0.63287615303186067766240485444366;
					rDOFs1D[4]  = -0.39953094096534893226434979156697;
					rDOFs1D[5]  = -0.13655293285492755486406185573969;
					rDOFs1D[6]  =  0.13655293285492755486406185573969;
					rDOFs1D[7]  =  0.39953094096534893226434979156697;
					rDOFs1D[8]  =  0.63287615303186067766240485444366;
					rDOFs1D[9]  =  0.81927932164400667834864158171690;
					rDOFs1D[10] =  0.94489927222288222340758013830322;
					rDOFs1D[11] =  1.0;

					break;
				}

				case 13:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.95330984664216391189690546475545;
					rDOFs1D[2]  = -0.84634756465187231686592560709875;
					rDOFs1D[3]  = -0.68618846908175742607275903956636;
					rDOFs1D[4]  = -0.48290982109133620174693723363693;
					rDOFs1D[5]  = -0.24928693010623999256867370037423;
					rDOFs1D[6]  =  0.0;
					rDOFs1D[7]  =  0.24928693010623999256867370037423;
					rDOFs1D[8]  =  0.48290982109133620174693723363693;
					rDOFs1D[9]  =  0.68618846908175742607275903956636;
					rDOFs1D[10] =  0.84634756465187231686592560709875;
					rDOFs1D[11] =  0.95330984664216391189690546475545;
					rDOFs1D[12] =  1.0;

					break;
				}

				case 14:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.95993504526726090135510016201542;
					rDOFs1D[2]  = -0.86780105383034725100022020290826;
					rDOFs1D[3]  = -0.72886859909132614058467240052088;
					rDOFs1D[4]  = -0.55063940292864705531662270585908;
					rDOFs1D[5]  = -0.34272401334271284504390340364167;
					rDOFs1D[6]  = -0.11633186888370386765877670973616;
					rDOFs1D[7]  =  0.11633186888370386765877670973616;
					rDOFs1D[8]  =  0.34272401334271284504390340364167;
					rDOFs1D[9]  =  0.55063940292864705531662270585908;
					rDOFs1D[10] =  0.72886859909132614058467240052088;
					rDOFs1D[11] =  0.86780105383034725100022020290826;
					rDOFs1D[12] =  0.95993504526726090135510016201542;
					rDOFs1D[13] =  1.0;

					break;
				}

				case 15:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.96524592650383857279585139206960;
					rDOFs1D[2]  = -0.88508204422297629882540163148223;
					rDOFs1D[3]  = -0.76351968995181520070411847597629;
					rDOFs1D[4]  = -0.60625320546984571112352993863673;
					rDOFs1D[5]  = -0.42063805471367248092189693873858;
					rDOFs1D[6]  = -0.21535395536379423822567944627292;
					rDOFs1D[7]  =  0.0;
					rDOFs1D[8]  =  0.21535395536379423822567944627292;
					rDOFs1D[9]  =  0.42063805471367248092189693873858;
					rDOFs1D[10] =  0.60625320546984571112352993863673;
					rDOFs1D[11] =  0.76351968995181520070411847597629;
					rDOFs1D[12] =  0.88508204422297629882540163148223;
					rDOFs1D[13] =  0.96524592650383857279585139206960;
					rDOFs1D[14] =  1.0;

					break;
				}

				case 16:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.96956804627021793295224273836746;
					rDOFs1D[2]  = -0.89920053309347209299462826151985;
					rDOFs1D[3]  = -0.79200829186181506393108827096315;
					rDOFs1D[4]  = -0.65238870288249308946788321964058;
					rDOFs1D[5]  = -0.48605942188713761178189078584687;
					rDOFs1D[6]  = -0.29983046890076320809835345472230;
					rDOFs1D[7]  = -0.10132627352194944784303300504592;
					rDOFs1D[8]  =  0.10132627352194944784303300504592;
					rDOFs1D[9]  =  0.29983046890076320809835345472230;
					rDOFs1D[10] =  0.48605942188713761178189078584687;
					rDOFs1D[11] =  0.65238870288249308946788321964058;
					rDOFs1D[12] =  0.79200829186181506393108827096315;
					rDOFs1D[13] =  0.89920053309347209299462826151985;
					rDOFs1D[14] =  0.96956804627021793295224273836746;
					rDOFs1D[15] =  1.0;

					break;
				}

				case 17:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.97313217663141831415697950187372;
					rDOFs1D[2]  = -0.91087999591557359562380250639773;
					rDOFs1D[3]  = -0.81569625122177030710675055323753;
					rDOFs1D[4]  = -0.69102898062768470539491935737245;
					rDOFs1D[5]  = -0.54138539933010153912373340750406;
					rDOFs1D[6]  = -0.37217443356547704190723468073526;
					rDOFs1D[7]  = -0.18951197351831738830426301475311;
					rDOFs1D[8]  =  0.0;
					rDOFs1D[9]  =  0.18951197351831738830426301475311;
					rDOFs1D[10] =  0.37217443356547704190723468073526;
					rDOFs1D[11] =  0.54138539933010153912373340750406;
					rDOFs1D[12] =  0.69102898062768470539491935737245;
					rDOFs1D[13] =  0.81569625122177030710675055323753;
					rDOFs1D[14] =  0.91087999591557359562380250639773;
					rDOFs1D[15] =  0.97313217663141831415697950187372;
					rDOFs1D[16] =  1.0;

					break;
				}

				case 18:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.976105557412198542864518924341700;
					rDOFs1D[2]  = -0.920649185347533873837854625431280;
					rDOFs1D[3]  = -0.835593535218090213713646362327940;
					rDOFs1D[4]  = -0.723679329283242681306210365302070;
					rDOFs1D[5]  = -0.588504834318661761173535893193560;
					rDOFs1D[6]  = -0.434415036912123975342287136740670;
					rDOFs1D[7]  = -0.266362652878280984167665332025600;
					rDOFs1D[8]  = -0.089749093484652111022645010088562;
					rDOFs1D[9]  =  0.089749093484652111022645010088562;
					rDOFs1D[10] =  0.266362652878280984167665332025600;
					rDOFs1D[11] =  0.434415036912123975342287136740670;
					rDOFs1D[12] =  0.588504834318661761173535893193560;
					rDOFs1D[13] =  0.723679329283242681306210365302070;
					rDOFs1D[14] =  0.835593535218090213713646362327940;
					rDOFs1D[15] =  0.920649185347533873837854625431280;
					rDOFs1D[16] =  0.976105557412198542864518924341700;
					rDOFs1D[17] =  1.0;

					break;
				}

				case 19:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.97861176622208009515263406311022;
					rDOFs1D[2]  = -0.92890152815258624371794025879655;
					rDOFs1D[3]  = -0.85246057779664609308595597004106;
					rDOFs1D[4]  = -0.75149420255261301416363748963394;
					rDOFs1D[5]  = -0.62890813726522049776683230622873;
					rDOFs1D[6]  = -0.48822928568071350277790963762492;
					rDOFs1D[7]  = -0.33350484782449861029850010384493;
					rDOFs1D[8]  = -0.16918602340928157137515415344488;
					rDOFs1D[9]  =  0.0;
					rDOFs1D[10] =  0.16918602340928157137515415344488;
					rDOFs1D[11] =  0.33350484782449861029850010384493;
					rDOFs1D[12] =  0.48822928568071350277790963762492;
					rDOFs1D[13] =  0.62890813726522049776683230622873;
					rDOFs1D[14] =  0.75149420255261301416363748963394;
					rDOFs1D[15] =  0.85246057779664609308595597004106;
					rDOFs1D[16] =  0.92890152815258624371794025879655;
					rDOFs1D[17] =  0.97861176622208009515263406311022;
					rDOFs1D[18] =  1.0;

					break;
				}

				case 20:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.980743704893914171925446438584230;
					rDOFs1D[2]  = -0.935934498812665435716181584930630;
					rDOFs1D[3]  = -0.866877978089950141309847214616290;
					rDOFs1D[4]  = -0.775368260952055870414317527594690;
					rDOFs1D[5]  = -0.663776402290311289846403322971160;
					rDOFs1D[6]  = -0.534992864031886261648135961828980;
					rDOFs1D[7]  = -0.392353183713909299386474703815820;
					rDOFs1D[8]  = -0.239551705922986495182401356927090;
					rDOFs1D[9]  = -0.080545937238821837975944518159554;
					rDOFs1D[10] =  0.080545937238821837975944518159554;
					rDOFs1D[11] =  0.239551705922986495182401356927090;
					rDOFs1D[12] =  0.392353183713909299386474703815820;
					rDOFs1D[13] =  0.534992864031886261648135961828980;
					rDOFs1D[14] =  0.663776402290311289846403322971160;
					rDOFs1D[15] =  0.775368260952055870414317527594690;
					rDOFs1D[16] =  0.866877978089950141309847214616290;
					rDOFs1D[17] =  0.935934498812665435716181584930630;
					rDOFs1D[18] =  0.980743704893914171925446438584230;
					rDOFs1D[19] =  1.0;

					break;
				}

				case 21:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.98257229660454802823448127655541;
					rDOFs1D[2]  = -0.94197629695974553429610265066144;
					rDOFs1D[3]  = -0.87929475532359046445115359630494;
					rDOFs1D[4]  = -0.79600192607771240474431258966036;
					rDOFs1D[5]  = -0.69405102606222323262731639319467;
					rDOFs1D[6]  = -0.57583196026183068692702187033809;
					rDOFs1D[7]  = -0.44411578327900210119451634960735;
					rDOFs1D[8]  = -0.30198985650876488727535186785875;
					rDOFs1D[9]  = -0.15278551580218546600635832848567;
					rDOFs1D[10] =  0.0;
					rDOFs1D[11] =  0.15278551580218546600635832848567;
					rDOFs1D[12] =  0.30198985650876488727535186785875;
					rDOFs1D[13] =  0.44411578327900210119451634960735;
					rDOFs1D[14] =  0.57583196026183068692702187033809;
					rDOFs1D[15] =  0.69405102606222323262731639319467;
					rDOFs1D[16] =  0.79600192607771240474431258966036;
					rDOFs1D[17] =  0.87929475532359046445115359630494;
					rDOFs1D[18] =  0.94197629695974553429610265066144;
					rDOFs1D[19] =  0.98257229660454802823448127655541;
					rDOFs1D[20] =  1.0;

					break;
				}

				default:
				Terminate("CElement::LocationDOFs1D", __FILE__, __LINE__,
									"LGL distribution for given polynomial order not (yet) implemented!");
			}

			break;
		}

		default:
			Terminate("CElement::LocationDOFs1D", __FILE__, __LINE__,
								"Unknown type of DOF distribution!");
	}

	// Make sure DOFs are good.
	as3double tmp = 0.0;
	for(unsigned short i=0; i<nDOFs1D; i++)
		tmp += rDOFs1D[i];
	assert( fabs(tmp) < TOL );
}

