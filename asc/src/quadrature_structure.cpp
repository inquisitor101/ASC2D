#include "quadrature_structure.hpp"




CGaussJacobiQuadrature::CGaussJacobiQuadrature
(
 void
)
 /*
	* Constructor, initializes CGaussJacobiQuadrature class.
	*/
{
	// Default configuration.
	maxIter = 10000;
	Err     = 1.0;
	TOL     = 1.0e-10;

	// Root-finding resolution.
	dx = 0.01;

	// Polynomial curve samples.
	nSamples = 50000;
}


void CGaussJacobiQuadrature::InitializeGLQuadrature
(
 as3double  nPoints,
 as3double *rQuad,
 as3double *wQuad
)
 /*
	* Function that computes the quadrature points and weights in 1D.
	* So far, this only uses the special case of alpha=beta=0, meaning
	* this is a Gauss-Legendre quadrature.
	*/
{
	// Temporary vector.
	std::vector<as3double> r(nSamples);
	as3double dx = 2.0/(nSamples-1);
	r[0] = -1.0;
	for(unsigned int i=1; i<nSamples; i++)
		r[i] = r[i-1]+dx;
	assert( fabs(r[nSamples-1]-1.0) < TOL );

	std::vector<as3double> Pn1p(nSamples);
	for(unsigned int i=0; i<nSamples; i++)
		Pn1p[i] = PolyLegendre(r[i], nPoints);

	// Find indices at which Pn changes sign.
	unsigned short nSize = 0;
	std::vector<unsigned int> rootIDX;
	for(unsigned short iSample=0; iSample<nSamples-1; iSample++){
		if(Pn1p[iSample]*Pn1p[iSample+1] < 0){
			// Book-keep 1st (left-most) value.
			rootIDX.push_back(iSample);
			nSize++;
		}
	}

	// Starting guess for Newton's method.
	std::vector<as3double> x0(nSize);
	for(unsigned short i=0; i<x0.size(); i++)
		x0[i] = r[rootIDX[i]];

	// Compute Gauss-Legendre points.
	as3double x = 0.0, xOld = 0.0, temp = 0.0, dPndx = 0.0;
	for(unsigned short n=0; n<nPoints; n++){
		x    = x0[n];
		xOld = x-dx;
		for(unsigned int k=0; k<maxIter; k++){
			dPndx = (PolyLegendre(x, nPoints)-PolyLegendre(xOld, nPoints))/(x-xOld);
			temp  = x;
			x    -= PolyLegendre(x, nPoints)/dPndx;
			Err   = fabs(x-temp);
			xOld  = temp;

			if(Err <= TOL)
				break;

		}
		rQuad[n] = x;
	}

	// Compute Gauss-Legendre weights.
	as3double dPn = 0.0;
	for(unsigned short i=0; i<nPoints; i++){

		dPn = dPolyLegendre(rQuad[i], nPoints);

		// Compute weight.
		wQuad[i] = 2.0/( (1.0-rQuad[i]*rQuad[i])*dPn*dPn);
	}
}


as3double CGaussJacobiQuadrature::PolyLegendre
(
 as3double      r,
 unsigned short iDegree
)
 /*
	* Function that computes the legendre polynomial.
	*/
{
	// P_{n+1}.
	as3double Pn1p;

	// P_{n-1} starting values.
	as3double Pn1m = 1.0;
	as3double Pn   = r;

	// Special cases.
	if(iDegree == 0){
		Pn1p = Pn1m;
		return Pn1p;
	}
	if(iDegree == 1){
		Pn1p = Pn;
		return Pn1p;
	}

	// Compute Pn_{n+1} via recurrence relation.
	for(unsigned short n=1; n<iDegree; n++){
		Pn1p = ((2.0*n+1.0)*r*Pn - n*Pn1m)/(n+1.0);

		// Update vector values.
		Pn1m = Pn;
		Pn   = Pn1p;
	}

	return Pn1p;
}


as3double CGaussJacobiQuadrature::dPolyLegendre
(
 as3double 			r,
 unsigned short iDegree
)
 /*
	* Function that computes the derivative of the legendre polynomial.
	*/
{
	as3double dPn = 0.0;

	as3double Pn   = PolyLegendre(r, iDegree);
	as3double Pn1m = PolyLegendre(r, iDegree-1);

	// Compute derivative.
	dPn = iDegree*( r*Pn-Pn1m)/(r*r-1.0);

	return dPn;
}
