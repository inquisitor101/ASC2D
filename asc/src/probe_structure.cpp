#include "probe_structure.hpp"





CProbe::CProbe
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 CElement                   *element_container,
 as3vector1d<as3double>      probe,
 as3vector1d<unsigned short> vars,
 unsigned short              iZone
)
 /*
  * Constructor used to initialize CProbe.
  */
{
  // Assign probe coordinates.
  ProbeLocation = probe;

	// Initialize the total number of variables required by this probe. 
	for(unsigned short i=0; i<vars.size(); i++)
	{
		switch( vars[i] )
		{
			case( PROBE_DENSITY     ): variables.push_back( GetValueDensity     ); break;
			case( PROBE_XMOMENTUM   ): variables.push_back( GetValueXMomentum   ); break;
			case( PROBE_YMOMENTUM   ): variables.push_back( GetValueYMomentum   ); break;
			case( PROBE_TOTALENERGY ): variables.push_back( GetValueTotalEnergy ); break;
			case( PROBE_XVELOCITY   ): variables.push_back( GetValueXVelocity   ); break;
			case( PROBE_YVELOCITY   ): variables.push_back( GetValueYVelocity   ); break;
			case( PROBE_PRESSURE    ): variables.push_back( GetValuePressure    ); break;
			default:
			  Terminate("CProbe::CProbe", __FILE__, __LINE__,
						      "Unknown probe variable specified.");
		}
	}
	

  // Extract current zone geometry.
  auto* geometry_zone = geometry_container->GetGeometryZone(iZone);

  // Number of DOFS of solution points on an element in 1D in this zone.
  unsigned short nDOFsSol1D = element_container->GetnDOFsSol1D();
  // Number of DOFS of solution points on an element in 2D in this zone.
  unsigned short nDOFsSol2D = element_container->GetnDOFsSol2D();

  // Search for the element containing the probe.
  const bool found = geometry_zone->SearchElementProbe(ProbeLocation, iElemProbe, true);

  // For now, this is only executed in the main zone. If element not found,
  // then throw an error and terminate.
  if( !found && (iZone == ZONE_MAIN) ){
    std::string message = "Probe: ("
                        + std::to_string(ProbeLocation[0]) + ", "
                        + std::to_string(ProbeLocation[1]) + ") "
                        + "could not be located, or is lying at an edge/corner.";
    Terminate("CProbe::CProbe", __FILE__, __LINE__, message);
  }

  // Reserve memory for the interpolation polynomial.
  Interpolation.resize(nDOFsSol2D);

  // Extract element.
  auto* geometry_element = geometry_zone->GetGeometryElem(iElemProbe);

  // Extract solution DOFs of element.
  auto& CoordSolDOFs = geometry_element->GetCoordSolDOFs();

  // Temporary basis storage in x and y.
  as3vector2d<as3double> ellx(nDOFsSol1D);
  as3vector2d<as3double> elly(nDOFsSol1D);
  for(unsigned short i=0; i<nDOFsSol1D; i++){
    ellx[i].resize(nDOFsSol1D);
    elly[i].resize(nDOFsSol1D);
  }

  // Extract tensor-product basis for interpolation.
  unsigned short idx = 0;
  for(unsigned short j=0; j<nDOFsSol1D; j++){
    for(unsigned short i=0; i<nDOFsSol1D; i++){
      // Extract basis per one-dimensional tensor-product.
      ellx[j][i] = CoordSolDOFs[0][idx]; // basis in x.
      elly[i][j] = CoordSolDOFs[1][idx]; // basis in y.

      // Update nodal index.
      idx++;
    }
  }

  // Form Lagrange interpolation polynomial.
  idx = 0;
  for(unsigned short j=0; j<nDOFsSol1D; j++)
    for(unsigned short i=0; i<nDOFsSol1D; i++)
      Interpolation[idx++] = EvaluateEntryLagrangePolynomial(ProbeLocation[0], ellx[j], i)
                            *EvaluateEntryLagrangePolynomial(ProbeLocation[1], elly[i], j);

  // Make sure Lagrange interpolation is computed correctly.
  as3double tmp = 0.0;
  for(unsigned short i=0; i<Interpolation.size(); i++) tmp += Interpolation[i];
  if( fabs(tmp - 1.0) > 1.0e-10 )
    Terminate("CProbe::CProbe", __FILE__, __LINE__,
              "Error in lagrange interpolation polynomial computation.");
}


CProbe::~CProbe
(
 void
)
 /*
  * Destructor for CProbe class, frees allocated memory.
  */
{

}


as3double CProbe::EvaluateEntryLagrangePolynomial
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
		if( j != iDegree ) ell *= (xPoint - rBasis[j])/(rBasis[iDegree]-rBasis[j]);

	// Return lagrange entry.
	return ell;
}


as3vector1d<as3double> CProbe::SampleProbeData
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSolver   *solver_container,
 CInitial  *initial_container
)
 /*
  * Function that samples the data in this probe.
  */
{
  // Some abbreviations.
  const as3double gm1 = GAMMA_MINUS_ONE;

	// Extract number of solution points in 2D in this zone.
  unsigned short nDOFsSol2D = solver_container->GetnDOFsSol2D();

  // Extract the required element solution on this probe.
  auto& sol = solver_container->GetDataContainer(iElemProbe)->GetDataDOFsSol();

  // Initialize interpolated value.
  as3vector1d<as3double> val(nVar, 0.0);

  // Interpolate the solution of the sensor element to target probe location.
  for(unsigned short iVar=0; iVar<nVar; iVar++)
	{
		as3double tmp = 0.0;
  	// Loop over all solution points and interpolate to the current location.
#pragma omp simd reduction(+:tmp)
    for(unsigned short l=0; l<nDOFsSol2D; l++) tmp += Interpolation[l]*sol[iVar][l];
		// Assign the interpolated value to val[iVar].
		val[iVar] = tmp;
	}

	// Reserve memory for the output probe data.
	as3vector1d<as3double> data( variables.size(), 0.0 );

	// Loop over all the variables required in the probe and compute them.
	for(unsigned int i=0; i<variables.size(); i++)
	{
		data[i] = variables[i](val);
	}

	// Return probed data variables.
	return data;
}


