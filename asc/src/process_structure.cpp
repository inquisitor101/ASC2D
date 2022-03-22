#include "process_structure.hpp"



CProcess::CProcess
(
  CConfig       *config_container,
  CGeometry     *geometry_container,
  COutput       *output_container,
  CElement     **element_container,
  CInitial     **initial_container,
  CSolver      **solver_container,
  CSpatial     **spatial_container,
  unsigned short iZone
)
 /*
  * Constructor, used to initialize CProcess in zone: iZone.
  */
{
  // Extract zone ID.
  zoneID = iZone;
}


CProcess::~CProcess
(
  void
)
 /*
  * Destructor for CProcess class, frees allocated memory.
  */
{

}


void CProcess::FilterSolution
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSolver   *solver_container,
 CSpatial  *spatial_container
)
 /*
  * Function that filters the solution.
  */
{
  // Extract total number of elements present in this zone.
  unsigned long nElem  = solver_container->GetnElem();

  // Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

  // Loop over all elements.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned long iElem=0; iElem<nElem; iElem++){

    // Compute filtered solution.
    spatial_container->ComputeFilteredSolution(config_container,
                                               geometry_container,
                                               element_container,
                                               data_container[iElem]);
  } // End of OpenMP parallel region.
}


CEEProcess::CEEProcess
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 COutput       *output_container,
 CElement     **element_container,
 CInitial     **initial_container,
 CSolver      **solver_container,
 CSpatial     **spatial_container,
 unsigned short iZone
)
  :
    CProcess
    (
      config_container,
      geometry_container,
      output_container,
      element_container,
      initial_container,
      solver_container,
      spatial_container,
      iZone
    )
 /*
  * Constructor, used to initialize CEEProcess in zone: iZone.
  */
{
  // Restrict processing to only main zone.
  if(iZone != ZONE_MAIN) return;

  // Extract the maximum number of temporal iterations.
  unsigned long MaxIter = config_container->GetMaxIter();

  // Reserve data for the size of the temporal iterations.
  TimeProcessed.reserve(MaxIter);

  // Obtain number of items to process.
  unsigned short nItem = 0;
  if( config_container->GetTypeProcessData() == PROCESS_IC ){
    // Determine which type of IC processing to use.
    switch( config_container->GetTypeIC(iZone) ){

      // Gaussian pressure pulse, entropy wave or isentropic vortex postprocessing.
      case(IC_GAUSSIAN_PRESSURE):
      case(IC_ENTROPY_WAVE):
      case(IC_ISENTROPIC_VORTEX):
      {
        // One item to process.
        nItem = 1;
        break;
      }

      // Vortex-rollup postprocessing.
      case(IC_VORTEX_ROLLUP):
      {
        // Obtain all probes, if specified.
        auto ProbeLocation = config_container->GetProbeLocation();

        // Consistency check.
        assert( ProbeLocation.size()%2 == 0 );

        // Total number of probes, division by nDim since these are pair of coordinates.
        unsigned short nProbe = ProbeLocation.size()/nDim;
        // Use 2 items to process (pressure and v-velocity), per probe.
        nItem = 2*nProbe;

        // Reserve memory for the sensors.
        sensor_container.resize(nProbe, nullptr);

        // Loop and initialize sensor for each probe.
        for(unsigned short iSensor=0; iSensor<nProbe; iSensor++){
          // Initialize current probe coordinates.
          as3vector1d<as3double> probe = { ProbeLocation[iSensor*nDim+0], ProbeLocation[iSensor*nDim+1] };
          // Create current sensor.
          sensor_container[iSensor] = new CSensorData(config_container,
                                                      geometry_container,
                                                      element_container[iZone],
                                                      probe, iZone);
        }

        break;
      }

      // Default, do nothing.
      default: break;
    }

    // Reserve data approximately the size of the temporal iterations.
    DataProcessed.resize(nItem);
    for(unsigned short i=0; i<nItem; i++) DataProcessed[i].reserve(MaxIter);
  }
}


CEEProcess::~CEEProcess
(
  void
)
 /*
  * Destructor for CEEProcess class, frees allocated memory.
  */
{
  for(unsigned short i=0; i<sensor_container.size(); i++)
    if( sensor_container[i] ) delete sensor_container[i];
}


void CEEProcess::WriteProcessedData
(
  CConfig   *config_container,
  CGeometry *geometry_container,
  COutput   *output_container
)
 /*
  * Function that writes the out the processed data to a file.
  */
{
  // Initialize the file information.
  std::string fileinfo = "\n";

  // Determine what processing was used.
  if( config_container->GetTypeProcessData() == PROCESS_IC ){

    // Determine which type of IC processing to use.
    switch( config_container->GetTypeIC(zoneID) ){

      // Pressure pulse IC.
      case(IC_GAUSSIAN_PRESSURE):
      {
        fileinfo += "# Reflection coefficient: R = max( p'(t) )/p'(t=0).";
        break;
      }

      // Entropy wave IC.
      case(IC_ENTROPY_WAVE):
      {
        fileinfo += "# Reflection coefficient: R = ( max( p'(t) )/pinf ) / ( max( T'(t=0) )/ Tinf ).";
        break;
      }

      // Isentropic vortex IC.
      case(IC_ISENTROPIC_VORTEX):
      {
        fileinfo += "# Reflection coefficient: R = ( max( p'(t) )/pinf ) / ( max( M'(t) )/ Minf ).";
        break;
      }

      // Vortex roll-up IC.
      case(IC_VORTEX_ROLLUP):
      {
        fileinfo += "# Probed data: pressure and v-velocity, for each sensor respectively.\n";
        fileinfo += "# Total number of sensors used: " + std::to_string(sensor_container.size()) + ".\n";
        fileinfo += "# As writing convention, probe locations are at:";
        for(unsigned short iProbe=0; iProbe<sensor_container.size(); iProbe++){
          fileinfo += "\n# Probe[" + std::to_string(iProbe) + "]: (";
          fileinfo += std::to_string(sensor_container[iProbe]->GetProbeLocation()[0]) + ", ";
          fileinfo += std::to_string(sensor_container[iProbe]->GetProbeLocation()[1]) + ")";
        }
        break;
      }

      // Do nothing.
      default: break;
    }
  }

  // Call the output write routine.
  output_container->WriteDataToFile(config_container,
                                    geometry_container,
                                    fileinfo.c_str(),
                                    TimeProcessed,
                                    DataProcessed);
}


void CEEProcess::ProcessData
(
  CConfig   *config_container,
  CGeometry *geometry_container,
  CElement  *element_container,
  CSpatial  *spatial_container,
  CSolver   *solver_container,
  CInitial  *initial_container,
  as3double  localTime
)
 /*
  * Function that selects what/how to process the data.
  */
{
  // Determine what processing to do.
  if( config_container->GetTypeProcessData() == PROCESS_IC ){

    // Determine which type of IC processing to use.
    switch( config_container->GetTypeIC(zoneID) ){

      // Pressure pulse IC.
      case(IC_GAUSSIAN_PRESSURE):
      {
        // Reflection coefficient: R = R(P).
        ComputeReflectionCoefficientFuncP(config_container,
                                          geometry_container,
                                          element_container,
                                          spatial_container,
                                          solver_container,
                                          initial_container,
                                          localTime);
        break;
      }

      // Entropy wave IC.
      case(IC_ENTROPY_WAVE):
      {
        // Reflection coefficient: R = R(P, T).
        ComputeReflectionCoefficientFuncPT(config_container,
                                           geometry_container,
                                           element_container,
                                           spatial_container,
                                           solver_container,
                                           initial_container,
                                           localTime);
        break;
      }

      // Isentropic vortex IC.
      case(IC_ISENTROPIC_VORTEX):
      {
        // Reflection coefficient: R = R(P, M).
        ComputeReflectionCoefficientFuncPM(config_container,
                                           geometry_container,
                                           element_container,
                                           spatial_container,
                                           solver_container,
                                           initial_container,
                                           localTime);
        break;
      }

      // Vortex roll-up IC.
      case(IC_VORTEX_ROLLUP):
      {
        // Sample data at predefined probes.
        SampleProbePoints(config_container,
                          geometry_container,
                          element_container,
                          spatial_container,
                          solver_container,
                          initial_container,
                          localTime);
        break;
      }

      // Do nothing.
      default: break;
    }
  }
}


void CEEProcess::ComputeReflectionCoefficientFuncP
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSpatial  *spatial_container,
 CSolver   *solver_container,
 CInitial  *initial_container,
 as3double  localTime
)
 /*
  * Function that computes the reflection coefficient as a function of (P).
  * That is: R = max( abs( p'(t)-p(t=0) ) )/p'(t=0).
  */
{
  // Extract total number of elements present in this zone.
  unsigned long nElem       = solver_container->GetnElem();
  // Extract number of solution points in 2D in this zone.
  unsigned short nDOFsSol2D = solver_container->GetnDOFsSol2D();

  // Extract mean data.
  const as3double pInf = initial_container->GetPinf();
  // Extract initial disturbance ratio.
  const as3double A0   = initial_container->GetA0();
  // Compute initial pressure disturbance.
  const as3double p0   = pInf*A0;
  // Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

  // Some abbreviations.
  const as3double gm1  = GAMMA_MINUS_ONE;

  // Maximum pressure ratio.
  as3double pratiomax = 0.0;

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {

    // Loop over all elements.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), reduction(max:pratiomax)
#endif
    for(unsigned long iElem=0; iElem<nElem; iElem++){

      // Extract current solution data.
      auto& sol = data_container[iElem]->GetDataDOFsSol();

      // Loop over all solution points and process the data.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsSol2D; l++){

        // Extract primitive data.
        const as3double rho   = sol[0][l];
        const as3double ovrho = 1.0/rho;
        const as3double u     = ovrho*sol[1][l];
        const as3double v     = ovrho*sol[2][l];
        const as3double p     = gm1*(sol[3][l]
                              - 0.5*(u*sol[1][l] + v*sol[2][l]) );

       // Compute ratio of pressure.
        const as3double pdiff = fabs(p - pInf);

        // Determine the max pressure quantity.
        pratiomax = std::max(pratiomax, pdiff);
      }
    }

  } // End of OpenMP parallel region.

  // Compute reflection coefficient.
  const as3double R = pratiomax/p0;

  // Register time.
  TimeProcessed.push_back(localTime);
  // Register data entry.
  DataProcessed[0].push_back(R);
}


void CEEProcess::ComputeReflectionCoefficientFuncPT
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSpatial  *spatial_container,
 CSolver   *solver_container,
 CInitial  *initial_container,
 as3double  localTime
)
 /*
  * Function that computes the reflection coefficient as a function of (P, T).
  * That is: R = ( max( p'(t) )/pinf ) / ( max( T'(t=0) )/ Tinf ).
  */
{
  // Extract total number of elements present in this zone.
  unsigned long nElem       = solver_container->GetnElem();
  // Extract number of solution points in 2D in this zone.
  unsigned short nDOFsSol2D = solver_container->GetnDOFsSol2D();

  // Extract max perturbation coefficient.
  const as3double A0   = initial_container->GetA0();
  // Extract mean data.
  const as3double pInf = initial_container->GetPinf();
  const as3double Tinf = initial_container->GetTinf();
  // Deduce max temperature perturbation as an IC.
  const as3double T0prime = Tinf*( A0 );

  // Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

  // Some abbreviations.
  const as3double gm1 = GAMMA_MINUS_ONE;

  // Maximum pressure.
  as3double pmax = 0.0;

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {

    // Loop over all elements.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), reduction(max:pmax)
#endif
    for(unsigned long iElem=0; iElem<nElem; iElem++){

      // Extract current solution data.
      auto& sol = data_container[iElem]->GetDataDOFsSol();

      // Loop over all solution points and process the data.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsSol2D; l++){

        // Extract primitive data.
        const as3double rho   = sol[0][l];
        const as3double ovrho = 1.0/rho;
        const as3double u     = ovrho*sol[1][l];
        const as3double v     = ovrho*sol[2][l];
        const as3double p     = gm1*(sol[3][l]
                              - 0.5*(u*sol[1][l] + v*sol[2][l]) );

        // Determine max pressure quantity.
        pmax = std::max(pmax, p);
      }
    }

  } // End of OpenMP parallel region.

  // Compute reflection coefficient.
  const as3double R = ( (pmax-pInf)/pInf )/( T0prime/Tinf );

  // Register time.
  TimeProcessed.push_back(localTime);
  // Register data entry.
  DataProcessed[0].push_back(R);
}


void CEEProcess::ComputeReflectionCoefficientFuncPM
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSpatial  *spatial_container,
 CSolver   *solver_container,
 CInitial  *initial_container,
 as3double  localTime
)
 /*
  * Function that computes the reflection coefficient as a function of (P, M).
  * That is: R = ( max( p(t)/pinf ) ) / ( max( M(t) )/ Minf ).
  */
{
  // Extract total number of elements present in this zone.
  unsigned long nElem       = solver_container->GetnElem();
  // Extract number of solution points in 2D in this zone.
  unsigned short nDOFsSol2D = solver_container->GetnDOFsSol2D();

  // Extract mean data.
  const as3double pInf = initial_container->GetPinf();
  const as3double Minf = initial_container->GetMach();

  // Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

  // Some abbreviations.
  const as3double gm1  = GAMMA_MINUS_ONE;

  // Maximum pressure and Mach number ratios.
  as3double pratiomax = 0.0, Mratiomax = 0.0;

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {

    // Loop over all elements.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), reduction(max:pratiomax), reduction(max:Mratiomax)
#endif
    for(unsigned long iElem=0; iElem<nElem; iElem++){

      // Extract current solution data.
      auto& sol = data_container[iElem]->GetDataDOFsSol();

      // Loop over all solution points and process the data.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsSol2D; l++){

        // Extract primitive data.
        const as3double rho   = sol[0][l];
        const as3double ovrho = 1.0/rho;
        const as3double u     = ovrho*sol[1][l];
        const as3double v     = ovrho*sol[2][l];
        const as3double p     = gm1*(sol[3][l]
                              - 0.5*(u*sol[1][l] + v*sol[2][l]) );

        // Magnitude of the velocity squared.
        const as3double umag2 = u*u + v*v;
        // Speed of sound squared.
        const as3double a2    = GAMMA*p*ovrho;

        // Compute the local Mach number.
        const as3double M = sqrt(umag2/a2);

        const as3double pratio = p/pInf;
        const as3double Mratio = M/Minf;

        // Determine the max pressure-ratio quantity.
        pratiomax = std::max(pratiomax, pratio);
        // Determine the max Mach-ratio quantity.
        Mratiomax = std::max(Mratiomax, Mratio);
      }
    }

  } // End of OpenMP parallel region.

  // Compute reflection coefficient.
  const as3double R = pratiomax/Mratiomax;

  // Register time.
  TimeProcessed.push_back(localTime);
  // Register data entry.
  DataProcessed[0].push_back(R);
}


void CEEProcess::SampleProbePoints
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSpatial  *spatial_container,
 CSolver   *solver_container,
 CInitial  *initial_container,
 as3double  localTime
)
 /*
  * Function that samples the pressure and v-velocity at different probe locations.
  */
{
  // Extract total number of elements present in this zone.
  unsigned long nElem       = solver_container->GetnElem();
  // Extract number of solution points in 2D in this zone.
  unsigned short nDOFsSol2D = solver_container->GetnDOFsSol2D();

  // Extract mean data.
  const as3double pInf = initial_container->GetPinf();
  const as3double Minf = initial_container->GetMach();

  // Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

  // Some abbreviations.
  const as3double gm1  = GAMMA_MINUS_ONE;

  // Loop over the number of probes specified.
  unsigned short idx = 0;
  for(unsigned short iProbe=0; iProbe<sensor_container.size(); iProbe++){

    // Make sure the current probe is in the same zone.
    if( sensor_container[iProbe]->GetIndexProbeZone() == zoneID ){

      // Extract element index carrying the probe sensor.
      unsigned long iElem = sensor_container[iProbe]->GetIndexProbeElem();

      // Extract interpolation polynomial.
      auto& ell = sensor_container[iProbe]->GetInterpolation();
      // Extract relevant solution data.
      auto& sol = data_container[iElem]->GetDataDOFsSol();

      // Consistency check.
      assert( ell.size() == nDOFsSol2D );

      // Initialize interpolated value.
      as3vector1d<as3double> val(nVar, 0.0);

      // Interpolate the solution of the sensor element to target probe location.
      for(unsigned short iVar=0; iVar<nVar; iVar++){
        // Loop over all solution points and process the data.
#pragma omp simd
        for(unsigned short l=0; l<nDOFsSol2D; l++) val[iVar] += ell[l]*sol[iVar][l];
      }

      // Convert the conservative data to primitive data.
      const as3double rho   = val[0];
      const as3double ovrho = 1.0/rho;
      const as3double u     = ovrho*val[1];
      const as3double v     = ovrho*val[2];
      const as3double p     = gm1*(val[3]
                            - 0.5*(u*val[1] + v*val[2]) );

      // Book-keep sensor data for pressure and v-velocity.
      DataProcessed[idx++].push_back(p);
      DataProcessed[idx++].push_back(v);
    }
  }

  // Register time.
  TimeProcessed.push_back(localTime);
}


CSensorData::CSensorData
(
 CConfig               *config_container,
 CGeometry             *geometry_container,
 CElement              *element_container,
 as3vector1d<as3double> probe,
 unsigned short         iZone
)
 /*
  * Constructor used to initialize CSensorData.
  */
{
  // Assign probe coordinates.
  ProbeLocation  = probe;
  // Assign current zone index.
  IndexProbeZone = iZone;

  // Extract current zone geometry.
  auto* geometry_zone = geometry_container->GetGeometryZone(iZone);

  // Number of DOFS of solution points on an element in 1D in this zone.
  unsigned short nDOFsSol1D = element_container->GetnDOFsSol1D();
  // Number of DOFS of solution points on an element in 2D in this zone.
  unsigned short nDOFsSol2D = element_container->GetnDOFsSol2D();

  // Search for the element containing the probe.
  const bool found = geometry_zone->SearchElementProbe(ProbeLocation, IndexProbeElem, true);

  // For now, this is only executed in the main zone. If element not found,
  // then throw an error and terminate.
  if( !found && (iZone == ZONE_MAIN) ){
    std::string message = "Probe: ("
                        + std::to_string(ProbeLocation[0]) + ", "
                        + std::to_string(ProbeLocation[1]) + ") "
                        + "could not be located, or is lying at an edge/corner.";
    Terminate("CSensorData::CSensorData", __FILE__, __LINE__, message);
  }

  // Reserve memory for the interpolation polynomial.
  Interpolation.resize(nDOFsSol2D);

  // Extract element.
  auto* geometry_element = geometry_zone->GetGeometryElem(IndexProbeElem);

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
    Terminate("CSensorData::CSensorData", __FILE__, __LINE__,
              "Error in lagrange interpolation polynomial computation.");
}


CSensorData::~CSensorData
(
 void
)
 /*
  * Destructor for CSensorData class, frees allocated memory.
  */
{

}


as3double CSensorData::EvaluateEntryLagrangePolynomial
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

