#include "spatial_structure.hpp"



CSpatial::CSpatial
(
 CConfig       *config_container,
 CGeometry 	   *geometry_container,
 CElement  	  **element_container,
 CInitial      *initial_container,
 unsigned short iZone
)
 /*
	* Constructor, used to initialize CSpatial in zone: iZone.
	*/
{
  // Assign zone ID.
  zoneID     = iZone;
  // Number of DOFS of solution points on an element in 1D in this zone.
  nDOFsSol1D = element_container[iZone]->GetnDOFsSol1D();
  // Number of DOFS of solution points on an element in 2D in this zone.
  nDOFsSol2D = element_container[iZone]->GetnDOFsSol2D();
  // Number of DOFS of integration points on an element in 1D.
	nDOFsInt1D = element_container[iZone]->GetnDOFsInt1D();
	// Number of DOFS of integration points on an element in 2D.
	nDOFsInt2D = element_container[iZone]->GetnDOFsInt2D();

  // Initialize the filtering matrix, if need be.
  if( config_container->GetTypeFilterSolution(iZone) != NO_FILTER )
    InitializeFilterMatrix(config_container,
                           geometry_container,
                           element_container[iZone]);
}


CSpatial::~CSpatial
(
 void
)
 /*
	* Destructor for CSpatial class, frees allocated memory.
	*/
{
  if( FilterMatrix      != nullptr ) delete [] FilterMatrix;
  if( riemann_container != nullptr ) delete riemann_container;
}


void CSpatial::InitializeFilterMatrix
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container
)
 /*
  * Function that initializes and computes the filtering matrix.
  */
{
  // Reserve memory for the filter matrix acting on the solution DOFs.
  FilterMatrix = new as3double[nDOFsSol2D*nDOFsSol2D]();

  // Check if allocation failed.
  if( FilterMatrix == nullptr )
    Terminate("CSpatial::InitializeFilterMatrix", __FILE__, __LINE__,
              "Allocation failed for FilterMatrix.");

  // Extract the Vandermonde matrix in 1D.
  auto* V1D    = element_container->GetVandermonde1D();
  // Extract the inverse of the Vandermonde matrix in 1D.
  auto* InvV1D = element_container->GetInvVandermonde1D();

  // Form temporary 2D versions of the Vandermonde matrix and its inversion.
  as3vector1d<as3double> V2D(nDOFsSol2D*nDOFsSol2D);
  as3vector1d<as3double> InvV2D(nDOFsSol2D*nDOFsSol2D);

  // Cast the filter matrix into a 2D array.
  as3double (*F2D)[nDOFsSol2D]    = (as3double (*)[nDOFsSol2D]) FilterMatrix;
  // Cast the Vandermonde matrices into 2D arrays.
  as3double (*v1D)[nDOFsSol1D]    = (as3double (*)[nDOFsSol1D]) V1D;
  as3double (*vInv1D)[nDOFsSol1D] = (as3double (*)[nDOFsSol1D]) InvV1D;
  as3double (*v2D)[nDOFsSol2D]    = (as3double (*)[nDOFsSol2D]) V2D.data();
  as3double (*vInv2D)[nDOFsSol2D] = (as3double (*)[nDOFsSol2D]) InvV2D.data();

  // Compute 2D Vandermonde matrices and inversion via 1D tensor-product.
  unsigned short iRow = 0, iCol = 0;
  for(unsigned short j2=0; j2<nDOFsSol1D; j2++){
    for(unsigned short i2=0; i2<nDOFsSol1D; i2++){
      // Reset column index.
      iCol = 0;
      for(unsigned short j1=0; j1<nDOFsSol1D; j1++){
        for(unsigned short i1=0; i1<nDOFsSol1D; i1++){
          v2D[iRow][iCol]    = v1D[j2][j1]*v1D[i2][i1];
          vInv2D[iRow][iCol] = vInv1D[j2][j1]*vInv1D[i2][i1];

          // Update column index.
          iCol++;
        }
      }
      // Update row index.
      iRow++;
    }
  }

  // Decide what type of filter to use.
  switch( config_container->GetTypeFilterSolution(zoneID) ){

    case(EXPONENTIAL_FILTER):
    {
      // Extract filter characteristics.
      auto& FilterCharacteristics = config_container->GetFilterCharacteristics();

      // Decompose filter characteristics.
      const as3double Nc    = (as3double) FilterCharacteristics[1];
      const as3double s     = (as3double) FilterCharacteristics[2];
      const as3double alpha = (as3double) FilterCharacteristics[3];

      // Damping function in 2D, initialize to one.
      as3vector1d<as3double> S2D(nDOFsSol2D, 1.0);

      // Step 1: compute the diagonal damping function.
      unsigned short sk = 0;
      for(int i=0; i<nDOFsSol1D; i++){
        for(int j=0; j<nDOFsSol1D; j++){
          if( (i+j) >= Nc ) S2D[sk] = exp( -alpha*pow( (i+j+1-Nc)/(nDOFsSol1D-Nc), s ) );
          // Update index.
          sk++;
        }
      }

      // Step 2: compute the C*Vinv and store it in vInv2D.
      for(unsigned short iRow=0; iRow<nDOFsSol2D; iRow++)
        for(unsigned short iCol=0; iCol<nDOFsSol2D; iCol++)
          vInv2D[iRow][iCol] = S2D[iRow]*vInv2D[iRow][iCol];

      // Step 3: compute the total filter V*C*Vinv and store it in F2D.
      for(unsigned short iRow=0; iRow<nDOFsSol2D; iRow++){
        for(unsigned short iCol=0; iCol<nDOFsSol2D; iCol++){
          as3double tmp = 0.0;
          for(unsigned short k=0; k<nDOFsSol2D; k++)
            // tmp += v1D[iRow][k]*F1D[k][iCol];
            tmp += v2D[iRow][k]*vInv2D[k][iCol];
          // Assign entry.
          F2D[iRow][iCol] = tmp;
        }
      }

      break;
    }

    default:
      Terminate("CSpatial::InitializeFilterMatrix", __FILE__, __LINE__,
                "Type of filter is not (yet) implemented!");
  }
}


void CSpatial::ComputeFilteredSolution
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CData     *data_container
)
 /*
  * Function that computes the filtered solution.
  */
{
  // Obtain current solution.
  auto& sol = data_container->GetDataDOFsSol();

  // Filter solution via matrix-matrix multiplication.
  for(unsigned short iVar=0; iVar<sol.size(); iVar++)
    gemv(nDOFsSol2D, nDOFsSol2D, FilterMatrix, sol[iVar], sol[iVar]);
}


CEESpatial::CEESpatial
(
 CConfig       *config_container,
 CGeometry 	   *geometry_container,
 CElement  	  **element_container,
 CInitial      *initial_container,
 unsigned short iZone
)
	:
		CSpatial
		(
		 config_container,
		 geometry_container,
		 element_container,
     initial_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CEESpatial in zone: iZone.
	*/
{
  // Preprocess riemann container.
  Riemann_Preprocessing(config_container, iZone);
}


CEESpatial::~CEESpatial
(
 void
)
 /*
	* Destructor for CEESpatial class, frees allocated memory.
	*/
{

}


void CEESpatial::InitializeSolution
(
  CConfig                *config_container,
  CInitial               *initial_container,
  CElement               *element_container,
  const CGeometryElement *grid_element,
  CData                  *data_element,
  as3double               time
)
 /*
  * Function that initializes the solution.
  */
{
  // If this is a restart, skip.
  if( config_container->GetRestartSolution() ) return;

  // Extract current solution DOFs.
  auto& data_nodes = data_element->GetDataDOFsSol();
  // Extract element solution nodes.
  auto& grid_nodes = grid_element->GetCoordSolDOFs();

  // Set initial condition specified in this element.
  initial_container->SetInitialCondition(grid_nodes,
                                         data_nodes,
                                         nDOFsSol2D,
                                         0.0);
}


void CEESpatial::Riemann_Preprocessing
(
  CConfig       *config_container,
  unsigned short iZone
)
 /*
  * Function that preprocesses the riemann container.
  */
{
  // Check which type of Riemann solver to use.
  switch( config_container->GetRiemannSolver(iZone) ){

    case(RIEMANN_ROE):       riemann_container = new CRoeRiemann(config_container);       break;
    case(RIEMANN_RUSANOV):   riemann_container = new CRusanovRiemann(config_container);   break;
    case(RIEMANN_ROEISMAIL): riemann_container = new CRoeIsmailRiemann(config_container); break;

    default:
      Terminate("CEESpatial::Riemann_Preprocessing", __FILE__, __LINE__,
                "Riemann solver is unknown!");
  }
}


void CEESpatial::ComputeVolumeResidual
(
  CConfig                *config_container,
  CGeometry              *geometry_container,
  CElement               *element_container,
  CInitial               *initial_container,
  CData                  *data_container,
  const CGeometryElement *geometry_element,
  as3data1d<as3double>   &work_array,
  as3double               localTime,
  as3vector1d<as3double> &MonitoringData
)
 /*
  * Function that computes the contribution of all volume terms to the residual.
  */
{
  // Some abbreviations.
  const as3double gm1 = GAMMA_MINUS_ONE;

	// Extract current element size.
  auto hElem = geometry_element->GetElemSize();

  // Step 0a: assign the pointers for the working solution and its gradients.
  as3double **Var    = work_array.data();
  as3double **dVarDx = Var    + nVar;
  as3double **dVarDy = dVarDx + nVar;

  // Step 0b: extract the solution and residual required.
  auto& sol = data_container->GetDataDOFsSol();
  auto& res = data_container->GetDataDOFsRes();

  // Step 0c: extract the required basis operators in 1D and weights in 2D.
  auto* ell  = element_container->GetLagrangeInt1D();
  auto* dell = element_container->GetDerLagrangeInt1D();
  auto* ellT = element_container->GetLagrangeInt1DTranspose();
  auto& wts  = element_container->GetwDOFsInt2D();

  // Step 1: interpolate the solution into the integration points.
  TensorProductSolAndGradVolume(nDOFsInt1D, nVar, nDOFsSol1D,
                                ellT, nullptr, sol.data(),
                                Var, nullptr, nullptr);

  // Step 2: compute the jacobian of the transformation.
  const as3double drdx = 2.0/hElem[XDIM];
  const as3double dsdy = 2.0/hElem[YDIM];

  // Step 3: loop over all integration points and compute the fluxes on them.
#pragma omp simd
  for(unsigned short l=0; l<nDOFsInt2D; l++){

    // Extract the integration weight multiplied with the jacobians.
    const as3double weightx = wts[l]*drdx;
    const as3double weighty = wts[l]*dsdy;

    // Determine primitive variables.
    const as3double rho   = Var[0][l];
    const as3double ovrho = 1.0/rho;
    const as3double u     = ovrho*Var[1][l];
    const as3double v     = ovrho*Var[2][l];
    const as3double p     = gm1*(Var[3][l]
                          - 0.5*(u*Var[1][l] + v*Var[2][l]) );

    // Compute the inviscid flux in the x-direction. Including the jacobian
    // and the integration weights.
    dVarDx[0][l] = weightx*(    Var[1][l]      ); // fx0
    dVarDx[1][l] = weightx*( u* Var[1][l] + p  ); // fx1
    dVarDx[2][l] = weightx*( v* Var[1][l]      ); // fx2
    dVarDx[3][l] = weightx*( u*(Var[3][l] + p) ); // fx3

    // Compute the inviscid flux in the y-direction. Including the jacobian
    // and the integration weights.
    dVarDy[0][l] = weighty*(    Var[2][l]      ); // fy0
    dVarDy[1][l] = weighty*( u* Var[2][l]      ); // fy1
    dVarDy[2][l] = weighty*( v* Var[2][l] + p  ); // fy2
    dVarDy[3][l] = weighty*( v*(Var[3][l] + p) ); // fy3
  }

  // Step 4a: compute flux contribution to the residual w.r.t. x-direction.
  TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                              dell, ell, dVarDx, res.data());

  // Step 4b: compute flux contribution to the residual w.r.t. y-direction.
  TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                              ell, dell, dVarDy, res.data());

  // Step 5: if there is a vortex-rollup source term, compute effect on residual.
  if( config_container->GetTypeIC(zoneID) == IC_VORTEX_ROLLUP ){

    // Extract pulse strength.
    const as3double A0   = initial_container->GetA0();
    // Extract pulse width.
    const as3double bb   = config_container->GetDisturbanceWidth();
    // Extract angular frequency.
    const as3double omg  = config_container->GetAngularFrequency();
    // Extract center of pulse.
    const as3double x0   = config_container->GetCenterX0()[0];
    const as3double y0   = config_container->GetCenterX0()[1];

    // Obtain boundary coordinates of current element.
    const as3double xmin = geometry_element->GetCoordBoundary(IDX_WEST);
    const as3double xmax = geometry_element->GetCoordBoundary(IDX_EAST);
    const as3double ymin = geometry_element->GetCoordBoundary(IDX_SOUTH);
    const as3double ymax = geometry_element->GetCoordBoundary(IDX_NORTH);

    // Flag to skip the current source term computations.
    bool skip = false;

    // Reduce computations by skipping elements who start past x = x0+50*bb.
    if( xmin > x0+50*bb ) skip = true;
    // Reduce computations by skipping elements who end before x = x0-50*bb.
    if( xmax < x0-50*bb ) skip = true;
    // Reduce computations by skipping elements who start past y = y0+50*bb.
    if( ymin > y0+50*bb ) skip = true;
    // Reduce computations by skipping elements who end before y = y0-50*bb.
    if( ymax < y0-50*bb ) skip = true;

    // Check if we may skip unnecessary computations.
    if( !skip ){
      // Extract reference speed of sound.
      const as3double aInf = initial_container->GetReferenceSpeedOfSound();
      // Extract reference pressure.
      const as3double pInf = initial_container->GetPinf();

      // Abbreviation.
      const as3double ovgm1   = 1.0/GAMMA_MINUS_ONE;
      const as3double ln2r2   = log(2.0)/(bb*bb);

      // Non-dimensionalize time, assuming reference length of one.
      const as3double tinf = localTime*aInf;
      // Dimensionalize the disturbance with the pressure.
      const as3double Ainf = ( initial_container->GetDimensionalProblem() ) ? A0*pInf*ovgm1 : A0;

      // Compute time-dependant disturbance coefficient.
      const as3double a0sinwt = Ainf*sin(omg*tinf);
      // Extract coordinates at integration nodes.
      auto& CoordInt = geometry_element->GetCoordIntDOFs();

      // Loop over all nodes and compute the source term.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){

        // Extract the integration weights.
        const as3double weight = -wts[l];

        // Extract coordinates.
        const as3double x = CoordInt[0][l];
        const as3double y = CoordInt[1][l];

        // Compute radial distance.
        const as3double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);

        // Store the source term temporarily in the first entry of Var.
        Var[0][l] = weight*a0sinwt*exp(-ln2r2*r2);
      }

      // Compute source term contribution to the residual. This is added only to
      // the energy equation entry in the residual.
      TensorProductVolumeResidual(nDOFsInt1D, 1, nDOFsSol1D,
                                  ell, ell, &Var[0], &res[3]);
    }
  }

  // Step 6: if there is a plane-wave source term, compute effect on residual.
  if( config_container->GetTypeIC(zoneID) == IC_ACOUSTIC_PLANE_WAVE ){

    // Extract pulse strength.
    const as3double A0   = initial_container->GetA0();
    // Extract spatial center of pulse.
    const as3double x0   = config_container->GetCenterX0()[0];
    // Extract temporal center of pulse.
    const as3double t0   = initial_container->Gett0();
    // Extract exponential coefficients in space and time.
    const as3double kx   = initial_container->GetKappax();
    const as3double kt   = initial_container->GetKappat();
    // Extract reference speed of sound.
    const as3double aInf = initial_container->GetReferenceSpeedOfSound();
    // Compute temporal radial distance.
    const as3double rt2  = (localTime-t0)*(localTime-t0);

    // Extract coordinates at integration nodes.
    auto& CoordInt = geometry_element->GetCoordIntDOFs();

    // Abbreviation.
    const as3double ovgm1  = 1.0/GAMMA_MINUS_ONE;
    const as3double ova    = 1.0/aInf;
    const as3double aovgm1 = aInf*ovgm1;
    const as3double A0expt = A0*exp(-kt*rt2);

    // Loop over all nodes and compute the source term.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt2D; l++){

      // Extract the integration weights.
      const as3double weight = wts[l];

      // Extract coordinates.
      const as3double x   = CoordInt[0][l];
      // Compute spatial radial distance.
      const as3double rx2 = (x-x0)*(x-x0);
      // Compute the source term.
      const as3double ss  = weight*A0expt*exp(-kx*rx2);

      // Store the source term in Var.
      Var[0][l] = ova*ss;
      Var[1][l] = ss;
      Var[2][l] = 0.0;
      Var[3][l] = aovgm1*ss;
    }

    // Compute source term contribution to the residual.
    TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                ell, ell, Var, res.data());
  }


  // If there is a periodic source-term, add its contribution on the residual.
  if( config_container->GetPeriodicPulse() ){

    // Physical time.
    const as3double time  = localTime;
		// Time step.
		const as3double dt    = config_container->GetTimeStep(); 
		// Starting angular frequency.
    const as3double omg0  = config_container->GetAngularFrequency();

    // Use a constant frequency definition, by default.
    as3double omg = omg0*time;
    // Check if the frequency is constant or varying. If it varies, supercede its value.
    if( !config_container->GetConstantFrequency() ){

			// Amplitude of pulse modification.
			const as3double Af  = config_container->GetSourceFrequencyParam()[0];
			// Period of pulse modification.
			const as3double tau = config_container->GetSourceFrequencyParam()[1];

			// Angular frequency function: 2*pi*fs, where fs = f0*Af^(t/tau).
      // What this means is the frequency is multiplied by Af (w.r.t. f0) every t=tau seconds.
      const as3double omgt = omg0*pow(Af, time/tau);
      // Actual angular frequency, taken by integrating omgt.
      omg  = omgt/( log(Af)/tau );
    }

    // Background pressure.
    const as3double p0  = initial_container->GetPinf();
    // Background speed of sound.
    const as3double a0  = initial_container->GetReferenceSpeedOfSound();

		// Disturbance center(s).
		auto C0Vec = config_container->GetCenterX0();
		
		// Deduce number of sources.
		unsigned short nSource = C0Vec.size()/nDim;
		// Extract coordinates of the source center(s).
		as3vector1d<as3double> C0x(nSource);
		as3vector1d<as3double> C0y(nSource);
		for(unsigned short iSource=0; iSource<nSource; iSource++){
		  C0x[iSource] = C0Vec[iSource*nDim+XDIM];
		  C0y[iSource] = C0Vec[iSource*nDim+YDIM];
		}

    // Disturbance amplitude.
    const as3double A0  = config_container->GetDisturbanceRatio();
    const as3double pp  = A0*p0;
    // Disturbance width.
    const as3double bb  = config_container->GetDisturbanceWidth();
		// Disturbance temporal-exponential coefficient.
		const as3double alpha = config_container->GetSourceFrequencyExponent();

    // Abbreviation.
    const as3double ova     = 1.0/a0;
		const as3double a2ovgm1 = a0*a0/gm1;


    // Time-dependant functions.
    const as3double sinwt = sin(omg);
		const as3double coswt = cos(omg);
    // Time-dependant source term expression.
		const as3double funct = ova*pp*pow( sinwt, alpha );
    // Disturbance spatial constant.
    const as3double kappa = log(2.0)/(bb*bb);

		// Initialize source-term modifications to zero.
		as3double ft = funct;
		
		// If source term is modified spatially, factor temporal-spatial coupling.
		if( !config_container->GetSourceTermCenterFixed() ){
		
		  // Initial time-transient part.
		  const as3double transient = exp(-100.0*dt/localTime);
		
		  // Obtain shift length.
		  const as3double dx = config_container->GetSourceTermCenterShift()[0];
		  const as3double dy = config_container->GetSourceTermCenterShift()[1];
		
		  // Modify source-term data, using the inputs specified.
		  ft *= transient;
		
		  // Modify source(s) centers.
		  for(unsigned short iSource=0; iSource<nSource; iSource++){
		    C0x[iSource] += dx*sinwt;
		    C0y[iSource] += dy*coswt;
		  }
		}

    // Extract coordinates at integration nodes.
    auto& CoordInt = geometry_element->GetCoordIntDOFs();

    // Check if source term is angled along the flow direction.
    if( config_container->GetAlignedPeriodicPulse() ){

      // Background velocity.
      const as3double uInf = initial_container->GetUinf();
      const as3double vInf = initial_container->GetVinf();

      // Abbreviation.
      const as3double umag = sqrt( uInf*uInf + vInf*vInf );
			// Flow direction, scaled with speed of sound.
			const as3double idir = a0*uInf/umag;
			const as3double jdir = a0*vInf/umag;

      // Loop over all nodes and compute the source term.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){

        // Extract the integration weights.
        const as3double weight = wts[l];

        // Extract coordinates.
        const as3double x  = CoordInt[0][l];
        const as3double y  = CoordInt[1][l];
				
				// Compute radial distance, per pulse center.
				as3double ss = 0.0;
				for(unsigned short iSource=0; iSource<nSource; iSource++){
				  // Extract pulse center.
				  const as3double x0 = C0x[iSource];
				  const as3double y0 = C0y[iSource];
				
				  // Compute spatial radial distance.
				  const as3double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);
				  // Compute the spatially-dependent part of the source term.
				  ss += exp(-kappa*r2);
				}
				// Compute overall source term(s).
				ss *= weight*ft;

        // Store the source term in Var.
				Var[0][l] = ss;
				Var[1][l] = ss*idir;
				Var[2][l] = ss*jdir;
				Var[3][l] = ss*a2ovgm1;
      }
      // Compute source term contribution to the residual.
      TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                  ell, ell, Var, res.data() );
    }
    else {

      // Loop over all nodes and compute the source term.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){

        // Extract the integration weights.
        const as3double weight = wts[l];

        // Extract coordinates.
        const as3double x  = CoordInt[0][l];
        const as3double y  = CoordInt[1][l];
				
				// Compute radial distance, per pulse center.
				as3double ss = 0.0;
				for(unsigned short iSource=0; iSource<nSource; iSource++){
				  // Extract pulse center.
				  const as3double x0 = C0x[iSource];
				  const as3double y0 = C0y[iSource];
				
				  // Compute spatial radial distance.
				  const as3double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);
				  // Compute the spatially-dependent part of the source term.
				  ss += exp(-kappa*r2);
				}
				// Compute overall source term(s).
				ss *= weight*ft;

        // Store the source term in Var.
				Var[0][l] = ss;
				Var[1][l] = 0.0;
				Var[2][l] = 0.0;
				Var[3][l] = ss*a2ovgm1;
      }
      // Compute source term contribution to the residual.
      TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                  ell, ell, Var, res.data() );
    }
  }

	// Compute max Mach number over the solution DOFs.
  as3double M2max = 0.0;
	for(unsigned short l=0; l<nDOFsSol2D; l++){

    // Determine primitive variables.
    const as3double rho   = sol[0][l];
    const as3double ovrho = 1.0/rho;
    const as3double u     = ovrho* sol[1][l];
    const as3double v     = ovrho* sol[2][l];
    const as3double p     = gm1*(  sol[3][l]
                          - 0.5*(u*sol[1][l] 
													     + v*sol[2][l]) );

    // Magnitude of the velocity squared.
    const as3double umag2 = u*u + v*v;
    // Speed of sound squared.
    const as3double a2    = GAMMA*p*ovrho;

    // Compute the local Mach number squared.
    const as3double M2 = umag2/a2;
    // Check if this value is the largest.
    M2max = std::max(M2max, M2);
	}


  // Assign the monitoring data.
  MonitoringData[0] = M2max;
}


void CEESpatial::ComputeSurfaceResidual
(
  CConfig              *config_container,
  CGeometry            *geometry_container,
  CElement             *element_container,
  as3element           &data_container,
  as3data1d<as3double> &work_array,
  as3double             localTime,
  unsigned long         iElem
)
 /*
  * Function that computes the contribution of all surfacce terms contribution
  * to the residual.
  */
{
  // Some abbreviations and initializations.
  auto* geometry_zone = geometry_container->GetGeometryZone(zoneID);
  auto& MatchingFace  = geometry_container->GetMatchingFace();

  // Extract current element size.
  auto hElem = geometry_zone->GetGeometryElem(iElem)->GetElemSize();

  // Step 0a: extract indices of internal elements and their neighbors.
  auto& InternalElemFace = geometry_zone->GetInternalElemFace()[iElem];
  auto& InternalNeighbor = geometry_zone->GetIndexNeighborInternalElement()[iElem];

  // Step 0b: assign the pointers for the working solution.
  as3double **VarI = work_array.data();
  as3double **VarJ = VarI + nVar;
  as3double **Flux = VarJ + nVar;

  // Step 0c: extract the current data container solution.
  auto& dataSolI = data_container[iElem]->GetDataDOFsSol();
  auto& res      = data_container[iElem]->GetDataDOFsRes();

  // Step 1: extract the required basis operators in 1D and integration weights in 1D.
  auto* ell  = element_container->GetLagrangeInt1D();
  auto* ellT = element_container->GetLagrangeInt1DTranspose();
  auto& wts  = element_container->GetwDOFsInt1D();

  // Step 2: loop over all the faces and interpolate the neighboring solution.
  for(unsigned short iFace=0; iFace<nFace; iFace++){

    // Check whether this is an internal or boundary face. In case this is a
    // boundary face, the data is already interpolated to the current face in
    // the preprocessing step, so dont do anything.
    if( InternalElemFace[iFace] ){

      // This is an internal face. Obtain the neighboring element index.
      unsigned long  jElem = InternalNeighbor[iFace];
      // Extract the matching face index.
      unsigned short jFace = MatchingFace[iFace];

      // Get the indices for the solution on the matcihng element face.
      auto& FaceIndexJ = element_container->GetIndexDOFsSol(jFace);

      // Extract the solution DOFs at the matching element.
      auto& dataSolJ = data_container[jElem]->GetDataDOFsSol();

      // Interpolate matching solution to integration points on the face.
      TensorProductSolAndGradFace(jFace, nDOFsInt1D, nVar, nDOFsSol1D,
                                  FaceIndexJ.data(), ellT,
                                  nullptr, nullptr,
                                  dataSolJ.data(), VarJ,
                                  nullptr, nullptr);
    }
    else {

      // This is a boundary face. Obtain the data that is already interpolated
      // to the integration points on this face. Note, we need to convert the
      // reference to a pointer in order to store the data in Var.
      VarJ = data_container[iElem]->GetDataDOFsIntFace(iFace).data();
    }

    // Get the indices for the solution on the current element face.
    auto& FaceIndexI = element_container->GetIndexDOFsSol(iFace);
    // Get unit-normal for this face.
    auto& UnitNormal = geometry_container->GetUnitNormal(iFace);

    // Step 3: Interpolate current solution to integration points on the face.
    TensorProductSolAndGradFace(iFace, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(), ellT,
                                nullptr, nullptr,
                                dataSolI.data(), VarI,
                                nullptr, nullptr);

    // Step 4: invoke a riemann solver to determine the flux at the face.
    riemann_container->ComputeFluxState(UnitNormal, wts, hElem, VarI, VarJ, Flux);

    // Step 5: compute flux contribution to the residual.
    TensorProductSurfaceResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                 FaceIndexI.data(), ell,
                                 Flux, res.data());
  }
}


CEESpongeSpatial::CEESpongeSpatial
(
 CConfig       *config_container,
 CGeometry 	   *geometry_container,
 CElement  	  **element_container,
 CInitial      *initial_container,
 unsigned short iZone
)
	:
		CEESpatial
		(
		 config_container,
		 geometry_container,
		 element_container,
     initial_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CEESpongeSpatial in zone: iZone.
	*/
{
  // Initialize the artificial-convection on each face.
  ArtificialConvectionFace.resize(nFace, false);
  if( config_container->GetArtificialConvection(iZone) ){

    // Check the type of zone used.
    switch( config_container->GetTypeZone(iZone) ){

      case(ZONE_EAST): case(ZONE_WEST):
      {
        ArtificialConvectionFace[IDX_EAST] = true;
        ArtificialConvectionFace[IDX_WEST] = true;
        break;
      }

      case(ZONE_SOUTH): case(ZONE_NORTH):
      {
        ArtificialConvectionFace[IDX_SOUTH] = true;
        ArtificialConvectionFace[IDX_NORTH] = true;
        break;
      }

      case(ZONE_CORNER_0): case(ZONE_CORNER_1): case(ZONE_CORNER_2): case(ZONE_CORNER_3):
      {
        for(unsigned short iFace=0; iFace<nFace; iFace++)
          ArtificialConvectionFace[iFace] = true;
        break;
      }

      default:
        Terminate("CEESpongeSpatial::CEESpongeSpatial", __FILE__, __LINE__,
                  "Wrong zone type inserted.");
    }
  }
}


CEESpongeSpatial::~CEESpongeSpatial
(
 void
)
 /*
	* Destructor for CEESpongeSpatial class, frees allocated memory.
	*/
{

}


void CEESpongeSpatial::InitializeSolution
(
  CConfig                *config_container,
  CInitial               *initial_container,
  CElement               *element_container,
  const CGeometryElement *grid_element,
  CData                  *data_element,
  as3double               time
)
 /*
  * Function that initializes the solution.
  */
{
  // Extract current solution DOFs.
  auto& DataSol  = data_element->GetDataDOFsSol();
  // Extract element solution nodes.
  auto& CoordSol = grid_element->GetCoordSolDOFs();

  // Extract current element size.
  auto hElem = grid_element->GetElemSize();

  // Set initial condition specified in this element. If this is a restart, skip.
  if( !config_container->GetRestartSolution() )
    initial_container->SetInitialCondition(CoordSol, DataSol, nDOFsSol2D, 0.0);

  // Extract target pseudo-mean solution.
  auto& DataIntMean = data_element->GetDataDOFsIntMean();
  // Extract element integration nodes.
  auto& CoordInt    = grid_element->GetCoordIntDOFs();

  // Set pseudo-mean flow in this element.
  initial_container->SetInitialCondition(CoordInt, DataIntMean, nDOFsInt2D, 0.0);

  // Check if this is a dimensional problem, if so dimensionalize the damping coefficients.
  if( initial_container->GetDimensionalProblem() ){

    // Extract the reference speed-of-sound.
    const as3double aInf = initial_container->GetReferenceSpeedOfSound();
    // Extract non-dimensional damping functions.
    auto& sigma          = data_element->GetDampDOFsInt();

    // Normalization w.r.t. x- and y-dimensions.
    const as3double aovhx = aInf/hElem[0];
    const as3double aovhy = aInf/hElem[1];

    // Loop over the integration points in 2D and dimensionalize the damping functions.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt2D; l++){

      // Normalize damping coefficients.
      sigma[0][l] *= aovhx;
      sigma[1][l] *= aovhy;
    }

    // Check if a supersonic velocity term is needed and determine the scaling ratio.
    if( config_container->GetArtificialConvection(zoneID) ){

      // Extract the non-dimensional artificial-convection at volume points.
      auto& vel_volume = data_element->GetDerLagrangeArtificialConvection1D();

      // Loop over the artificial-velocity and dimensionalize the velocity functions.
      for(unsigned short i=0; i<vel_volume.size(); i++){
        if( vel_volume[i] ){
#pragma omp simd
          for(unsigned short l=0; l<nDOFsSol1D*nDOFsInt1D; l++){
            // Dimensionalize velocity coefficients.
            vel_volume[i][l] *= aInf;
          }
        }
      }

      // Extract the non-dimensional artificial velocities at the surfaces.
      auto& vel_face = data_element->GetDerLagrangeArtificialConvectionFace();

      // Dimensionalize the surface boundary velocity coefficients on each face.
      for(unsigned short i=0; i<vel_face.size(); i++) vel_face[i] *= aInf;
    }
  }

}


void CEESpongeSpatial::ComputeVolumeResidual
(
  CConfig                *config_container,
  CGeometry              *geometry_container,
  CElement               *element_container,
  CInitial               *initial_container,
  CData                  *data_container,
  const CGeometryElement *geometry_element,
  as3data1d<as3double>   &work_array,
  as3double               localTime,
  as3vector1d<as3double> &MonitoringData
)
 /*
  * Function that computes the contribution of all volume terms to the residual.
  */
{
  // Some abbreviations.
  const as3double gm1 = GAMMA_MINUS_ONE;

  // Extract current element size.
  auto hElem = geometry_element->GetElemSize();

  // Step 0a: assign the pointers for the working solution and its gradients.
  as3double **Var    = work_array.data();
  as3double **dVarDx = Var    + nVar;
  as3double **dVarDy = dVarDx + nVar;
  as3double **VarSrc = dVarDy + nVar;

  // Step 0b: extract the solution and residual required.
  auto& sol = data_container->GetDataDOFsSol();
  auto& res = data_container->GetDataDOFsRes();

  // Step 0c: extract damping functions and pseudo-mean flow.
  auto& sigma  = data_container->GetDampDOFsInt();
  auto& VarInf = data_container->GetDataDOFsIntMean();

  // Step 0d: extract the required basis operators in 1D and weights in 2D.
  auto* ell  = element_container->GetLagrangeInt1D();
  auto* dell = element_container->GetDerLagrangeInt1D();
  auto* ellT = element_container->GetLagrangeInt1DTranspose();
  auto& wts  = element_container->GetwDOFsInt2D();

  // Grid-stretching lagrange differential operators.
  auto& derLagGridStretching = data_container->GetDerLagrangeGridStretching1D();
  // Select lagrange differential operators, explicitly.
  const as3double *dellx = ( derLagGridStretching[0] ) ? derLagGridStretching[0] : dell;
  const as3double *delly = ( derLagGridStretching[1] ) ? derLagGridStretching[1] : dell;

  // Artificial-convection lagrange differential operators.
  auto& u0dell = data_container->GetDerLagrangeArtificialConvection1D();
  // Select the lagrange differential operators for artificial-convection, explicitly.
  const as3double *u0dellx = u0dell[0];
  const as3double *u0delly = u0dell[1];

  // Step 1: interpolate the solution into the integration points.
  TensorProductSolAndGradVolume(nDOFsInt1D, nVar, nDOFsSol1D,
                                ellT, nullptr, sol.data(),
                                Var, nullptr, nullptr);

  // Step 2: compute the jacobian of the transformation.
  const as3double drdx = 2.0/hElem[XDIM];
  const as3double dsdy = 2.0/hElem[YDIM];

  // Step 3: loop over all integration points and compute the fluxes on them.
#pragma omp simd
  for(unsigned short l=0; l<nDOFsInt2D; l++){

		// NOTE, the residual is defined on the right-hand side, opposite to the definition
		// in VCP3D (which is on the LHS).
    // Extract the integration weight multiplied with the jacobians.
    const as3double weightx = wts[l]*drdx;
    const as3double weighty = wts[l]*dsdy;

    // Extract damping function coefficients.
    const as3double sx = sigma[0][l];
    const as3double sy = sigma[1][l];

    // Weights associated with the damping terms.
    const as3double weights = -wts[l]*(sx+sy);

    // Determine primitive variables.
    const as3double rho   = Var[0][l];
    const as3double ovrho = 1.0/rho;
    const as3double u     = ovrho*Var[1][l];
    const as3double v     = ovrho*Var[2][l];
    const as3double p     = gm1*(Var[3][l]
                          - 0.5*(u*Var[1][l] + v*Var[2][l]) );

    // Compute the inviscid flux in the x-direction. Including the jacobian
    // and the integration weights.
    dVarDx[0][l] = weightx*(    Var[1][l]      ); // fx0
    dVarDx[1][l] = weightx*( u* Var[1][l] + p  ); // fx1
    dVarDx[2][l] = weightx*( v* Var[1][l]      ); // fx2
    dVarDx[3][l] = weightx*( u*(Var[3][l] + p) ); // fx3

    // Compute the inviscid flux in the y-direction. Including the jacobian
    // and the integration weights.
    dVarDy[0][l] = weighty*(    Var[2][l]      ); // fy0
    dVarDy[1][l] = weighty*( u* Var[2][l]      ); // fy1
    dVarDy[2][l] = weighty*( v* Var[2][l] + p  ); // fy2
    dVarDy[3][l] = weighty*( v*(Var[3][l] + p) ); // fy3

    // Compute the damping terms and store them in Var.
    for(unsigned short i=0; i<nVar; i++)
      VarSrc[i][l] = weights*( Var[i][l] - VarInf[i][l] );
  }


  // Step 6a: compute flux contribution to the residual w.r.t. x-direction.
  TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                              dellx, ell, dVarDx, res.data());

  // Step 6b: compute flux contribution to the residual w.r.t. y-direction.
  TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                              ell, delly, dVarDy, res.data());

  // Step 6c: compute source contribution to the residual.
  TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                              ell, ell, VarSrc, res.data());


	// Check if a characteristic layer is needed.
	if( config_container->GetCharacteristicMatching(zoneID) )
	{
		// Abbreviation involving gamma.
		const as3double ovgm1 = 1.0/gm1;

		// Extract the characteristic-matching profile.
		auto& match = data_container->GetMatchDOFsInt();

		// Extract transpose of 1D derivative Lagrange polynomial.
		auto* dellT = element_container->GetDerLagrangeInt1DTranspose();

		// Compute the data at the integation nodes.
  	TensorProductSolAndGradVolume(nDOFsInt1D, nVar, nDOFsSol1D,
  	                              ellT, dellT, sol.data(),
  	                              Var, dVarDx, dVarDy);

    // Convert the parametric gradients into physical space.
    for(unsigned short iVar=0; iVar<nVar; iVar++){
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){
        dVarDx[iVar][l] *= drdx;
        dVarDy[iVar][l] *= dsdy;
      }
    }

		// Deduce which indices are the outgoing indices in each direction. Then, assign
		// a value of zero for each of the dxI and dyI values. This way, only the incoming
		// indices are modified, as expected.

		// Obtain the layer orientation, w.r.t. the main zone.
		const as3double nx = geometry_container->GetGeometryZone(zoneID)->GetLayerOrientation()[0];
		const as3double ny = geometry_container->GetGeometryZone(zoneID)->GetLayerOrientation()[1];


  	// Loop over all integration points and compute the additional characteristic
		// matching term on them.
#pragma omp simd
  	for(unsigned short l=0; l<nDOFsInt2D; l++){
    	
	   	// Weights associated source terms.
    	const as3double wx = wts[l]*match[0][l];
			const as3double wy = wts[l]*match[1][l];

			// Determine primitive variables.
    	const as3double rho   = Var[0][l];
			const as3double rhou  = Var[1][l];
			const as3double rhov  = Var[2][l];
			const as3double ovrho = 1.0/rho;
    	const as3double u     = ovrho*rhou;
    	const as3double v     = ovrho*rhov;
    	const as3double p     = gm1*(Var[3][l]
    	                      - 0.5*(u*rhou + v*rhov) );

			// Compute the speed of sound and its square.
			const as3double a2    = GAMMA*p*ovrho;
			const as3double a     = sqrt( a2 );

			// Some abbreviations commonly used.
			const as3double ova2  = 1.0/a2;
			const as3double rova  = rho/a;
			const as3double rhoa  = rho*a;
			const as3double ovra  = 1.0/rhoa;
			const as3double ek    = 0.5*( u*u + v*v );

			// Compute the eigenvalues in both directions. Note, the if statements are
			// there in order to identify which eigenvalues are outgoing, thus completely
			// eliminate them -- since we only want to match the incoming characteristics,
			// with respect to the main, physical zone.
			
			// First in the x-direction.
			as3double lmbx1 = (u - a);
			as3double lmbx2 = (u    );
			as3double lmbx3 = (u + a);

			// Then in the y-direction.
			as3double lmby1 = (v - a);
			as3double lmby2 = (v    );
			as3double lmby3 = (v + a);

			// Eliminate any outgoing eigenvalues in x-direction.
			if( nx*lmbx1 > 0.0 ) lmbx1 = 0.0;
			if( nx*lmbx2 > 0.0 ) lmbx2 = 0.0;
			if( nx*lmbx3 > 0.0 ) lmbx3 = 0.0;

			// Eliminate any outgoing eigenvalues in y-direction.
			if( ny*lmby1 > 0.0 ) lmby1 = 0.0;
			if( ny*lmby2 > 0.0 ) lmby2 = 0.0;
			if( ny*lmby3 > 0.0 ) lmby3 = 0.0;

			// Extract the derivative of the primitive variables in the x-direction.
      const as3double drdx = dVarDx[0][l];
      const as3double dudx = ovrho*(  dVarDx[1][l] - u*drdx );
    	const as3double dvdx = ovrho*(  dVarDx[2][l] - v*drdx );
    	const as3double dpdx = gm1*( ek*dVarDx[0][l]
      										 -  		  u*dVarDx[1][l]
      										 -        v*dVarDx[2][l]
      										 +          dVarDx[3][l] );


			// Extract the derivative of the primitive variables in the y-direction.
      const as3double drdy = dVarDy[0][l];
      const as3double dudy = ovrho*(  dVarDy[1][l] - u*drdy );
    	const as3double dvdy = ovrho*(  dVarDy[2][l] - v*drdy );
    	const as3double dpdy = gm1*( ek*dVarDy[0][l]
      										 -  		  u*dVarDy[1][l]
      										 -        v*dVarDy[2][l]
      										 +          dVarDy[3][l] );


			// Compute the directional wave amplitudes and factor in the indices that are 
			// only incoming, zero all outgoing information.
			
			// First in the x-direction.
			const as3double Lx1 = 0.5*lmbx1*( dudx - ovra*dpdx );
			const as3double Lx2 =     lmbx2*( drdx - ova2*dpdx );
			const as3double Lx3 =     lmbx2*( dvdx             );
			const as3double Lx4 = 0.5*lmbx3*( dudx + ovra*dpdx );

			// Then in the y-direction.
			const as3double Ly1 = 0.5*lmby1*( dvdy - ovra*dpdy );
			const as3double Ly2 =     lmby2*( dudy             );
			const as3double Ly3 =     lmby2*( drdy - ova2*dpdy );
			const as3double Ly4 = 0.5*lmby3*( dvdy + ovra*dpdy );


			// Some abbreviations involving the wave amplitudes: x-direction.
			const as3double Lx14p     = Lx1 + Lx4;
			const as3double Lx14m     = Lx1 - Lx4;
			const as3double rhoLx14m  = rho*Lx14m;
			const as3double termx     = Lx2 + rova*Lx14p;

			// Some abbreviations involving the wave amplitudes: y-direction.
			const as3double Ly14p     = Ly1 + Ly4;
			const as3double Ly14m     = Ly1 - Ly4;
			const as3double rhoLy14m  = rho*Ly14m;
			const as3double termy     = Ly3 + rova*Ly14p;
			
			// Subtract the modified/incoming quasi-linear flux in the x-direction: dFdx.
			VarSrc[0][l] = wx*(    termx                                            );
			VarSrc[1][l] = wx*(  u*termx - rhoLx14m                                 );
			VarSrc[2][l] = wx*(  v*termx +  rho*Lx3                                 ); 
			VarSrc[3][l] = wx*( ek*termx + rhov*Lx3 - rhou*Lx14m + ovgm1*rhoa*Lx14p );
	
			// Subtract the modified/incoming quasi-linear flux in the y-direction: dGdy.
			VarSrc[0][l] += wy*(    termy                                            );
			VarSrc[1][l] += wy*(  u*termy +  rho*Ly2                                 );
			VarSrc[2][l] += wy*(  v*termy - rhoLy14m                                 ); 
			VarSrc[3][l] += wy*( ek*termy + rhou*Ly2 - rhov*Ly14m + ovgm1*rhoa*Ly14p );
		}

  	// Compute source contribution to the residual.
  	TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
  	                            ell, ell, VarSrc, res.data());
	}



  // Check if an artificial-convection term is needed in the x-direction.
  if( u0dellx ){

    // Loop over all integration points and compute the artificial convection.
    for(unsigned short i=0; i<nVar; i++){
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){
        // Add artificial convection contributions in x-direction.
        dVarDx[i][l] = drdx*wts[l]*Var[i][l];
      }
    }

    // Compute artificial-convection flux contribution to the residual w.r.t. x-direction.
    TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                u0dellx, ell, dVarDx, res.data());
  }

  // Check if an artificial-convection term is needed in the y-direction.
  if( u0delly ){

    // Loop over all integration points and compute the artificial convection.
    for(unsigned short i=0; i<nVar; i++){
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){
        // Add artificial convection contributions in y-direction.
        dVarDy[i][l] = dsdy*wts[l]*Var[i][l];
      }
    }

    // Compute artificial-convection flux contribution to the residual w.r.t. y-direction.
    TensorProductVolumeResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                ell, u0delly, dVarDy, res.data());
  }


	// Compute max Mach number over the solution DOFs.
  as3double M2max = 0.0;
	for(unsigned short l=0; l<nDOFsSol2D; l++){

    // Determine primitive variables.
    const as3double rho   = sol[0][l];
    const as3double ovrho = 1.0/rho;
    const as3double u     = ovrho* sol[1][l];
    const as3double v     = ovrho* sol[2][l];
    const as3double p     = gm1*(  sol[3][l]
                          - 0.5*(u*sol[1][l] 
													     + v*sol[2][l]) );

    // Magnitude of the velocity squared.
    const as3double umag2 = u*u + v*v;
    // Speed of sound squared.
    const as3double a2    = GAMMA*p*ovrho;

    // Compute the local Mach number squared.
    const as3double M2 = umag2/a2;
    // Check if this value is the largest.
    M2max = std::max(M2max, M2);
	}

  // Assign the monitoring data.
  MonitoringData[0] = M2max;
}


void CEESpongeSpatial::ComputeSurfaceResidual
(
  CConfig              *config_container,
  CGeometry            *geometry_container,
  CElement             *element_container,
  as3element           &data_container,
  as3data1d<as3double> &work_array,
  as3double             localTime,
  unsigned long         iElem
)
 /*
  * Function that computes the contribution of all surfacce terms contribution
  * to the residual.
  */
{
  // Some abbreviations and initializations.
  auto* geometry_zone = geometry_container->GetGeometryZone(zoneID);
  auto& MatchingFace  = geometry_container->GetMatchingFace();

  // Extract current element size.
  auto hElem = geometry_zone->GetGeometryElem(iElem)->GetElemSize();

  // Step 0a: extract indices of internal elements and their neighbors.
  auto& InternalElemFace = geometry_zone->GetInternalElemFace()[iElem];
  auto& InternalNeighbor = geometry_zone->GetIndexNeighborInternalElement()[iElem];

  // Step 0b: assign the pointers for the working solution.
  as3double **VarI = work_array.data();
  as3double **VarJ = VarI + nVar;
  as3double **Flux = VarJ + nVar;

  // Step 0c: extract the current data container solution.
  auto& dataSolI = data_container[iElem]->GetDataDOFsSol();
  auto& res      = data_container[iElem]->GetDataDOFsRes();

  // Step 1: extract the required basis operators in 1D and integration weights in 1D.
  auto* ell  = element_container->GetLagrangeInt1D();
  auto* ellT = element_container->GetLagrangeInt1DTranspose();
  auto& wts  = element_container->GetwDOFsInt1D();	


  // Step 2: loop over all the faces and interpolate the neighboring solution.
  for(unsigned short iFace=0; iFace<nFace; iFace++){

    // Initialize element sizes.
    as3vector1d<as3double> hh = { hElem[0], hElem[1] };

    // If there is grid-stretching, account for that, depending on the current face.
    // Note, the stretching coefficients are included in the element sizes.
    if( config_container->GetGridStretching(zoneID) ){

      // Extract the grid-streching factor on this face.
      const as3double aa = data_container[iElem]->GetDerLagrangeGridStretchingFace(iFace);

      // Factor the grid-stretching coefficient in the element sizes.
      hh[0] /= aa; hh[1] /= aa;
    }

    // Check whether this is an internal or boundary face. In case this is a
    // boundary face, the data is already interpolated to the current face in
    // the preprocessing step, so dont do anything.
    if( InternalElemFace[iFace] ){

      // This is an internal face. Obtain the neighboring element index.
      unsigned long  jElem = InternalNeighbor[iFace];
      // Extract the matching face index.
      unsigned short jFace = MatchingFace[iFace];

      // Get the indices for the solution on the matcihng element face.
      auto& FaceIndexJ = element_container->GetIndexDOFsSol(jFace);

      // Extract the solution DOFs at the matching element.
      auto& dataSolJ = data_container[jElem]->GetDataDOFsSol();

      // Interpolate matching solution to integration points on the face.
      TensorProductSolAndGradFace(jFace, nDOFsInt1D, nVar, nDOFsSol1D,
                                  FaceIndexJ.data(), ellT,
                                  nullptr, nullptr,
                                  dataSolJ.data(), VarJ,
                                  nullptr, nullptr);
    }
    else {

      // This is a boundary face. Obtain the data that is already interpolated
      // to the integration points on this face. Note, we need to convert the
      // reference to a pointer in order to store the data in Var.
      VarJ = data_container[iElem]->GetDataDOFsIntFace(iFace).data();
    }

    // Get the indices for the solution on the current element face.
    auto& FaceIndexI = element_container->GetIndexDOFsSol(iFace);
    // Get unit-normal for this face.
    auto& UnitNormal = geometry_container->GetUnitNormal(iFace);

    // Step 3: Interpolate current solution to integration points on the face.
    TensorProductSolAndGradFace(iFace, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(), ellT,
                                nullptr, nullptr,
                                dataSolI.data(), VarI,
                                nullptr, nullptr);

    // Step 4: invoke a riemann solver to determine the flux at the face.
		riemann_container->ComputeFluxState(UnitNormal, wts, hh, VarI, VarJ, Flux);

    // Step 5: compute flux contribution to the residual.
    TensorProductSurfaceResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                 FaceIndexI.data(), ell,
                                 Flux, res.data());


    // Step 6: check if there need be any artificial convection needed.
    if( ArtificialConvectionFace[iFace] ){

      // Extract artificial velocity.
      const as3double u0 = data_container[iElem]->GetDerLagrangeArtificialConvectionFace(iFace);

      // Unit normals.
      const as3double nx = UnitNormal[0];
      const as3double ny = UnitNormal[1];
      // Element sizes.
      const as3double hx = hh[0];
      const as3double hy = hh[1];

      // Some abbreviations.
      const as3double unovh = u0*( nx/hx + ny/hy );

      // Normal convective speed.
      const as3double u0n   = u0*nx  + u0*ny;
      // Sign of the normal artificial dissipation component.
      const as3double sign  = ( fabs(u0n) > 1.0e-15 ) ? u0n/fabs(u0n) : 0.0;

      // Use full-upwinding to determine the boundary state. Note, this is
      // temporarily stored in Flux.
      for(unsigned short i=0; i<nVar; i++){
#pragma omp simd
        for(unsigned short l=0; l<nDOFsInt1D; l++){

          // Compute the weights and Jacobian. Note, the factor 2.0 in the
          // Jacobian is ignored, since it is being multiplied by a 0.5 in the
          // upwinding/central-differencing step.
          const as3double scale = -unovh*wts[l];

          // Compute upwinding. Note, the half is removed because it is
          // implicitly taken into account by the 2.0 in the Jacobian.
          Flux[i][l] = scale*(       ( VarJ[i][l] + VarI[i][l] )
                     -          sign*( VarJ[i][l] - VarI[i][l] ) );
        }
      }

      // Compute flux contribution to the residual.
      TensorProductSurfaceResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                   FaceIndexI.data(), ell,
                                   Flux, res.data());
    }
  }
}


CEEPMLSpatial::CEEPMLSpatial
(
 CConfig       *config_container,
 CGeometry 	   *geometry_container,
 CElement  	  **element_container,
 CInitial      *initial_container,
 unsigned short iZone
)
	:
		CEESpatial
		(
		 config_container,
		 geometry_container,
		 element_container,
     initial_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CEEPMLSpatial in zone: iZone.
	*/
{
  // Set the dispersion-relation correction coefficient.
  DispersionCorrection = initial_container->GetBetaPML();
  // Set the transverse velocity component. Note, in this work
	// this is assumed to be the v-velocity component.
  VelocityTransverse   = initial_container->GetVinf();
}


CEEPMLSpatial::~CEEPMLSpatial
(
 void
)
 /*
	* Destructor for CEEPMLSpatial class, frees allocated memory.
	*/
{

}


void CEEPMLSpatial::InitializeSolution
(
  CConfig                *config_container,
  CInitial               *initial_container,
  CElement               *element_container,
  const CGeometryElement *grid_element,
  CData                  *data_element,
  as3double               time
)
 /*
  * Function that initializes the solution.
  */
{
  // Extract current solution DOFs.
  auto& DataSol  = data_element->GetDataDOFsSol();
  // Extract element solution nodes.
  auto& CoordSol = grid_element->GetCoordSolDOFs();

  // Extract current element size.
  auto hElem = grid_element->GetElemSize();

  // Set initial condition specified in this element. If this is a restart, skip.
  if( !config_container->GetRestartSolution() )
    initial_container->SetInitialCondition(CoordSol, DataSol, nDOFsSol2D, 0.0);

  // Extract target pseudo-mean solution.
  auto& DataIntMean = data_element->GetDataDOFsIntMean();
  // Extract element integration nodes.
  auto& CoordInt    = grid_element->GetCoordIntDOFs();
  // Extract transpose of lagrange interpolation polynomial in 1D.
  auto* ellT        = element_container->GetLagrangeInt1DTranspose();

  // Set pseudo-mean flow in this element.
  initial_container->SetInitialCondition(CoordInt, DataIntMean, nDOFsInt2D, 0.0);

  // Loop over all faces and determine the pseudo-mean data on that face.
  for(unsigned short iFace=0; iFace<nFace; iFace++){

    // Get the indices for the solution on the current element face.
    auto& FaceIndexI = element_container->GetIndexDOFsSol(iFace);
    // Extract pseudo-mean solution at this face.
    auto& DataIntMeanFace = data_element->GetDataDOFsIntMeanFace(iFace);

    // Interpolate the grid coordinates from 2D solution points into the
    // 1D integration points on the current face. Note, these are temporarily
    // stored in DataIntMeanFace for now.
    TensorProductSolAndGradFace(iFace, nDOFsInt1D, nDim, nDOFsSol1D,
                                FaceIndexI.data(), ellT,
                                nullptr, nullptr,
                                CoordSol.data(),
                                DataIntMeanFace.data(),
                                nullptr, nullptr);

    // Set the pseudo-mean flow data of this face.
    initial_container->SetInitialCondition(DataIntMeanFace, DataIntMeanFace, nDOFsInt1D, 0.0);
  }

  // Check if this is a dimensional problem, if so the dimensionalize the damping coefficient.
  if( initial_container->GetDimensionalProblem() ){

    // Extract non-dimensional damping functions.
    auto& sigma = data_element->GetDampDOFsInt();
    // Determine the reference speed-of-sound.
    const as3double aInf = initial_container->GetReferenceSpeedOfSound();

    // Normalization w.r.t. x- and y-dimensions.
    const as3double aovhx = aInf/hElem[XDIM];
    const as3double aovhy = aInf/hElem[YDIM];

    // Loop over the integration points in 2D and dimensionalize the damping functions.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt2D; l++){

      // Normalize damping coefficients.
      sigma[0][l] *= aovhx;
      sigma[1][l] *= aovhy;
    }
  }
}


void CEEPMLSpatial::ComputeVolumeResidual
(
  CConfig                *config_container,
  CGeometry              *geometry_container,
  CElement               *element_container,
  CInitial               *initial_container,
  CData                  *data_container,
  const CGeometryElement *geometry_element,
  as3data1d<as3double>   &work_array,
  as3double               localTime,
  as3vector1d<as3double> &MonitoringData
)
 /*
  * Function that computes the contribution of all volume terms to the residual.
  */
{
  // Some abbreviations.
  const as3double gm1 = GAMMA_MINUS_ONE;
  // Dispersion-relation correction coefficient.
  const as3double bb = DispersionCorrection;
  // Number of variables for the cross terms contributoin, default nVar (no cross terms).
  unsigned short nVarCross = nVar;

  // Extract current element size.
  auto hElem = geometry_element->GetElemSize();

  // Step 0a: assign the pointers for the working solution and auxiliary data.
  as3double **VarU    = work_array.data();
  as3double **VarQ    = VarU    + nVar;
  as3double **VarSrcU = VarQ    + nVar;
  as3double **VarSrcQ = VarSrcU + nVar;
  as3double **dVarDxU = VarSrcQ + nVar;
  as3double **dVarDxQ = dVarDxU + nVar;
  as3double **dVarDyU = dVarDxQ + nVar;
  as3double **dVarDyQ = dVarDyU + nVar;

  // Abbreviation for the total variables.
  as3double **VarTotal    = VarU;
  as3double **VarSrcTotal = VarSrcU;
  as3double **dVarDxTotal = dVarDxU;
  as3double **dVarDyTotal = dVarDyU;

  // Step 0b: extract total solution.
  auto& sol = data_container->GetDataDOFsSol();
  // Step 0c: extract total residual.
  auto& res = data_container->GetDataDOFsRes();

  // Step 0d: extract damping functions and pseudo-mean flow.
  auto& sigma  = data_container->GetDampDOFsInt();
  auto& VarInf = data_container->GetDataDOFsIntMean();

  // Step 0e: extract the required basis operators in 1D and weights in 2D.
  auto* ell  = element_container->GetLagrangeInt1D();
  auto* dell = element_container->GetDerLagrangeInt1D();
  auto* ellT = element_container->GetLagrangeInt1DTranspose();
  auto& wts  = element_container->GetwDOFsInt2D();

  // Grid-stretching lagrange differential operators.
  auto& derLagGridStretching = data_container->GetDerLagrangeGridStretching1D();

  // Select lagrange differential operators.
  const as3double *dellx = ( derLagGridStretching[0] ) ? derLagGridStretching[0] : dell;
  const as3double *delly = ( derLagGridStretching[1] ) ? derLagGridStretching[1] : dell;


  // Step 1: interpolate the total solution into the integration points.
  // Note, VarTotal already has VarU and VarQ in it. No need to explicitly update them.
  TensorProductSolAndGradVolume(nDOFsInt1D, 2*nVar, nDOFsSol1D,
                                ellT, nullptr, sol.data(),
                                VarTotal, nullptr, nullptr);

  // Step 2: compute the jacobian of the transformation.
  const as3double drdx = 2.0/hElem[XDIM];
  const as3double dsdy = 2.0/hElem[YDIM];

  // Step 3: loop over all integration points and compute the fluxes on them.
#pragma omp simd
  for(unsigned short l=0; l<nDOFsInt2D; l++){

    // Extract integration weights only. Note, the residual is defined on the
    // right-hand-side of the equality. Hence, this weight is used on source
    // terms, which are negative since they are on the right or the equality.
    const as3double weight = -wts[l];

    // Extract the integration weight multiplied with the jacobians.
    const as3double weightJacx = wts[l]*drdx;
    const as3double weightJacy = wts[l]*dsdy;

    // Extract damping function coefficients.
    const as3double sx = sigma[0][l];
    const as3double sy = sigma[1][l];

    // Determine the primitive variables based on the pseudo-mean data.
    const as3double rhoInf   = VarInf[0][l];
    const as3double ovrhoInf = 1.0/rhoInf;
    const as3double uInf     = ovrhoInf*VarInf[1][l];
    const as3double vInf     = ovrhoInf*VarInf[2][l];
    const as3double pInf     = gm1*(VarInf[3][l]
                             - 0.5*(uInf*VarInf[1][l] + vInf*VarInf[2][l]) );

    // Determine primitive variables.
    const as3double rho   = VarU[0][l];
    const as3double ovrho = 1.0/rho;
    const as3double u     = ovrho*VarU[1][l];
    const as3double v     = ovrho*VarU[2][l];
    const as3double p     = gm1*(VarU[3][l]
                          - 0.5*(u*VarU[1][l] + v*VarU[2][l]) );

    // Compute the inviscid flux in the x-direction. Including the jacobian
    // and the integration weights.
    dVarDxU[0][l] =    VarU[1][l];      // fx0
    dVarDxU[1][l] = u* VarU[1][l] + p;  // fx1
    dVarDxU[2][l] = v* VarU[1][l];      // fx2
    dVarDxU[3][l] = u*(VarU[3][l] + p); // fx3

    // Compute the inviscid flux in the y-direction.
    // Note, the integration weights and Jacobian are not included.
    dVarDyU[0][l] =    VarU[2][l];      // fy0
    dVarDyU[1][l] = u* VarU[2][l];      // fy1
    dVarDyU[2][l] = v* VarU[2][l] + p;  // fy2
    dVarDyU[3][l] = v*(VarU[3][l] + p); // fy3

    // Compute the x-flux based on the pseudo-mean flow. Note, these are
    // stored temporarily in dVarDxQ for now.
    dVarDxQ[0][l] =       VarInf[1][l];         // fxInf0
    dVarDxQ[1][l] = uInf* VarInf[1][l] + pInf;  // fxInf1
    dVarDxQ[2][l] = vInf* VarInf[1][l];         // fxInf2
    dVarDxQ[3][l] = uInf*(VarInf[3][l] + pInf); // fxInf3

    for(unsigned short i=0; i<nVar; i++){

      // Step a: compute: fprime = fx - fxInf.
      dVarDxQ[i][l] = dVarDxU[i][l] - dVarDxQ[i][l];

      // Step b: compute: sx*b*fprime + sx*q.
      VarSrcQ[i][l] = sx*( bb*dVarDxQ[i][l] + VarQ[i][l] );

      // Step c: compute: sx*b*fprime + sx*q + sy*( u - uInf - q ).
      VarSrcU[i][l] = VarSrcQ[i][l]
                    + sy*( VarU[i][l] - VarInf[i][l] - VarQ[i][l] );

      // Step d: factor jacobians and weights for physical residual
      VarSrcU[i][l] *= weight;
      dVarDxU[i][l] *= weightJacx;
      dVarDyU[i][l] *= weightJacy;

      // Step e: factor jacobians and weights for auxiliary residual
      VarSrcQ[i][l] *= weight;
      dVarDxQ[i][l] *= weightJacx;
    }

  }

  // Check if there is a cross-flow, if so, add the cross-term contribution
  // to the auxiliary residual.
  if( config_container->GetCrossFlow() ){

    // Extract transverse velocity component.
    const as3double V0 = VelocityTransverse;

    for(unsigned short iVar=0; iVar<nVar; iVar++){
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){
        dVarDyQ[iVar][l] = wts[l]*dsdy*V0*VarQ[iVar][l];
      }
    }

    // Correct the cross-term variable contribution to the y-flux.
    nVarCross = 2*nVar;
  }


  // Step 4b: compute totol flux contribution to the residual w.r.t. x-direction.
  TensorProductVolumeResidual(nDOFsInt1D, 2*nVar, nDOFsSol1D,
                              dellx, ell, dVarDxTotal, res.data());

  // Step 4a: compute flux contribution to the residual w.r.t. y-direction.
  TensorProductVolumeResidual(nDOFsInt1D, nVarCross, nDOFsSol1D,
                              ell, delly, dVarDyTotal, res.data());

  // Step 4c: compute total source contribution to the physical residual.
  TensorProductVolumeResidual(nDOFsInt1D, 2*nVar, nDOFsSol1D,
                              ell, ell, VarSrcTotal, res.data());


	// Compute max Mach number over the solution DOFs.
  as3double M2max = 0.0;
	for(unsigned short l=0; l<nDOFsSol2D; l++){

    // Determine primitive variables.
    const as3double rho   = sol[0][l];
    const as3double ovrho = 1.0/rho;
    const as3double u     = ovrho* sol[1][l];
    const as3double v     = ovrho* sol[2][l];
    const as3double p     = gm1*(  sol[3][l]
                          - 0.5*(u*sol[1][l] 
														   + v*sol[2][l]) );

    // Magnitude of the velocity squared.
    const as3double umag2 = u*u + v*v;
    // Speed of sound squared.
    const as3double a2    = GAMMA*p*ovrho;

    // Compute the local Mach number squared.
    const as3double M2 = umag2/a2;
    // Check if this value is the largest.
    M2max = std::max(M2max, M2);
	}

  // Assign the monitoring data.
  MonitoringData[0] = M2max;
}


void CEEPMLSpatial::ComputeSurfaceResidual
(
  CConfig              *config_container,
  CGeometry            *geometry_container,
  CElement             *element_container,
  as3element           &data_container,
  as3data1d<as3double> &work_array,
  as3double             localTime,
  unsigned long         iElem
)
 /*
  * Function that computes the contribution of all surfacce terms contribution
  * to the residual.
  */
{
  // Some abbreviations and initializations.
  auto* geometry_zone = geometry_container->GetGeometryZone(zoneID);
  auto& MatchingFace  = geometry_container->GetMatchingFace();

  // Extract current element size.
  auto hElem = geometry_zone->GetGeometryElem(iElem)->GetElemSize();

  // Step 0a: extract indices of internal elements and their neighbors.
  auto& InternalElemFace = geometry_zone->GetInternalElemFace()[iElem];
  auto& InternalNeighbor = geometry_zone->GetIndexNeighborInternalElement()[iElem];

  // Step 0b: assign the pointers for the working solution.
  as3double **VarUI = work_array.data();
  as3double **VarUJ = VarUI + nVar;
  as3double **Flux  = VarUJ + nVar;
  as3double **VarQI = Flux  + nVar;
  as3double **VarQJ = VarQI + nVar;

  // Step 0c: extract the solution and residual required.
  auto& solI = data_container[iElem]->GetDataDOFsSol();

  // Step 0d: extract the residual required.
  auto& res  = data_container[iElem]->GetDataDOFsRes();

  // Step 0e: explicitly separate physical and auxiliary solution.
  as3double **sol_phy = solI.data();
  as3double **sol_aux = sol_phy + nVar;

  // Step 0f: explicitly separate physical and auxiliary residual.
  as3double **res_phy = res.data();
  as3double **res_aux = res_phy + nVar;

  // Step 1: extract the required basis operators in 1D and integration weights in 1D.
  auto* ell  = element_container->GetLagrangeInt1D();
  auto* ellT = element_container->GetLagrangeInt1DTranspose();
  auto& wts  = element_container->GetwDOFsInt1D();


  // Step 2: loop over all the faces and interpolate the neighboring solution.
  for(unsigned short iFace=0; iFace<nFace; iFace++){

    // Check if this is a face whose normal is in only the x-direction.
    const bool VerticalFace = ( (iFace==IDX_EAST) || (iFace==IDX_WEST) ) ? true : false;

    // Initialize element sizes.
    as3vector1d<as3double> hh = { hElem[0], hElem[1] };

    // If there is grid-stretching, account for that, depending on the current face.
    // Note, the stretching coefficients are included in the element sizes.
    if( config_container->GetGridStretching(zoneID) ){

      // Extract the grid-streching factor on this face.
      const as3double aa = data_container[iElem]->GetDerLagrangeGridStretchingFace(iFace);

      // Factor the grid-stretching coefficient in the element sizes.
      hh[0] /= aa; hh[1] /= aa;
    }


    // Check whether this is an internal or boundary face. In case this is a
    // boundary face, the data is already interpolated to the current face in
    // the preprocessing step, so dont do anything.
    if( InternalElemFace[iFace] ){

      // This is an internal face. Obtain the neighboring element index.
      unsigned long  jElem = InternalNeighbor[iFace];
      // Extract the matching face index.
      unsigned short jFace = MatchingFace[iFace];

      // Get the indices for the solution on the matcihng element face.
      auto& FaceIndexJ = element_container->GetIndexDOFsSol(jFace);
      // Extract the solution DOFs at the matching element.
      auto& solJ = data_container[jElem]->GetDataDOFsSol();

      // Interpolate matching solution to integration points on the face.
      TensorProductSolAndGradFace(jFace, nDOFsInt1D, nVar, nDOFsSol1D,
                                  FaceIndexJ.data(), ellT,
                                  nullptr, nullptr,
                                  solJ.data(), VarUJ,
                                  nullptr, nullptr);

      // Interpolate the auxiliary solution on the surface if there is some
      // cross-flow and this is a south or north face.
      if( config_container->GetCrossFlow() && !VerticalFace ){

        // Extract the auxiliary solution DOFs at the matching element.
        auto** auxJ = &data_container[jElem]->GetDataDOFsSol()[nVar];

        // Interpolate the matching auxiliary solution to the integration points on this face.
        TensorProductSolAndGradFace(jFace, nDOFsInt1D, nVar, nDOFsSol1D,
                                    FaceIndexJ.data(), ellT,
                                    nullptr, nullptr,
                                    auxJ, VarQJ,
                                    nullptr, nullptr);
      }
    }
    else {

      // This is a boundary face. Obtain the data that is already interpolated
      // to the integration points on this face.
      VarUJ = data_container[iElem]->GetDataDOFsIntFace(iFace).data();

      // Extract the boundary state of the auxiliary solution, if it is needed.
      if( config_container->GetCrossFlow() && !VerticalFace )
        VarQJ = data_container[iElem]->GetDataDOFsIntAuxFace(iFace).data();
    }

    // Get the indices for the solution on the current element face.
    auto& FaceIndexI = element_container->GetIndexDOFsSol(iFace);
    // Get unit-normal for this face.
    auto& UnitNormal = geometry_container->GetUnitNormal(iFace);

    // Step 3: Interpolate current solution to integration points on the face.
    TensorProductSolAndGradFace(iFace, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(), ellT,
                                nullptr, nullptr,
                                sol_phy, VarUI,
                                nullptr, nullptr);

    // Step 4: invoke a riemann solver to determine the flux at the face.
    riemann_container->ComputeFluxState(UnitNormal, wts, hh, 
																				VarUI, VarUJ, Flux);

    // Step 5: compute flux contribution to the physical residual.
    TensorProductSurfaceResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                 FaceIndexI.data(), ell,
                                 Flux, res_phy);

    // Add the pseudo-mean flux contribution, in case this is a face
    // that has its normal in the x-direction.
    if( VerticalFace ){

      // Extract pseudo-mean flow.
      auto& VarInf = data_container[iElem]->GetDataDOFsIntMeanFace(iFace);

      // Step 6: invoke a riemann solver to determine the pseudo-mean flux at the face.
      // Note, the result is stored in VarUI temporarily.
      riemann_container->ComputeFluxState(UnitNormal, wts, hh,
                                          VarInf.data(), VarInf.data(), VarUI);

      // Step 7: compute the difference between the x-fluxes.
      // Note, VarUI contains the pseudo-mean flux state at the boundary.
      for(unsigned short i=0; i<nVar; i++){
#pragma omp simd
        for(unsigned short l=0; l<nDOFsInt1D; l++)
          Flux[i][l] -= VarUI[i][l];
      }

      // Step 8: compute flux contribution to the physical residual.
      TensorProductSurfaceResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                   FaceIndexI.data(), ell,
                                   Flux, res_aux);
    }
    else {
      // This is either a south or a north face. If there is a cross-flow present,
      // then compute the auxiliary variable at this boundary and then add the
      // contribution to the residual accordingly.
      if( config_container->GetCrossFlow() ){

        // Extract the transverse velocity.
        const as3double V0 = VelocityTransverse;

        // Interpolate the current auxiliary solution at the current face.
        TensorProductSolAndGradFace(iFace, nDOFsInt1D, nVar, nDOFsSol1D,
                                    FaceIndexI.data(), ellT,
                                    nullptr, nullptr,
                                    sol_aux, VarQI,
                                    nullptr, nullptr);


				// Compute the unique state of the auxiliary variable via upwinding.
				riemann_container->ComputeVariableStateUpwinding(UnitNormal, wts, hh, V0,
				                                                 VarQI, VarQJ, Flux);

        // Compute flux contribution to the auxiliary residual.
        TensorProductSurfaceResidual(nDOFsInt1D, nVar, nDOFsSol1D,
                                     FaceIndexI.data(), ell,
                                     Flux, res_aux);
      }
    }
  }
}




