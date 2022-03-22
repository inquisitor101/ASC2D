#include "boundary_structure.hpp"



CBoundary::CBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
 /*
	* Constructor, used to initialize CBoundary per zone: iZone.
	*/
{
	// Assign zone ID.
	zoneID     = iZone;
	// Assign boundary ID that is unique to this zone only.
	boundaryID = iBoundary;
  // Assign type of zone.
  typeZone   = config_container->GetTypeZone(iZone);
  // Extract number of solution nodes per each element's boundary face in this zone.
  nDOFsSol1D = element_container[iZone]->GetnDOFsSol1D();
	// Extract number of integration nodes per each element's boundary face in this zone.
	nDOFsInt1D = element_container[iZone]->GetnDOFsInt1D();
  // Extract the current unit-normal for this boundary.
  UnitNormal = geometry_container->GetUnitNormal(iBoundary);

  // Kronecker-delta to facilitate some computations.
  KronDelta11 = fabs( UnitNormal[0] );
  KronDelta22 = fabs( UnitNormal[1] );

  // Element indices that share this boundary.
  ElemIndexI = geometry_container->GetElemBoundaryIndex(iZone, iBoundary);
}


CBoundary::~CBoundary
(
 void
)
 /*
	* Destructor for CBoundary class, frees allocated memory.
	*/
{

}


CEEBoundary::CEEBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEBoundary per zone: iZone.
	*/
{

}


void CEEBoundary::InitializePrescribedState
(
  CConfig       *config_container,
  CGeometry     *geometry_container,
  CElement      *element_container,
  CInitial      *initial_container,
  unsigned short iZone,
  unsigned short iBoundary
)
 /*
  * Function that initializes and computes the prescribed data at the boundary.
  */
{
  // Get current grid zone.
  auto* grid_zone = geometry_container->GetGeometryZone(iZone);

  // Get the transposed lagrange interpolation function in 1D.
  auto* lagrangeInt1DTranspose = element_container->GetLagrangeInt1DTranspose();
  // Get the indices for the solution nodes on an element in the current boundary.
  auto& FaceIndexI = element_container->GetIndexDOFsSol(iBoundary);

  // Number of elements on this boundary.
  unsigned long nsElem = ElemIndexI.size();

  // Reserve memory for the target state.
  DataDOFsIntBoundary.resize(nsElem);

  // Allocate the memory per element.
  for(unsigned long iElem=0; iElem<nsElem; iElem++){

    // Reserve memory for each variable.
    DataDOFsIntBoundary[iElem].resize(nVar, nullptr);

    // Allocate memory per each variable.
    for(unsigned short iVar=0; iVar<nVar; iVar++){

      // Allocate actual memory.
      DataDOFsIntBoundary[iElem][iVar] = new as3double[nDOFsInt1D]();

      // Check if allocation failed.
      if( !DataDOFsIntBoundary[iElem][iVar] )
        Terminate("CEEBoundary::InitializePrescribedState", __FILE__, __LINE__,
                   "Allocation failed for DataDOFsIntBoundary.");
    }
  }

  // Loop over all the elements and determine the target state.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract element grid nodes at solution points.
    auto& CoordSol = grid_zone->GetGeometryElem(iElem)->GetCoordSolDOFs();

    // Interpolate the coordinates to the integration points at the current face.
    // Note, the interpolated coordinates are temporarily stored in DataDOFsIntBoundary[i].
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nDim, nDOFsSol1D,
                                FaceIndexI.data(),
                                lagrangeInt1DTranspose,
                                nullptr, nullptr,
                                CoordSol.data(),
                                DataDOFsIntBoundary[i].data(),
                                nullptr, nullptr);

    // Compute target state, given the initial condition.
    initial_container->ComputeTargetStatePrimitive(DataDOFsIntBoundary[i],
                                                   DataDOFsIntBoundary[i],
                                                   nDOFsInt1D);
  }
}


void CEEBoundary::InitializeModifiedBC
(
  CConfig       *config_container,
  CGeometry     *geometry_container,
  CElement      *element_container,
  CInitial      *initial_container,
  unsigned short iZone,
  unsigned short iBoundary
)
 /*
  * Function that initializes the modified boundary data needed.
  */
{
  // Get current grid zone.
  auto* grid_zone = geometry_container->GetGeometryZone(iZone);

  // Get the transposed lagrange interpolation function in 1D.
  auto* lagrangeInt1DTranspose = element_container->GetLagrangeInt1DTranspose();
  // Get the indices for the solution nodes on an element in the current boundary.
  auto& FaceIndexI = element_container->GetIndexDOFsSol(iBoundary);

  // Number of elements on this boundary.
  unsigned long nsElem = ElemIndexI.size();

  // Reserve memory for the original target state.
  OrigDOFsIntBoundary.resize(nsElem);

  // Allocate the memory per element.
  for(unsigned long iElem=0; iElem<nsElem; iElem++){

    // Reserve memory for each variable.
    OrigDOFsIntBoundary[iElem].resize(nVar, nullptr);

    // Allocate memory per each variable.
    for(unsigned short iVar=0; iVar<nVar; iVar++){

      // Allocate actual memory.
      OrigDOFsIntBoundary[iElem][iVar] = new as3double[nDOFsInt1D]();

      // Check if allocation failed.
      if( !OrigDOFsIntBoundary[iElem][iVar] )
        Terminate("CEEBoundary::InitializeModifiedBC", __FILE__, __LINE__,
                   "Allocation failed for OrigDOFsIntBoundary.");

      // Assign data to initial conditions defined.
      for(unsigned short iNode=0; iNode<nDOFsInt1D; iNode++)
        OrigDOFsIntBoundary[iElem][iVar][iNode] = DataDOFsIntBoundary[iElem][iVar][iNode];
    }
  }


  // Reserve memory for the boundary coordinates.
  GridDOFsIntBoundary.resize(nsElem);

  // Allocate the memory per element.
  for(unsigned long iElem=0; iElem<nsElem; iElem++){

    // Reserve memory for each dimension.
    GridDOFsIntBoundary[iElem].resize(nDim, nullptr);

    // Allocate memory per each dimension.
    for(unsigned short iDim=0; iDim<nDim; iDim++){

      // Allocate actual memory.
      GridDOFsIntBoundary[iElem][iDim] = new as3double[nDOFsInt1D]();

      // Check if allocation failed.
      if( !GridDOFsIntBoundary[iElem][iDim] )
        Terminate("CEEBoundary::InitializeModifiedBC", __FILE__, __LINE__,
                   "Allocation failed for GridDOFsIntBoundary.");
    }
  }


  // Loop over all the elements and determine the boundary coordinates at integration points.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract element grid nodes at solution points.
    auto& CoordSol = grid_zone->GetGeometryElem(iElem)->GetCoordSolDOFs();

    // Interpolate the coordinates to the integration points at the current face.
    // Note, the interpolated coordinates are stored in GridDOFsIntBoundary[i].
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nDim, nDOFsSol1D,
                                FaceIndexI.data(),
                                lagrangeInt1DTranspose,
                                nullptr, nullptr,
                                CoordSol.data(),
                                GridDOFsIntBoundary[i].data(),
                                nullptr, nullptr);
  }
}


void CEEBoundary::ModifyBoundaryCondition
(
  CConfig    *config_container,
  CGeometry  *geometry_container,
  CSolver    *solver_container,
  CElement   *element_container,
  CSpatial   *spatial_container,
  as3double   localTime
)
 /*
  * Function that modifies the prescribed boundary data.
  */
{
  // Decide what modification to make.
  switch( config_container->GetTypeModifyBC(boundaryID) ){

    // Do nothing.
    case(NO_BC_MODIFICATION): break;

    // Gaussian pressure profile.
    case(GAUSSIAN_PRESSURE_BC_MODIFICATION):
    {
      // Abbreviations.
      const as3double gm1   = GAMMA_MINUS_ONE;
      const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;

      // Extract angular frequency specified.
      auto W0Vec = config_container->GetModifyFreqBC();
      // Extract pulse strength ratio w.r.t. background flow.
      auto A0Vec = config_container->GetModifyStrengthBC();
      // Extract pulse width w.r.t. background flow.
      auto bbVec = config_container->GetModifyWidthBC();
      // Extract pulse(s) center(s).
      auto C0Vec = config_container->GetModifyCenterBC();

      // Deduce number of pulses.
      unsigned short nPulse = C0Vec.size()/nDim;

      // Extract data pertaining to only the pressure variable.
      const as3double omg = W0Vec[3];
      const as3double A0  = A0Vec[3];
      const as3double bb  = bbVec[3];

      // Extract coordinates of the pulse center(s).
      as3vector1d<as3double> C0x(nPulse);
      as3vector1d<as3double> C0y(nPulse);
      for(unsigned short iPulse=0; iPulse<nPulse; iPulse++){
        C0x[iPulse] = C0Vec[iPulse*nDim+XDIM];
        C0y[iPulse] = C0Vec[iPulse*nDim+YDIM];
      }

      // Compute time-dependent factor.
      const as3double a0sinwt = A0*sin(omg*localTime);
      // Compute spatial exponential factor.
      const as3double kappa   = log(2.0)/(bb*bb);

      // Loop over all elements sharing this boundary.
      for(unsigned long i=0; i<ElemIndexI.size(); i++){

        // Abbreviation for the primitive data.
        as3double **Var0 = OrigDOFsIntBoundary[i].data(); // original data.
        as3double **Var1 = DataDOFsIntBoundary[i].data(); // modified data.

        // Abbreviation for the grid data.
        as3double **Coord = GridDOFsIntBoundary[i].data();

        // Loop over integration points and compute the modified boundary state.
#pragma omp simd
        for(unsigned short l=0; l<nDOFsInt1D; l++){

          // Extract grid coordinates.
          const as3double x = Coord[0][l];
          const as3double y = Coord[1][l];

          // Reset modified data to original data.
          Var1[0][l] = Var0[0][l];
          Var1[1][l] = Var0[1][l];
          Var1[2][l] = Var0[2][l];
          Var1[3][l] = Var0[3][l];

          // Compute radial distance, per pulse center.
          as3double fexp = 0.0;
          for(unsigned short iPulse=0; iPulse<nPulse; iPulse++){

            // Extract pulse center.
            const as3double x0 = C0x[iPulse];
            const as3double y0 = C0y[iPulse];

            // Radial distance squared.
            const as3double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);

            // Spatial exponential term.
            fexp += exp(-kappa*r2 );
          }

          // Modify current boundary state.
          Var1[3][l] *= ( 1.0 + fexp*a0sinwt );
        }
      }

      break;
    }

    // Sinusoidal spatial velocity profile.
    case(SINUSOIDAL_VELOCITY_BC_MODIFICATION):
    {
      // Abbreviations.
      const as3double gm1   = GAMMA_MINUS_ONE;
      const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;
      const as3double d11   = KronDelta11;
      const as3double d22   = KronDelta22;

      // Extract angular frequency specified.
      auto W0Vec = config_container->GetModifyFreqBC();
      // Extract pulse strength ratio w.r.t. background flow.
      auto A0Vec = config_container->GetModifyStrengthBC();
      // Extract pulse width w.r.t. background flow.
      auto bbVec = config_container->GetModifyWidthBC();
      // Extract pulse(s) center(s).
      auto C0Vec = config_container->GetModifyCenterBC();
      // Extract center shift.
      auto S0Vec = config_container->GetModifyShiftCenterBC();
      // Extract time step.
      auto dt    = config_container->GetTimeStep();

      // Deduce number of pulses.
      unsigned short nPulse = C0Vec.size()/nDim;

      // Extract data pertaining to only the u-velocity variable.
      const as3double omg = W0Vec[1];
      const as3double A0  = A0Vec[1];
      const as3double bb  = bbVec[1];
      const as3double s0  = S0Vec[1];

      // Extract coordinates of the pulse center(s).
      as3vector1d<as3double> C0x(nPulse);
      as3vector1d<as3double> C0y(nPulse);
      for(unsigned short iPulse=0; iPulse<nPulse; iPulse++){
        C0x[iPulse] = C0Vec[iPulse*nDim+XDIM];
        C0y[iPulse] = C0Vec[iPulse*nDim+YDIM];
      }

      // Compute spatial exponential factor.
      const as3double kappa  = log(2.0)/(bb*bb);
      // Spatial-temporal modification term for pulse center.
      const as3double factor = cos(localTime*omg);
      const as3double c0xfunc = d22*s0*factor;
      const as3double c0yfunc = d11*s0*factor;

      // Initial transient is smooth.
      const as3double texp = fabs(factor)*A0*exp(-1000.0*dt/localTime);

      // Loop over all elements sharing this boundary.
      for(unsigned long i=0; i<ElemIndexI.size(); i++){

        // Abbreviation for the primitive data.
        as3double **Var0 = OrigDOFsIntBoundary[i].data(); // original data.
        as3double **Var1 = DataDOFsIntBoundary[i].data(); // modified data.

        // Abbreviation for the grid data.
        as3double **Coord = GridDOFsIntBoundary[i].data();

        // Loop over integration points and compute the modified boundary state.
#pragma omp simd
        for(unsigned short l=0; l<nDOFsInt1D; l++){

          // Extract grid coordinates.
          const as3double x = Coord[0][l];
          const as3double y = Coord[1][l];

          // Reset modified data to original data.
          Var1[0][l] = Var0[0][l];
          Var1[1][l] = Var0[1][l];
          Var1[2][l] = Var0[2][l];
          Var1[3][l] = Var0[3][l];

          // Compute radial distance, per pulse center.
          as3double fexp = 0.0;
          for(unsigned short iPulse=0; iPulse<nPulse; iPulse++){

            // Extract pulse center.
            const as3double x0 = C0x[iPulse] + c0xfunc;
            const as3double y0 = C0y[iPulse] + c0yfunc;

            // Radial distance squared.
            const as3double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);

            // Spatial exponential term.
            fexp += exp(-kappa*r2 );
          }

          // Modify current boundary state.
          Var1[1][l] *= ( 1.0 + fexp*texp );
        }
      }

      break;
    }

    default:
      Terminate("CEEBoundary::ModifyBoundaryCondition", __FILE__, __LINE__,
                "Unknown TypeModifyBC used.");
  }
}


CEEBoundary::~CEEBoundary
(
 void
)
 /*
	* Destructor for CEEBoundary class, frees allocated memory.
	*/
{
  for(unsigned short i=0; i<DataDOFsIntBoundary.size(); i++)
    for(unsigned short j=0; j<DataDOFsIntBoundary[i].size(); j++)
      if( DataDOFsIntBoundary[i][j] ) delete [] DataDOFsIntBoundary[i][j];

  for(unsigned short i=0; i<OrigDOFsIntBoundary.size(); i++)
    for(unsigned short j=0; j<OrigDOFsIntBoundary[i].size(); j++)
      if( OrigDOFsIntBoundary[i][j] ) delete [] OrigDOFsIntBoundary[i][j];

  for(unsigned short i=0; i<GridDOFsIntBoundary.size(); i++)
    for(unsigned short j=0; j<GridDOFsIntBoundary[i].size(); j++)
      if( GridDOFsIntBoundary[i][j] ) delete [] GridDOFsIntBoundary[i][j];
}


void CEEInterfaceBoundary::ReportOutput
(
  void
)
 /*
  * Function that reports output for monitoring.
  */
{
  // Report output.
  std::cout << "  "
            << "iBoundary(" << boundaryID      << "): "
            << DisplayBoundarySide(boundaryID) << " "
            << "jZone("     << zoneMatchID     << "): "
            << DisplayTypeZone(typeZoneMatch)  << " "
            << "jBoundary(" << boundaryMatchID << "): "
            << DisplayBoundarySide(boundaryMatchID)
            << std::endl;
}


CEEInterfaceBoundary::CEEInterfaceBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEInterfaceBoundary per zone: iZone.
	*/
{
  // Determine matching zone ID.
  zoneMatchID     = config_container->GetInterfaceID(iZone, iBoundary)[0];
  // Determine matching boundary ID.
  boundaryMatchID = config_container->GetInterfaceID(iZone, iBoundary)[1];
  // Determine type of matching zone.
  typeZoneMatch   = config_container->GetTypeZone(zoneMatchID);
  // Determine number of DOFs of the solution points in 1D in the matching zone.
  nDOFsSol1DMatch = element_container[zoneMatchID]->GetnDOFsSol1D();

  // Element indices that share the matching boundary.
  ElemIndexJ = geometry_container->GetElemBoundaryIndex(zoneMatchID, boundaryMatchID);

  // Consistency check.
  assert( ElemIndexI.size() == ElemIndexJ.size() && "Interface elements not of same size!" );

  // Get the basis of the solution used in the matching zone.
  auto& rDOFsSol1D_J = element_container[zoneMatchID]->GetrDOFsSol1D();
  // Get the nodes we are interpolating to, based on this current zone.
  auto& rDOFsInt1D_I = element_container[iZone]->GetrDOFsInt1D();

  // Initialize lagrange interpolation matrix (transposed) in 1D:
  // SolDOFs in the matching zone to IntDOFs in this zone.
  lagrangeIntExt1DTranspose = new as3double[nDOFsSol1DMatch*nDOFsInt1D]();

  // Compute lagrange polynomial that interpolates from the nDOFsSol1D of the
  // matching zone to the nDOFsInt1D of this current zone. Note, it doesn't matter
  // which zone this element container operation is being called from.
  element_container[iZone]->LagrangeBasisFunctions(rDOFsSol1D_J, rDOFsInt1D_I,
                                                   nullptr, nullptr,
  											                           lagrangeIntExt1DTranspose,
                                                   nullptr, nullptr, nullptr);

  // Report output.
  ReportOutput();
}


CEEInterfaceBoundary::~CEEInterfaceBoundary
(
 void
)
 /*
	* Destructor for CEEInterfaceBoundary class, frees allocated memory.
	*/
{
	if( lagrangeIntExt1DTranspose != nullptr ) delete [] lagrangeIntExt1DTranspose;
}


void CEEInterfaceBoundary::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current interface boundary condition.
	*/
{
  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get data of matching zone.
  auto& dataJ = solver_container[zoneMatchID]->GetDataContainer();
  // Get the indicial number for the solution on the matching zone and face.
  auto& FaceIndexJ = element_container[zoneMatchID]->GetIndexDOFsSol(boundaryMatchID);

  // Loop over all elements sharing this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element index.
    unsigned long iElem = ElemIndexI[i];
    // Deduce matching element index.
    unsigned long jElem = ElemIndexJ[i];

    // Extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);
    // Extract matching element's solution.
    auto& dataSolJ = dataJ[jElem]->GetDataDOFsSol();

    // Interpolate matching solution to current element's integration
    // points of the current boundary.
    TensorProductSolAndGradFace(boundaryMatchID, nDOFsInt1D, nVar, nDOFsSol1DMatch,
                                FaceIndexJ.data(), lagrangeIntExt1DTranspose,
                                nullptr, nullptr,
                                dataSolJ.data(), dataIntI.data(),
                                nullptr, nullptr);
  }
}


CEESymmetryBoundary::CEESymmetryBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEESymmetryBoundary per zone: iZone.
	*/
{

}


CEESymmetryBoundary::~CEESymmetryBoundary
(
 void
)
 /*
	* Destructor for CEESymmetryBoundary class, frees allocated memory.
	*/
{

}


void CEESymmetryBoundary::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current symmetry boundary condition.
	*/
{
  // Abbreviations.
  const as3double nxabs = fabs( UnitNormal[XDIM] );
  const as3double nyabs = fabs( UnitNormal[YDIM] );
  const as3double gm1   = GAMMA_MINUS_ONE;
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;

  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get the indicial number for the solution on the face.
  auto& FaceIndexI = element_container[zoneID]->GetIndexDOFsSol(boundaryID);
  // Get the transposed lagrange interpolation function in 1D.
  auto *lagrangeInt1DTranspose = element_container[zoneID]->GetLagrangeInt1DTranspose();

  // Loop over all elements sharing this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);
    // Extract current element's solution.
    auto& dataSolI = dataI[iElem]->GetDataDOFsSol();

    // Interpolate solution to integration points of the current boundary.
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(), lagrangeInt1DTranspose,
                                nullptr, nullptr,
                                dataSolI.data(), dataIntI.data(),
                                nullptr, nullptr);

    // Abbreviation for consiceness.
    as3double **Var = dataIntI.data();

    // Loop over integration points and compute the symmetry variables.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Extract primitive data.
      as3double rho   = Var[0][l];
      as3double ovrho = 1.0/rho;
      as3double u     = ovrho*Var[1][l];
      as3double v     = ovrho*Var[2][l];
      as3double p     = gm1*(Var[3][l]
                      - 0.5*(u*Var[1][l] + v*Var[2][l]) );

      // Compute symmetry/slip condition on primitive variables.
      u *= ( 1.0 - 2.0*nxabs );
      v *= ( 1.0 - 2.0*nyabs );

      // Assemble conservative data.
      Var[1][l] = rho*u;
      Var[2][l] = rho*v;
      Var[3][l] = p*ovgm1 + 0.5*rho*( u*u + v*v );
    }
  }
}


CEEStaticOutletBoundary::CEEStaticOutletBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEStaticOutletBoundary per zone: iZone.
	*/
{
  // Initialize and compute the prescribed boundary data.
  InitializePrescribedState(config_container,
                            geometry_container,
                            element_container[iZone],
                            initial_container,
                            iZone, iBoundary);
}


CEEStaticOutletBoundary::~CEEStaticOutletBoundary
(
 void
)
 /*
	* Destructor for CEEStaticOutletBoundary class, frees allocated memory.
	*/
{

}


void CEEStaticOutletBoundary::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current subsonic static outlet boundary condition.
	*/
{
  // Abbreviations.
  const as3double nx    = UnitNormal[XDIM];
  const as3double ny    = UnitNormal[YDIM];
  const as3double gm1   = GAMMA_MINUS_ONE;
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;
  const as3double ovg   = 1.0/GAMMA;

  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get the indicial number for the solution on the face.
  auto& FaceIndexI = element_container[zoneID]->GetIndexDOFsSol(boundaryID);
  // Get the transposed lagrange interpolation function in 1D.
  auto *lagrangeInt1DTranspose = element_container[zoneID]->GetLagrangeInt1DTranspose();

  // Loop over all elements sharing this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);
    // Extract current element's solution.
    auto& dataSolI = dataI[iElem]->GetDataDOFsSol();

    // Interpolate solution to integration points of the current boundary.
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(), lagrangeInt1DTranspose,
                                nullptr, nullptr,
                                dataSolI.data(), dataIntI.data(),
                                nullptr, nullptr);

    // Abbreviation for consiceness.
    as3double **Var = dataIntI.data();

    // Loop over integration points and compute the symmetry variables.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Extract primitive data.
      const as3double rho   = Var[0][l];
      const as3double ovrho = 1.0/rho;
      const as3double u     = ovrho*Var[1][l];
      const as3double v     = ovrho*Var[2][l];
      const as3double p     = gm1*(Var[3][l]
                            - 0.5*(u*Var[1][l] + v*Var[2][l]) );

      // Compute normal velocity component.
      const as3double un = u*nx + v*ny;

      // Compute speed of sound and its square.
      const as3double a2 = GAMMA*p*ovrho;
      const as3double a  = sqrt(a2);

      // Compute the entropy.
      const as3double s = p*pow(ovrho, GAMMA);
      // Compute the Riemann invariant to be extrapolated.
      const as3double riemann = 2.0*a*ovgm1 + un;

      // Extract reference data needed.
      const as3double pInf = DataDOFsIntBoundary[i][3][l];

      // Compute the extrapolated state variables.
      const as3double rhoInf = pow(pInf/s, ovg);
      const as3double aInf   = sqrt(GAMMA*pInf/rhoInf);
      const as3double unInf  = riemann - 2.0*aInf*ovgm1;

      // Compute the velocity of the extrapolated state.
      const as3double uInf = u + (unInf - un)*nx;
      const as3double vInf = v + (unInf - un)*ny;

      // Compute kinetic energy, based on prescribed conditions.
      const as3double ekInf = 0.5*( uInf*uInf + vInf*vInf );

      // Assemble conservative data.
      Var[0][l] = rhoInf;
      Var[1][l] = rhoInf*uInf;
      Var[2][l] = rhoInf*vInf;
      Var[3][l] = pInf*ovgm1 + rhoInf*ekInf;
    }
  }
}


CEESupersonicOutletBoundary::CEESupersonicOutletBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEESupersonicOutletBoundary per zone: iZone.
	*/
{
  // Initialize and compute the prescribed boundary data.
  InitializePrescribedState(config_container,
                            geometry_container,
                            element_container[iZone],
                            initial_container,
                            iZone, iBoundary);
}


CEESupersonicOutletBoundary::~CEESupersonicOutletBoundary
(
 void
)
 /*
	* Destructor for CEESupersonicOutletBoundary class, frees allocated memory.
	*/
{

}


void CEESupersonicOutletBoundary::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current supersonic outlet boundary condition.
	*/
{
  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get the indicial number for the solution on the face.
  auto& FaceIndexI = element_container[zoneID]->GetIndexDOFsSol(boundaryID);
  // Get the transposed lagrange interpolation function in 1D.
  auto *lagrangeInt1DTranspose = element_container[zoneID]->GetLagrangeInt1DTranspose();

  // Loop over all elements sharing this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);
    // Extract current element's solution.
    auto& dataSolI = dataI[iElem]->GetDataDOFsSol();

    // Interpolate solution to integration points of the current boundary. Note,
    // this is the boundary condition, since we are taking the interiot-state as
    // the boundary -- which is what we need for a supersonic outlet.
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(), lagrangeInt1DTranspose,
                                nullptr, nullptr,
                                dataSolI.data(), dataIntI.data(),
                                nullptr, nullptr);
  }
}


CEEStaticInletBoundary::CEEStaticInletBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEStaticInletBoundary per zone: iZone.
	*/
{
  // Initialize and compute the prescribed boundary data.
  InitializePrescribedState(config_container,
                            geometry_container,
                            element_container[iZone],
                            initial_container,
                            iZone, iBoundary);

  // If a varying BC profile is specified, initialize required coordinates and extra data.
  if( config_container->GetTypeModifyBC(iBoundary) != NO_BC_MODIFICATION )
    InitializeModifiedBC(config_container,
                         geometry_container,
                         element_container[iZone],
                         initial_container,
                         iZone, iBoundary);
}


CEEStaticInletBoundary::~CEEStaticInletBoundary
(
 void
)
 /*
	* Destructor for CEEStaticInletBoundary class, frees allocated memory.
	*/
{

}


void CEEStaticInletBoundary::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current subsonic static inlet boundary condition.
	*/
{
  // Abbreviations.
  const as3double nx    = UnitNormal[XDIM];
  const as3double ny    = UnitNormal[YDIM];
  const as3double gm1   = GAMMA_MINUS_ONE;
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;
  const as3double ovg   = 1.0/GAMMA;

  // If this boundary condition needs modification, do so.
  if( OrigDOFsIntBoundary.size() )
    ModifyBoundaryCondition(config_container,
                            geometry_container,
                            solver_container[zoneID],
                            element_container[zoneID],
                            spatial_container[zoneID],
                            localTime);

  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get the indicial number for the solution on the face.
  auto& FaceIndexI = element_container[zoneID]->GetIndexDOFsSol(boundaryID);
  // Get the transposed lagrange interpolation function in 1D.
  auto *lagrangeInt1DTranspose = element_container[zoneID]->GetLagrangeInt1DTranspose();

  // Loop over all elements sharing this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);
    // Extract current element's solution.
    auto& dataSolI = dataI[iElem]->GetDataDOFsSol();

    // Interpolate solution to integration points of the current boundary.
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(), lagrangeInt1DTranspose,
                                nullptr, nullptr,
                                dataSolI.data(), dataIntI.data(),
                                nullptr, nullptr);

    // Abbreviation for consiceness.
    as3double **Var = dataIntI.data();

    // Loop over integration points and compute the symmetry variables.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Extract primitive data.
      const as3double rho   = Var[0][l];
      const as3double ovrho = 1.0/rho;
      const as3double u     = ovrho*Var[1][l];
      const as3double v     = ovrho*Var[2][l];
      const as3double p     = gm1*(Var[3][l]
                            - 0.5*(u*Var[1][l] + v*Var[2][l]) );

      // Compute normal velocity component.
      const as3double un = u*nx + v*ny;

      // Compute speed of sound and its square.
      const as3double a2 = GAMMA*p*ovrho;
      const as3double a  = sqrt(a2);

      // Compute the Riemann invariant to be extrapolated.
      const as3double riemann = 2.0*a*ovgm1 + un;

      // Extract reference data needed.
      const as3double rhoInf = DataDOFsIntBoundary[i][0][l];
      const as3double uInf   = DataDOFsIntBoundary[i][1][l];
      const as3double vInf   = DataDOFsIntBoundary[i][2][l];

      // Compute kinetic energy, based on prescribed conditions.
      const as3double ekInf = 0.5*( uInf*uInf + vInf*vInf );

      // Compute the normal velocity based on the prescribed conditions.
      const as3double unInf = uInf*nx + vInf*ny;

      // Compute the speed of sound, based on the riemann extract.
      const as3double aInf = 0.5*gm1*( riemann - unInf );

      // Compute the extrapolated pressure.
      const as3double pInf = ovg*rhoInf*aInf*aInf;

      // Assemble conservative data.
      Var[0][l] = rhoInf;
      Var[1][l] = rhoInf*uInf;
      Var[2][l] = rhoInf*vInf;
      Var[3][l] = pInf*ovgm1 + rhoInf*ekInf;
    }
  }
}


CEETotalInletBoundary::CEETotalInletBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEETotalInletBoundary per zone: iZone.
	*/
{
  // Initialize and compute the prescribed boundary data.
  InitializePrescribedState(config_container,
                            geometry_container,
                            element_container[iZone],
                            initial_container,
                            iZone, iBoundary);

  // If a varying BC profile is specified, initialize required coordinates and extra data.
  if( config_container->GetTypeModifyBC(iBoundary) != NO_BC_MODIFICATION )
    InitializeModifiedBC(config_container,
                         geometry_container,
                         element_container[iZone],
                         initial_container,
                         iZone, iBoundary);

  // Extract reference data needed.
  const as3double uInf = initial_container->GetVelocityNormal();
  const as3double vInf = initial_container->GetVelocityTransverse();

  // Determine velocity magnitude and its inverse.
  const as3double umag   = sqrt( uInf*uInf + vInf*vInf );
  const as3double ovumag = 1.0/umag;

  // Direction of the prescribed velocities normalized.
  udir = uInf*ovumag;
  vdir = vInf*ovumag;

  // Abbreviation for the normals.
  const as3double nx = UnitNormal[XDIM];
  const as3double ny = UnitNormal[YDIM];

  // Compute the dot product between the normal and the
  // velocity direction. This value should be negative.
  alpha = nx*udir + ny*vdir;
}


CEETotalInletBoundary::~CEETotalInletBoundary
(
 void
)
 /*
	* Destructor for CEETotalInletBoundary class, frees allocated memory.
	*/
{

}


void CEETotalInletBoundary::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current subsonic total-condition inlet boundary condition.
	*/
{
  // Abbreviations.
  const as3double nx     = UnitNormal[XDIM];
  const as3double ny     = UnitNormal[YDIM];
  const as3double gm1    = GAMMA_MINUS_ONE;
  const as3double ovgm1  = 1.0/GAMMA_MINUS_ONE;
  const as3double govgm1 = GAMMA*ovgm1;
  const as3double ovrg   = 1.0/GAS_CONSTANT;

  // If this boundary condition needs modification, do so.
  if( OrigDOFsIntBoundary.size() )
    ModifyBoundaryCondition(config_container,
                            geometry_container,
                            solver_container[zoneID],
                            element_container[zoneID],
                            spatial_container[zoneID],
                            localTime);

  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get the indicial number for the solution on the face.
  auto& FaceIndexI = element_container[zoneID]->GetIndexDOFsSol(boundaryID);
  // Get the transposed lagrange interpolation function in 1D.
  auto *lagrangeInt1DTranspose = element_container[zoneID]->GetLagrangeInt1DTranspose();

  // Loop over all elements sharing this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);
    // Extract current element's solution.
    auto& dataSolI = dataI[iElem]->GetDataDOFsSol();

    // Interpolate solution to integration points of the current boundary.
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(), lagrangeInt1DTranspose,
                                nullptr, nullptr,
                                dataSolI.data(), dataIntI.data(),
                                nullptr, nullptr);

    // Abbreviation for consiceness.
    as3double **Var = dataIntI.data();

    // Loop over integration points and compute the symmetry variables.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Extract reference data needed.
      const as3double R0   = DataDOFsIntBoundary[i][0][l];
      const as3double U0   = DataDOFsIntBoundary[i][1][l];
      const as3double V0   = DataDOFsIntBoundary[i][2][l];
      const as3double P0   = DataDOFsIntBoundary[i][3][l];

      // Abbreviation: 1/rho.
      const as3double OVR0 = 1.0/R0;
      // Deduce reference temperature.
      const as3double T0   = P0*ovrg/R0;
      // Deduce reference magnitude of the velocity squared.
      const as3double UM02 = U0*U0 + V0*V0;;
      // Deduce reference speed-of-sound squared.
      const as3double C02  = GAMMA*P0*OVR0;
      // Deduce reference Mach number squared.
      const as3double M02  = UM02/C02;

      // Abbreviation.
      const as3double gg   = 1.0 + 0.5*gm1*M02;

      // Compute total pressure, temperature and enthalpy needed.
      const as3double Ptot = P0*pow(gg, govgm1);
      const as3double Ttot = T0*gg;
      const as3double Htot = CP_CONSTANT*Ttot;

      // Extract primitive data.
      const as3double rho   = Var[0][l];
      const as3double ovrho = 1.0/rho;
      const as3double u     = ovrho*Var[1][l];
      const as3double v     = ovrho*Var[2][l];
      const as3double p     = gm1*(Var[3][l]
                            - 0.5*(u*Var[1][l] + v*Var[2][l]) );

      // Compute normal velocity component.
      const as3double un = u*nx + v*ny;

      // Compute speed of sound and its square.
      const as3double a2 = GAMMA*p*ovrho;
      const as3double a  = sqrt(a2);

      // Compute the Riemann invariant to be extrapolated.
      as3double riemann = 2.0*a*ovgm1 + un;

      // Apply total enthalpy scaling to increase stability.
      // If this is not desired, comment the following lines.
      const as3double H = ovrho*(Var[3][l] + p);
      riemann          *= sqrt(Htot/H);

      // Coefficients in the quadratic equation for the magnitude of the velocity.
      const as3double aa =  1.0 + 0.5*gm1*alpha*alpha;
      const as3double bb = -gm1*alpha*riemann;
      const as3double cc =  0.5*gm1*riemann*riemann - 2.0*Htot;

      // Solve the equation for the magnitude of the velocity. As this value
      // must be positive and both aa and bb are positive (alpha is negative and
      // riemann is positive up till Mach = 5.0 or so, which is not really subsonic
      // anymore), it is clear which of the two possible solutions must be taken.
      // Some clipping is present, but this is normally not active.
      as3double dd = bb*bb - 4.0*aa*cc;   dd = sqrt(std::max(0.0, dd));
      as3double qR = (-bb + dd)/(2.0*aa); qR = std::max(0.0, qR);

      // Compute the square of the Mach number and clip it between 0 and 1,
      // because this is a subsonic inflow boundary.
      as3double qR2 = qR*qR;
      as3double aR2 = gm1*(Htot - 0.5*qR2);
      as3double MR2 = qR2/aR2; MR2 = std::min(1.0, MR2);

      // Compute the final value of the magnitude of the velocity.
      const as3double tmp  = 1.0/(1.0 + 0.5*gm1*MR2);
      aR2 = gm1*Htot*tmp;
      qR2 = MR2*aR2;
      qR  = sqrt(qR2);

      // Compute the pressure from the prescribed total pressure
      // and the temperature from the pressure and total temperature.
      const as3double pR   = Ptot*pow(tmp, govgm1);
      const as3double rhoR = GAMMA*pR/aR2;

      // Assemble conservative data.
      Var[0][l] = rhoR;
      Var[1][l] = rhoR*qR*udir;
      Var[2][l] = rhoR*qR*vdir;
      Var[3][l] = pR*ovgm1 + 0.5*rhoR*qR2;
    }
  }
}


CEESupersonicInletBoundary::CEESupersonicInletBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEESupersonicInletBoundary per zone: iZone.
	*/
{
  // Initialize and compute the prescribed boundary data.
  InitializePrescribedState(config_container,
                            geometry_container,
                            element_container[iZone],
                            initial_container,
                            iZone, iBoundary);

  // If a varying BC profile is specified, initialize required coordinates and extra data.
  if( config_container->GetTypeModifyBC(iBoundary) != NO_BC_MODIFICATION )
    InitializeModifiedBC(config_container,
                         geometry_container,
                         element_container[iZone],
                         initial_container,
                         iZone, iBoundary);
}


CEESupersonicInletBoundary::~CEESupersonicInletBoundary
(
 void
)
 /*
	* Destructor for CEESupersonicInletBoundary class, frees allocated memory.
	*/
{

}


void CEESupersonicInletBoundary::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current supersonic inlet boundary condition.
	*/
{
  // Abbreviations.
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;

  // If this boundary condition needs modification, do so.
  if( OrigDOFsIntBoundary.size() )
    ModifyBoundaryCondition(config_container,
                            geometry_container,
                            solver_container[zoneID],
                            element_container[zoneID],
                            spatial_container[zoneID],
                            localTime);

  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();

  // Loop over all elements sharing this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element integration points of the current boundary.
    auto& Var = dataI[iElem]->GetDataDOFsIntFace(boundaryID);

    // Loop over integration points and compute the symmetry variables.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Extract reference data needed.
      const as3double rhoInf = DataDOFsIntBoundary[i][0][l];
      const as3double uInf   = DataDOFsIntBoundary[i][1][l];
      const as3double vInf   = DataDOFsIntBoundary[i][2][l];
      const as3double pInf   = DataDOFsIntBoundary[i][3][l];

      // Kinetic energy.
      const as3double ekInf  = 0.5*( uInf*uInf + vInf*vInf );

      // Assemble conservative data.
      Var[0][l] = rhoInf;
      Var[1][l] = rhoInf*uInf;
      Var[2][l] = rhoInf*vInf;
      Var[3][l] = pInf*ovgm1 + rhoInf*ekInf;
    }
  }
}


CEECharacteristicBoundary::CEECharacteristicBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEECharacteristicBoundary per zone: iZone.
	*/
{
  // Determine the index of the incoming wave.
  if( (iBoundary == IDX_EAST) || (iBoundary == IDX_NORTH) ) PsiIndex = 0;
  else                                                      PsiIndex = 3;

  // Make distinction between the normal and transverse direction.
  if( (iBoundary == IDX_WEST) || (iBoundary == IDX_EAST) ){
    // This is a face whose normal is aligned with the x-direction (i.e. ny=0).
    IndexNormal = 1;   IndexTransverse = 2;
  }
  else {
    // This is a face whose normal is aligned with the y-direction (i.e. nx=0).
    IndexNormal = 2;   IndexTransverse = 1;
  }

  // Reserve memory for the working array.
  WorkingDataInt1D.resize(nData);

  // Reserve memory per data of the working array.
  for(unsigned short iData=0; iData<nData; iData++){

    // Reserve memory for each variable.
    WorkingDataInt1D[iData].resize(nVar, nullptr);

    // Reserve memory per variable.
    for(unsigned short iVar=0; iVar<nVar; iVar++){

      // Allocate actual memory.
      WorkingDataInt1D[iData][iVar] = new as3double[nDOFsInt1D]();

      // Check if allocation failed.
      if( !WorkingDataInt1D[iData][iVar] )
        Terminate("CEECharacteristicBoundary::CEECharacteristicBoundary", __FILE__, __LINE__,
                  "Allocation failed for WorkingDataInt1D.");
    }
  }

  // Get the relevant differentiation expansion acting on the current boudnary.
  auto* dellFace = element_container[iZone]->GetDerLagrangeDOFsSol1DFace(iBoundary);
  // Get the lagrange interpolation function in 1D.
  auto *lagrangeInt1D  = element_container[iZone]->GetLagrangeInt1D();
  // Get the transposed lagrange interpolation function in 1D.
  auto *lagrangeInt1DTranspose = element_container[iZone]->GetLagrangeInt1DTranspose();

  // Determine the normal-gradient-based coefficient.
  Coef_dell = ( PsiIndex == 0 ) ? dellFace[nDOFsSol1D-1] : dellFace[0];

  // Reserve memory for the least-squares matrix used in the PC-formulation.
  MatrixC = new as3double[nDOFsInt1D*nDOFsInt1D]();

  // Compute the least-square matrix.
  ComputeLeastSquaresMatrix(lagrangeInt1D,
                            lagrangeInt1DTranspose,
                            Coef_dell, MatrixC);

  // Initialize and compute the prescribed boundary data, aka target state.
  InitializePrescribedState(config_container,
                            geometry_container,
                            element_container[iZone],
                            initial_container,
                            iZone, iBoundary);
}


CEECharacteristicBoundary::~CEECharacteristicBoundary
(
 void
)
 /*
	* Destructor for CEECharacteristicBoundary class, frees allocated memory.
	*/
{
  if( MatrixC != nullptr ) delete [] MatrixC;

  for(unsigned short i=0; i<WorkingDataInt1D.size(); i++)
    for(unsigned short j=0; j<WorkingDataInt1D[i].size(); j++)
      if( WorkingDataInt1D[i][j] ) delete [] WorkingDataInt1D[i][j];
}


void CEECharacteristicBoundary::ComputeLeastSquaresMatrix
(
  const as3double *lagrangeInt1D,
  const as3double *lagrangeInt1DTranspose,
  const as3double  dell,
  as3double       *MatrixLeastSquares
)
 /*
  * Function that computes the least-squares matrix.
  */
{
  // Temporary matrices.
  as3vector1d<as3double> tmp1(nDOFsSol1D*nDOFsSol1D, 0.0);
  as3vector1d<as3double> tmp2(nDOFsSol1D*nDOFsInt1D, 0.0);

  // Abbreviation. Note, the actual metric of the normal Jacobian is included
  // per element multiplication, as it can vary depending on the element size.
  const as3double ovdell = 1.0/dell;

  // Step 1: compute ell^T * ell.
  lgemm(nDOFsSol1D, nDOFsSol1D, nDOFsInt1D,
        lagrangeInt1DTranspose, lagrangeInt1D,
        tmp1.data());

  // Step 2: invert the matrix.
  ComputeInverseMatrix(nDOFsSol1D, tmp1.data());

  // Step 3: multiply from the right by ell^T.
  lgemm(nDOFsSol1D, nDOFsInt1D, nDOFsSol1D,
        tmp1.data(), lagrangeInt1DTranspose,
        tmp2.data());

  // Step 4: multiply from the left by ell and scale by ovdell.
  lgemm(nDOFsInt1D, nDOFsInt1D, nDOFsSol1D, ovdell,
        lagrangeInt1D, tmp2.data(),
        MatrixLeastSquares);
}


as3double CEECharacteristicBoundary::ComputeAverageMachLocal
(
  const as3vector1d<as3double> &weights,
  as3double                   **Var
)
 /*
  * Function that computes the average of the Mach number on the boundary of an
  * element lying on the current boundary.
  */
{
  // Initialize Mach number.
  as3double Mavg = 0.0;

  // Loop over all integration points on the surface.
#pragma omp simd
  for(unsigned short l=0; l<nDOFsInt1D; l++){

    // Extract primitive data.
    const as3double rho   = Var[0][l];
    const as3double ovrho = 1.0/rho;
    const as3double u     = ovrho*Var[1][l];
    const as3double v     = ovrho*Var[2][l];
    const as3double p     = GAMMA_MINUS_ONE*(Var[3][l]
                          - 0.5*(u*Var[1][l] + v*Var[2][l]) );

    // Magnitude of the velocity.
    const as3double umag = sqrt( u*u + v*v );
    // Speed of sound.
    const as3double a = sqrt(GAMMA*p*ovrho);

    // Accumulate the averaged Mach number.
    Mavg += weights[l]*umag/a;
  }

  // Normalize by the sum of the weights (i.e. 0.5) to obtain the average.
  Mavg *= 0.5;

  // Return averaged value.
  return Mavg;
}


as3double CEECharacteristicBoundary::ComputeAverageMachGlobal
(
  const CGeometryZone               *geometry_zone,
  const as3element                  &data_container,
  const as3vector1d<as3double>      &weights,
  const as3vector1d<unsigned short> &FaceIndexI,
  const as3double                   *ellT,
  as3double                        **Var
)
 /*
  * Function that computes the average of the Mach number on the entire boundary.
  */
{
  // Initialize averaged Mach number.
  as3double Mavg = 0.0;
  // Initialize surface area of boundary face.
  as3double area = 0.0;

  // Abbreviations.
  const as3double d11 = KronDelta11;
  const as3double d22 = KronDelta22;

  // Loop over all elements to compute the average.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element size.
    auto ElemSize = geometry_zone->GetGeometryElem(iElem)->GetElemSize();
    // Extract element sizes explicitly.
    const as3double hx = ElemSize[0];
    const as3double hy = ElemSize[1];
    // Deduce the boundary surface area of the current element.
    const as3double ds = d11*hy + d22*hx;

    // Extract current element integration points of the current boundary.
    auto& dataSolI = data_container[iElem]->GetDataDOFsSol();

    // Interpolate the data at the integration points.
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(),
                                ellT, nullptr, nullptr,
                                dataSolI.data(), Var,
                                nullptr, nullptr);

    // Compute local-averaged Mach number.
    const as3double Mach = ComputeAverageMachLocal(weights, Var);

    // Add contribution to global Mach number average.
    Mavg += Mach*ds;
    // Accumulate total boundary area.
    area += ds;
  }

  // Normalize by entire surface area of boundary.
  Mavg /= area;

  // Return averaged value.
  return Mavg;
}


CEEOutletCBC::CEEOutletCBC
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEECharacteristicBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEOutletCBC per zone: iZone.
	*/
{
  // Extract tuning paramers.
  auto& ParamTuning = config_container->GetParamOutletNSCBC();

  // Extract the normalized length scale.
  const as3double lenInv = 1.0/ParamTuning[3];
  // Transverse safety-factor coefficient.
  Coef_eta       = ParamTuning[4];
  // Check if the coupled transverse terms are adaptive or not.
  AdaptiveBeta_l = ( ParamTuning[1] < 0 ) ? true : false;
  // Check if the uncoupled transverse terms are adaptive or not.
  AdaptiveBeta_t = ( ParamTuning[2] < 0 ) ? true : false;

  // Coupled transverse term relaxation coefficient.
  beta_l  = Coef_eta*ParamTuning[1];
  // Uncoupled transverse term relaxation coefficient.
  beta_t  = Coef_eta*ParamTuning[2];
  // Normal term relaxation coefficient.
  coefK   = 0.5*lenInv*ParamTuning[0];
}


CEEOutletCBC::~CEEOutletCBC
(
 void
)
 /*
	* Destructor for CEEOutletCBC class, frees allocated memory.
	*/
{

}


void CEEOutletCBC::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current outlet CBC boundary condition.
	*/
{
  // Some abbreviations.
  const as3double gm1   = GAMMA_MINUS_ONE;
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;

  // Geometry in current zone.
  auto* geometry_zone = geometry_container->GetGeometryZone(zoneID);

  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get the indicial number for the solution in each element sharing this face.
  auto& FaceIndexI = element_container[zoneID]->GetIndexDOFsSol(boundaryID);
  // Get the 1D lagrange interpolation operator.
  auto* ellT    = element_container[zoneID]->GetLagrangeInt1DTranspose();
  // Get the 1D differentiation lagrange operator.
  auto* dellT   = element_container[zoneID]->GetDerLagrangeInt1DTranspose();
  // Get the 1D surface differentiation operator.
  auto* dellS   = element_container[zoneID]->GetDerLagrangeDOFsSol1DFace(boundaryID);
  // Get the integration weights in 1D.
  auto& weights = element_container[zoneID]->GetwDOFsInt1D();

  // Abbreviations.
  const as3double d11  = KronDelta11;
  const as3double d22  = KronDelta22;

  // Assign the correct indices of the working array w.r.t. the boundary orientation.
  as3double **Var    = WorkingDataInt1D[0].data();
  as3double **dVarDn = WorkingDataInt1D[IndexNormal].data();
  as3double **dVarDt = WorkingDataInt1D[IndexTransverse].data();

  // Compute the average mach number on all element faces on this boundary.
  const as3double Mavg = ComputeAverageMachGlobal(geometry_zone, dataI,
                                                  weights, FaceIndexI, ellT, Var);
  // Abbreviation involving Mavg.
  const as3double omM2 = 1.0 - Mavg*Mavg;
  // Override relaxation term specification, if specified.
  if( AdaptiveBeta_l ) beta_l = Coef_eta*Mavg;
  if( AdaptiveBeta_t ) beta_t = Coef_eta*Mavg;

  // Loop over all elements on this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element size.
    auto ElemSize = geometry_zone->GetGeometryElem(iElem)->GetElemSize();
    // Extract element sizes explicitly.
    const as3double hx = ElemSize[0];
    const as3double hy = ElemSize[1];
    // Compute Jacobians on this element in normal and transverse terms.
    const as3double Jn = d11*(2.0/hx) + d22*(2.0/hy);
    const as3double Jt = d11*(2.0/hy) + d22*(2.0/hx);

    // Abbreviation for inverse of normal Jacobian component.
    const as3double ovJn = 1.0/Jn;
    // Include the normal gradient with metric.
    const as3double dell = Coef_dell*Jn;

    // Step 0a: extract current element integration points of the current boundary.
    auto& dataSolI = dataI[iElem]->GetDataDOFsSol();
    // Step 0b: extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);

    // Step 1a: compute the data and the gradients at the integration points.
    // Note, these are in parametric coordinates.
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(),
                                ellT, dellT, dellS,
                                dataSolI.data(), WorkingDataInt1D[0].data(),
                                WorkingDataInt1D[1].data(),
                                WorkingDataInt1D[2].data());

    // Step 1b: convert the parametric gradients into physical space.
    for(unsigned short iVar=0; iVar<nVar; iVar++){
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt1D; l++){
        dVarDn[iVar][l] *= Jn;
        dVarDt[iVar][l] *= Jt;
      }
    }

    // Loop over all the integration points living on this element face.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Step 2: compute the primitive variables.
      const as3double rho   = Var[0][l];
      const as3double ovrho = 1.0/rho;
      const as3double u     = ovrho*Var[1][l];
      const as3double v     = ovrho*Var[2][l];
      const as3double p     = gm1*(Var[3][l]
                            - 0.5*(u*Var[1][l] + v*Var[2][l]) );

    	// Kinetic energy.
    	const as3double ek = 0.5*(u*u + v*v);

      // Step 3a: compute the primitive gradient from conservative in normal direction.
      const as3double drdn = dVarDn[0][l];
      const as3double dudn = ovrho*( dVarDn[1][l] - u*drdn );
    	const as3double dvdn = ovrho*( dVarDn[2][l] - v*drdn );
    	const as3double dpdn = gm1*( ek*dVarDn[0][l]
      										 -  		  u*dVarDn[1][l]
      										 -        v*dVarDn[2][l]
      										 +          dVarDn[3][l] );

      // Step 3b: compute the primitive gradient from conservative in transverse direction.
      const as3double drdt = dVarDt[0][l];
      const as3double dudt = ovrho*( dVarDt[1][l] - u*drdt );
      const as3double dvdt = ovrho*( dVarDt[2][l] - v*drdt );
      const as3double dpdt = gm1*( ek*dVarDt[0][l]
                           -  		  u*dVarDt[1][l]
                           -        v*dVarDt[2][l]
                           +          dVarDt[3][l] );

      // Compute the local speed of sound and its square.
      const as3double a2 = GAMMA*p*ovrho;
      const as3double a  = sqrt(a2);

      // Some abbreviations.
      const as3double ova    = 1.0/a;
      const as3double ova2   = ova*ova;
      const as3double ovrhoa = ovrho*ova;

      // Reserve memory for the wave amplitudes.
      as3double vecL[nVar], vecLt[nVar], vecT[nVar];

      // Select the velocity component in the normal direction.
      const as3double un = d11*u + d22*v;
      // Select the velocity component in the transverse direction.
      const as3double ut = d22*u + d11*v;
      // Select the derivative w.r.t. normal direction of the normal velocity component.
      const as3double dundn = d11*dudn + d22*dvdn;
      // Select the derivative w.r.t. transverse direction of the normal velocity component.
      const as3double dundt = d11*dudt + d22*dvdt;
      // Select the derivative w.r.t. normal direction of the transverse velocity component.
      const as3double dutdn = d11*dvdn + d22*dudn;
      // Select the derivative w.r.t. transverse direction of the transverse velocity component.
      const as3double dutdt = d11*dvdt + d22*dudt;
      // Compute entropy change in the normal direction.
      const as3double dsdn  = drdn - ova2*dpdn;
      // Compute entropy change in the transverse direction.
      const as3double dsdt  = drdt - ova2*dpdt;

      // Compute the eigenvalues in the normal direction.
      const as3double lmb1 = un - a;
      const as3double lmb3 = un + a;

      // Step 4a: compute the normal wave-amplitudes.
      // Note, the lmb2 is removed here and ovlmb2 in the reconstruction
      // in order to avoid problems with division by zero when un = 0.
      vecL[0] = 0.5*lmb1*( ovrhoa*dpdn - dundn  );
      vecL[1] =          ( d11*dsdn + d22*dutdn );
      vecL[2] =          ( d11*dutdn + d22*dsdn );
      vecL[3] = 0.5*lmb3*( ovrhoa*dpdn + dundn  );

      // Step 4b: compute the coupled transverse wave-amplitudes.
      vecLt[0] = 0.5*ut*( ovrhoa*dpdt - dundt  );
      vecLt[1] =     ut*( d11*dsdt + d22*dutdt );
      vecLt[2] =     ut*( d11*dutdt + d22*dsdt );
      vecLt[3] = 0.5*ut*( ovrhoa*dpdt + dundt  );

      // Step 4c: compute the uncoupled transverse wave-amplitudes.
      vecT[0] =    -0.5*a*dutdt;
      vecT[1] = -d22*ovrho*dpdt;
      vecT[2] = -d11*ovrho*dpdt;
      vecT[3] =         vecT[0];

      // Extract reference data needed.
      const as3double pInf = DataDOFsIntBoundary[i][3][l];

      // Step 5: apply boundary condition through incoming acoustic wave-amplitude.
      vecL[PsiIndex] = coefK*ovrho*omM2*( p - pInf )
                     - beta_l*vecLt[PsiIndex]
                     + beta_t*vecT[PsiIndex];

      // More abbvreviations.
      const as3double rhoa   = rho*a;
      const as3double rhoova = rho*ova;
      const as3double ovlmb1 = 1.0/lmb1;
      const as3double ovlmb3 = 1.0/lmb3;

      // Step 6: compute the primitive gradient from the normal wave-amplitude.
      // Store this value inside the Var entry in the working array.

      // Note, multiplication by ovlmb2 is completely avoided, since it is
      // implicitly done by avoiding a multiplication of lmb2 in the first place.

      // First entry of the BC-imposed primitive gradient: drdn.
      Var[0][l] =  ovlmb1*rhoova*vecL[0]
                +         (  d11*vecL[1]
                +            d22*vecL[2] )
                +  ovlmb3*rhoova*vecL[3];

      // Second entry of the BC-imposed primitive gradient: dudn.
      Var[1][l] = -ovlmb1*d11*vecL[0]
                +         d22*vecL[1]
                +  ovlmb3*d11*vecL[3];

      // Third entry of the BC-imposed primitive gradient: dvdn.
      Var[2][l] = -ovlmb1*d22*vecL[0]
                +         d11*vecL[2]
                +  ovlmb3*d22*vecL[3];

      // Fourth entry of the BC-imposed primitive gradient: dpdn.
      Var[3][l] =  ovlmb1*rhoa*vecL[0]
                +  ovlmb3*rhoa*vecL[3];

      // Step 7: compute the B matrix entry.
      Var[0][l] -= ( drdn - rho*dell );
      Var[1][l] -= ( dudn -   u*dell );
      Var[2][l] -= ( dvdn -   v*dell );
      Var[3][l] -= ( dpdn -   p*dell );
    }

    // Step 8: compute the boundary state from the B matrix using a least-squares approach.
    for(unsigned short iVar=0; iVar<nVar; iVar++)
      lgemv(nDOFsInt1D, nDOFsInt1D, ovJn, MatrixC, Var[iVar], dataIntI[iVar]);

    // Step 9: convert back into the conservative variables.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Extract primitive variables.
      const as3double rho = dataIntI[0][l];
      const as3double u   = dataIntI[1][l];
      const as3double v   = dataIntI[2][l];
      const as3double p   = dataIntI[3][l];

      // Compute the conservative variables.
      dataIntI[1][l] = rho*u;
      dataIntI[2][l] = rho*v;
      dataIntI[3][l] = p*ovgm1 + 0.5*rho*( u*u + v*v );
    }
  }
}


CEEInletCBC::CEEInletCBC
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEECharacteristicBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEInletCBC per zone: iZone.
	*/
{
  // If a varying BC profile is specified, initialize required coordinates and extra data.
  if( config_container->GetTypeModifyBC(iBoundary) != NO_BC_MODIFICATION )
    InitializeModifiedBC(config_container,
                         geometry_container,
                         element_container[iZone],
                         initial_container,
                         iZone, iBoundary);

  // Extract tuning paramers.
  auto& ParamTuning = config_container->GetParamInletNSCBC();

  // Extract the normalized length scale.
  const as3double lenInv = 1.0/ParamTuning[2];
  // Transverse safety-factor coefficient.
  Coef_eta       = ParamTuning[3];
  // Check if the uncoupled transverse terms are adaptive or not.
  AdaptiveBeta_t = ( ParamTuning[1] < 0 ) ? true : false;

  // Transverse term relaxation coefficient.
  beta_t = Coef_eta*ParamTuning[1];
  // Normal acoustic term relaxation coefficient.
  coefS  = 0.5*lenInv*ParamTuning[0];
  // Normal vorticity/entropy term relaxation coefficient.
  coefE  = lenInv*ParamTuning[0];

  // Account for -ve sign, in case this is a max boundary face. The reason is
  // due to (ovrhoa*dp +/- du), if it is -ve, then account for that.
  if( PsiIndex == 0 ) coefS *= -1.0;
}


CEEInletCBC::~CEEInletCBC
(
 void
)
 /*
	* Destructor for CEEInletCBC class, frees allocated memory.
	*/
{

}


void CEEInletCBC::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current inlet CBC boundary condition.
	*/
{
  // Some abbreviations.
  const as3double gm1   = GAMMA_MINUS_ONE;
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;
  // const as3double ovrg  = 1.0/GAS_CONSTANT;

  // If this boundary condition needs modification, do so.
  if( OrigDOFsIntBoundary.size() )
    ModifyBoundaryCondition(config_container,
                            geometry_container,
                            solver_container[zoneID],
                            element_container[zoneID],
                            spatial_container[zoneID],
                            localTime);

  // Geometry in current zone.
  auto* geometry_zone = geometry_container->GetGeometryZone(zoneID);

  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get the indicial number for the solution in each element sharing this face.
  auto& FaceIndexI = element_container[zoneID]->GetIndexDOFsSol(boundaryID);
  // Get the 1D lagrange interpolation operator.
  auto* ellT    = element_container[zoneID]->GetLagrangeInt1DTranspose();
  // Get the 1D differentiation lagrange operator.
  auto* dellT   = element_container[zoneID]->GetDerLagrangeInt1DTranspose();
  // Get the 1D surface differentiation operator.
  auto* dellS   = element_container[zoneID]->GetDerLagrangeDOFsSol1DFace(boundaryID);
  // Get the integration weights in 1D.
  auto& weights = element_container[zoneID]->GetwDOFsInt1D();

  // Abbreviations.
  const as3double d11  = KronDelta11;
  const as3double d22  = KronDelta22;

  // Assign the correct indices of the working array w.r.t. the boundary orientation.
  as3double **Var    = WorkingDataInt1D[0].data();
  as3double **dVarDn = WorkingDataInt1D[IndexNormal].data();
  as3double **dVarDt = WorkingDataInt1D[IndexTransverse].data();

  // Compute the average mach number on all element faces on this boundary.
  const as3double Mavg = ComputeAverageMachGlobal(geometry_zone, dataI,
                                                  weights, FaceIndexI, ellT, Var);
  // Abbreviation involving Mavg.
  const as3double omM2 = 1.0 - Mavg*Mavg;
  // Override relaxation term specification, if specified.
  if( AdaptiveBeta_t ) beta_t = Coef_eta*Mavg;

  // Loop over all elements on this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element global index.
    unsigned long iElem = ElemIndexI[i];

    // Extract current element size.
    auto ElemSize = geometry_zone->GetGeometryElem(iElem)->GetElemSize();
    // Extract element sizes explicitly.
    const as3double hx = ElemSize[0];
    const as3double hy = ElemSize[1];
    // Compute Jacobians on this element in normal and transverse terms.
    const as3double Jn = d11*(2.0/hx) + d22*(2.0/hy);
    const as3double Jt = d11*(2.0/hy) + d22*(2.0/hx);

    // Abbreviation for inverse of normal Jacobian component.
    const as3double ovJn = 1.0/Jn;
    // Include the normal gradient with metric.
    const as3double dell = Coef_dell*Jn;

    // Step 0a: extract current element integration points of the current boundary.
    auto& dataSolI = dataI[iElem]->GetDataDOFsSol();
    // Step 0b: extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);

    // Step 1a: compute the data and the gradients at the integration points.
    // Note, these are in parametric coordinates.
    TensorProductSolAndGradFace(boundaryID, nDOFsInt1D, nVar, nDOFsSol1D,
                                FaceIndexI.data(),
                                ellT, dellT, dellS,
                                dataSolI.data(), WorkingDataInt1D[0].data(),
                                WorkingDataInt1D[1].data(),
                                WorkingDataInt1D[2].data());

    // Step 1b: convert the parametric gradients into physical space.
    for(unsigned short iVar=0; iVar<nVar; iVar++){
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt1D; l++){
        dVarDn[iVar][l] *= Jn;
        dVarDt[iVar][l] *= Jt;
      }
    }

    // Loop over all the integration points living on this element face.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Step 2: compute the primitive variables.
      const as3double rho   = Var[0][l];
      const as3double ovrho = 1.0/rho;
      const as3double u     = ovrho*Var[1][l];
      const as3double v     = ovrho*Var[2][l];
      const as3double p     = gm1*(Var[3][l]
                            - 0.5*(u*Var[1][l] + v*Var[2][l]) );

    	// Kinetic energy.
    	const as3double ek = 0.5*(u*u + v*v);
      // Temperature.
      // const as3double T  = p*ovrho*ovrg;

      // Step 3a: compute the primitive gradient from conservative in normal direction.
      const as3double drdn = dVarDn[0][l];
      const as3double dudn = ovrho*( dVarDn[1][l] - u*drdn);
    	const as3double dvdn = ovrho*( dVarDn[2][l] - v*drdn );
    	const as3double dpdn = gm1*( ek*dVarDn[0][l]
      										 -  		  u*dVarDn[1][l]
      										 -        v*dVarDn[2][l]
      										 +          dVarDn[3][l] );

      // Step 3b: compute the primitive gradient from conservative in transverse direction.
      const as3double drdt = dVarDt[0][l];
      const as3double dudt = ovrho*( dVarDt[1][l] - u*drdt );
      const as3double dvdt = ovrho*( dVarDt[2][l] - v*drdt );
      const as3double dpdt = gm1*( ek*dVarDt[0][l]
                           -  		  u*dVarDt[1][l]
                           -        v*dVarDt[2][l]
                           +          dVarDt[3][l] );

      // Compute the local speed of sound and its square.
      const as3double a2 = GAMMA*p*ovrho;
      const as3double a  = sqrt(a2);

      // Some abbreviations.
      const as3double ova    = 1.0/a;
      const as3double ova2   = ova*ova;
      const as3double ovrhoa = ovrho*ova;
      // const as3double rrova2 = GAS_CONSTANT*rho*ova2;

      // Reserve memory for the wave amplitudes.
      as3double vecL[nVar], vecT[nVar];

      // Select the velocity component in the normal direction.
      const as3double un = d11*u + d22*v;
      // Select the velocity component in the transverse direction.
      const as3double ut = d22*u + d11*v;
      // Select the derivative w.r.t. normal direction of the normal velocity component.
      const as3double dundn = d11*dudn + d22*dvdn;
      // Select the derivative w.r.t. transverse direction of the normal velocity component.
      const as3double dundt = d11*dudt + d22*dvdt;
      // Select the derivative w.r.t. normal direction of the transverse velocity component.
      const as3double dutdn = d11*dvdn + d22*dudn;
      // Select the derivative w.r.t. transverse direction of the transverse velocity component.
      const as3double dutdt = d11*dvdt + d22*dudt;
      // Compute entropy change in the normal direction.
      const as3double dsdn  = drdn - ova2*dpdn;
      // Compute entropy change in the transverse direction.
      const as3double dsdt  = drdt - ova2*dpdt;

      // Compute the eigenvalues in the normal direction.
      const as3double lmb1 = un - a;
      const as3double lmb2 = un;
      const as3double lmb3 = un + a;

      // Step 4a: compute the normal wave-amplitudes.
      vecL[0] = 0.5*lmb1*( ovrhoa*dpdn - dundn  );
      vecL[1] =     lmb2*( d11*dsdn + d22*dutdn );
      vecL[2] =     lmb2*( d11*dutdn + d22*dsdn );
      vecL[3] = 0.5*lmb3*( ovrhoa*dpdn + dundn  );

      // Step 4b: compute the coupled transverse wave-amplitudes.
      vecT[0] = 0.5*ut*( ovrhoa*dpdt - dundt  ) +    0.5*a*dutdt;
      vecT[1] =     ut*( d11*dsdt + d22*dutdt ) + d22*ovrho*dpdt;
      vecT[2] =     ut*( d11*dutdt + d22*dsdt ) + d11*ovrho*dpdt;
      vecT[3] = 0.5*ut*( ovrhoa*dpdt + dundt  ) +    0.5*a*dutdt;

      // Extract reference data needed.
      const as3double rhoInf = DataDOFsIntBoundary[i][0][l];
      const as3double uInf   = DataDOFsIntBoundary[i][1][l];
      const as3double vInf   = DataDOFsIntBoundary[i][2][l];
      // const as3double pInf   = DataDOFsIntBoundary[i][3][l];
      // Deduce reference temperature.
      // const as3double Tinf   = pInf*ovrg/rhoInf;

      // Step 5a: apply boundary condition through incoming acoustic wave-amplitude.
      vecL[PsiIndex] = coefS*ovrho*a*omM2*( d11*(u - uInf) + d22*(v - vInf) )
                     - beta_t*vecT[PsiIndex];

      // Step 5b: apply boundary condition through incoming entropy/vorticity wave-amplitude.
      // vecL[1] = coefE*a*( -d11*rrova2*(T - Tinf) + d22*(u - uInf) )
      //         - beta_t*vecT[1];
      vecL[1] = coefE*a*( d11*(rho - rhoInf) + d22*(u - uInf) )
              - beta_t*vecT[1];

      // Step 5c: apply boundary condition through incoming entropy/vorticity wave-amplitude.
      // vecL[2] = coefE*a*( d11*(v - vInf) - d22*rrova2*(T - Tinf) )
      //         - beta_t*vecT[2];
      vecL[2] = coefE*a*( d11*(v - vInf) + d22*(rho - rhoInf) )
              - beta_t*vecT[2];

      // More abbvreviations.
      const as3double rhoa   = rho*a;
      const as3double rhoova = rho*ova;
      const as3double ovlmb1 = 1.0/lmb1;
      const as3double ovlmb2 = ( fabs(lmb2) > EPS_VALUE ) ? 1.0/(lmb2 + EPS_VALUE) : 0.0;
      const as3double ovlmb3 = 1.0/lmb3;

      // Step 6: compute the primitive gradient from the normal wave-amplitude.
      // Store this value inside the Var entry in the working array.

      // First entry of the BC-imposed primitive gradient: drdn.
      Var[0][l] =  ovlmb1*rhoova*vecL[0]
                +  ovlmb2*(  d11*vecL[1]
                +            d22*vecL[2] )
                +  ovlmb3*rhoova*vecL[3];

      // Second entry of the BC-imposed primitive gradient: dudn.
      Var[1][l] = -ovlmb1*d11*vecL[0]
                +  ovlmb2*d22*vecL[1]
                +  ovlmb3*d11*vecL[3];

      // Third entry of the BC-imposed primitive gradient: dvdn.
      Var[2][l] = -ovlmb1*d22*vecL[0]
                +  ovlmb2*d11*vecL[2]
                +  ovlmb3*d22*vecL[3];

      // Fourth entry of the BC-imposed primitive gradient: dpdn.
      Var[3][l] =  ovlmb1*rhoa*vecL[0]
                +  ovlmb3*rhoa*vecL[3];

      // Step 7: compute the B matrix entry.
      Var[0][l] -= ( drdn - rho*dell );
      Var[1][l] -= ( dudn -   u*dell );
      Var[2][l] -= ( dvdn -   v*dell );
      Var[3][l] -= ( dpdn -   p*dell );
    }

    // Step 8: compute the boundary state from the B matrix using a least-squares approach.
    for(unsigned short iVar=0; iVar<nVar; iVar++)
      lgemv(nDOFsInt1D, nDOFsInt1D, ovJn, MatrixC, Var[iVar], dataIntI[iVar]);

    // Step 9: convert back into the conservative variables.
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Extract primitive variables.
      const as3double rho = dataIntI[0][l];
      const as3double u   = dataIntI[1][l];
      const as3double v   = dataIntI[2][l];
      const as3double p   = dataIntI[3][l];

      // Compute the conservative variables.
      dataIntI[1][l] = rho*u;
      dataIntI[2][l] = rho*v;
      dataIntI[3][l] = p*ovgm1 + 0.5*rho*( u*u + v*v );
    }
  }
}


CEEPMLInterfaceBoundary::CEEPMLInterfaceBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CInitial      *initial_container,
 CElement     **element_container,
 unsigned short iZone,
 unsigned short iBoundary
)
	:
		CEEInterfaceBoundary
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEPMLInterfaceBoundary per zone: iZone.
	*/
{

}


CEEPMLInterfaceBoundary::~CEEPMLInterfaceBoundary
(
 void
)
 /*
	* Destructor for CEEPMLInterfaceBoundary class, frees allocated memory.
	*/
{

}


void CEEPMLInterfaceBoundary::ImposeBoundaryCondition
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CSolver   **solver_container,
 CElement  **element_container,
 CSpatial  **spatial_container,
 as3double   localTime
)
 /*
	* Function that imposes the current PML interface boundary condition.
	*/
{
  // Get data of current zone.
  auto& dataI = solver_container[zoneID]->GetDataContainer();
  // Get data of matching zone.
  auto& dataJ = solver_container[zoneMatchID]->GetDataContainer();
  // Get the indicial number for the solution on the matching zone and face.
  auto& FaceIndexJ = element_container[zoneMatchID]->GetIndexDOFsSol(boundaryMatchID);

  // Loop over all elements sharing this boundary.
  for(unsigned long i=0; i<ElemIndexI.size(); i++){

    // Deduce current element index.
    unsigned long iElem = ElemIndexI[i];
    // Deduce matching element index.
    unsigned long jElem = ElemIndexJ[i];

    // Extract current element integration points of the current boundary.
    auto& dataIntI = dataI[iElem]->GetDataDOFsIntFace(boundaryID);
    // Extract matching element's solution.
    auto& dataSolJ = dataJ[jElem]->GetDataDOFsSol();

    // Interpolate matching solution to current element's integration
    // points of the current boundary.
    TensorProductSolAndGradFace(boundaryMatchID, nDOFsInt1D, nVar, nDOFsSol1DMatch,
                                FaceIndexJ.data(), lagrangeIntExt1DTranspose,
                                nullptr, nullptr,
                                dataSolJ.data(), dataIntI.data(),
                                nullptr, nullptr);
  }
}




