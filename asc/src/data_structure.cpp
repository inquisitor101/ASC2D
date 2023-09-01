#include "data_structure.hpp"



CData::CData
(
 CConfig   	   *config_container,
 CGeometry 	   *geometry_container,
 CInitial      *initial_container,
 CElement  	   *element_container,
 unsigned short iZone,
 unsigned long  iElem
)
 /*
	* Constructor, used to initialize CData per element.
	*/
{
  // Number of DOFS of solution points on an element in 2D.
  nDOFsSol2D = element_container->GetnDOFsSol2D();
}


CData::~CData
(
 void
)
 /*
	* Destructor for CData class, frees allocated memory.
	*/
{
	for(unsigned short i=0; i<MatchDOFsInt.size(); i++)
		if( MatchDOFsInt[i] ) delete [] MatchDOFsInt[i];

  for(unsigned short i=0; i<DataDOFsSol.size(); i++)
    if( DataDOFsSol[i] ) delete [] DataDOFsSol[i];

  for(unsigned short i=0; i<DataDOFsRes.size(); i++)
    if( DataDOFsRes[i] ) delete [] DataDOFsRes[i];

  for(unsigned short i=0; i<DampDOFsInt.size(); i++)
    if( DampDOFsInt[i] ) delete [] DampDOFsInt[i];

  for(unsigned short i=0; i<derLagrangeGridStretching1D.size(); i++)
    if( derLagrangeGridStretching1D[i] ) delete [] derLagrangeGridStretching1D[i];

  for(unsigned short i=0; i<derLagrangeArtificialConvection1D.size(); i++)
    if( derLagrangeArtificialConvection1D[i] ) delete [] derLagrangeArtificialConvection1D[i];

  for(unsigned short i=0; i<DataDOFsIntMean.size(); i++)
    if( DataDOFsIntMean[i] ) delete [] DataDOFsIntMean[i];

  for(unsigned short i=0; i<DataDOFsIntFace.size(); i++)
    if( !DataDOFsIntFace[i].empty() )
      for(unsigned short j=0; j<DataDOFsIntFace[i].size(); j++)
        if( DataDOFsIntFace[i][j] ) delete [] DataDOFsIntFace[i][j];

  for(unsigned short i=0; i<DataDOFsIntMeanFace.size(); i++)
    if( !DataDOFsIntMeanFace[i].empty() )
      for(unsigned short j=0; j<DataDOFsIntMeanFace[i].size(); j++)
        if( DataDOFsIntMeanFace[i][j] ) delete [] DataDOFsIntMeanFace[i][j];

  for(unsigned short i=0; i<DataDOFsIntAuxFace.size(); i++)
    if( !DataDOFsIntAuxFace[i].empty() )
      for(unsigned short j=0; j<DataDOFsIntAuxFace[i].size(); j++)
        if( DataDOFsIntAuxFace[i][j] ) delete [] DataDOFsIntAuxFace[i][j];
}


CEEData::CEEData
(
 CConfig   	   *config_container,
 CGeometry 	   *geometry_container,
 CInitial      *initial_container,
 CElement  	   *element_container,
 unsigned short iZone,
 unsigned long  iElem
)
	:
		CData
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone, iElem
		)
 /*
	* Constructor, used to initialize CEEData per element.
	*/
{
  // Number of DOFS of integration points on an element in 1D.
  unsigned short nDOFsInt1D = element_container->GetnDOFsInt1D();

  // Check whether or not this is a PML simulation.
  const bool SimPML = (config_container->GetTypeBufferLayer(iZone) == PML_LAYER);

  // Check how many variables needed for the solution and residual DOFs.
  unsigned short nVarTotal = ( SimPML ) ? 2*nVar : nVar;

	// Allocate memory for all needed volume (surface) terms in 2D.
	DataDOFsSol.resize(nVarTotal, nullptr);
  DataDOFsRes.resize(nVarTotal, nullptr);

	// Allocate in every variable.
	for(unsigned short iVar=0; iVar<nVarTotal; iVar++){

		// Allocate actual memory.
		DataDOFsSol[iVar] = new as3double[nDOFsSol2D]();
    DataDOFsRes[iVar] = new as3double[nDOFsSol2D]();

		// Check if allocation failed.
		if( !DataDOFsSol[iVar] || !DataDOFsRes[iVar] )
			Terminate("CEEData::CEEData", __FILE__, __LINE__,
								"Allocation failed for DataDOFsSol and DataDOFsRes.");
	}

  // Extract element face types (internal or boundary).
  auto& ElemType = geometry_container->GetGeometryZone(iZone)->GetInternalElemFace();

  // Reserve memory for number of faces in this element.
  DataDOFsIntFace.resize(nFace);

  // Check if the current face is internal or external. If internal, do nothing.
  // Otherwise, if external, that means it is a boundary face and so allocate
  // the needed amount of memory, respectively.
  for(unsigned short iFace=0; iFace<nFace; iFace++){

    // Check if this is an external/boundary face.
    if( !ElemType[iElem][iFace] ){

      // Reserve memory for number of variables.
      DataDOFsIntFace[iFace].resize(nVar, nullptr);

      for(unsigned short iVar=0; iVar<nVar; iVar++){
        // Allocate actual memory.
        DataDOFsIntFace[iFace][iVar] = new as3double[nDOFsInt1D]();

        // Check if allocation failed.
        if( !DataDOFsIntFace[iFace][iVar] )
          Terminate("CEEData::CEEData", __FILE__, __LINE__,
                    "Allocation failed for DataDOFsIntFace");
      }
    }
  }
}


CEEData::~CEEData
(
 void
)
 /*
	* Destructor for CEEData class, frees allocated memory.
	*/
{

}


CEESpongeData::CEESpongeData
(
 CConfig   	   *config_container,
 CGeometry 	   *geometry_container,
 CInitial      *initial_container,
 CElement  	   *element_container,
 unsigned short iZone,
 unsigned long  iElem
)
	:
		CEEData
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone, iElem
		)
 /*
	* Constructor, used to initialize CEESpongeData per element.
	*/
{
  // Number of DOFS of integration points on an element in 2D.
  unsigned short nDOFsInt2D = element_container->GetnDOFsInt2D();

	// Allocate memory for the pseudo-mean flow.
  DataDOFsIntMean.resize(nVar, nullptr);

	for(unsigned short iVar=0; iVar<nVar; iVar++){

    DataDOFsIntMean[iVar] = new as3double[nDOFsInt2D]();

		// Check if allocation failed.
		if( !DataDOFsIntMean[iVar] )
			Terminate("CEESpongeData::CEESpongeData", __FILE__, __LINE__,
								"Allocation failed for DataDOFsIntMean.");
	}

  // Initialize needed grid-stretching dimensions.
  derLagrangeGridStretching1D.resize(2, nullptr);
  // Initialize needed artificial-velocity dimensions.
  derLagrangeArtificialConvection1D.resize(2, nullptr);

  // Initialize the sponge damping functions.
  InitializeSpongeDamping(config_container,
                          geometry_container,
                          element_container,
                          iZone, iElem);

  // If grid-stretching is specified, then factor it in the differential operators.
  if( config_container->GetGridStretching(iZone) )
    InitializeGridStretching(config_container,
                             geometry_container,
                             initial_container,
                             element_container,
                             iZone, iElem);

  // If an artificial-convection term is specified, initialize velocity functions.
  // Note, this must be executed after the grid-stretching.
  if( config_container->GetArtificialConvection(iZone) )
    InitializeArtificialConvection(config_container,
                                   geometry_container,
                                   initial_container,
                                   element_container,
                                   iZone, iElem);


	// If a characteristic-layer term is specified, initialize the matching profile.
	// Note, this must be executed after the grid-stretching.
	if( config_container->GetCharacteristicMatching(iZone) )
		InitializeCharacteristicMatching(config_container,
				                             geometry_container,
																		 element_container,
																		 iZone, iElem);
}


CEESpongeData::~CEESpongeData
(
 void
)
 /*
	* Destructor for CEESpongeData class, frees allocated memory.
	*/
{

}


void CEESpongeData::InitializeCharacteristicMatching
(
  CConfig       *config_container,
  CGeometry     *geometry_container,
  CElement      *element_container,
  unsigned short iZone,
  unsigned long  iElem
)
 /*
  * Function that initializes and defines the characteristic matching coefficients.
  */
{
  // Number of integration DOFs in 2D.
  const unsigned short nDOFsInt2D = element_container->GetnDOFsInt2D();

  // Get grid zone.
  auto* grid_zone  = geometry_container->GetGeometryZone(iZone);
  // Get grid element.
  auto* grid_elem  = grid_zone->GetGeometryElem(iElem);
  // Get zone dimensions/bounds.
  auto& zone_bound = grid_zone->GetZoneSize();
  // Main domain bounding box.
  auto& MainBox    = config_container->GetDomainBound();

  // To avoid confusion, explicitly extract the main-zone bounding coordinates.
  const as3double xMin = MainBox[0]; // west.
  const as3double xMax = MainBox[1]; // east.
  const as3double yMin = MainBox[2]; // south.
  const as3double yMax = MainBox[3]; // north.

  // Type of zone.
  unsigned short TypeZone = config_container->GetTypeZone(iZone);

  // Extract characteristic-matching constant.
  const as3double ms = config_container->GetCharacteristicConstant(iZone);
  // Extract characteristic-matching exponential.
  const as3double mb = config_container->GetCharacteristicExponent(iZone);

  // Extract grid-stretching constant.
  const as3double gs = config_container->GetGridStretchingConstant(iZone);
  // Extract grid-stretching exponential.
  const as3double gb = config_container->GetGridStretchingExponent(iZone);

	// Extract the coordinates of the volume integration nodes.
	auto& q = grid_elem->GetCoordIntDOFs();

  // Obtain min x-coordinate of face.
  const as3double xmin = grid_elem->GetCoordBoundary(IDX_WEST);
  // Obtain max x-coordinate of face.
  const as3double xmax = grid_elem->GetCoordBoundary(IDX_EAST);
  // Obtain min y-coordinate of face.
  const as3double ymin = grid_elem->GetCoordBoundary(IDX_SOUTH);
  // Obtain max y-coordinate of face.
  const as3double ymax = grid_elem->GetCoordBoundary(IDX_NORTH);

	// Initialize characteristic matching on the volume integration nodes in 2D.
	MatchDOFsInt.resize( nDim, nullptr );
	MatchDOFsInt[0] = new as3double[nDOFsInt2D]();
	MatchDOFsInt[1] = new as3double[nDOFsInt2D]();

	// Check if allocation failed.
	if( !MatchDOFsInt[0] || !MatchDOFsInt[1] )
		Terminate("InitializeCharacteristicMatching", __FILE__, __LINE__,
				      "Allocation failed for MatchDOFsInt.");


  // Determine which type of zone we are dealing with.
  switch( TypeZone ){

    // In these zones: alpha(x) is non-zero and alpha(y) = 1.
    case(ZONE_EAST): case(ZONE_WEST):
    {
      // Extract inverse of width of zone in x-direction.
      const as3double ovDx = 1.0/zone_bound[0];
      // Interface location in x-direction.
      const as3double x0   = ( TypeZone == ZONE_WEST ) ? xMin : xMax;

			// Loop over volume integration nodes and populate the matching profile.
			for(unsigned short l=0; l<nDOFsInt2D; l++)
			{
				// Compute the necessary functions.
				const as3double fx = gs*pow( ovDx*fabs( q[0][l] - x0 ), gb );
				const as3double gx = ms*pow( ovDx*fabs( q[0][l] - x0 ), mb );
				// Assign the matching profiles.
				MatchDOFsInt[0][l] = gx/( 1.0 + fx ); 
				MatchDOFsInt[1][l] = 0.0;
			}

      break;
    }

    // In these zones: alpha(y) is non-zero and alpha(x) = 1.
    case(ZONE_SOUTH): case(ZONE_NORTH):
    {
      // Extract inverse of width of zone in y-direction.
      const as3double ovDy = 1.0/zone_bound[1];
      // Interface location in y-direction.
      const as3double y0   = ( TypeZone == ZONE_SOUTH ) ? yMin : yMax;
     
			// Loop over volume integration nodes and populate the matching profile.
			for(unsigned short l=0; l<nDOFsInt2D; l++)
			{
				// Compute the necessary functions.
				const as3double fy = gs*pow( ovDy*fabs( q[1][l] - y0 ), gb );
				const as3double gy = ms*pow( ovDy*fabs( q[1][l] - y0 ), mb );
				// Assign the matching profiles.
				MatchDOFsInt[0][l] = 0.0; 
				MatchDOFsInt[1][l] = gy/( 1.0 + fy );
			}

      break;
    }

    // In these zones: alpha(x) and alpha(y) are non-zero.
    case(ZONE_CORNER_0): case(ZONE_CORNER_1): case(ZONE_CORNER_2): case(ZONE_CORNER_3):
    {
      // Extract inverse of width of zone in x-direction.
      const as3double ovDx = 1.0/zone_bound[0];
      // Extract inverse of width of zone in y-direction.
      const as3double ovDy = 1.0/zone_bound[1];
      // Interface location in x- and y-directions.
      const as3double x0   = ( (TypeZone == ZONE_CORNER_0) || (TypeZone == ZONE_CORNER_2) ) ? xMin : xMax;
      const as3double y0   = ( (TypeZone == ZONE_CORNER_0) || (TypeZone == ZONE_CORNER_1) ) ? yMin : yMax;

			// Loop over volume integration nodes and populate the matching profile.
			for(unsigned short l=0; l<nDOFsInt2D; l++)
			{
				// Compute the necessary functions.
				const as3double fx = gs*pow( ovDx*fabs( q[0][l] - x0 ), gb );
				const as3double gx = ms*pow( ovDx*fabs( q[0][l] - x0 ), mb );
				const as3double fy = gs*pow( ovDy*fabs( q[1][l] - y0 ), gb );
				const as3double gy = ms*pow( ovDy*fabs( q[1][l] - y0 ), mb );
				// Assign the matching profiles.
				MatchDOFsInt[0][l] = gx/( 1.0 + fx ); 
				MatchDOFsInt[1][l] = gy/( 1.0 + fy );
			}

      break;
    }

    // If this is the main zone, exit immediately.
    default:
      Terminate("CEESpongeData::InitializeCharacteristicMatching", __FILE__, __LINE__,
                "Wrong zone type specified.");
  }
}


void CEESpongeData::InitializeGridStretching
(
  CConfig       *config_container,
  CGeometry     *geometry_container,
  CInitial      *initial_container,
  CElement      *element_container,
  unsigned short iZone,
  unsigned long  iElem
)
 /*
  * Function that initializes and defines the grid-stretching coefficients.
  */
{
  // Number of solution DOFs in 1D.
  const unsigned short nDOFsSol1D = element_container->GetnDOFsSol1D();
  // Number of integration DOFs in 1D.
  const unsigned short nDOFsInt1D = element_container->GetnDOFsInt1D();

  // Get grid zone.
  auto* grid_zone  = geometry_container->GetGeometryZone(iZone);
  // Get grid element.
  auto* grid_elem  = grid_zone->GetGeometryElem(iElem);
  // Get zone dimensions/bounds.
  auto& zone_bound = grid_zone->GetZoneSize();
  // Main domain bounding box.
  auto& MainBox    = config_container->GetDomainBound();
  // Lagrange differential operator.
  auto* dell       = element_container->GetDerLagrangeInt1D();

  // To avoid confusion, explicitly extract the main-zone bounding coordinates.
  const as3double xMin = MainBox[0]; // west.
  const as3double xMax = MainBox[1]; // east.
  const as3double yMin = MainBox[2]; // south.
  const as3double yMax = MainBox[3]; // north.

  // Type of zone.
  unsigned short TypeZone = config_container->GetTypeZone(iZone);

  // Extract grid-stretching constant.
  const as3double sigma = config_container->GetGridStretchingConstant(iZone);
  // Extract grid-stretching exponential.
  const as3double beta  = config_container->GetGridStretchingExponent(iZone);

  // Flag used for error detection.
  bool ErrorDetected = false;

  // Determine which type of zone we are dealing with.
  switch( TypeZone ){

    // In these zones: alpha(x) is non-zero and alpha(y) = 1.
    case(ZONE_EAST): case(ZONE_WEST):
    {
      // Allocate actual data.
      derLagrangeGridStretching1D[0] = new as3double[nDOFsSol1D*nDOFsInt1D]();

      // Extract inverse of width of zone in x-direction.
      const as3double ovDx = 1.0/zone_bound[0];
      // Interface location in x-direction.
      const as3double x0   = ( TypeZone == ZONE_WEST ) ? xMin : xMax;
      // Obtain min x-coordinate of face.
      const as3double xmin = grid_elem->GetCoordBoundary(IDX_WEST);
      // Obtain max x-coordinate of face.
      const as3double xmax = grid_elem->GetCoordBoundary(IDX_EAST);

      // Initialize grid-stretching on each face.
      derLagrangeGridStretchingFace.resize(nFace, 1.0);

      // Change the required faces, depending on the gradients.
      derLagrangeGridStretchingFace[IDX_WEST] = 1.0/( 1.0 + sigma*pow(ovDx*fabs(xmin - x0), beta) );
      derLagrangeGridStretchingFace[IDX_EAST] = 1.0/( 1.0 + sigma*pow(ovDx*fabs(xmax - x0), beta) );

      // Extract solution points.
      auto& coordSol = grid_elem->GetCoordSolDOFs();

      // Extract the varying x-coordinate.
      as3vector1d<as3double> xCoord(nDOFsSol1D);
      for(unsigned short i=0; i<nDOFsSol1D; i++) xCoord[i] = coordSol[0][i];
      // Consistency check.
      assert( fabs(xmin-xCoord[0]) < 1e-10 && fabs(xmax-xCoord[nDOFsSol1D-1]) < 1e-10 );

      // Populate grid-stretching data.
      unsigned short idx = 0;
      for(unsigned short i=0; i<nDOFsInt1D; i++){
        for(unsigned short l=0; l<nDOFsSol1D; l++){

          // Extract coordinates.
          const as3double x = xCoord[l];

          // Compute the actual grid-stretching function.
          derLagrangeGridStretching1D[0][idx] = dell[idx]/( 1.0 + sigma*pow(ovDx*fabs(x - x0), beta) );

          // Update index.
          idx++;
        }
      }

      // Check if allocation failed.
      if( !derLagrangeGridStretching1D[0] ) ErrorDetected = true;

      break;
    }

    // In these zones: alpha(y) is non-zero and alpha(x) = 1.
    case(ZONE_SOUTH): case(ZONE_NORTH):
    {
      // Allocate actual data.
      derLagrangeGridStretching1D[1] = new as3double[nDOFsSol1D*nDOFsInt1D]();

      // Extract inverse of width of zone in y-direction.
      const as3double ovDy = 1.0/zone_bound[1];
      // Interface location in y-direction.
      const as3double y0   = ( TypeZone == ZONE_SOUTH ) ? yMin : yMax;
      // Obtain min y-coordinate of face.
      const as3double ymin = grid_elem->GetCoordBoundary(IDX_SOUTH);
      // Obtain max y-coordinate of face.
      const as3double ymax = grid_elem->GetCoordBoundary(IDX_NORTH);

      // Initialize Grid-stretching on each face.
      derLagrangeGridStretchingFace.resize(nFace, 1.0);

      // Change the required faces, depending on the gradients.
      derLagrangeGridStretchingFace[IDX_SOUTH] = 1.0/( 1.0 + sigma*pow(ovDy*fabs(ymin - y0), beta) );
      derLagrangeGridStretchingFace[IDX_NORTH] = 1.0/( 1.0 + sigma*pow(ovDy*fabs(ymax - y0), beta) );

      // Extract solution points.
      auto& coordSol = grid_elem->GetCoordSolDOFs();

      // Extract the varying y-coordinate.
      as3vector1d<as3double> yCoord(nDOFsSol1D);
      for(unsigned short i=0; i<nDOFsSol1D; i++) yCoord[i] = coordSol[1][i*nDOFsSol1D];
      // Consistency check.
      assert( fabs(ymin-yCoord[0]) < 1e-10 && fabs(ymax-yCoord[nDOFsSol1D-1]) < 1e-10 );

      // Populate grid-stretching data.
      unsigned short idx = 0;
      for(unsigned short i=0; i<nDOFsInt1D; i++){
        for(unsigned short l=0; l<nDOFsSol1D; l++){

          // Extract coordinates.
          const as3double y = yCoord[l];

          // Compute the actual grid-stretching function.
          derLagrangeGridStretching1D[1][idx] = dell[idx]/( 1.0 + sigma*pow(ovDy*fabs(y - y0), beta) );

          // Update index.
          idx++;
        }
      }

      // Check if allocation failed.
      if( !derLagrangeGridStretching1D[1] ) ErrorDetected = true;

      break;
    }

    // In these zones: alpha(x) and alpha(y) are non-zero.
    case(ZONE_CORNER_0): case(ZONE_CORNER_1): case(ZONE_CORNER_2): case(ZONE_CORNER_3):
    {
      // Allocate actual data.
      derLagrangeGridStretching1D[0] = new as3double[nDOFsSol1D*nDOFsInt1D]();
      derLagrangeGridStretching1D[1] = new as3double[nDOFsSol1D*nDOFsInt1D]();

      // Extract inverse of width of zone in x-direction.
      const as3double ovDx = 1.0/zone_bound[0];
      // Extract inverse of width of zone in y-direction.
      const as3double ovDy = 1.0/zone_bound[1];
      // Interface location in x- and y-directions.
      const as3double x0   = ( (TypeZone == ZONE_CORNER_0) || (TypeZone == ZONE_CORNER_2) ) ? xMin : xMax;
      const as3double y0   = ( (TypeZone == ZONE_CORNER_0) || (TypeZone == ZONE_CORNER_1) ) ? yMin : yMax;
      // Obtain min x-coordinate of face.
      const as3double xmin = grid_elem->GetCoordBoundary(IDX_WEST);
      // Obtain max x-coordinate of face.
      const as3double xmax = grid_elem->GetCoordBoundary(IDX_EAST);
      // Obtain min y-coordinate of face.
      const as3double ymin = grid_elem->GetCoordBoundary(IDX_SOUTH);
      // Obtain max y-coordinate of face.
      const as3double ymax = grid_elem->GetCoordBoundary(IDX_NORTH);

      // Initialize grid-stretching on each face.
      derLagrangeGridStretchingFace.resize(nFace, 1.0);

      // Change the required faces, depending on the gradients.
      derLagrangeGridStretchingFace[IDX_WEST]  = 1.0/( 1.0 + sigma*pow(ovDx*fabs(xmin - x0), beta) );
      derLagrangeGridStretchingFace[IDX_EAST]  = 1.0/( 1.0 + sigma*pow(ovDx*fabs(xmax - x0), beta) );
      derLagrangeGridStretchingFace[IDX_SOUTH] = 1.0/( 1.0 + sigma*pow(ovDy*fabs(ymin - y0), beta) );
      derLagrangeGridStretchingFace[IDX_NORTH] = 1.0/( 1.0 + sigma*pow(ovDy*fabs(ymax - y0), beta) );

      // Extract solution points.
      auto& coordSol = grid_elem->GetCoordSolDOFs();

      // Extract the varying x- and y-coordinates.
      as3vector1d<as3double> xCoord(nDOFsSol1D);
      as3vector1d<as3double> yCoord(nDOFsSol1D);
      for(unsigned short i=0; i<nDOFsSol1D; i++) xCoord[i] = coordSol[0][i];
      for(unsigned short i=0; i<nDOFsSol1D; i++) yCoord[i] = coordSol[1][i*nDOFsSol1D];
      // Consistency check.
      assert( fabs(xmin-xCoord[0]) < 1e-10 && fabs(xmax-xCoord[nDOFsSol1D-1]) < 1e-10 );
      assert( fabs(ymin-yCoord[0]) < 1e-10 && fabs(ymax-yCoord[nDOFsSol1D-1]) < 1e-10 );

      // Populate grid-stretching data.
      unsigned short idx = 0;
      for(unsigned short i=0; i<nDOFsInt1D; i++){
        for(unsigned short l=0; l<nDOFsSol1D; l++){

          // Extract coordinates.
          const as3double x = xCoord[l];
          const as3double y = yCoord[l];

          // Compute the actual grid-stretching function.
          derLagrangeGridStretching1D[0][idx] = dell[idx]/( 1.0 + sigma*pow(ovDx*fabs(x - x0), beta) );
          derLagrangeGridStretching1D[1][idx] = dell[idx]/( 1.0 + sigma*pow(ovDy*fabs(y - y0), beta) );

          // Update index.
          idx++;
        }
      }

      // Check if allocation failed.
      if( !derLagrangeGridStretching1D[0] || !derLagrangeGridStretching1D[1] ) ErrorDetected = true;

      break;
    }

    // If this is the main zone, exit immediately.
    default:
      Terminate("CEESpongeData::InitializeGridStretching", __FILE__, __LINE__,
                "Wrong zone type specified.");
  }

  // Check if allocation failed for grid-stretching, if specified.
  if( ErrorDetected )
    Terminate("CEESpongeData::InitializeGridStretching", __FILE__, __LINE__,
              "Allocation failed for derLagrangeGridStretching1D.");
}


void CEESpongeData::InitializeArtificialConvection
(
  CConfig       *config_container,
  CGeometry     *geometry_container,
  CInitial      *initial_container,
  CElement      *element_container,
  unsigned short iZone,
  unsigned long  iElem
)
 /*
  * Function that initializes and defines the artificial-velocity coefficients.
  */
{
  // Consistency check, this cannot be a PML layer.
  if( config_container->GetTypeBufferLayer(iZone) == PML_LAYER )
    Terminate("CEESpongeData::InitializeArtificialConvection", __FILE__, __LINE__,
              "Cannot use artificial convection with a PML layer.");

  // Number of solution DOFs in 1D.
  const unsigned short nDOFsSol1D = element_container->GetnDOFsSol1D();
  // Number of integration DOFs in 1D.
  const unsigned short nDOFsInt1D = element_container->GetnDOFsInt1D();

  // Get grid zone.
  auto* grid_zone  = geometry_container->GetGeometryZone(iZone);
  // Get grid element.
  auto* grid_elem  = grid_zone->GetGeometryElem(iElem);
  // Get zone dimensions/bounds.
  auto& zone_bound = grid_zone->GetZoneSize();
  // Main domain bounding box.
  auto& MainBox    = config_container->GetDomainBound();
  // Lagrange differential operator.
  auto* dell       = element_container->GetDerLagrangeInt1D();

  // To avoid confusion, explicitly extract the main-zone bounding coordinates.
  const as3double xMin = MainBox[0]; // west.
  const as3double xMax = MainBox[1]; // east.
  const as3double yMin = MainBox[2]; // south.
  const as3double yMax = MainBox[3]; // north.

  // Type of zone.
  unsigned short TypeZone = config_container->GetTypeZone(iZone);

  // Extract artificial-convection constant.
  const as3double maxmach = config_container->GetArtificialConvectionConstant(iZone);
  // Extract artificial-convection exponential.
  const as3double beta    = config_container->GetArtificialConvectionExponent(iZone);

  // Flag used for error detection.
  bool ErrorDetected = false;

  // Determine which type of zone we are dealing with.
  switch( TypeZone ){

    // In these zones: U(x) is non-zero and U(y) = 0.
    case(ZONE_EAST): case(ZONE_WEST):
    {
      // Allocate actual data.
      derLagrangeArtificialConvection1D[0] = new as3double[nDOFsSol1D*nDOFsInt1D]();

      // Extract inverse of width of zone in x-direction.
      const as3double ovDx = 1.0/zone_bound[0];
      // Interface location in x-direction.
      const as3double x0   = ( TypeZone == ZONE_WEST ) ? xMin : xMax;
      // Obtain min x-coordinate of face.
      const as3double xmin = grid_elem->GetCoordBoundary(IDX_WEST);
      // Obtain max x-coordinate of face.
      const as3double xmax = grid_elem->GetCoordBoundary(IDX_EAST);

      // Check sign and value of max artificial max number at boundary.
      as3double maxmachx = maxmach;
      if( TypeZone == ZONE_WEST )
        if( config_container->GetTypeBC(iZone, IDX_WEST) == BC_SUPERSONIC_OUTLET )
          maxmachx *= -1.0;

      // Initialize artificial-velocity on each face.
      derLagrangeArtificialConvectionFace.resize(nFace, 0.0);

      // Change the required faces, depending on the gradients.
      derLagrangeArtificialConvectionFace[IDX_WEST] = maxmachx*pow( ovDx*fabs(xmin - x0), beta );
      derLagrangeArtificialConvectionFace[IDX_EAST] = maxmachx*pow( ovDx*fabs(xmax - x0), beta );

      // Extract solution points.
      auto& coordSol = grid_elem->GetCoordSolDOFs();

      // Extract the varying x-coordinate.
      as3vector1d<as3double> xCoord(nDOFsSol1D);
      for(unsigned short i=0; i<nDOFsSol1D; i++) xCoord[i] = coordSol[0][i];
      // Consistency check.
      assert( fabs(xmin-xCoord[0]) < 1e-10 && fabs(xmax-xCoord[nDOFsSol1D-1]) < 1e-10 );

      // Populate artificial-velocity data.
      unsigned short idx = 0;
      for(unsigned short i=0; i<nDOFsInt1D; i++){
        for(unsigned short l=0; l<nDOFsSol1D; l++){

          // Extract coordinates.
          const as3double x = xCoord[l];

          // Compute the actual artificial-velocity function.
          derLagrangeArtificialConvection1D[0][idx] = maxmachx*pow( ovDx*fabs(x - x0), beta );

          // Update index.
          idx++;
        }
      }

      // Check if allocation failed.
      if( !derLagrangeArtificialConvection1D[0] ) ErrorDetected = true;

      // If grid-stretching is used, factor that in too.
      if( derLagrangeGridStretching1D[0] ){
        for(unsigned short i=0; i<nDOFsSol1D*nDOFsInt1D; i++)
          derLagrangeArtificialConvection1D[0][i] *= derLagrangeGridStretching1D[0][i];
      }
      else {
        // Include the standard lagrange differential operator.
        for(unsigned short i=0; i<nDOFsSol1D*nDOFsInt1D; i++)
          derLagrangeArtificialConvection1D[0][i] *= dell[i];
      }

      break;
    }

    // In these zones: U(y) is non-zero and U(x) = 0.
    case(ZONE_SOUTH): case(ZONE_NORTH):
    {
      // Allocate actual data.
      derLagrangeArtificialConvection1D[1] = new as3double[nDOFsSol1D*nDOFsInt1D]();

      // Extract inverse of width of zone in y-direction.
      const as3double ovDy = 1.0/zone_bound[1];
      // Interface location in y-direction.
      const as3double y0   = ( TypeZone == ZONE_SOUTH ) ? yMin : yMax;
      // Obtain min y-coordinate of face.
      const as3double ymin = grid_elem->GetCoordBoundary(IDX_SOUTH);
      // Obtain max y-coordinate of face.
      const as3double ymax = grid_elem->GetCoordBoundary(IDX_NORTH);

      // Check sign and value of max artificial max number at boundary.
      as3double maxmachy = maxmach;
      if( TypeZone == ZONE_SOUTH )
        if( config_container->GetTypeBC(iZone, IDX_SOUTH) == BC_SUPERSONIC_OUTLET )
          maxmachy *= -1.0;

      // Initialize Grid-stretching on each face.
      derLagrangeArtificialConvectionFace.resize(nFace, 0.0);

      // Change the required faces, depending on the gradients.
      derLagrangeArtificialConvectionFace[IDX_SOUTH] = maxmachy*pow( ovDy*fabs(ymin - y0), beta );
      derLagrangeArtificialConvectionFace[IDX_NORTH] = maxmachy*pow( ovDy*fabs(ymax - y0), beta );

      // Extract solution points.
      auto& coordSol = grid_elem->GetCoordSolDOFs();

      // Extract the varying y-coordinate.
      as3vector1d<as3double> yCoord(nDOFsSol1D);
      for(unsigned short i=0; i<nDOFsSol1D; i++) yCoord[i] = coordSol[1][i*nDOFsSol1D];
      // Consistency check.
      assert( fabs(ymin-yCoord[0]) < 1e-10 && fabs(ymax-yCoord[nDOFsSol1D-1]) < 1e-10 );

      // Populate artificial-velocity data.
      unsigned short idx = 0;
      for(unsigned short i=0; i<nDOFsInt1D; i++){
        for(unsigned short l=0; l<nDOFsSol1D; l++){

          // Extract coordinates.
          const as3double y = yCoord[l];

          // Compute the actual artificial-velocity function.
          derLagrangeArtificialConvection1D[1][idx] = maxmachy*pow( ovDy*fabs(y - y0), beta );

          // Update index.
          idx++;
        }
      }

      // Check if allocation failed.
      if( !derLagrangeArtificialConvection1D[1] ) ErrorDetected = true;

      // If grid-stretching is used, factor that in too.
      if( derLagrangeGridStretching1D[1] ){
        for(unsigned short i=0; i<nDOFsSol1D*nDOFsInt1D; i++)
          derLagrangeArtificialConvection1D[1][i] *= derLagrangeGridStretching1D[1][i];
      }
      else {
        // Include the standard lagrange differential operator.
        for(unsigned short i=0; i<nDOFsSol1D*nDOFsInt1D; i++)
          derLagrangeArtificialConvection1D[1][i] *= dell[i];
      }

      break;
    }

    // In these zones: U(x) and U(y) are non-zero.
    case(ZONE_CORNER_0): case(ZONE_CORNER_1): case(ZONE_CORNER_2): case(ZONE_CORNER_3):
    {
      // Allocate actual data.
      derLagrangeArtificialConvection1D[0] = new as3double[nDOFsSol1D*nDOFsInt1D]();
      derLagrangeArtificialConvection1D[1] = new as3double[nDOFsSol1D*nDOFsInt1D]();

      // Extract inverse of width of zone in x-direction.
      const as3double ovDx = 1.0/zone_bound[0];
      // Extract inverse of width of zone in y-direction.
      const as3double ovDy = 1.0/zone_bound[1];
      // Interface location in x- and y-directions.
      const as3double x0   = ( (TypeZone == ZONE_CORNER_0) || (TypeZone == ZONE_CORNER_2) ) ? xMin : xMax;
      const as3double y0   = ( (TypeZone == ZONE_CORNER_0) || (TypeZone == ZONE_CORNER_1) ) ? yMin : yMax;
      // Obtain min x-coordinate of face.
      const as3double xmin = grid_elem->GetCoordBoundary(IDX_WEST);
      // Obtain max x-coordinate of face.
      const as3double xmax = grid_elem->GetCoordBoundary(IDX_EAST);
      // Obtain min y-coordinate of face.
      const as3double ymin = grid_elem->GetCoordBoundary(IDX_SOUTH);
      // Obtain max y-coordinate of face.
      const as3double ymax = grid_elem->GetCoordBoundary(IDX_NORTH);

      // Check sign and value of max artificial max number at boundary.
      as3double maxmachx = maxmach, maxmachy = maxmach;
      if( TypeZone == ZONE_CORNER_0 ){
        if( config_container->GetTypeBC(iZone, IDX_WEST)  == BC_SUPERSONIC_OUTLET )
          maxmachx *= -1.0;
        if( config_container->GetTypeBC(iZone, IDX_SOUTH) == BC_SUPERSONIC_OUTLET )
          maxmachy *= -1.0;
      }
      if( TypeZone == ZONE_CORNER_1 ){
        if( config_container->GetTypeBC(iZone, IDX_SOUTH) == BC_SUPERSONIC_OUTLET )
          maxmachy *= -1.0;
      }
      if( TypeZone == ZONE_CORNER_2 ){
        if( config_container->GetTypeBC(iZone, IDX_WEST)  == BC_SUPERSONIC_OUTLET )
          maxmachx *= -1.0;
      }

      // Correct convection term in case there is no convection on one side of a south-west corner zone.
      if( TypeZone == ZONE_CORNER_0 ){
        if( !config_container->GetArtificialConvection(ZONE_WEST)  ) maxmachx = 0.0;
        if( !config_container->GetArtificialConvection(ZONE_SOUTH) ) maxmachy = 0.0;
      }

      // Correct convection term in case there is no convection on one side of a south-east corner zone
      if( TypeZone == ZONE_CORNER_1 ){
        if( !config_container->GetArtificialConvection(ZONE_EAST)  ) maxmachx = 0.0;
        if( !config_container->GetArtificialConvection(ZONE_SOUTH) ) maxmachy = 0.0;
      }

      // Correct convection term in case there is no convection on one side of a north-west corner zone
      if( TypeZone == ZONE_CORNER_2 ){
        if( !config_container->GetArtificialConvection(ZONE_WEST)  ) maxmachx = 0.0;
        if( !config_container->GetArtificialConvection(ZONE_NORTH) ) maxmachy = 0.0;
      }

      // Correct convection term in case there is no convection on one side of a north-east corner zone
      if( TypeZone == ZONE_CORNER_3 ){
        if( !config_container->GetArtificialConvection(ZONE_EAST)  ) maxmachx = 0.0;
        if( !config_container->GetArtificialConvection(ZONE_NORTH) ) maxmachy = 0.0;
      }


      // Initialize artificial-velocity on each face.
      derLagrangeArtificialConvectionFace.resize(nFace, 0.0);

      // Change the required faces, depending on the gradients.
      derLagrangeArtificialConvectionFace[IDX_WEST]  = maxmachx*pow( ovDx*fabs(xmin - x0), beta );
      derLagrangeArtificialConvectionFace[IDX_EAST]  = maxmachx*pow( ovDx*fabs(xmax - x0), beta );
      derLagrangeArtificialConvectionFace[IDX_SOUTH] = maxmachy*pow( ovDy*fabs(ymin - y0), beta );
      derLagrangeArtificialConvectionFace[IDX_NORTH] = maxmachy*pow( ovDy*fabs(ymax - y0), beta );

      // Extract solution points.
      auto& coordSol = grid_elem->GetCoordSolDOFs();

      // Extract the varying x- and y-coordinates.
      as3vector1d<as3double> xCoord(nDOFsSol1D);
      as3vector1d<as3double> yCoord(nDOFsSol1D);
      for(unsigned short i=0; i<nDOFsSol1D; i++) xCoord[i] = coordSol[0][i];
      for(unsigned short i=0; i<nDOFsSol1D; i++) yCoord[i] = coordSol[1][i*nDOFsSol1D];
      // Consistency check.
      assert( fabs(xmin-xCoord[0]) < 1e-10 && fabs(xmax-xCoord[nDOFsSol1D-1]) < 1e-10 );
      assert( fabs(ymin-yCoord[0]) < 1e-10 && fabs(ymax-yCoord[nDOFsSol1D-1]) < 1e-10 );

      // Populate artificial-velocity data.
      unsigned short idx = 0;
      for(unsigned short i=0; i<nDOFsInt1D; i++){
        for(unsigned short l=0; l<nDOFsSol1D; l++){

          // Extract coordinates.
          const as3double x = xCoord[l];
          const as3double y = yCoord[l];

          // Compute the actual artificial-velocity function.
          derLagrangeArtificialConvection1D[0][idx] = maxmachx*pow( ovDx*fabs(x - x0), beta );
          derLagrangeArtificialConvection1D[1][idx] = maxmachy*pow( ovDy*fabs(y - y0), beta );

          // Update index.
          idx++;
        }
      }

      // Check if allocation failed.
      if( !derLagrangeArtificialConvection1D[0] || !derLagrangeArtificialConvection1D[1] ) ErrorDetected = true;

      // If grid-stretching is used, factor that in too.
      for(unsigned short k=0; k<nDim; k++){
        if( derLagrangeGridStretching1D[k] ){
          for(unsigned short i=0; i<nDOFsSol1D*nDOFsInt1D; i++)
            derLagrangeArtificialConvection1D[k][i] *= derLagrangeGridStretching1D[k][i];
        }
        else {
          // Include the standard lagrange differential operator.
          for(unsigned short i=0; i<nDOFsSol1D*nDOFsInt1D; i++)
            derLagrangeArtificialConvection1D[k][i] *= dell[i];
        }
      }

      break;
    }

    // If this is the main zone, exit immediately.
    default:
      Terminate("CEESpongeData::InitializeArtificialConvection", __FILE__, __LINE__,
                "Wrong zone type specified.");
  }

  // Check if allocation failed for grid-stretching, if specified.
  if( ErrorDetected )
    Terminate("CEESpongeData::InitializeArtificialConvection", __FILE__, __LINE__,
              "Allocation failed for derLagrangeArtificialConvection1D.");
}


void CEESpongeData::InitializeSpongeDamping
(
  CConfig       *config_container,
  CGeometry     *geometry_container,
  CElement      *element_container,
  unsigned short iZone,
  unsigned long  iElem
)
 /*
  * Function that initializes and defines the sponge damping coefficients.
  */
{
  // Number of DOFs of integration points on an element in 2D.
  unsigned short nDOFsInt2D = element_container->GetnDOFsInt2D();
  // Get grid zone.
  auto* grid_zone  = geometry_container->GetGeometryZone(iZone);
  // Get grid element.
  auto* grid_elem  = grid_zone->GetGeometryElem(iElem);
  // Get zone dimensions/bounds.
  auto& zone_bound = grid_zone->GetZoneSize();
  // Get grid element coordinates at integration points.
  auto& coord      = grid_elem->GetCoordIntDOFs();
  // Main domain bounding box.
  auto& MainBox    = config_container->GetDomainBound();

  // To avoid confusion, explicitly extract the main-zone bounding coordinates.
  const as3double xMin = MainBox[0]; // west.
  const as3double xMax = MainBox[1]; // east.
  const as3double yMin = MainBox[2]; // south.
  const as3double yMax = MainBox[3]; // north.

  // Type of zone.
  unsigned short TypeZone = config_container->GetTypeZone(iZone);

  // Get damping constant.
  const as3double sigma   = config_container->GetDampingConstant(iZone);
  // Get damping exponential.
  const as3double beta    = config_container->GetDampingExponent(iZone);

  // Determine which type of zone we are dealing with.
  switch( TypeZone ){

    // In these zones: sigma(x) is non-zero and sigma(y) = 0.
    case(ZONE_EAST): case(ZONE_WEST):
    {
      // Initialize needed damping dimensions.
      DampDOFsInt.resize(2, nullptr);
      // Allocate actual data.
      DampDOFsInt[0] = new as3double[nDOFsInt2D]();
      DampDOFsInt[1] = new as3double[nDOFsInt2D]();

      // Extract inverse of width of zone in x-direction.
      const as3double ovDx = 1.0/zone_bound[0];

      // Interface location in x-direction.
      const as3double x0 = ( TypeZone == ZONE_WEST ) ? xMin : xMax;

      // Populate damping data.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){

        // Extract coordinates.
        const as3double x = coord[0][l];

        // Compute the actual damping function.
        DampDOFsInt[0][l] = sigma*pow( ovDx*fabs(x - x0), beta );
        DampDOFsInt[1][l] = 0.0;
      }

      break;
    }

    // In these zones: sigma(y) is non-zero and sigma(x) = 0.
    case(ZONE_SOUTH): case(ZONE_NORTH):
    {
      // Initialize needed damping dimensions.
      DampDOFsInt.resize(2, nullptr);
      // Allocate actual data.
      DampDOFsInt[0] = new as3double[nDOFsInt2D]();
      DampDOFsInt[1] = new as3double[nDOFsInt2D]();

      // Extract inverse of width of zone in y-direction.
      const as3double ovDy = 1.0/zone_bound[1];

      // Interface location in y-direction.
      const as3double y0 = ( TypeZone == ZONE_SOUTH ) ? yMin : yMax;

      // Populate damping data.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){

        // Extract coordinates.
        const as3double y = coord[1][l];

        // Compute the actual damping function.
        DampDOFsInt[0][l] = 0.0;
        DampDOFsInt[1][l] = sigma*pow( ovDy*fabs(y - y0), beta );
      }

      break;
    }

    // In these zones: sigma(x) and sigma(y) are non-zero.
    case(ZONE_CORNER_0): case(ZONE_CORNER_1): case(ZONE_CORNER_2): case(ZONE_CORNER_3):
    {
      // Initialize needed damping dimensions.
      DampDOFsInt.resize(2, nullptr);
      // Allocate actual data.
      DampDOFsInt[0] = new as3double[nDOFsInt2D]();
      DampDOFsInt[1] = new as3double[nDOFsInt2D]();

      // Extract inverse of width of zone in x-direction.
      const as3double ovDx = 1.0/zone_bound[0];
      // Extract inverse of width of zone in y-direction.
      const as3double ovDy = 1.0/zone_bound[1];

      // Interface location in x- and y-directions.
      as3double x0, y0;
      if( (TypeZone == ZONE_CORNER_0) || (TypeZone == ZONE_CORNER_2) ) x0 = xMin; else x0 = xMax;
      if( (TypeZone == ZONE_CORNER_0) || (TypeZone == ZONE_CORNER_1) ) y0 = yMin; else y0 = yMax;

      // Populate damping data.
#pragma omp simd
      for(unsigned short l=0; l<nDOFsInt2D; l++){

        // Extract coordinates.
        const as3double x = coord[0][l];
        const as3double y = coord[1][l];

        // Compute the actual damping function.
        DampDOFsInt[0][l] = sigma*pow( ovDx*fabs(x - x0), beta );
        DampDOFsInt[1][l] = sigma*pow( ovDy*fabs(y - y0), beta );
      }

      break;
    }

    // If this is the main zone, exit immediately.
    default:
      Terminate("CEESpongeData::InitializeSpongeDamping", __FILE__, __LINE__,
                "Wrong zone type specified.");
  }

  // Check if allocation failed.
  for(unsigned short i=0; i<DampDOFsInt.size(); i++)
    if( !DampDOFsInt[i] )
      Terminate("CEESpongeData::InitializeSpongeDamping", __FILE__, __LINE__,
                "Allocation failed for DampDOFsInt.");
}


CEEPMLData::CEEPMLData
(
 CConfig   	   *config_container,
 CGeometry 	   *geometry_container,
 CInitial      *initial_container,
 CElement  	   *element_container,
 unsigned short iZone,
 unsigned long  iElem
)
	:
		CEESpongeData
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 iZone, iElem
		)
 /*
	* Constructor, used to initialize CEEPMLData per element.
	*/
{
  // Number of DOFS of integration points on an element in 1D.
  unsigned short nDOFsInt1D = element_container->GetnDOFsInt1D();

  // Initialize memory for the pseudo-mean flow on the surface integration points.
  DataDOFsIntMeanFace.resize(nFace);

  // Loop over all faces and reserve memory.
  for(unsigned short iFace=0; iFace<nFace; iFace++){

    // Reserve memory for each variable.
    DataDOFsIntMeanFace[iFace].resize(nVar, nullptr);

    // Initialize memory per each variable.
    for(unsigned short iVar=0; iVar<nVar; iVar++){

      // Allocate actual memory.
      DataDOFsIntMeanFace[iFace][iVar] = new as3double[nDOFsInt1D]();

      // Check if allocation failed.
      if( !DataDOFsIntMeanFace[iFace][iVar] )
        Terminate("CEEPMLData::CEEPMLData", __FILE__, __LINE__,
                  "Allocation failed for DataDOFsIntMeanFace.");
    }
  }

  // Initialize memory for the auxiliary variable on the faces, in case the
  // flow has a cross-flow direction. According to the expressions, this quantity
  // is differentiable w.r.t. y-direction, hence only the south and north faces
  // are needed for storage.
  if( config_container->GetCrossFlow() ){

    // Initialize memory for the auxiliary state at the surface.
    DataDOFsIntAuxFace.resize(nFace);

    // Loop over all faces, but allocate only for the ones with a y-parallel normal.
    for(unsigned short iFace=0; iFace<nFace; iFace++){

      // Allocate only if this is a south or north face, since the term we are
      // dealing with is dQdy, hence nx = 0, always.
      if( (iFace == IDX_SOUTH) || (iFace == IDX_NORTH) ){

        // Reserve memory for each variable.
        DataDOFsIntAuxFace[iFace].resize(nVar, nullptr);

        // Initialize memory per each variable.
        for(unsigned short iVar=0; iVar<nVar; iVar++){

          // Allocate actual memory.
          DataDOFsIntAuxFace[iFace][iVar] = new as3double[nDOFsInt1D]();

          // Check if allocation failed.
          if( !DataDOFsIntAuxFace[iFace][iVar] )
            Terminate("CEEPMLData::CEEPMLData", __FILE__, __LINE__,
                      "Allocation failed for DataDOFsIntAuxFace.");
        }
      }
    }
  }
}


CEEPMLData::~CEEPMLData
(
 void
)
 /*
	* Destructor for CEEPMLData class, frees allocated memory.
	*/
{

}






