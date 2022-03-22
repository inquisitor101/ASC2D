#include "geometry_structure.hpp"



CGeometry::CGeometry
(
 CConfig   *config_container,
 CElement **element_container
)
 /*
  * Constructor, Reads input grid.
  */
{
	// Number of zones.
	nZone      = config_container->GetnZone();
  // Extract solution polynomial, per zone.
  nPolyZone  = config_container->GetnPolySolZone();
  // Extract number of elements in x, per zone.
  nxElemZone = config_container->GetnxElemZone();
  // Extract number of elements in y, per zone.
  nyElemZone = config_container->GetnyElemZone();

  // Deduce total number of elements in each zone.
  nElemZone.resize(nZone);
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    nElemZone[iZone] = nxElemZone[iZone]*nyElemZone[iZone];

	// Reserve memory for geometry zone container.
	geometry_zone.resize(nZone);

  // Create each zone.
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    geometry_zone[iZone] = new CGeometryZone(config_container,
                                             element_container[iZone],
                                             iZone);

  // Compute the total number of nPointSubElemP1.
  ComputePointsSubElementsP1();

  // Identify boundary elements indices, per zone, per boundary.
  IdentifyBoundaryElementIndices();

  // Maps every face to its matching face.
  IdentifyMatchingFace();

  // Assigns the unit-normal to each face of any element.
  ComputeUnitNormal();
}


CGeometry::~CGeometry
(
  void
)
 /*
 * Destructor for geometry class, frees allocated memory.
  */
{
  for(unsigned short i=0; i<geometry_zone.size(); i++)
    if( geometry_zone[i] ) delete geometry_zone[i];
}


void CGeometry::ComputeUnitNormal
(
  void
)
 /*
  * Function that computes the unit-normal for a generic structured
  * Cartesian-based element.
  */
{
  // Initialize unit-normals.
  UnitNormal.resize(nFace);
  for(unsigned short iFace=0; iFace<nFace; iFace++)
    UnitNormal[iFace].resize(nDim);

  // Compute unit-normals following design convention.
  UnitNormal[IDX_SOUTH][XDIM] =  0.0; UnitNormal[IDX_SOUTH][YDIM] = -1.0;
  UnitNormal[IDX_NORTH][XDIM] =  0.0; UnitNormal[IDX_NORTH][YDIM] =  1.0;
  UnitNormal[IDX_EAST][XDIM]  =  1.0; UnitNormal[IDX_EAST][YDIM]  =  0.0;
  UnitNormal[IDX_WEST][XDIM]  = -1.0; UnitNormal[IDX_WEST][YDIM]  =  0.0;
}


void CGeometry::IdentifyMatchingFace
(
  void
)
 /*
  * Function that matches each face to its opposite matching counterpart.
  */
{
  // Initialize matching faces.
  MatchingFace.resize(nFace);
  // Link faces.
  MatchingFace[IDX_SOUTH] = IDX_NORTH;
  MatchingFace[IDX_NORTH] = IDX_SOUTH;
  MatchingFace[IDX_WEST]  = IDX_EAST;
  MatchingFace[IDX_EAST]  = IDX_WEST;
}


void CGeometry::IdentifyBoundaryElementIndices
(
  void
)
 /*
  * Function that determines the indices of the elements that share a boundary.
  */
{
  // Reserve memory for the number of zones.
  ElemBoundaryIndex.resize(nZone);

  // Local index for boundary-only elements.
  unsigned long idx;

  for(unsigned short iZone=0; iZone<nZone; iZone++){

    // Extract vector containing all number of elements per each boundary.
    auto& nsElem = geometry_zone[iZone]->GetnsElem();

    // Extract total number of elements in this zone.
    unsigned long nElem  = nElemZone[iZone];
    // Extract number of elements in x-direction.
    unsigned long nxElem = nxElemZone[iZone];
    // Extract number of elements in y-direction.
    unsigned long nyElem = nyElemZone[iZone];

    // Extract number of boundaries in this zone.
    unsigned short nBoundary = nsElem.size();

    // Reserve memory for the number of boundaries.
    ElemBoundaryIndex[iZone].resize(nBoundary);

    // Identify elements living on south-most boundary faces.
    ElemBoundaryIndex[iZone][IDX_SOUTH].resize(nsElem[IDX_SOUTH]);
    idx = 0; // Reset local index.
    for(unsigned long i=0; i<nxElem; i++)
      ElemBoundaryIndex[iZone][IDX_SOUTH][idx++] = i;

    // Identify elements living on north-most boundary faces.
    ElemBoundaryIndex[iZone][IDX_NORTH].resize(nsElem[IDX_NORTH]);
    idx = 0; // Reset local index.
    for(unsigned long i=nxElem*(nyElem-1); i<nElem; i++)
      ElemBoundaryIndex[iZone][IDX_NORTH][idx++] = i;

    // Identify elements living on west-most boundary faces.
    ElemBoundaryIndex[iZone][IDX_WEST].resize(nsElem[IDX_WEST]);
    idx = 0; // Reset local index.
    for(unsigned long i=0; i<nElem; i+=nxElem)
      ElemBoundaryIndex[iZone][IDX_WEST][idx++] = i;

    // Identify elements living on east-most boundary faces.
    ElemBoundaryIndex[iZone][IDX_EAST].resize(nsElem[IDX_EAST]);
    idx = 0; // Reset local index.
    for(unsigned long i=nxElem-1; i<nElem; i+=nxElem)
      ElemBoundaryIndex[iZone][IDX_EAST][idx++] = i;
  }
}


void CGeometry::ComputePointsSubElementsP1
(
  void
)
 /*
  * Function that computes the total number of points based on a nPoly=1
  * sub-elements per each zone. This is used in the VTK output.
  */
{
  // Compute total number of points in all zones, based on nPoly=1 sub-elements.
  nPointSubElemP1 = 0;
  for(unsigned short iZone=0; iZone<nZone; iZone++){

    // Compute number of nPoly=1 sub-elements.
    unsigned long nSubElemP1 = nPolyZone[iZone]*nPolyZone[iZone];

    // Total number of points, under the assumption of nPoly=1 sub-elements.
    nPointSubElemP1 += nSubElemP1*N_POINTS_QUADRILATERAL*nElemZone[iZone];
  }
}


CGeometryZone::CGeometryZone
(
 CConfig 			 *config_container,
 CElement      *element_container,
 unsigned short iZone
)
 /*
	* Constructor, initializes grid in each zone.
	*/
{
	// Set current zone ID.
	zoneID    = iZone;
	// Solution polynomial.
	nPolySol  = config_container->GetnPolySolZone()[iZone];

	// Number of solution nodes in 1D.
	nDOFsSol1D = nPolySol+1;
	// Number of solution nodes in 2D.
	nDOFsSol2D = nDOFsSol1D*nDOFsSol1D;

  // Extract number of elements in x.
  nxElem = config_container->GetnxElemZone()[iZone];
  // Extract number of elements in y.
  nyElem = config_container->GetnyElemZone()[iZone];
  // Deduce total number of elements.
  nElem  = nxElem*nyElem;

  // Initialize number of elements sharing the external boundary.
  nsElem.resize(nFace);
  nsElem[IDX_SOUTH] = nxElem; nsElem[IDX_NORTH] = nxElem;
  nsElem[IDX_EAST]  = nyElem; nsElem[IDX_WEST]  = nyElem;

  // Type of zone.
  TypeZone = config_container->GetTypeZone()[iZone];

  // Create grid in this zone.
  GenerateGridZone(config_container,
                   element_container);

  // Identify the internal/boundary element faces.
  IdentifyElementFaceType();

  // Identify the internal elements neighbors that share the same face.
  IdentifyNeighborInternalElement();
}


CGeometryZone::~CGeometryZone
(
 void
)
 /*
	* Destructor for CGeometryZone class, frees allocated memory.
	*/
{
  for(unsigned long i=0; i<geometry_element.size(); i++)
    if( geometry_element[i] ) delete geometry_element[i];
}


bool CGeometryZone::SearchElementProbe
(
 as3vector1d<as3double> probe,
 unsigned long         &index,
 bool                   unique
) const
 /*
  * Function that searches the grid zone to locate the element containing the probe.
  */
{
  // For now, only a unique value must be located.
  if( !unique )
    Terminate("CGeometryZone::SearchElementProbe", __FILE__, __LINE__,
              "Search must be unique to resume (for now).");

  // Number of element instances found.
  unsigned short nFound = 0;

  // Explicitly extract probing coordinates.
  const as3double xprobe = probe[0];
  const as3double yprobe = probe[1];

  // Loop over all elements.
  for(unsigned long iElem=0; iElem<nElem; iElem++){

    // Extract dimensions of current element.
    const as3double ymin = geometry_element[iElem]->GetCoordBoundary(IDX_SOUTH);
    const as3double ymax = geometry_element[iElem]->GetCoordBoundary(IDX_NORTH);
    const as3double xmin = geometry_element[iElem]->GetCoordBoundary(IDX_WEST);
    const as3double xmax = geometry_element[iElem]->GetCoordBoundary(IDX_EAST);

    // Check if element fits each dimension.
    if( xprobe > xmin )
      if( xprobe < xmax )
        if( yprobe > ymin )
          if( yprobe < ymax ) { index = iElem; nFound++;  }
  }

  // Check if this is a non-unique element ownership.
  if( unique && (nFound > 1) ){
    std::string message = "Probe: ("
                        + std::to_string(xprobe) + ", "
                        + std::to_string(yprobe) + ") "
                        + "exists in more than one element.";
    Terminate("CGeometryZone::SearchElementProbe", __FILE__, __LINE__, message);
  }

  // Specify whether or not the owner element has been found.
  bool found = (nFound == 1) ? true : false;

  // Return whether the element with the probe has been identified.
  return found;
}


void CGeometryZone::IdentifyNeighborInternalElement
(
  void
)
 /*
  * Function that determines the indices of internal element neighbors that
  * share the same face.
  */
{
  // Initialize total elements and their faces.
  IndexNeighborInternalElement.resize(nElem);

  // Upon initialization, assign index out of nElem size that way if things
  // go wrong, the solver should throw a segmentation fault and exit.
  for(unsigned long iElem=0; iElem<nElem; iElem++)
    IndexNeighborInternalElement[iElem].resize(nFace, 99*nElem);

  // Identify element whose face is matched with internal elements.
  for(unsigned long iElem=0; iElem<nElem; iElem++){

    // Match east-face, if internal.
    if( InternalElemFace[iElem][IDX_EAST] )
      IndexNeighborInternalElement[iElem][IDX_EAST] = iElem+1;

    // Match west-face, if internal.
    if( InternalElemFace[iElem][IDX_WEST] )
      IndexNeighborInternalElement[iElem][IDX_WEST] = iElem-1;

    // Match north-face, if internal.
    if( InternalElemFace[iElem][IDX_NORTH] )
      IndexNeighborInternalElement[iElem][IDX_NORTH] = iElem+nxElem;

    // Match south-face, if internal.
    if( InternalElemFace[iElem][IDX_SOUTH] )
      IndexNeighborInternalElement[iElem][IDX_SOUTH] = iElem-nxElem;
  }
}


void CGeometryZone::IdentifyElementFaceType
(
  void
)
 /*
  * Function that flags whether an element's face is internal or boundary.
  */
{
  // Initialize element face type.
  InternalElemFace.resize(nElem);
  for(unsigned long iElem=0; iElem<nElem; iElem++)
    InternalElemFace[iElem].resize(nFace, true);

  // Flag south-most boundary faces.
  for(unsigned long i=0; i<nxElem; i++)
    InternalElemFace[i][IDX_SOUTH] = false;
  // Flag north-most boundary faces.
  for(unsigned long i=nxElem*(nyElem-1); i<nElem; i++)
    InternalElemFace[i][IDX_NORTH] = false;
  // Flag west-most boundary faces.
  for(unsigned long i=0; i<nElem; i+=nxElem)
    InternalElemFace[i][IDX_WEST]  = false;
  // Flag east-most boundary faces.
  for(unsigned long i=nxElem-1; i<nElem; i+=nxElem)
    InternalElemFace[i][IDX_EAST]  = false;
}


void CGeometryZone::GenerateGridZone
(
  CConfig  *config_container,
  CElement *element_container
)
 /*
  * Function that generates the grid zone.
  */
{
  // Main/physical domain bounding box.
  auto MainBox = config_container->GetDomainBound();

  // Partition the main domain explicitly for readability.
  const as3double xmin0 = MainBox[0];
  const as3double xmax0 = MainBox[1];
  const as3double ymin0 = MainBox[2];
  const as3double ymax0 = MainBox[3];

  // Number of elements in x-direction in main zone.
  const unsigned long nx0 = config_container->GetnxElemZone()[0];
  // Number of elements in y-direction in main zone.
  const unsigned long ny0 = config_container->GetnyElemZone()[0];

  // Domain size in x-direction of main zone.
  const as3double lx0 = xmax0 - xmin0;
  // Domain size in y-direction of main zone.
  const as3double ly0 = ymax0 - ymin0;

  // Element size in x-direction in the main zone.
  const as3double hx0 = lx0 / (as3double) nx0;
  // Element size in y-direction in the main zone.
  const as3double hy0 = ly0 / (as3double) ny0;


  // Bounding grid dimensions: (xmin, xmax, ymin, ymax).
  as3vector1d<as3double> DomainBox(4);
  // Element sizes per dimension. Default: (hx0, hy0).
  as3vector2d<as3double> ElementSizeZone(nDim);
  ElementSizeZone[0].resize(nxElem);
  ElementSizeZone[1].resize(nyElem);

  // // Extract expansion ratio in main zone in x- and y-directions.
  // const as3double rx = 0.95;//config_container->GetExpansionRatioMainZone()[0];
  // const as3double ry = 1.0;//config_container->GetExpansionRatioMainZone()[1];
  //
  // // Compute geometric sum of ratios in x and y directions.
  // as3double wx = 0.0; for(unsigned long i=0; i<nx0; i++) wx += pow(rx, i);
  // as3double wy = 0.0; for(unsigned long j=0; j<ny0; j++) wy += pow(ry, j);
  // // Determine the starting element sizes.
  // const as3double dx0 = lx0/wx;
  // const as3double dy0 = ly0/wy;
  //
  // // Compute step size in x-direction in the main zone.
  // as3vector1d<as3double> hx(nx0); hx[0] = dx0;
  // for(unsigned long i=0; i<nx0-1; i++) hx[i+1] = hx[i]*rx;
  //
  // // Compute step size in y-direction in the main zone.
  // as3vector1d<as3double> hy(ny0); hy[0] = dy0;
  // for(unsigned long j=0; j<ny0-1; j++) hy[j+1] = hy[j]*ry;



  // Compute step size in x-direction in the main zone.
  as3vector1d<as3double> hx(nx0);
  // Compute step size in y-direction in the main zone.
  as3vector1d<as3double> hy(ny0);

  // Check whether this is a uniform grid or not.
  if( config_container->GetUniformGridResolution() ){

    // Specify constant resolution.
    for(unsigned long i=0; i<nx0; i++) hx[i] = hx0;
    for(unsigned long j=0; j<ny0; j++) hy[j] = hy0;
  }
  else {

    // Number of elements in first block in x-direction.
    const unsigned long nxb1 = config_container->GetnxBlockElem()[0];
    const unsigned long nxb2 = config_container->GetnxBlockElem()[1];
    // Number of elements in first block in y-direction.
    const unsigned long nyb1 = config_container->GetnyBlockElem()[0];
    const unsigned long nyb2 = config_container->GetnyBlockElem()[1];

    // Consistency check.
    assert( (nxb1+nxb2) == nx0 );
    assert( (nyb1+nyb2) == ny0 );

    // Expansion ratios.
    const as3double rxb1 = config_container->GetDomainExpansionRatio()[0];
    const as3double rxb2 = config_container->GetDomainExpansionRatio()[1];
    const as3double ryb1 = config_container->GetDomainExpansionRatio()[2];
    const as3double ryb2 = config_container->GetDomainExpansionRatio()[3];

    // Interface location.
    const as3double xi = config_container->GetBlockInterfaceLocation()[0];
    const as3double yi = config_container->GetBlockInterfaceLocation()[1];

    // Block domain sizes.
    const as3double lxb1 = xi - xmin0, lyb1 = yi - ymin0;
    const as3double lxb2 = xmax0 - xi, lyb2 = ymax0 - yi;

    // Assemble block information in vector format in x-direction.
    as3vector1d<as3double>     rxb = { rxb1, rxb2 };
    as3vector1d<as3double>     lxb = { lxb1, lxb2 };
    as3vector1d<unsigned long> nxb = { nxb1, nxb2 };
    as3vector1d<unsigned long> ixb = {    0, nxb1 };

    // Assemble block information in vector format in y-direction.
    as3vector1d<as3double>     ryb = { ryb1, ryb2 };
    as3vector1d<as3double>     lyb = { lyb1, lyb2 };
    as3vector1d<unsigned long> nyb = { nyb1, nyb2 };
    as3vector1d<unsigned long> iyb = {    0, nyb1 };

    // Loop over all blocks and compute element sizes.
    for(unsigned short iBlock=0; iBlock<2; iBlock++){

      // Explicitly extract information in this block domain.
      const as3double     rx = rxb[iBlock], ry = ryb[iBlock];
      const as3double     lx = lxb[iBlock], ly = lyb[iBlock];
      const unsigned long nx = nxb[iBlock], ny = nyb[iBlock];
      const unsigned long ix = ixb[iBlock], iy = iyb[iBlock];

      // Compute geometric sum of ratios in  both direction.
      as3double wx = 0.0; for(unsigned long i=0; i<nx; i++) wx += pow(rx, i);
      as3double wy = 0.0; for(unsigned long j=0; j<ny; j++) wy += pow(ry, j);
      // Determine the starting element sizes.
      const as3double dx0 = lx/wx, dy0 = ly/wy;

      // Compute step size in x-direction in the main zone.
      hx[ix] = dx0; for(unsigned long i=ix; i<ix+nx-1; i++) hx[i+1] = hx[i]*rx;
      hy[iy] = dy0; for(unsigned long i=iy; i<iy+ny-1; i++) hy[i+1] = hy[i]*ry;
    }
  }

  // Determine what zone we are dealing with.
  switch( TypeZone ){

    // Main zone.
    case(ZONE_MAIN):
    {
      // Use main-zone as domain box.
      DomainBox = MainBox;

      // Initialize the element size in the x-direction in the main zone.
      for(unsigned long iElem=0; iElem<nxElem; iElem++) ElementSizeZone[0][iElem] = hx[iElem];
      // Initialize the element size in the y-direction in the main zone.
      for(unsigned long jElem=0; jElem<nyElem; jElem++) ElementSizeZone[1][jElem] = hy[jElem];

      break;
    }

    // West zone.
    case(ZONE_WEST):
    {
      // Expansion ratio in x-dimension for the west.
      const as3double rx = config_container->GethElemRatioZone()[0];

      // Compute each element size in x-dimension.
      for(unsigned long iElem=0; iElem<nxElem; iElem++)
        ElementSizeZone[0][iElem] = hx[0]*pow(rx, nxElem-iElem);

      // Initialize the element size in the y-direction in the main zone.
      for(unsigned long jElem=0; jElem<nyElem; jElem++)
        ElementSizeZone[1][jElem] = hy[jElem];


      // Total zone width in x-dimension.
      as3double ww = 0.0; for(auto h : ElementSizeZone[0]) ww += h;

      // Bounding coordinates.
      DomainBox[0] = xmin0 - ww; // xmin.
      DomainBox[1] = xmin0;      // xmax.
      DomainBox[2] = ymin0;      // ymin.
      DomainBox[3] = ymax0;      // ymax.

      break;
    }

    // East zone.
    case(ZONE_EAST): {

      // Expansion ratio in x-dimension for the east.
      const as3double rx = config_container->GethElemRatioZone()[1];

      // Compute each element size in x-dimension.
      for(unsigned long iElem=0; iElem<nxElem; iElem++)
        ElementSizeZone[0][iElem] = hx[nx0-1]*pow(rx, iElem+1);

      // Initialize the element size in the y-direction in the main zone.
      for(unsigned long jElem=0; jElem<nyElem; jElem++)
        ElementSizeZone[1][jElem] = hy[jElem];


      // Total zone width in x-dimension.
      as3double ww = 0.0; for(auto h : ElementSizeZone[0]) ww += h;

      // Bounding coordinates.
      DomainBox[0] = xmax0;      // xmin.
      DomainBox[1] = xmax0 + ww; // xmax.
      DomainBox[2] = ymin0;      // ymin.
      DomainBox[3] = ymax0;      // ymax.

      break;
    }

    // South zone.
    case(ZONE_SOUTH): {

      // Expansion ratio in y-dimension for the south.
      const as3double ry = config_container->GethElemRatioZone()[2];

      // Initialize the element size in the x-direction in the main zone.
      for(unsigned long iElem=0; iElem<nxElem; iElem++)
        ElementSizeZone[0][iElem] = hx[iElem];

      // Compute each element size in y-dimension.
      for(unsigned long iElem=0; iElem<nyElem; iElem++)
        ElementSizeZone[1][iElem] = hy[0]*pow(ry, nyElem-iElem);

      // Total zone width in y-dimension.
      as3double ww = 0.0; for(auto h : ElementSizeZone[1]) ww += h;

      // Bounding coordinates.
      DomainBox[0] = xmin0;      // xmin.
      DomainBox[1] = xmax0;      // xmax.
      DomainBox[2] = ymin0 - ww; // ymin.
      DomainBox[3] = ymin0;      // ymax.

      break;
    }

    // North zone.
    case(ZONE_NORTH): {

      // Expansion ratio in y-dimension for the north.
      const as3double ry = config_container->GethElemRatioZone()[3];

      // Initialize the element size in the x-direction in the main zone.
      for(unsigned long iElem=0; iElem<nxElem; iElem++)
        ElementSizeZone[0][iElem] = hx[iElem];

      // Compute each element size in y-dimension.
      for(unsigned long iElem=0; iElem<nyElem; iElem++)
        ElementSizeZone[1][iElem] = hy[ny0-1]*pow(ry, iElem+1);

      // Total zone width in y-dimension.
      as3double ww = 0.0; for(auto h : ElementSizeZone[1]) ww += h;

      // Bounding coordinates.
      DomainBox[0] = xmin0;      // xmin.
      DomainBox[1] = xmax0;      // xmax.
      DomainBox[2] = ymax0;      // ymin.
      DomainBox[3] = ymax0 + ww; // ymax.

      break;
    }

    // Corner0 zone.
    case(ZONE_CORNER_0): {

      // Expansion ratio in x-dimension for the west.
      const as3double rx = config_container->GethElemRatioZone()[0];
      // Expansion ratio in y-dimension for the south.
      const as3double ry = config_container->GethElemRatioZone()[2];

      // Compute each element size in x-dimension.
      for(unsigned long iElem=0; iElem<nxElem; iElem++)
        ElementSizeZone[0][iElem] = hx[0]*pow(rx, nxElem-iElem);
      // Compute each element size in y-dimension.
      for(unsigned long iElem=0; iElem<nyElem; iElem++)
        ElementSizeZone[1][iElem] = hy[0]*pow(ry, nyElem-iElem);

      // Total zone width in x-dimension.
      as3double wwx = 0.0; for(auto hx : ElementSizeZone[0]) wwx += hx;
      // Total zone width in y-dimension.
      as3double wwy = 0.0; for(auto hy : ElementSizeZone[1]) wwy += hy;

      // Bounding coordinates.
      DomainBox[0] = xmin0 - wwx; // xmin.
      DomainBox[1] = xmin0;       // xmax.
      DomainBox[2] = ymin0 - wwy; // ymin.
      DomainBox[3] = ymin0;       // ymax.

      break;
    }

    // Corner1 zone.
    case(ZONE_CORNER_1): {

      // Expansion ratio in x-dimension for the east.
      const as3double rx = config_container->GethElemRatioZone()[1];
      // Expansion ratio in y-dimension for the south.
      const as3double ry = config_container->GethElemRatioZone()[2];

      // Compute each element size in x-dimension.
      for(unsigned long iElem=0; iElem<nxElem; iElem++)
        ElementSizeZone[0][iElem] = hx[nx0-1]*pow(rx, iElem+1);
      // Compute each element size in y-dimension.
      for(unsigned long iElem=0; iElem<nyElem; iElem++)
        ElementSizeZone[1][iElem] = hy[0]*pow(ry, nyElem-iElem);

      // Total zone width in x-dimension.
      as3double wwx = 0.0; for(auto hx : ElementSizeZone[0]) wwx += hx;
      // Total zone width in y-dimension.
      as3double wwy = 0.0; for(auto hy : ElementSizeZone[1]) wwy += hy;

      // Bounding coordinates.
      DomainBox[0] = xmax0;       // xmin.
      DomainBox[1] = xmax0 + wwx; // xmax.
      DomainBox[2] = ymin0 - wwy; // ymin.
      DomainBox[3] = ymin0;       // ymax.

      break;
    }

    // Corner2 zone.
    case(ZONE_CORNER_2): {

      // Expansion ratio in x-dimension for the west.
      const as3double rx = config_container->GethElemRatioZone()[0];
      // Expansion ratio in y-dimension for the north.
      const as3double ry = config_container->GethElemRatioZone()[3];

      // Compute each element size in x-dimension.
      for(unsigned long iElem=0; iElem<nxElem; iElem++)
        ElementSizeZone[0][iElem] = hx[0]*pow(rx, nxElem-iElem);
      // Compute each element size in y-dimension.
      for(unsigned long iElem=0; iElem<nyElem; iElem++)
        ElementSizeZone[1][iElem] = hy[ny0-1]*pow(ry, iElem+1);

      // Total zone width in x-dimension.
      as3double wwx = 0.0; for(auto hx : ElementSizeZone[0]) wwx += hx;
      // Total zone width in y-dimension.
      as3double wwy = 0.0; for(auto hy : ElementSizeZone[1]) wwy += hy;

      // Bounding coordinates.
      DomainBox[0] = xmin0 - wwx; // xmin.
      DomainBox[1] = xmin0;       // xmax.
      DomainBox[2] = ymax0;       // ymin.
      DomainBox[3] = ymax0 + wwy; // ymax.

      break;
    }

    // Corner3 zone.
    case(ZONE_CORNER_3): {

      // Expansion ratio in x-dimension for the east.
      const as3double rx = config_container->GethElemRatioZone()[1];
      // Expansion ratio in y-dimension for the north.
      const as3double ry = config_container->GethElemRatioZone()[3];

      // Compute each element size in x-dimension.
      for(unsigned long iElem=0; iElem<nxElem; iElem++)
        ElementSizeZone[0][iElem] = hx[nx0-1]*pow(rx, iElem+1);
      // Compute each element size in y-dimension.
      for(unsigned long iElem=0; iElem<nyElem; iElem++)
        ElementSizeZone[1][iElem] = hy[ny0-1]*pow(ry, iElem+1);

      // Total zone width in x-dimension.
      as3double wwx = 0.0; for(auto hx : ElementSizeZone[0]) wwx += hx;
      // Total zone width in y-dimension.
      as3double wwy = 0.0; for(auto hy : ElementSizeZone[1]) wwy += hy;

      // Bounding coordinates.
      DomainBox[0] = xmax0;       // xmin.
      DomainBox[1] = xmax0 + wwx; // xmax.
      DomainBox[2] = ymax0;       // ymin.
      DomainBox[3] = ymax0 + wwy; // ymax.

      break;
    }

    default:
      Terminate("CGeometryZone::GenerateGridZone", __FILE__, __LINE__,
                "Unknown/wrong zone type detected");
  }

  // Assemble grid from elements, given the bounding box computed.
  AssembleGridZoneElements(config_container,
                           element_container,
                           DomainBox, ElementSizeZone);
}


void CGeometryZone::AssembleGridZoneElements
(
  CConfig                      *config_container,
  CElement                     *element_container,
  const as3vector1d<as3double> &BoundingBox,
  const as3vector2d<as3double> &ElementSizeZone
)
 /*
  * Function that assembles the elements in this grid zone.
  */
{
  // Obtain the degrees-of-freedom.
  const auto& rBasis = element_container->GetrDOFsSol1D();

  // Coordinates of the zone bounding box.
  const as3double xmin = BoundingBox[0]; // xmin
  const as3double xmax = BoundingBox[1]; // xmax
  const as3double ymin = BoundingBox[2]; // ymin
  const as3double ymax = BoundingBox[3]; // ymax

  // Record the width of the zone.
  ZoneSize = { xmax-xmin, ymax-ymin };

  // Local bounding box, per element: (xmin, xmax, ymin, ymax).
  as3vector1d<as3double> LocalBox(4);

  // Reserve memory for all the elements needed.
  geometry_element.resize(nElem, nullptr);

  // Loop over the grid zone and create each element.
  unsigned long idx = 0; as3double y0 = ymin;
  for(unsigned long jElem=0; jElem<nyElem; jElem++){
    as3double x0 = xmin;
    for(unsigned long iElem=0; iElem<nxElem; iElem++){

      // Extract current element sizes explicitly.
      const as3double hx = ElementSizeZone[0][iElem];
      const as3double hy = ElementSizeZone[1][jElem];

      // Determine current element bounding box.
      LocalBox[0] = x0;      // xmin.
      LocalBox[1] = x0 + hx; // xmax.
      LocalBox[2] = y0;      // ymin.
      LocalBox[3] = y0 + hy; // ymax.

      // Create geometry element.
      geometry_element[idx] = new CGeometryElement(config_container,
                                                   element_container,
                                                   rBasis, LocalBox);

      // Update global element index.
      idx++;

      // Update starting value of x-coordinate.
      x0 = LocalBox[1];
    }
    // Update starting value of y-coordinate.
    y0 = LocalBox[3];
  }

  // Consitency check, if grid matches specified input on its borders.
  bool ErrorDetected  = false;
  // Tolerance value.
  const as3double TOL = 1.0e-10;

  // Bottom-left node.
  unsigned short idx0 = 0;
  auto& coord0 = geometry_element[0]->GetCoordSolDOFs();
  if( (fabs(coord0[0][idx0]-xmin) > TOL) || (fabs(coord0[1][idx0]-ymin) > TOL) )
    ErrorDetected = true;

  // Top-right node.
  unsigned short idx1 = nDOFsSol2D-1;
  auto& coord1 = geometry_element[nElem-1]->GetCoordSolDOFs();
  if( (fabs(coord1[0][idx1]-xmax) > TOL) || (fabs(coord1[1][idx1]-ymax) > TOL) )
    ErrorDetected = true;

  // If an inconsistency is detected, terminate program.
  if( ErrorDetected )
    Terminate("CGeometryZone::AssembleGridZoneElements", __FILE__, __LINE__,
              "Grid does not abide by input bounding coordinates.");
}


CGeometryElement::CGeometryElement
(
  CConfig                      *config_container,
  CElement                     *element_container,
  const as3vector1d<as3double> &rBasis,
  const as3vector1d<as3double> &LocalBox
)
 /*
	* Constructor, initializes surface element.
	*/
{
  // Extract number of solution nodes in 1D.
  unsigned short nDOFsSol1D = rBasis.size();
  // Deduce number of solution nodes in 2D.
  unsigned short nDOFsSol2D = nDOFsSol1D*nDOFsSol1D;

  // Extract number of integration nodes in 1D.
  unsigned short nDOFsInt1D = element_container->GetnDOFsInt1D();
  // Extract number of integration nodes in 2D.
  unsigned short nDOFsInt2D = element_container->GetnDOFsInt2D();

  // Extract needed lagrange operators.
  auto* ellT = element_container->GetLagrangeInt1DTranspose();

  // Physical coordinates that bound the element (xmin, xmax, ymin, ymax).
  const as3double xmin = LocalBox[0];
  const as3double xmax = LocalBox[1];
  const as3double ymin = LocalBox[2];
  const as3double ymax = LocalBox[3];

  // Deduce element size.
  const as3double hx = xmax-xmin;
  const as3double hy = ymax-ymin;

  // Assign current element dimension.
  ElemSize = {hx, hy};

  // Resize local boundary coordinates.
  CoordBoundary.resize(nFace);
  // Populate bounding coordinates.
  CoordBoundary[IDX_SOUTH] = ymin;
  CoordBoundary[IDX_NORTH] = ymax;
  CoordBoundary[IDX_WEST]  = xmin;
  CoordBoundary[IDX_EAST]  = xmax;

  // Reserve memory for coordinates at solution points.
  CoordSolDOFs.resize(nDim, nullptr);
  for(unsigned short iDim=0; iDim<nDim; iDim++){

    // Allocate actual memory.
    CoordSolDOFs[iDim] = new as3double[nDOFsSol2D]();

    // Check if allocation failed.
    if( !CoordSolDOFs[iDim] )
      Terminate("CGeometryElement::CGeometryElement", __FILE__, __LINE__,
                 "Allocation failed for CoordSolDOFs.");
  }

  // Populate nodes.
  unsigned short idx = 0;
  for(unsigned short jNode=0; jNode<nDOFsSol1D; jNode++){
    for(unsigned short iNode=0; iNode<nDOFsSol1D; iNode++){
      CoordSolDOFs[0][idx] = xmin + 0.5*(1.0 + rBasis[iNode])*hx; // x-coord.
      CoordSolDOFs[1][idx] = ymin + 0.5*(1.0 + rBasis[jNode])*hy; // y-coord.
      // Update local index.
      idx++;
    }
  }

  // Reserve memory for coordinates at integration points.
  CoordIntDOFs.resize(nDim, nullptr);
  for(unsigned short iDim=0; iDim<nDim; iDim++){

    // Allocate actual memory.
    CoordIntDOFs[iDim] = new as3double[nDOFsInt2D]();

    // Check if allocation failed.
    if( !CoordIntDOFs[iDim] )
      Terminate("CGeometryElement::CGeometryElement", __FILE__, __LINE__,
                 "Allocation failed for CoordIntDOFs.");

  }

  // Interpolate the coordinates at the integration points.
  TensorProductSolAndGradVolume(nDOFsInt1D, nDim, nDOFsSol1D,
                                ellT, nullptr, CoordSolDOFs.data(),
                                CoordIntDOFs.data(),
                                nullptr, nullptr);
}


CGeometryElement::~CGeometryElement
(
 void
)
 /*
	* Destructor for CGeometryElement class, frees allocated memory.
	*/
{
  for(unsigned short i=0; i<CoordSolDOFs.size(); i++)
    if( CoordSolDOFs[i] ) delete [] CoordSolDOFs[i];

  for(unsigned short i=0; i<CoordIntDOFs.size(); i++)
    if( CoordIntDOFs[i] ) delete [] CoordIntDOFs[i];
}





