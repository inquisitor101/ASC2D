#include "solver_structure.hpp"




CSolver::CSolver
(
 CConfig   		  *config_container,
 CGeometry 		  *geometry_container,
 CInitial       *initial_container,
 CElement      **element_container,
 CSpatial      **spatial_container,
 unsigned short  iZone
)
 /*
	* Constructor, used to initialize CSolver in zone: iZone.
	*/
{
	// Zone ID.
	zoneID     = iZone;
	// Polynomial order.
	nPoly      = config_container->GetnPolySolZone(iZone);
	// Number of nodes for solution DOFs in 2D.
	nDOFsSol2D = element_container[iZone]->GetnDOFsSol2D();
	// Number of elements.
	nElem      = geometry_container->GetGeometryZone(iZone)->GetnElem();
  // Number of boundaries in this zone. Note, this is always fixed as 4.
	nBoundary  = nFace;

  // Initialize element data container.
  Data_Preprocessing(config_container,
                     geometry_container,
                     initial_container,
                     element_container[iZone],
                     iZone);

  // Initialize boundary data container.
  Boundary_Preprocessing(config_container,
                         geometry_container,
                         initial_container,
                         element_container,
                         iZone);
}


CSolver::~CSolver
(
 void
)
 /*
	* Destructor for CSolver class, frees allocated memory.
	*/
{
  if( InvMassMatrix != nullptr ) delete [] InvMassMatrix;

	for(unsigned long i=0; i<data_container.size(); i++)
		if( data_container[i] ) delete data_container[i];

	for(unsigned short i=0; i<boundary_container.size(); i++)
		if( boundary_container[i] ) delete boundary_container[i];
}


void CSolver::Data_Preprocessing
(
	CConfig       *config_container,
	CGeometry     *geometry_container,
  CInitial      *initial_container,
	CElement      *element_container,
	unsigned short iZone
)
 /*
	* Function that preprocess the data container.
	*/
{
	// Initialize element data structure container.
	data_container.resize(nElem, nullptr);

  // Flag for error detection.
  bool ErrorDetected = false;

  // Check what type of solver we are dealing with.
  switch( config_container->GetTypeSolver(iZone) ){

    // Pure Euler solver.
    case(SOLVER_EE):
    {
      // Check what type of buffer layer this is, if any.
      switch( config_container->GetTypeBufferLayer(iZone) ){

        // EE-type data container.
        case(NO_LAYER):
        {
          // Initialize the data structure.
          for(unsigned long iElem=0; iElem<nElem; iElem++){
        		// Data container of type Euler.
        		data_container[iElem] = new CEEData(config_container,
        																			  geometry_container,
                                                initial_container,
        																			  element_container,
        																			  iZone, iElem);

        		// Check if allocation failed.
        		if( !data_container[iElem] )
        			Terminate("CSolver::Data_Preprocessing", __FILE__, __LINE__,
        								"Allocation failed for data_container.");
        	}
          break;
        }

        // Sponge EE-type data container.
        case(SPONGE_LAYER):
        {
          // Initialize the data structure.
          for(unsigned long iElem=0; iElem<nElem; iElem++){
        		// Data container of type sponge Euler.
        		data_container[iElem] = new CEESpongeData(config_container,
              																			  geometry_container,
                                                      initial_container,
              																			  element_container,
              																			  iZone, iElem);

        		// Check if allocation failed.
        		if( !data_container[iElem] )
        			Terminate("CSolver::Data_Preprocessing", __FILE__, __LINE__,
        								"Allocation failed for data_container.");
        	}
          break;
        }

        // PML EE-type data container.
        case(PML_LAYER):
        {
          // Initialize the data structure.
          for(unsigned long iElem=0; iElem<nElem; iElem++){
        		// Data container of type PML Euler.
        		data_container[iElem] = new CEEPMLData(config_container,
            																			 geometry_container,
                                                   initial_container,
            																			 element_container,
            																			 iZone, iElem);

        		// Check if allocation failed.
        		if( !data_container[iElem] )
        			Terminate("CSolver::Data_Preprocessing", __FILE__, __LINE__,
        								"Allocation failed for data_container.");
        	}
          break;
        }

        // Otherwise, flag for an error.
        default: ErrorDetected = true;
      }

      break;
    }

    default:
      Terminate("CEESolver::Data_Preprocessing", __FILE__, __LINE__,
                "Type of solver specified is unknown.");
  }

  if( ErrorDetected )
    Terminate("CEESolver::Data_Preprocessing", __FILE__, __LINE__,
              "Type of buffer layer specified is unknown.");
}


void CSolver::Boundary_Preprocessing
(
 CConfig  	   *config_container,
 CGeometry	   *geometry_container,
 CInitial      *initial_container,
 CElement 	  **element_container,
 unsigned short iZone
)
 /*
	* Function that preprocess the boundary container and initialize it.
	*/
{
  // Header for output.
  std::cout << "----------------------------------------------"
               "----------------------------------------------\n";
  std::cout << "Processing boundary of type: "
            << DisplaySolverType( config_container->GetTypeSolver(iZone) )
            << " in: "
            << "iZone(" << iZone << "): "
            << DisplayTypeZone(config_container->GetTypeZone(iZone))
            << std::endl;

	// Initialize the needed number of boundaries in this zone.
	boundary_container.resize(nBoundary, nullptr);

	// Loop over every boundary and specify input condition.
	for(unsigned short iBoundary=0; iBoundary<nBoundary; iBoundary++){

		// Check which boundary condition to apply.
		switch( config_container->GetTypeBC(iZone, iBoundary) ){

			case(BC_INTERFACE):
			{
				// Initialize interface/periodic boundary.
        switch( config_container->GetTypeSolver(iZone) ){

          // Pure EE interface boundary.
          case( SOLVER_EE ):
          {
            // Check what type of buffer layer this is, if any.
            switch( config_container->GetTypeBufferLayer(iZone) ){

              // Physical or sponge EE-type spatial container.
              case(NO_LAYER): case(SPONGE_LAYER):
              {
                // Initialize interface/periodic boundary.
                boundary_container[iBoundary] = new CEEInterfaceBoundary(config_container,
                                                                         geometry_container,
                                                                         initial_container,
                                                                         element_container,
                                                                         iZone, iBoundary);
                break;
              }

              // PML EE-type spatial container.
              case(PML_LAYER):
              {
                // Initialize interface/periodic boundary.
                boundary_container[iBoundary] = new CEEPMLInterfaceBoundary(config_container,
                                                                            geometry_container,
                                                                            initial_container,
                                                                            element_container,
                                                                            iZone, iBoundary);
                break;
              }

              default:
                Terminate("CSolver::Boundary_Preprocessing", __FILE__, __LINE__,
                          "Unknown buffer layer specified for the interface boundary.");
            }

            break;
          }

          default:
            Terminate("CSolver::Boundary_Preprocessing", __FILE__, __LINE__,
                      "Unknown solver specified for the interface boundary.");
        }

        // Break out of the interface case.
        break;
			}

      case(BC_SYMMETRY):
      {
        // Initialize symmetry boundary.
        boundary_container[iBoundary] = new CEESymmetryBoundary(config_container,
                                                                geometry_container,
                                                                initial_container,
                                                                element_container,
                                                                iZone, iBoundary);
        break;
      }

      case(BC_CBC_OUTLET):
      {
        // Initialize characteristic outlet boundary.
				boundary_container[iBoundary] = new CEEOutletCBC(config_container,
																												 geometry_container,
                                                         initial_container,
																												 element_container,
																												 iZone, iBoundary);
        break;
      }

      case(BC_CBC_INLET):
      {
        // Initialize characteristic inlet boundary.
				boundary_container[iBoundary] = new CEEInletCBC(config_container,
																												geometry_container,
                                                        initial_container,
																												element_container,
																												iZone, iBoundary);
        break;
      }

      case(BC_STATIC_OUTLET):
      {
        // Initialize subsonic static-condition outlet boundary.
				boundary_container[iBoundary] = new CEEStaticOutletBoundary(config_container,
            																											  geometry_container,
                                                                    initial_container,
            																											  element_container,
            																											  iZone, iBoundary);
        break;
      }

      case(BC_STATIC_INLET):
      {
        // Initialize subsonic static-condition inlet boundary.
				boundary_container[iBoundary] = new CEEStaticInletBoundary(config_container,
            																											 geometry_container,
                                                                   initial_container,
            																											 element_container,
            																											 iZone, iBoundary);
        break;
      }

      case(BC_TOTAL_INLET):
      {
        // Initialize subsonic static-condition inlet boundary.
				boundary_container[iBoundary] = new CEETotalInletBoundary(config_container,
            																										  geometry_container,
                                                                  initial_container,
            																											element_container,
            																											iZone, iBoundary);
        break;
      }

      case(BC_SUPERSONIC_OUTLET):
      {
        // Initialize supersonic outlet boundary.
				boundary_container[iBoundary] = new CEESupersonicOutletBoundary(config_container,
                  																										  geometry_container,
                                                                        initial_container,
                  																											element_container,
                  																											iZone, iBoundary);
        break;
      }

      case(BC_SUPERSONIC_INLET):
      {
        // Initialize supersonic inlet boundary.
				boundary_container[iBoundary] = new CEESupersonicInletBoundary(config_container,
                																										   geometry_container,
                                                                       initial_container,
                																											 element_container,
                																											 iZone, iBoundary);
        break;
      }

			default:
				Terminate("CSolver::Boundary_Preprocessing", __FILE__, __LINE__,
									"Boundary condition is unknown!");
		}
	}

  // Report progress.
  std::cout << "Done." << std::endl;
}


CEESolver::CEESolver
(
 CConfig        *config_container,
 CGeometry      *geometry_container,
 CInitial       *initial_container,
 CElement      **element_container,
 CSpatial      **spatial_container,
 unsigned short  iZone
)
	:
		CSolver
		(
		 config_container,
		 geometry_container,
     initial_container,
		 element_container,
		 spatial_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CEESolver in zone: iZone.
	*/
{
  // Compute the inverse mass matrix in this zone.
  ComputeInvMassMatrix(config_container,
                       geometry_container,
                       element_container[iZone]);
}


CEESolver::~CEESolver
(
 void
)
 /*
	* Destructor for CEESolver class, frees allocated memory.
	*/
{

}


as3double CEESolver::ComputeTimeStep
(
  unsigned long                 iElem,
  const as3vector1d<as3double> &hElem
)
 /*
  * Function that computes the time step per input element.
  */
{
  // CFL-dependant coefficient.
  // const as3double f1 = 1.0/( nDim*(2.0*nPoly+1.0) );
  const as3double f1 = 1.0/(nPoly*nPoly);

  // Abbreviation.
  const as3double gm1  = GAMMA_MINUS_ONE;
  const as3double ovhx = 1.0/hElem[XDIM];
  const as3double ovhy = 1.0/hElem[YDIM];

  // Extract current element solution.
  auto& sol = data_container[iElem]->GetDataDOFsSol();

  // Inverse of spectral radius and time step.
  as3double radInv = 1.0, dt = 0.0;

  // Loop across all nodes and determine the largest eigenvalue.
  for(unsigned short l=0; l<nDOFsSol2D; l++){

    // Abbreviation.
    const as3double ovrho = 1.0/sol[0][l];
    // Compute primitive variables.
    const as3double u     = ovrho*sol[1][l];
    const as3double v     = ovrho*sol[2][l];
    const as3double p     = gm1*(sol[3][l]
                          - 0.5*(u*sol[1][l] + v*sol[2][l]) );

    // Compute the local speed of sound.
    const as3double a = sqrt(GAMMA*p*ovrho);

    // Compute largest eigenvalues.
    const as3double lmbx = fabs(u)+a;
    const as3double lmby = fabs(v)+a;

    // Compute the inverse of the spectral radius.
    radInv = std::max( radInv, ovhx*lmbx );
    radInv = std::max( radInv, ovhy*lmby );
    // radInv = std::max( radInv, ovhx*lmbx + ovhy*lmby );
  }

  // Compute actual time step by including the stability coefficient.
  dt = f1/radInv;

  // Return computed value.
  return dt;
}


void CEESolver::ComputeInvMassMatrix
(
  CConfig   *config_container,
  CGeometry *geometry_container,
  CElement  *element_container
)
 /*
  * Function that computes the inverse mass matrix.
  */
{
  // Obtain number of solution DOFs in 2D.
  unsigned short nDOFsSol2D = element_container->GetnDOFsSol2D();

  // Reserve memory needed for the inverse mass matrix.
  InvMassMatrix = new as3double[nDOFsSol2D*nDOFsSol2D]();

  // Check if allocation failed.
  if( InvMassMatrix == nullptr )
    Terminate("CEESolver::ComputeInvMassMatrix", __FILE__, __LINE__,
              "Allocation failed for InvMassMatrix.");

  // Compute mass matrix.
  element_container->ComputeMassMatrix(nDOFsSol2D, InvMassMatrix);

  // Compute inverse of mass matrix.
  ComputeInverseMatrix(nDOFsSol2D, InvMassMatrix);
}




