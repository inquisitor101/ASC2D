#include "driver_structure.hpp"



CDriver::CDriver
(
 CConfig *config
)
 /*
	* Constructor, used to set-up, run and (post-)process data.
	*/
{
	// Set number of zones.
	nZone = config->GetnZone();

	// Initialize all containers to nullptr.
	SetContainers_Null();

  // Choose active config container as input container.
  config_container = config;

	// Preprocess element container.
  Element_Preprocessing(config_container);

  // Preprocess geometry container.
  Geometry_Preprocessing(config_container,
                         element_container);

  // Preprocess the parallelization. Note, this must be set before the
  // temporal container.
  Parallelization_Preprocessing();


  // Preprocess input container.
  Input_Preprocessing(config_container,
  										geometry_container);

  // Preprocess output container.
  Output_Preprocessing(config_container,
  										 geometry_container);

  // Preprocess initial container.
  Initial_Preprocessing(config_container,
  											geometry_container,
  											element_container);

  // Preprocess spatial container.
  Spatial_Preprocessing(config_container,
  											geometry_container,
  											element_container,
                        initial_container);

  // Preprocess solver container.
  Solver_Preprocessing(config_container,
  										 geometry_container,
                       initial_container,
  										 element_container,
  										 spatial_container);

  // Preprocess iteration container.
  Iteration_Preprocessing(config_container,
  												geometry_container,
  												solver_container,
  												element_container,
  												spatial_container);

	// Preprocess temporal container.
  Temporal_Preprocessing(config_container,
  											 geometry_container,
  											 iteration_container,
  											 solver_container,
  											 element_container,
  											 spatial_container);

   // Preprocess the process container.
   Process_Preprocessing(config_container,
                         geometry_container,
                         output_container,
                         element_container,
                         initial_container,
                         solver_container,
                         spatial_container);


  // Simulation start time.
  SimTimeStart = config_container->GetSimulationStartTime();
  // Simulation end time.
  SimTimeFinal = config_container->GetSimulationFinalTime();

  // Max temporal iterations specified.
  MaxTimeIter  = config_container->GetMaxIter();
}


CDriver::~CDriver
(
 void
)
 /*
	* Destructor for CDriver class, frees allocated memory.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Freeing memory..." << std::endl;


	std::cout << "  deleting CConfig........ ";
	if(config_container != nullptr) delete config_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CGeometry...... ";
	if(geometry_container != nullptr) delete geometry_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CTemporal...... ";
	if(temporal_container != nullptr) delete temporal_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CInput......... ";
	if(input_container != nullptr) delete input_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting COutput........ ";
	if(output_container != nullptr) delete output_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CIteration..... ";
	if(iteration_container != nullptr){
		for(unsigned short iZone=0; iZone<nZone; iZone++){
			if(iteration_container[iZone] != nullptr) delete iteration_container[iZone];
		}
		delete [] iteration_container;
	}
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CSolver........ ";
	if(solver_container != nullptr){
		for(unsigned short iZone=0; iZone<nZone; iZone++){
			if(solver_container[iZone] != nullptr) delete solver_container[iZone];
		}
		delete [] solver_container;
	}
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CInitial....... ";
	if(initial_container != nullptr){
		for(unsigned short iZone=0; iZone<nZone; iZone++){
			if(initial_container[iZone] != nullptr) delete initial_container[iZone];
		}
		delete [] initial_container;
	}
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CSpatial....... ";
	if(spatial_container != nullptr){
		for(unsigned short iZone=0; iZone<nZone; iZone++){
			if(spatial_container[iZone] != nullptr) delete spatial_container[iZone];
		}
		delete [] spatial_container;
	}
	std::cout << "Done." << std::endl;

  std::cout << "  deleting CElement....... ";
	if(element_container != nullptr){
		for(unsigned short iZone=0; iZone<nZone; iZone++){
			if(element_container[iZone] != nullptr) delete element_container[iZone];
		}
		delete [] element_container;
	}
	std::cout << "Done." << std::endl;

  std::cout << "  deleting CProcess....... ";
	if(process_container != nullptr){
		for(unsigned short iZone=0; iZone<nZone; iZone++){
			if(process_container[iZone] != nullptr) delete process_container[iZone];
		}
		delete [] process_container;
	}
	std::cout << "Done." << std::endl;

	// Set all containers deleted to nullptr.
	config_container    = nullptr;
	geometry_container  = nullptr;
	input_container     = nullptr;
	output_container    = nullptr;
	solver_container    = nullptr;
	element_container   = nullptr;
	temporal_container  = nullptr;
	iteration_container = nullptr;
	spatial_container   = nullptr;
	initial_container   = nullptr;
  process_container   = nullptr;

	std::cout << "Done." << std::endl;
}


void CDriver::Parallelization_Preprocessing
(
  void
)
 /*
  * Function that sets the necessary steps needed for a parallel implementation.
  */
{
  // Extract the number of elements per each zone individually.
  auto& nElemZone = geometry_container->GetnElemZone();

  // Compute total number of elements in all the grid zones combined.
  nElemTotal = 0;
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    nElemTotal += nElemZone[iZone];

  // Initialize global element ID mapper.
  MapGlobalToLocal.resize(nElemTotal);
  // For each global element index, reserve memory for two entries corresponding
  // to the zone and the local (zone) element indices.
  for(unsigned long i=0; i<nElemTotal; i++)
    MapGlobalToLocal[i].resize(2);

  // Assign unique indices that map in between global-to-local ID.
  unsigned long idx = 0;
  for(unsigned short iZone=0; iZone<nZone; iZone++){
    for(unsigned long iElem=0; iElem<nElemZone[iZone]; iElem++){
      MapGlobalToLocal[idx][0] = iZone;
      MapGlobalToLocal[idx][1] = iElem;

      // Update global index.
      idx++;
    }
  }

  // TODO
  // Add a checkup that makes sure the number of threads provided and the
  // number of elements in total is efficient, otherwise return a warning.
  // Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
#ifdef HAVE_OPENMP
  // Get max number of threads specified.
  unsigned short nThreads = omp_get_max_threads();
  std::cout << "This is a parallel implementation using: "
            << nThreads << " threads." << std::endl;

  // Check if the work load is divisible by the number of threads.
  if( nElemTotal%nThreads != 0 )
    std::cout << "\n**********************************************"
  						<< "**********************************************\n"
              << "Warning, inefficient implementation!\n"
              << "... Number of total elements is not a multiple of the number of threads."
              << "\n**********************************************"
              << "**********************************************\n" << std::endl;

  // Estimate computational work load of each thread in terms of elements.
  unsigned long WorkLoadElem = nElemTotal/omp_get_max_threads();

  // Report work load per element.
  std::cout << "Workload of elements shared among each thread is: "
            << WorkLoadElem << " [Elem/Thread]. " << std::endl;

  // Estimate computational work load of each thread in terms of DOFs.
  as3vector1d<unsigned long> WorkLoadDOFs(nThreads);

#pragma omp parallel for
  for(unsigned long i=0; i<nElemTotal; i++){

    // Extract local zone number.
    unsigned short iZone = MapGlobalToLocal[i][0];
    // Extract local zone element index.
    unsigned long iElem  = MapGlobalToLocal[i][1];

    // Thread index.
    unsigned short iThread = omp_get_thread_num();

    // Accumulate the number of solution DOFs per thread.
    WorkLoadDOFs[iThread] += geometry_container->GetGeometryZone(iZone)->GetnDOFsSol2D();
  }

  // Compute number of max digits needed for output.
  unsigned long nDigits = 0;
  for(unsigned short i=0; i<nThreads; i++)
    nDigits = std::max( nDigits, WorkLoadDOFs[i] );
  // Deduce max value needed for the digits width.
  nDigits = std::to_string(nDigits).size();

  // Total number of work load.
  std::cout << "Workload of solution DOFs shared among each thread is: " << std::endl;
  for(unsigned short i=0; i<nThreads; i++)
    std::cout << "  thread(" << i << ") has: " << std::setw(nDigits)
              << WorkLoadDOFs[i] << " [DOFsSol/Thread]. " << std::endl;
#else
  std::cout << "This is a serial implementation." << std::endl;
#endif
}


void CDriver::SetContainers_Null
(
 void
)
 /*
	* Function that initializes all containers to nullptr.
	*/
{
	// Containers with first dimension.
	config_container    = nullptr;
	geometry_container  = nullptr;
	temporal_container  = nullptr;
	input_container     = nullptr;
	output_container    = nullptr;
	solver_container    = nullptr;
	element_container   = nullptr;
	iteration_container = nullptr;
	spatial_container   = nullptr;
	initial_container   = nullptr;
  process_container   = nullptr;

	// Containers that are zone dependant.
	iteration_container = new CIteration*[nZone];
	solver_container  	= new CSolver*[nZone];
	element_container 	= new CElement*[nZone];
	spatial_container   = new CSpatial*[nZone];
	initial_container   = new CInitial*[nZone];
  process_container   = new CProcess*[nZone];

	// Containers per zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		iteration_container[iZone] = nullptr;
		solver_container[iZone]    = nullptr;
		element_container[iZone]   = nullptr;
		spatial_container[iZone]   = nullptr;
		initial_container[iZone]   = nullptr;
    process_container[iZone]   = nullptr;
	}
}


void CDriver::Geometry_Preprocessing
(
 CConfig   *config_container,
 CElement **element_container
)
 /*
	* Function that preprocesses the geometry container.
	*/
{
	// Assign a geometry container.
	geometry_container = new CGeometry(config_container,
																		 element_container);
}


void CDriver::Input_Preprocessing
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that preprocesses the input container.
	*/
{
	// Assign an input container.
  input_container = new CInput(config_container,
									 						 geometry_container);
}


void CDriver::Output_Preprocessing
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that preprocesses the output container.
	*/
{
	// Assign an output container.
	output_container = new COutput(config_container,
																 geometry_container,
																 MapGlobalToLocal);
}


void CDriver::Process_Preprocessing
(
 CConfig     *config_container,
 CGeometry   *geometry_container,
 COutput     *output_container,
 CElement   **element_container,
 CInitial   **initial_container,
 CSolver    **solver_container,
 CSpatial   **spatial_container
)
 /*
  * Function that preprocesses the process container.
  */
{
  // Assign process per each zone.
  for(unsigned short iZone=0; iZone<nZone; iZone++){

    // Check which process to specify.
    switch( config_container->GetTypeSolver(iZone) ){

      // Pure EE solver.
      case(SOLVER_EE):
      {
        // Assign process container.
        process_container[iZone] = new CEEProcess(config_container,
                                                  geometry_container,
                                                  output_container,
                                                  element_container,
                                                  initial_container,
                                                  solver_container,
                                                  spatial_container,
                                                  iZone);
        break;
      }

      // Otherwise, exit immediately.
      default:
        Terminate("CDriver::Process_Preprocessing", __FILE__, __LINE__,
                  "Process container for specified solver is not (yet) implemented!");

    }
  }
}


void CDriver::Initial_Preprocessing
(
 CConfig      *config_container,
 CGeometry    *geometry_container,
 CElement    **element_container
)
 /*
	* Function that preprocesses the initial condition container.
	*/
{
	// Assign initial condition per each zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){

		// Check which initial condition to specify.
		switch( config_container->GetTypeIC(iZone) ){

			// Gaussian IC.
			case(IC_GAUSSIAN_PRESSURE): case(IC_GAUSSIAN_PRESSURE_1D_X): case(IC_GAUSSIAN_PRESSURE_1D_Y):
			{
        // Assign initial container.
				initial_container[iZone] = new CGaussianInitial(config_container,
																								  			geometry_container,
																								  			element_container,
																								  			iZone);
				break;
			}

      // Isentropic vortex IC.
      case(IC_ISENTROPIC_VORTEX):
      {
        // Assign initial container.
        initial_container[iZone] = new CIsentropicVortexInitial(config_container,
                                                                geometry_container,
                                                                element_container,
                                                                iZone);
        break;
      }

      // Entropy wave IC.
			case(IC_ENTROPY_WAVE):
			{
        // Assign initial container.
				initial_container[iZone] = new CEntropyWave(config_container,
																						  			geometry_container,
																						  			element_container,
																						  			iZone);
				break;
			}

      // Vortex roll-up IC.
			case(IC_VORTEX_ROLLUP):
			{
        // Assign initial container.
				initial_container[iZone] = new CVortexRollup(config_container,
																						  			 geometry_container,
																						  			 element_container,
																						  			 iZone);
				break;
			}

      // Acoustic Gaussian plane wave IC.
      case(IC_ACOUSTIC_PLANE_WAVE):
      {
        // Assign initial container.
        initial_container[iZone] = new CAcousticPlane(config_container,
                                                      geometry_container,
                                                      element_container,
                                                      iZone);
        break;
      }

			// Otherwise, exit immediately.
			default:
				Terminate("CDriver::Initial_Preprocessing", __FILE__, __LINE__,
									"Initial condition prescribed is not (yet) implemented!");
		}
	}
}


void CDriver::Element_Preprocessing
(
 CConfig  *config_container
)
 /*
	* Function that preprocesses the element container.
	*/
{
	// Assign element per each zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++)
		element_container[iZone] = new CElement(config_container, iZone);
}


void CDriver::Temporal_Preprocessing
(
 CConfig      *config_container,
 CGeometry    *geometry_container,
 CIteration  **iteration_container,
 CSolver     **solver_container,
 CElement    **element_container,
 CSpatial    **spatial_container
)
 /*
	* Function that preprocesses the temporal container.
	*/
{
	// Assign temporal container.
	switch( config_container->GetTypeTemporalScheme() ){

		// Low-storage 4th-order 5-stage Runge-Kutta (explicit).
		case(TEMPORAL_SCHEME_LSRK4):
		{
      // Assign temporal container.
			temporal_container = new CLSRK4Temporal(config_container,
																							geometry_container,
																							iteration_container,
																							solver_container,
																							element_container,
																							spatial_container,
                                              MapGlobalToLocal);
			break;
		}

    // Strong-stability-preserving 3-stage Runge-Kutta (explicit).
		case(TEMPORAL_SCHEME_SSPRK3):
		{
      // Assign temporal container.
			temporal_container = new CSSPRK3Temporal(config_container,
																							 geometry_container,
																							 iteration_container,
																							 solver_container,
																							 element_container,
																							 spatial_container,
                                               MapGlobalToLocal);
			break;
		}

    // Otherwise, exit immediately.
		default:
			Terminate("CDriver::Temporal_Preprocessing", __FILE__, __LINE__,
								"Temporal scheme specified is not (yet) implemented!");
	}
}


void CDriver::Iteration_Preprocessing
(
 CConfig      *config_container,
 CGeometry    *geometry_container,
 CSolver     **solver_container,
 CElement    **element_container,
 CSpatial    **spatial_container
)
 /*
	* Function that preprocesses the iteration container.
	*/
{
	// Assign iteration per each zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){

		// Check which iteration to specify.
		switch( config_container->GetTypeSolver(iZone) ){

			// Euler solver.
			case(SOLVER_EE):
			{
        // Assign iteration container.
				iteration_container[iZone] = new CEEIteration(config_container,
																											geometry_container,
																											solver_container,
																											element_container,
																											iZone);
				break;
			}

			// Otherwise, exit immediately.
			default:
				Terminate("CDriver::Iteration_Preprocessing", __FILE__, __LINE__,
									"Iteration strategy for solver specified is not (yet) implemented!");
		}
	}
}


void CDriver::Solver_Preprocessing
(
 CConfig   	  *config_container,
 CGeometry 	  *geometry_container,
 CInitial    **initial_container,
 CElement  	 **element_container,
 CSpatial    **spatial_container
)
 /*
	* Function that preprocesses the solver container.
	*/
{
	// Assign solver per each zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){

		// Check which solver to specify.
		switch( config_container->GetTypeSolver(iZone) ){

			// Euler solver.
			case(SOLVER_EE):
			{
        // Assign solver container.
				solver_container[iZone] = new CEESolver(config_container,
																								geometry_container,
                                                initial_container[iZone],
																								element_container,
																								spatial_container,
																								iZone);
				break;
			}

			// Otherwise, exit immediately.
			default:
				Terminate("CDriver::Solver_Preprocessing", __FILE__, __LINE__,
									"Solver container specified is not (yet) implemented!");
		}
	}
}


void CDriver::Spatial_Preprocessing
(
 CConfig   	  *config_container,
 CGeometry 	  *geometry_container,
 CElement  	 **element_container,
 CInitial    **initial_container
)
 /*
	* Function that preprocesses the spatial discretization container.
	*/
{
  // Flag for error detection.
  bool ErrorDetected = false;

	// Assign spatial discretization per each zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){

		// Check which spatial discretization to specify.
		switch( config_container->GetTypeSolver(iZone) ){

			// Euler solver.
			case(SOLVER_EE):
			{

        // Check what type of buffer layer this is, if any.
        switch( config_container->GetTypeBufferLayer(iZone) ){

          // EE-type spatial container.
          case(NO_LAYER):
          {
            spatial_container[iZone] = new CEESpatial(config_container,
    																								  geometry_container,
    																								  element_container,
                                                      initial_container[iZone],
    																								  iZone);
            break;
          }

          // Sponge EE-type spatial container.
          case(SPONGE_LAYER):
          {
            spatial_container[iZone] = new CEESpongeSpatial(config_container,
          																								  geometry_container,
          																								  element_container,
                                                            initial_container[iZone],
          																								  iZone);
            break;
          }

          // PML EE-type spatial container.
          case(PML_LAYER):
          {
            spatial_container[iZone] = new CEEPMLSpatial(config_container,
        																								 geometry_container,
        																								 element_container,
                                                         initial_container[iZone],
        																								 iZone);
            break;
          }

          // Otherwise, flag for an error.
          default: ErrorDetected = true;
        }

        break;
      }

      // Otherwise, exit immediately.
  		default: ErrorDetected = true;
		}

	}

  if( ErrorDetected )
    Terminate("CDriver::Spatial_Preprocessing", __FILE__, __LINE__,
              "Spatial container for solver specified is not (yet) implemented!");
}


void CDriver::MonitorOutput
(
 unsigned long           iIter,
 as3double               time,
 as3double               dt,
 as3vector1d<as3double> &MonitoringData,
 bool                    MonitorData = true
)
 /*
	* Function that outputs the header of the information being displayed.
	*/
{
	// Number of output reports for monitoring progress.
	unsigned long nOutput = std::max(1ul, MaxTimeIter/100);
	// Compute number of max digits needed for output.
	unsigned long nDigits = std::to_string(MaxTimeIter).size();
  // Monitoring output frequency.
  unsigned long OutputFreq = config_container->GetOutputFreq();

	// Display header.
	if( iIter%(50*OutputFreq) == 0 ){
		std::cout << "**********************************************"
							<< "**********************************************" << std::endl;
		std::cout << " Iteration\tPhysical Time \t Time step \t Max(Mach) \t Res[RMS(rho)]" << std::endl;
		std::cout << "**********************************************"
							<< "**********************************************" << std::endl;
	}

	// Display data in this iteration.
	if( MonitorData ){
		if( (iIter%OutputFreq==0) || (MaxTimeIter<OutputFreq) ){

    	// Extract the maximum Mach number.
    	const as3double Mmax = MonitoringData[0];
			// Extract the RMS of the density residual.
			const as3double RMSr = MonitoringData[1];

			// Display progress.
			std::cout << std::scientific << "   "
								<< std::setw(nDigits) << iIter
								<< " \t "  << time
								<< " \t "  << dt
    	          << " \t "  << Mmax
								<< " \t "  << RMSr << std::endl;
		}

		// Check if there need be a gnuplot file written and write one.
		if( config_container->GetWriteGNUplot() )
			output_container->WriteGNUplot(config_container, iIter, time, MonitoringData);
	}


  // Check for floating-point errors at run-time.
#ifdef ENABLE_NAN_CHECK
    CheckFloatingError();
#endif
}


as3double CDriver::ComputeTimeStep
(
 void
)
 /*
	* Function that estimates the needed time step.
	*/
{
  // Check whether or not a fixed time step is specified.
  if( config_container->GetCFL() < 0.0 ){

    // Make sure no adaptive time-stepping is used.
    if( config_container->GetAdaptTime() )
      Terminate("CDriver::ComputeTimeStep", __FILE__, __LINE__,
                "This is a fixed-time-step simulation, there should be no adaptive time-stepping specified.");

    // Specify the time step input.
    const as3double deltaT = config_container->GetTimeStep();

    // Return data.
    return deltaT;
  }

	// For now, use a fixed time step input by the user.
	// as3double dt = config_container->GetTimeStep();

  // Initialize time step.
  as3double deltaT = 1.0e6;

  // Loop over all the elements in all the zones and compute time step.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(min:deltaT)
#endif
  for(unsigned long i=0; i<nElemTotal; i++){

    // Extract local zone number.
    unsigned short iZone = MapGlobalToLocal[i][0];
    // Extract local zone element index.
    unsigned long iElem  = MapGlobalToLocal[i][1];

    // Extract current element size.
    auto hElem = geometry_container->GetGeometryZone(iZone)->GetGeometryElem(iElem)->GetElemSize();

    // Compute time step each each element.
    const as3double dtElem = solver_container[iZone]->ComputeTimeStep(iElem, hElem);

    // Determine the largest possible stable time step.
    deltaT = std::min(dtElem, deltaT);

  } // End of parallel loop.

  // Account for the CFL number.
  deltaT *= config_container->GetCFL();

	// Check if user has specified a fixed time step, otherwise exit.
	if( deltaT < 0 )
		Terminate("CDriver::ComputeTimeStep", __FILE__, __LINE__,
							"Must expilictly specify a time step.");

  // Return estimated time step.
	return deltaT;
}


void CDriver::StartSolver
(
 void
)
 /*
	* Function that initiates the solver driver.
	*/
{
	// Gauge start time.
	as3double startTime = as3double(clock())/as3double(CLOCKS_PER_SEC);

	// Run a preprocessing step to initialize the solution and condition the data
	// in case there need be. Note, this in only executed once and before marching
	// in time.
	Preprocess();

	// Begin actual solver.
	Run();

	// Gauge end time used by the solver.
	as3double stopTime = as3double(clock())/as3double(CLOCKS_PER_SEC);

	// Lapse time used by the entire solver.
	as3double lapsedTime = stopTime - startTime;


	// Report lapsed time.
	std::cout << "\n% % % % % % % % % % % % % % % % %" << std::endl;
	std::cout << std::scientific << "lapsed time [sec]: " << lapsedTime << " %" << std::endl;
}


void CDriver::Preprocess
(
 void
)
 /*
	* Function that preprocesses the solver, before commencing with the iterative
	* solution.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Initializing solution... ";

	// Check whether a restart is issued, if so then read the solution from restart file first.
	if( config_container->GetRestartSolution() )
    input_container->ReadSolutionRestartFile(config_container,
                                             geometry_container,
                                             element_container,
                                             solver_container,
                                             SimTimeStart);

  // Loop over all the elements in all the zones and compute time step.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned long i=0; i<nElemTotal; i++){

    // Extract local zone number.
    unsigned short iZone = MapGlobalToLocal[i][0];
    // Extract local zone element index.
    unsigned long iElem  = MapGlobalToLocal[i][1];

    // Extract current grid zone.
  	auto* grid_zone = geometry_container->GetGeometryZone(iZone);
		// Extract data belonging element of interest.
		auto* data_element = solver_container[iZone]->GetDataContainer(iElem);
		// Extract element solution nodes.
    const auto* grid_element = grid_zone->GetGeometryElem(iElem);

    spatial_container[iZone]->InitializeSolution(config_container,
                                                 initial_container[iZone],
                                                 element_container[iZone],
                                                 grid_element,
                                                 data_element,
                                                 0.0);
  } // End of parallel loop.


	// Report progress.
	std::cout << "Done." << std::endl;

	// Write output VTK for initial condition.
	output_container->WriteFileVTK(config_container,
																 geometry_container,
																 element_container,
																 solver_container);

	// Write solution restart data to file.
  output_container->WriteSolutionToFile(config_container,
                                        geometry_container,
                                        element_container,
                                        solver_container,
                                        SimTimeStart);

	// Write the processed data, if specified. Note, this is only done in the main zone.
	if( config_container->GetTypeProcessData() != PROCESS_NOTHING || config_container->GetProbeSpecified() )
  	process_container[ZONE_MAIN]->ProcessData(config_container,
  	                                          geometry_container,
  	                                          element_container[ZONE_MAIN],
  	                                          spatial_container[ZONE_MAIN],
  	                                          solver_container[ZONE_MAIN],
  	                                          initial_container[ZONE_MAIN],
																							output_container,
  	                                          SimTimeStart);
}


void CDriver::Run
(
 void
)
 /*
	* Function that runs the entire solver.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
	for(unsigned short iZone=0; iZone<nZone; iZone++)
		switch( config_container->GetTypeSolver(iZone) ){
			case(SOLVER_EE):
        std::cout << "Beginning EE Solver in iZone(" << iZone << "): "
                  << DisplayTypeZone(config_container->GetTypeZone(iZone))
                  << std::endl;
      break;
		}
	std::cout << std::endl;

  // Output writing VTK file frequency.
  const unsigned long WriteVTKFreq     = config_container->GetWriteVTKFreq();
  // Output writing restart file frequency.
  const unsigned long WriteRestartFreq = config_container->GetWriteRestartFreq();
  // Filtering frequency.
  const unsigned long FilterFreq       = (unsigned long) config_container->GetFilterCharacteristics()[0];
  // Sample zone writing frequency.
  const unsigned long WriteFreqZone    = config_container->GetWriteFreqZoneData();

  // Monitoring data.
  // Thus far, use only [0]: max(Mach) and [1]: RMS(res[RHO]).
  as3vector1d<as3double> MonitoringData(2);

	// Estimate time step needed.
	as3double dt = ComputeTimeStep();

	// Current time.
	as3double SimTime       = SimTimeStart;
	// Current iteration.
	unsigned long IterCount = 0;
  // Marker to check if final time step is written or not.
  bool FinalStep          = false;


  // Check if there needs to be any sampling of entire zone data.
  if( config_container->GetSampleZoneData() )
    output_container->WriteZoneDataToFile(config_container,
                                          geometry_container,
                                          solver_container,
                                          SimTime);

	// Check if there needs to be any sampling for the boundary data. 
	// Note, this only applies to the main/internal zone.
	if( config_container->GetSampleSurfaceData() )
		for(unsigned short iSample=0; iSample<config_container->GetSampleDataBoundaryID().size(); iSample++)
			output_container->WriteBoundaryDataToFile(config_container,
					                                      geometry_container,
																								element_container[ZONE_MAIN],
																								initial_container[ZONE_MAIN],
																								solver_container[ZONE_MAIN],
																								iSample, SimTime, IterCount);


	// Display header for output format.
	MonitorOutput(IterCount, SimTime, dt, MonitoringData, false);

	// March in time, until target time is reached.
	while( (SimTime < SimTimeFinal) && (IterCount < MaxTimeIter) )
	{

		// Execute a single temporal update.
		temporal_container->TimeMarch(config_container,
																	geometry_container,
																	iteration_container,
																	solver_container,
																	element_container,
																	spatial_container,
                                  initial_container,
																	SimTime, dt,
                                  MonitoringData);

		// Update (physical) time.
		SimTime += dt;

		// Update iteration count.
		IterCount++;

    // Check if there needs be any filtering done as a processing step.
    if( IterCount%FilterFreq == 0 )
      for(unsigned short iZone=0; iZone<nZone; iZone++)
        if( config_container->GetTypeFilterSolution(iZone) != NO_FILTER )
          process_container[iZone]->FilterSolution(config_container,
                                                   geometry_container,
                                                   element_container[iZone],
                                                   solver_container[iZone],
                                                   spatial_container[iZone]);

    // Check if there needs to be any processing done. For now, limit the
    // processing to only the main physical zone.
		if( config_container->GetTypeProcessData() != PROCESS_NOTHING || config_container->GetProbeSpecified() )
      process_container[ZONE_MAIN]->ProcessData(config_container,
                                                geometry_container,
                                                element_container[ZONE_MAIN],
                                                spatial_container[ZONE_MAIN],
                                                solver_container[ZONE_MAIN],
                                                initial_container[ZONE_MAIN],
																								output_container,
                                                SimTime);

    // Check if there needs to be any sampling of entire zone data.
    if( config_container->GetSampleZoneData() && (IterCount%WriteFreqZone) == 0 )
      output_container->WriteZoneDataToFile(config_container,
                                            geometry_container,
                                            solver_container,
                                            SimTime);

	  // Check if there needs to be any sampling for the boundary data. 
	  // Note, this only applies to the main/internal zone.
	  if( config_container->GetSampleSurfaceData() )
	  	for(unsigned short iSample=0; iSample<config_container->GetSampleDataBoundaryID().size(); iSample++)
	  		output_container->WriteBoundaryDataToFile(config_container,
	  				                                      geometry_container,
	  																							element_container[ZONE_MAIN],
	  																							initial_container[ZONE_MAIN],
	  																							solver_container[ZONE_MAIN],
	  																							iSample, SimTime, IterCount);

    // Check if this is the final time step.
    if( (SimTime >= SimTimeFinal) || (IterCount >= MaxTimeIter) )
      FinalStep = true;

    // Output VTK data every input iterations.
    if( IterCount%WriteVTKFreq == 0 || FinalStep )
		{
      // Write output VTK for initial condition.
      output_container->WriteFileVTK(config_container,
    																 geometry_container,
																		 element_container,
    																 solver_container);		
		}

		// Output restart data every input iterations.
		if( IterCount%WriteRestartFreq == 0 || FinalStep )
		{
			// Write solution restart data to file.
      output_container->WriteSolutionToFile(config_container,
                                            geometry_container,
                                            element_container,
                                            solver_container,
                                            SimTime);
		}
	
		// If this is an adaptive time-stepping, then compute new time-step.
		if( config_container->GetAdaptTime() ) dt = ComputeTimeStep();

		// Display output for progress monitoring.
		MonitorOutput(IterCount, SimTime, dt, MonitoringData);
	}
}




