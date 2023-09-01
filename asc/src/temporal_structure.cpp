#include "temporal_structure.hpp"



CTemporal::CTemporal
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 CIteration                **iteration_container,
 CSolver                   **solver_container,
 CElement                  **element_container,
 CSpatial                  **spatial_container,
 as3vector2d<unsigned long> &input_MapGlobalToLocal
)
 /*
	* Constructor, used to initialize CTemporal.
	*/
{
	// Extract number of zones.
	nZone = config_container->GetnZone();

	// Initialize number of elements in zones.
	nElemZone.resize(nZone);

	// Extract number of elements per zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++)
		nElemZone[iZone] = solver_container[iZone]->GetnElem();

	// Initialize number of nodes in every element in each zone.
	nNodeZone.resize(nZone);

	// Extract number of nodes per element per zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++)
		nNodeZone[iZone] = solver_container[iZone]->GetnDOFsSol2D();

  // Total number of elements in all the grid zones combined.
  nElemTotal = input_MapGlobalToLocal.size();
  // Initialize global element ID mapper that converts from global to local
  // indices that have the local zone and its element index in it.
  MapGlobalToLocal = input_MapGlobalToLocal;

  // Initialize working array parameter dimensions.
  InitializeWorkArrayDimension(config_container, element_container);
}


CTemporal::~CTemporal
(
 void
)
 /*
	* Destructor for CTemporal class, frees allocated memory.
	*/
{

}


void CTemporal::InitializeWorkArrayDimension
(
  CConfig   *config_container,
  CElement **element_container
)
 /*
  * Function that defines the dimension parameters used in the working array.
  */
{
  // Determine largest integration DOFs in all zones and assign as the inner
  // dimension of the working array.
  nWorkingArrayDOFs = 0;
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    nWorkingArrayDOFs = std::max(nWorkingArrayDOFs,
                                 element_container[iZone]->GetnDOFsInt2D());

  // Determine max number of variables needed across all zones, used this as a
  // multiplier in the outer dimension of the working array.
  nWorkingArrayVar = nVar;

  // Determine max number of data entries needed, used as a multiplier in the
  // outer dimension of the working array. Initially, the below is that of a
  // physical domain:
  // [0]: volume integration,
  // [1]: derivative w.r.t. x for F-flux,
  // [2]: derivative w.r.t. y for G-flux.
  nWorkingArrayEntry = 3;


  // Loop over all the extra zones and check if this is a PML layer. If so, then
  // increment the number of entries in the working array as need be.
  as3vector1d<unsigned short> nEntries(nZone);
  for(unsigned short iZone=0; iZone<nZone; iZone++){
    switch( config_container->GetTypeBufferLayer(iZone) ){
      case(PML_LAYER):
      {
        // Account for the extra terms used in the auxiliary variables, such as:
        // [3]: volume integration for auxiliary variable,
        // [4]: derivative w.r.t. x for auxiliary terms,
        // [5]: derivative w.r.t. y for auxiliary terms,
        // [6]: source terms for physical residual,
        // [7]: source terms for auxiliary residual.
        nEntries[iZone] = nWorkingArrayEntry + 5;

        break;
      }

      case(SPONGE_LAYER):
      {
        // Add an entry for storing the source term, such as:
        // [3]: volume integration for damping terms.
        nEntries[iZone] = nWorkingArrayEntry + 1;

        break;
      }
    }
  }

  // Take the maximum number of entries as that of the working variables.
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    if( nEntries[iZone] > nWorkingArrayEntry )
      nWorkingArrayEntry = nEntries[iZone];
}


void CTemporal::InitializeWorkArray
(
  as3data1d<as3double> &work_array
)
 /*
  * Function that reserves and initializes the required memory for a work array.
  */
{
  // Total number of data needed for 1st/outer index in the working array.
  const unsigned short nDataOuter = nWorkingArrayEntry*nWorkingArrayVar;
  // Total number of data needed for 2nd/inner index in the working array.
  const unsigned short nDataInner = nWorkingArrayDOFs;

  // Initialize the work array needed to carry out a residual update over an
  // entire grid sweep iteration.
  work_array.resize(nDataOuter, nullptr);

  // Allocate data per every entry of the working array.
  for(unsigned short i=0; i<work_array.size(); i++){

    // Allocate the actual memory.
    work_array[i] = new as3double[nDataInner]();

    // Check if allocation failed.
    if( !work_array[i] )
      Terminate("CTemporal::InitializeWorkArray", __FILE__, __LINE__,
                "Allocation failed for work_array.");
  }
}


CLSRK4Temporal::CLSRK4Temporal
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 CIteration                **iteration_container,
 CSolver                   **solver_container,
 CElement                  **element_container,
 CSpatial                  **spatial_container,
 as3vector2d<unsigned long> &input_MapGlobalToLocal
)
	:
		CTemporal
		(
		 config_container,
		 geometry_container,
		 iteration_container,
		 solver_container,
		 element_container,
		 spatial_container,
     input_MapGlobalToLocal
		)
 /*
	* Constructor, used to initialize CLSRK4Temporal.
	*/
{
	// Initialize the LSRK4 coefficients.
	rk4a.resize(nStageRK);
	rk4b.resize(nStageRK);
	rk4c.resize(nStageRK);

	// Low-storage 4th-order Runge-Kutta coefficients.
	rk4a[0] =  0.0;
	rk4a[1] = -567301805773.0/1357537059087.0;
	rk4a[2] = -2404267990393.0/2016746695238.0;
	rk4a[3] = -3550918686646.0/2091501179385.0;
	rk4a[4] = -1275806237668.0/842570457699.0;

	rk4b[0] =  1432997174477.0/9575080441755.0;
	rk4b[1] =  5161836677717.0/13612068292357.0;
	rk4b[2] =  1720146321549.0/2090206949498.0;
	rk4b[3] =  3134564353537.0/4481467310338.0;
	rk4b[4] =  2277821191437.0/14882151754819.0;

	rk4c[0] =  0.0;
	rk4c[1] =  1432997174477.0/9575080441755.0;
	rk4c[2] =  2526269341429.0/6820363962896.0;
	rk4c[3] =  2006345519317.0/3224310063776.0;
	rk4c[4] =  2802321613138.0/2924317926251.0;


	// Reserve tentative data.
	DataDOFsSolTentative.resize(nZone);

	// Initialize every zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){

		// Reserve data in zone.
		DataDOFsSolTentative[iZone].resize(nElemZone[iZone]);

    // Get generic data container.
    auto& data_container = solver_container[iZone]->GetDataContainer();

		// Initialize every element.
		for(unsigned long iElem=0; iElem<nElemZone[iZone]; iElem++){

      // Extract the number of variables in total (physical + auxiliary).
      unsigned short nVarTotal = data_container[iElem]->GetnVarTotal();

			// Reserve data in element.
			DataDOFsSolTentative[iZone][iElem].resize(nVarTotal, nullptr);

			// Initialize every variable.
			for(unsigned short iVar=0; iVar<nVarTotal; iVar++)
				DataDOFsSolTentative[iZone][iElem][iVar] = new as3double[nNodeZone[iZone]]();
		}
	}
}


CLSRK4Temporal::~CLSRK4Temporal
(
 void
)
 /*
	* Destructor for CLSRK4Temporal class, frees allocated memory.
	*/
{
	for(unsigned short iZone=0; iZone<DataDOFsSolTentative.size(); iZone++)
		for(unsigned long iElem=0; iElem<DataDOFsSolTentative[iZone].size(); iElem++)
			for(unsigned short iVar=0; iVar<DataDOFsSolTentative[iZone][iElem].size(); iVar++)
				if( DataDOFsSolTentative[iZone][iElem][iVar] ) delete [] DataDOFsSolTentative[iZone][iElem][iVar];
}


void CLSRK4Temporal::TimeMarch
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CIteration            **iteration_container,
 CSolver               **solver_container,
 CElement              **element_container,
 CSpatial              **spatial_container,
 CInitial              **initial_container,
 as3double               physicalTime,
 as3double               dtTime,
 as3vector1d<as3double> &MonitoringData
)
 /*
	* Function that updates the simulation by a single LSRK4 step.
	*/
{
	// Loop over all RK stages.
	for(unsigned short iStageRK=0; iStageRK<nStageRK; iStageRK++){

		// Local physical time.
		const as3double localTime = physicalTime + rk4c[iStageRK]*dtTime;

		// Local RK coefficients.
		const as3double alpha = rk4a[iStageRK];
		const as3double beta  = rk4b[iStageRK];

    // Perform an entire LSRK4 sweep and update the residual.
    UpdateTime(config_container,
               geometry_container,
               iteration_container,
               solver_container,
               element_container,
               spatial_container,
               initial_container,
               localTime, dtTime,
               alpha, beta,
               MonitoringData);
	}
}


void CLSRK4Temporal::UpdateTime
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CIteration            **iteration_container,
 CSolver               **solver_container,
 CElement              **element_container,
 CSpatial              **spatial_container,
 CInitial              **initial_container,
 as3double               localTime,
 as3double               dtTime,
 as3double               alpha,
 as3double               beta,
 as3vector1d<as3double> &MonitoringData

)
 /*
	* Function that performs a single LSRK4 stage sweep and update the residual.
	*/
{
  // Initialize max of the Mach number squared.
  as3double M2max      = 0.0;
	// Initialize the RMS of the density.
	as3double RMS_ResRho = 0.0;

	// Compute the total number of DOFs in the main zone only.
	const as3double nDOFsZoneMain = nElemZone[ZONE_MAIN]*nNodeZone[ZONE_MAIN];


  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Work array needed for performing a single grid sweep.
    as3data1d<as3double> work_array;
    // Reserve needed memory for the working array.
    InitializeWorkArray(work_array);


  // Impose the boundary conditions in all boundary elements across all zones.
#ifdef HAVE_OPENMP
#pragma omp for schedule(dynamic), collapse(2)
#endif
    for(unsigned short iZone=0; iZone<nZone; iZone++){
      for(unsigned short iBoundary=0; iBoundary<nFace; iBoundary++){

        // Extract the relevant boundary container.
        auto* boundary_container = solver_container[iZone]->GetBoundaryContainer(iBoundary);

        // Impose boundary condition.
        boundary_container->ImposeBoundaryCondition(config_container,
                                                    geometry_container,
                                                    solver_container,
                                                    element_container,
                                                    spatial_container,
                                                    localTime);
      }
    }


    // Loop over all the elements in all the zones and compute their residual.
#ifdef HAVE_OPENMP
#pragma omp for schedule(dynamic), reduction(max:M2max)
#endif
    for(unsigned long i=0; i<nElemTotal; i++){

      // Extract local zone number.
      unsigned short iZone = MapGlobalToLocal[i][0];
      // Extract local zone element index.
      unsigned long  iElem = MapGlobalToLocal[i][1];

			// Temporary storage for the Monitoring values.
			as3vector1d<as3double> monitordata( MonitoringData.size(), 0.0 ); 

      // Preprocess the data and apply the boundary conditions.
      iteration_container[iZone]->Preprocess(config_container,
  																					 geometry_container,
  																					 solver_container,
  																					 element_container,
  																					 spatial_container,
                                             work_array,
  																					 localTime, iElem);

      // Perform a single grid sweep over a single unique element.
      iteration_container[iZone]->ComputeResidual(config_container,
      																						geometry_container,
      																						solver_container[iZone],
      																						element_container[iZone],
      																						spatial_container[iZone],
                                                  initial_container[iZone],
                                                  work_array,
      																						localTime, iElem,
                                                  monitordata);

      // Set the max of the Mach squared in this element.
      M2max = std::max(M2max, monitordata[0]);
    }

		// Loop over all the elements in the main zone only and compute the RMS
		// of the residual of the density.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), reduction(+:RMS_ResRho)
#endif
		for(unsigned long iElem=0; iElem<nElemZone[ZONE_MAIN]; iElem++){

			// Extract current data container.
			auto* data = solver_container[ZONE_MAIN]->GetDataContainer(iElem);

      // Extract current total residual.
      auto& resu = data->GetDataDOFsRes();

			// Loop over the nodes and compute the RMS of the variable.
			for(unsigned short iNode=0; iNode<nNodeZone[ZONE_MAIN]; iNode++){
				
				// Compute the contribution onto the overall RMS of error.
				RMS_ResRho += resu[0][iNode]*resu[0][iNode]; 
			}
		}

    // Loop over all elements in all zones and update the solution.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned long i=0; i<nElemTotal; i++){

      // Extract local zone number.
      unsigned short iZone = MapGlobalToLocal[i][0];
      // Extract local zone element index.
      unsigned long  iElem = MapGlobalToLocal[i][1];

      // Extract current data container.
      auto* data_container = solver_container[iZone]->GetDataContainer(iElem);

      // Number of nodes per zone.
      unsigned short nNode = nNodeZone[iZone];

      // Extract current total solution.
      auto& sol = data_container->GetDataDOFsSol();
      // Extract current total residual.
      auto& res = data_container->GetDataDOFsRes();
      // Extract tentative solution.
      auto& tmp = DataDOFsSolTentative[iZone][iElem];

      // Loop over every variable.
      for(unsigned short iVar=0; iVar<sol.size(); iVar++){

        // Loop over nodes and update the residual.
#pragma omp simd
        for(unsigned short iNode=0; iNode<nNode; iNode++){
          // Predictor step, update the tentative solution.
          tmp[iVar][iNode]  = alpha*tmp[iVar][iNode] + dtTime*res[iVar][iNode];
        }
#pragma omp simd
        for(unsigned short iNode=0; iNode<nNode; iNode++){
          // Corrector step, update the true solution.
          sol[iVar][iNode] += beta*tmp[iVar][iNode];
        }
      }
    }

    // Free the local work array.
    for(unsigned short i=0; i<work_array.size(); i++)
      if( work_array[i] ) delete [] work_array[i];

  } // End of OpenMP parallel region.

  // Assign the actual max of the Mach number.
  MonitoringData[0] = sqrt(M2max);
	// Assign the RMS of the density normalized.
	MonitoringData[1] = sqrt( RMS_ResRho/nDOFsZoneMain );
}


CSSPRK3Temporal::CSSPRK3Temporal
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 CIteration                **iteration_container,
 CSolver                   **solver_container,
 CElement                  **element_container,
 CSpatial                  **spatial_container,
 as3vector2d<unsigned long> &input_MapGlobalToLocal
)
	:
		CTemporal
		(
		 config_container,
		 geometry_container,
		 iteration_container,
		 solver_container,
		 element_container,
		 spatial_container,
     input_MapGlobalToLocal
		)
 /*
	* Constructor, used to initialize CSSPRK3Temporal.
	*/
{
	// Initialize the SSPRK3 coefficients.
	rk4a.resize(nStageRK);
	rk4b.resize(nStageRK);
	rk4c.resize(nStageRK);

	// Strong-stability-preserving 3rd-order Runge-Kutta coefficients.
	rk4a[0] = 1.0;
	rk4a[1] = 3.0/4.0;
	rk4a[2] = 1.0/3.0;

	rk4b[0] = 1.0;
	rk4b[1] = 1.0/4.0;
	rk4b[2] = 2.0/3.0;

	rk4c[0] = 0.0;
	rk4c[1] = 1.0;
	rk4c[2] = 1.0/2.0;


	// Reserve tentative data.
	DataDOFsSolTentative.resize(nZone);

	// Initialize every zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){

		// Reserve data in zone.
		DataDOFsSolTentative[iZone].resize(nElemZone[iZone]);

    // Get generic data container.
    auto& data_container = solver_container[iZone]->GetDataContainer();

		// Initialize every element.
		for(unsigned long iElem=0; iElem<nElemZone[iZone]; iElem++){

      // Extract the number of variables in total (physical + auxiliary).
      unsigned short nVarTotal = data_container[iElem]->GetnVarTotal();

			// Reserve data in element.
			DataDOFsSolTentative[iZone][iElem].resize(nVarTotal, nullptr);

			// Initialize every variable.
			for(unsigned short iVar=0; iVar<nVarTotal; iVar++)
				DataDOFsSolTentative[iZone][iElem][iVar] = new as3double[nNodeZone[iZone]]();
		}
	}
}


CSSPRK3Temporal::~CSSPRK3Temporal
(
 void
)
 /*
	* Destructor for CSSPRK3Temporal class, frees allocated memory.
	*/
{
	for(unsigned short iZone=0; iZone<DataDOFsSolTentative.size(); iZone++)
		for(unsigned long iElem=0; iElem<DataDOFsSolTentative[iZone].size(); iElem++)
			for(unsigned short iVar=0; iVar<DataDOFsSolTentative[iZone][iElem].size(); iVar++)
				if( DataDOFsSolTentative[iZone][iElem][iVar] ) delete [] DataDOFsSolTentative[iZone][iElem][iVar];
}


void CSSPRK3Temporal::TimeMarch
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CIteration            **iteration_container,
 CSolver               **solver_container,
 CElement              **element_container,
 CSpatial              **spatial_container,
 CInitial              **initial_container,
 as3double               physicalTime,
 as3double               dtTime,
 as3vector1d<as3double> &MonitoringData
)
 /*
	* Function that updates the simulation by a single SSPRK3 step.
	*/
{
	// Loop over all RK stages.
	for(unsigned short iStageRK=0; iStageRK<nStageRK; iStageRK++){

		// Local physical time.
		const as3double localTime = physicalTime + rk4c[iStageRK]*dtTime;

		// Local RK coefficients.
		const as3double alpha = rk4a[iStageRK];
		const as3double beta  = rk4b[iStageRK];

    // Perform an entire SSPRK3 sweep and update the residual.
    UpdateTime(config_container,
               geometry_container,
               iteration_container,
               solver_container,
               element_container,
               spatial_container,
               initial_container,
               localTime, dtTime,
               alpha, beta,
               MonitoringData);
	}
}


void CSSPRK3Temporal::UpdateTime
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CIteration            **iteration_container,
 CSolver               **solver_container,
 CElement              **element_container,
 CSpatial              **spatial_container,
 CInitial              **initial_container,
 as3double               localTime,
 as3double               dtTime,
 as3double               alpha,
 as3double               beta,
 as3vector1d<as3double> &MonitoringData

)
 /*
	* Function that performs a single SSPRK3 stage sweep and update the residual.
	*/
{
  // Initialize max of the Mach number squared.
  as3double M2max      = 0.0;
	// Initialize the RMS of the density.
	as3double RMS_ResRho = 0.0;

	// Compute the total number of DOFs in the main zone only.
	const as3double nDOFsZoneMain = nElemZone[ZONE_MAIN]*nNodeZone[ZONE_MAIN];

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Work array needed for performing a single grid sweep.
    as3data1d<as3double> work_array;
    // Reserve needed memory for the working array.
    InitializeWorkArray(work_array);


  // Impose the boundary conditions in all boundary elements across all zones.
#ifdef HAVE_OPENMP
#pragma omp for schedule(dynamic), collapse(2)
#endif
    for(unsigned short iZone=0; iZone<nZone; iZone++){
      for(unsigned short iBoundary=0; iBoundary<nFace; iBoundary++){

        // Extract the relevant boundary container.
        auto* boundary_container = solver_container[iZone]->GetBoundaryContainer(iBoundary);

        // Impose boundary condition.
        boundary_container->ImposeBoundaryCondition(config_container,
                                                    geometry_container,
                                                    solver_container,
                                                    element_container,
                                                    spatial_container,
                                                    localTime);
      }
    }


    // Loop over all the elements in all the zones and compute their residual.
#ifdef HAVE_OPENMP
#pragma omp for schedule(dynamic), reduction(max:M2max)
#endif
    for(unsigned long i=0; i<nElemTotal; i++){

      // Extract local zone number.
      unsigned short iZone = MapGlobalToLocal[i][0];
      // Extract local zone element index.
      unsigned long  iElem = MapGlobalToLocal[i][1];

			// Temporary storage for the Monitoring values.
			as3vector1d<as3double> monitordata( MonitoringData.size(), 0.0 ); 

      // Preprocess the data and apply the boundary conditions.
      iteration_container[iZone]->Preprocess(config_container,
  																					 geometry_container,
  																					 solver_container,
  																					 element_container,
  																					 spatial_container,
                                             work_array,
  																					 localTime, iElem);

      // Perform a single grid sweep over a single unique element.
      iteration_container[iZone]->ComputeResidual(config_container,
      																						geometry_container,
      																						solver_container[iZone],
      																						element_container[iZone],
      																						spatial_container[iZone],
                                                  initial_container[iZone],
                                                  work_array,
      																						localTime, iElem,
                                                  monitordata);

      // Set the max of the Mach squared in this element.
      M2max = std::max(M2max, monitordata[0]);
    }

		// Loop over all the elements in the main zone only and compute the RMS
		// of the residual of the density.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), reduction(+:RMS_ResRho)
#endif
		for(unsigned long iElem=0; iElem<nElemZone[ZONE_MAIN]; iElem++){

			// Extract current data container.
			auto* data = solver_container[ZONE_MAIN]->GetDataContainer(iElem);

      // Extract current total residual.
      auto& resu = data->GetDataDOFsRes();

			// Loop over the nodes and compute the RMS of the variable.
			for(unsigned short iNode=0; iNode<nNodeZone[ZONE_MAIN]; iNode++){
				
				// Compute the contribution onto the overall RMS of error.
				RMS_ResRho += resu[0][iNode]*resu[0][iNode]; 
			}
		}


    // Abbreviations for updating the residual according to a SSPRK3 scheme.
    const as3double oma = 1.0 - alpha;
    const as3double bdt = beta*dtTime;

    // Distinguish between the first stage evaluation and the rest.
    if( fabs( oma ) < 1.0e-10 ){

      // Loop over all elements in all zones and update the solution.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
      for(unsigned long i=0; i<nElemTotal; i++){

        // Extract local zone number.
        unsigned short iZone = MapGlobalToLocal[i][0];
        // Extract local zone element index.
        unsigned long  iElem = MapGlobalToLocal[i][1];

        // Extract current data container.
        auto* data_container = solver_container[iZone]->GetDataContainer(iElem);

        // Number of nodes per zone.
        unsigned short nNode = nNodeZone[iZone];

        // Extract current total solution.
        auto& sol = data_container->GetDataDOFsSol();
        // Extract current total residual.
        auto& res = data_container->GetDataDOFsRes();
        // Extract tentative solution.
        auto& tmp = DataDOFsSolTentative[iZone][iElem];

        // Loop over every variable.
        for(unsigned short iVar=0; iVar<sol.size(); iVar++){

          // Loop over nodes and update the residual.
#pragma omp simd
          for(unsigned short iNode=0; iNode<nNode; iNode++){
            // Storage for initial solution at the first stage.
            tmp[iVar][iNode]  = sol[iVar][iNode];
            // Predictor-corrector step, update the true solution.
            sol[iVar][iNode] += bdt*res[iVar][iNode];
          }
        }
      }
    }
    else {
      // This is the second or third stage of the SSPRK3.

      // Loop over all elements in all zones and update the solution.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
      for(unsigned long i=0; i<nElemTotal; i++){

        // Extract local zone number.
        unsigned short iZone = MapGlobalToLocal[i][0];
        // Extract local zone element index.
        unsigned long  iElem = MapGlobalToLocal[i][1];

        // Extract current data container.
        auto* data_container = solver_container[iZone]->GetDataContainer(iElem);

        // Number of nodes per zone.
        unsigned short nNode = nNodeZone[iZone];

        // Extract current total solution.
        auto& sol = data_container->GetDataDOFsSol();
        // Extract current total residual.
        auto& res = data_container->GetDataDOFsRes();
        // Extract tentative solution.
        auto& tmp = DataDOFsSolTentative[iZone][iElem];

        // Loop over every variable.
        for(unsigned short iVar=0; iVar<sol.size(); iVar++){

          // Loop over nodes and update the residual.
#pragma omp simd
          for(unsigned short iNode=0; iNode<nNode; iNode++){
            // Predictor-corrector step, update the true solution.
            sol[iVar][iNode] =   oma*sol[iVar][iNode]
                             + alpha*tmp[iVar][iNode]
                             +   bdt*res[iVar][iNode];
          }
        }
      }
    }


    // Free the local work array.
    for(unsigned short i=0; i<work_array.size(); i++)
      if( work_array[i] ) delete [] work_array[i];

  } // End of OpenMP parallel region.

  // Assign the actual max of the Mach number.
  MonitoringData[0] = sqrt(M2max);
	// Assign the RMS of the density normalized.
	MonitoringData[1] = sqrt( RMS_ResRho/nDOFsZoneMain );
}
