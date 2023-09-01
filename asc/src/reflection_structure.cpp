#include "reflection_structure.hpp"





CReflection::CReflection
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CElement      *element_container,
 CInitial      *initial_container,
 CSolver       *solver_container,
 unsigned short iZone
)
 /*
  * Constructor, used to initialize CReflection in zone: iZone.
  */
{
  // Extract zone ID.
  zoneID = iZone;

	// Issue warning if the mean flow of the IC used is spatially-dependent.
	if( config_container->GetTypeIC(iZone) == IC_VORTEX_ROLLUP )
		Terminate("CReflection::CReflection", __FILE__, __LINE__,
				      "Reflection coefficient cannot be used with spatially-varying mean flow.");
}


CReflection::~CReflection
(
 void
)
 /*
  * Destructor for CReflection class, frees allocated memory.
  */
{

}



CRatioR1Reflection::CRatioR1Reflection
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CElement      *element_container,
 CInitial      *initial_container,
 CSolver       *solver_container,
 unsigned short iZone
)
	:
		CReflection
		(
		 config_container,
		 geometry_container,
		 element_container,
		 initial_container,
		 solver_container,
		 iZone
		)
 /*
  * Constructor, used to initialize CRatioR1Reflection in zone: iZone.
  */
{
  // Extract mean data.
	uInf = initial_container->GetUinf();
  pInf = initial_container->GetPinf();
}


CRatioR1Reflection::~CRatioR1Reflection
(
 void
)
 /*
  * Destructor for CRatioR1Reflection class, frees allocated memory.
  */
{

}


as3double CRatioR1Reflection::ComputeCoefficient
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSolver   *solver_container
)
 /*
	* Function which computes the reflection expression of R1. See header for definition.
	*/
{
	// Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

	// Abbreviation involving gamma.
	const as3double gm1  = GAMMA_MINUS_ONE;

	// Maximum pressure.
	as3double pmax = 0.0;
	// Maximum u-velocity.
	as3double umax = 0.0;

	// Determine the location of the processed data.
	switch( config_container->GetProcessLocation() )
	{

		// XMIN boundary.
		case( PROCESS_LOCATION_XMIN ):
		{
			// Extract node indices on xmin boundary.
			auto& NodeList = element_container->GetIndexDOFsSol(IDX_WEST); 
			// Extract element indices on xmin boundary.
			auto& ElemList = solver_container->GetBoundaryContainer(IDX_WEST)->GetElemIndexI();

			// Total number of elements.
			unsigned long nElem = ElemList.size();

			// Loop over the points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:umax)
#endif		
			for(unsigned long i=0; i<nElem; i++)
			{
				// Explicitly extract element index.
				unsigned long iElem = ElemList[i];

				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over the points and process the data.
				for(unsigned short l : NodeList)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max u-velocity.
					umax = std::max(umax, u);
				}
			}

			break; 
		}
		
		// XMAX boundary.
		case( PROCESS_LOCATION_XMAX ): 
		{
			// Extract node indices on xmax boundary.
			auto& NodeList = element_container->GetIndexDOFsSol(IDX_EAST); 
			// Extract element indices on xmax boundary.
			auto& ElemList = solver_container->GetBoundaryContainer(IDX_EAST)->GetElemIndexI();
	
			// Total number of elements.
			unsigned long nElem = ElemList.size();
		
			// Loop over the points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:umax)
#endif		
			for(unsigned long i=0; i<nElem; i++)
			{
				// Explicitly extract element index.
				unsigned long iElem = ElemList[i];

				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over the points and process the data.
				for(unsigned short l : NodeList)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

    			// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max u-velocity.
					umax = std::max(umax, u);
				}
			}

			break;
		}

		// YMIN boundary.
		case( PROCESS_LOCATION_YMIN ):
		{
			// Extract node indices on ymin boundary.
			auto& NodeList = element_container->GetIndexDOFsSol(IDX_SOUTH); 
			// Extract element indices on ymin boundary.
			auto& ElemList = solver_container->GetBoundaryContainer(IDX_SOUTH)->GetElemIndexI();
	
			// Total number of elements.
			unsigned long nElem = ElemList.size();
				
			// Loop over the points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:umax)
#endif
			for(unsigned long i=0; i<nElem; i++)
			{
				// Explicitly extract element index.
				unsigned long iElem = ElemList[i];

				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over the points and process the data.
				for(unsigned short l : NodeList)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max u-velocity.
					umax = std::max(umax, u);
				}
			}

			break;
		}

		// YMAX boundary.
		case( PROCESS_LOCATION_YMAX ):
		{
			// Extract node indices on ymax boundary.
			auto& NodeList = element_container->GetIndexDOFsSol(IDX_NORTH); 
			// Extract element indices on ymax boundary.
			auto& ElemList = solver_container->GetBoundaryContainer(IDX_NORTH)->GetElemIndexI();
	
			// Total number of elements.
			unsigned long nElem = ElemList.size();
						
			// Loop over the points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:umax)
#endif
			for(unsigned long i=0; i<nElem; i++)
			{
				// Explicitly extract element index.
				unsigned long iElem = ElemList[i];

				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over the points and process the data.
				for(unsigned short l : NodeList)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max u-velocity.
					umax = std::max(umax, u);
				}
			}

			break;
		}

		// ZONE_MAIN.
		case( PROCESS_LOCATION_DOMAIN ):
		{
			// Extract the total number of solution DOFs in this zone.
			unsigned short nNode = element_container->GetnDOFsSol2D();
			// Extract the total number of elements in this zone.
			unsigned long  nElem = solver_container->GetnElem();

			// Loop over ther points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:umax)
#endif
			for(unsigned long iElem=0; iElem<nElem; iElem++)
			{
				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over all the points of interest.
				for(unsigned short l=0; l<nNode; l++)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max u-velocity.
					umax = std::max(umax, u);
				}
			}
			break;
		}

		default:
			Terminate("CRatioR1Reflection::ComputeCoefficient", __FILE__, __LINE__,
					      "Unknown process location.");
	}

	// Compute the reflection coefficient: R_1(t).
	const as3double R1 = (pmax/pInf)/(umax/uInf);

	// Return the computed reflection coefficient.
	return R1; 
}



CRatioR2Reflection::CRatioR2Reflection
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CElement      *element_container,
 CInitial      *initial_container,
 CSolver       *solver_container,
 unsigned short iZone
)
	:
		CReflection
		(
		 config_container,
		 geometry_container,
		 element_container,
		 initial_container,
		 solver_container,
		 iZone
		)
 /*
  * Constructor, used to initialize CRatioR2Reflection in zone: iZone.
  */
{
  // Extract mean data.
	Minf = initial_container->GetMach();
  pInf = initial_container->GetPinf();
}


CRatioR2Reflection::~CRatioR2Reflection
(
 void
)
 /*
  * Destructor for CRatioR2Reflection class, frees allocated memory.
  */
{

}


as3double CRatioR2Reflection::ComputeCoefficient
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSolver   *solver_container
)
 /*
	* Function which computes the reflection expression of R2. See header for definition.
	*/
{
	// Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

	// Abbreviation involving gamma.
	const as3double gm1  = GAMMA_MINUS_ONE;

	// Maximum pressure.
	as3double pmax = 0.0;
	// Maximum Mach number.
	as3double mmax = 0.0;

	// Determine the location of the processed data.
	switch( config_container->GetProcessLocation() )
	{

		// XMIN boundary.
		case( PROCESS_LOCATION_XMIN ):
		{
			// Extract node indices on xmin boundary.
			auto& NodeList = element_container->GetIndexDOFsSol(IDX_WEST); 
			// Extract element indices on xmin boundary.
			auto& ElemList = solver_container->GetBoundaryContainer(IDX_WEST)->GetElemIndexI();
	
			// Total number of elements.
			unsigned long nElem = ElemList.size();
			
			// Loop over the points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:mmax)
#endif		
			for(unsigned long i=0; i<nElem; i++)
			{
				// Explicitly extract element index.
				unsigned long iElem = ElemList[i];

				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over the points and process the data.
				for(unsigned short l : NodeList)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Magnitude of the velocity squared.
					const as3double magu2 = u*u + v*v;
					// Speed of sound squared.
					const as3double a2    = GAMMA*p*ovrho;

					// Compute the Mach number.
					const as3double M     = sqrt( magu2/a2 );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max Mach number.
					mmax = std::max(mmax, M);
				}
			}

			break; 
		}
		
		// XMAX boundary.
		case( PROCESS_LOCATION_XMAX ): 
		{
			// Extract node indices on xmax boundary.
			auto& NodeList = element_container->GetIndexDOFsSol(IDX_EAST); 
			// Extract element indices on xmax boundary.
			auto& ElemList = solver_container->GetBoundaryContainer(IDX_EAST)->GetElemIndexI();
	
			// Total number of elements.
			unsigned long nElem = ElemList.size();
					
			// Loop over the points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:mmax)
#endif		
			for(unsigned long i=0; i<nElem; i++)
			{		
				// Explicitly extract element index.
				unsigned long iElem = ElemList[i];
				
				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over the points and process the data.
				for(unsigned short l : NodeList)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Magnitude of the velocity squared.
					const as3double magu2 = u*u + v*v;
					// Speed of sound squared.
					const as3double a2    = GAMMA*p*ovrho;

					// Compute the Mach number.
					const as3double M     = sqrt( magu2/a2 );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max Mach number.
					mmax = std::max(mmax, M);
				}
			}

			break;
		}

		// YMIN boundary.
		case( PROCESS_LOCATION_YMIN ):
		{
			// Extract node indices on ymin boundary.
			auto& NodeList = element_container->GetIndexDOFsSol(IDX_SOUTH); 
			// Extract element indices on ymin boundary.
			auto& ElemList = solver_container->GetBoundaryContainer(IDX_SOUTH)->GetElemIndexI();
			
			// Total number of elements.
			unsigned long nElem = ElemList.size();
			
			// Loop over the points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:mmax)
#endif
			for(unsigned long i=0; i<nElem; i++)
			{
				// Explicitly extract element index.
				unsigned long iElem = ElemList[i];

				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over the points and process the data.
				for(unsigned short l : NodeList)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Magnitude of the velocity squared.
					const as3double magu2 = u*u + v*v;
					// Speed of sound squared.
					const as3double a2    = GAMMA*p*ovrho;

					// Compute the Mach number.
					const as3double M     = sqrt( magu2/a2 );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max Mach number.
					mmax = std::max(mmax, M);	
				}
			}

			break;
		}

		// YMAX boundary.
		case( PROCESS_LOCATION_YMAX ):
		{
			// Extract node indices on ymax boundary.
			auto& NodeList = element_container->GetIndexDOFsSol(IDX_NORTH); 
			// Extract element indices on ymax boundary.
			auto& ElemList = solver_container->GetBoundaryContainer(IDX_NORTH)->GetElemIndexI();
				
			// Total number of elements.
			unsigned long nElem = ElemList.size();
					
			// Loop over the points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:mmax)
#endif
			for(unsigned long i=0; i<nElem; i++)
			{
				// Explicitly extract element index.
				unsigned long iElem = ElemList[i];

				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over the points and process the data.
				for(unsigned short l : NodeList)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Magnitude of the velocity squared.
					const as3double magu2 = u*u + v*v;
					// Speed of sound squared.
					const as3double a2    = GAMMA*p*ovrho;

					// Compute the Mach number.
					const as3double M     = sqrt( magu2/a2 );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max Mach number.
					mmax = std::max(mmax, M);	
				}
			}

			break;
		}

		// ZONE_MAIN.
		case( PROCESS_LOCATION_DOMAIN ):
		{
			// Extract the total number of solution DOFs in this zone.
			unsigned short nNode = element_container->GetnDOFsSol2D();
			// Extract the total number of elements in this zone.
			unsigned long  nElem = solver_container->GetnElem();

			// Loop over ther points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:pmax), reduction(max:mmax)
#endif
			for(unsigned long iElem=0; iElem<nElem; iElem++)
			{
				// Extract current solution data.
      	auto& sol = data_container[iElem]->GetDataDOFsSol();

				// Loop over all the points of interest.
				for(unsigned short l=0; l<nNode; l++)
				{
					// Extract primitive data.
        	const as3double rho   = sol[0][l];
        	const as3double ovrho = 1.0/rho;
        	const as3double u     = ovrho* sol[1][l];
        	const as3double v     = ovrho* sol[2][l];
        	const as3double p     = gm1*(  sol[3][l]
        	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

					// Magnitude of the velocity squared.
					const as3double magu2 = u*u + v*v;
					// Speed of sound squared.
					const as3double a2    = GAMMA*p*ovrho;

					// Compute the Mach number.
					const as3double M     = sqrt( magu2/a2 );

					// Determine the max pressure.
        	pmax = std::max(pmax, p);
					// Determine the max Mach number.
					mmax = std::max(mmax, M);	
				}
			}
			break;
		}

		default:
			Terminate("CRatioR2Reflection::ComputeCoefficient", __FILE__, __LINE__,
					      "Unknown process location.");
	}

	// Compute the reflection coefficient: R_2(t).
	const as3double R2 = (pmax/pInf)/(mmax/Minf);

	// Return the computed reflection coefficient.
	return R2;
}




CRatioR3Reflection::CRatioR3Reflection
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CElement      *element_container,
 CInitial      *initial_container,
 CSolver       *solver_container,
 unsigned short iZone
)
	:
		CReflection
		(
		 config_container,
		 geometry_container,
		 element_container,
		 initial_container,
		 solver_container,
		 iZone
		)
 /*
  * Constructor, used to initialize CRatioR3Reflection in zone: iZone.
  */
{
  // Extract mean data.
	uInf = initial_container->GetUinf();
	vInf = initial_container->GetVinf();
  pInf = initial_container->GetPinf();

	// Determine boundary ID needed for the current coefficient.
	switch( config_container->GetProcessLocation() )
	{
		case( PROCESS_LOCATION_XMIN ): { iBoundary = IDX_WEST;  break; }
		case( PROCESS_LOCATION_XMAX ): { iBoundary = IDX_EAST;  break; }
		case( PROCESS_LOCATION_YMIN ): { iBoundary = IDX_SOUTH; break; }
		case( PROCESS_LOCATION_YMAX ): { iBoundary = IDX_NORTH; break; }

		// ZONE_MAIN.
		case( PROCESS_LOCATION_DOMAIN ):
		{
			Terminate("CRatioR3Reflection::CRatioR3Reflection", __FILE__, __LINE__,
					      "Cannot process this condition over entire domain.");
		}

		default:
			Terminate("CRatioR3Reflection::CRatioR3Reflection", __FILE__, __LINE__,
					      "Unknown process location.");
	}
}


CRatioR3Reflection::~CRatioR3Reflection
(
 void
)
 /*
  * Destructor for CRatioR3Reflection class, frees allocated memory.
  */
{

}


as3double CRatioR3Reflection::ComputeCoefficient
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSolver   *solver_container
)
 /*
	* Function which computes the reflection expression of R3. See header for definition.
	*/
{
	// Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

	// Abbreviation involving gamma.
	const as3double gm1  = GAMMA_MINUS_ONE;

	// Maximum reflection coefficient.
	as3double rmax = 0.0;

	// Extract node indices on this boundary.
	auto& NodeList = element_container->GetIndexDOFsSol(iBoundary); 
	// Extract element indices on this boundary.
	auto& ElemList = solver_container->GetBoundaryContainer(iBoundary)->GetElemIndexI();
  // Extract the current unit-normal for this boundary.
  auto& unitnorm = geometry_container->GetUnitNormal(iBoundary);

	// Extract unit-normals explicitly.
	const as3double nx = unitnorm[0];
	const as3double ny = unitnorm[1];

	// Extract mean normal-velocity.
	const as3double unInf = nx*uInf + ny*vInf;
				
	// Total number of elements.
	unsigned long nElem = ElemList.size();
	
	// Loop over ther points and process the data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:rmax)
#endif
	for(unsigned long i=0; i<nElem; i++)
	{
		// Explicitly extract element index.
		unsigned long iElem = ElemList[i];

		// Extract current solution data.
  	auto& sol = data_container[iElem]->GetDataDOFsSol();

		// Loop over the points and process the data.
		for(unsigned short l : NodeList)
		{	
			// Extract primitive data.
    	const as3double rho   = sol[0][l];
    	const as3double ovrho = 1.0/rho;
    	const as3double u     = ovrho* sol[1][l];
    	const as3double v     = ovrho* sol[2][l];
    	const as3double p     = gm1*(  sol[3][l]
    	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

			// Compute the local speed of sound.
			const as3double a  = sqrt( GAMMA*p*ovrho );
			// Compute the normal velocity.
			const as3double un = nx*u + ny*v; 

			// Compute the fluctuations of the normal velocity and pressure.
			const as3double dun = un - unInf;
			const as3double dp  = p  - pInf;

			// Assemble the incoming and outgoing characteristic.
			const as3double wp = dp + rho*a*dun;
			const as3double wm = dp - rho*a*dun;

			// Acoustic reflection coefficient.
			const as3double r  = fabs( wm/wp );

			// Determine the max acoustic reflection.
    	rmax = std::max(rmax, r);
		}
	}


	// Compute the reflection coefficient: R_3(t).
	const as3double R3 = rmax;

	// Return the computed reflection coefficient.
	return R3;
}



CRatioR4Reflection::CRatioR4Reflection
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CElement      *element_container,
 CInitial      *initial_container,
 CSolver       *solver_container,
 unsigned short iZone
)
	:
		CReflection
		(
		 config_container,
		 geometry_container,
		 element_container,
		 initial_container,
		 solver_container,
		 iZone
		)
 /*
  * Constructor, used to initialize CRatioR4Reflection in zone: iZone.
  */
{
	// Determine boundary ID needed for the current coefficient.
	switch( config_container->GetProcessLocation() )
	{
		case( PROCESS_LOCATION_XMIN ): { iBoundary = IDX_WEST;  break; }
		case( PROCESS_LOCATION_XMAX ): { iBoundary = IDX_EAST;  break; }
		case( PROCESS_LOCATION_YMIN ): { iBoundary = IDX_SOUTH; break; }
		case( PROCESS_LOCATION_YMAX ): { iBoundary = IDX_NORTH; break; }

		// ZONE_MAIN.
		case( PROCESS_LOCATION_DOMAIN ):
		{
			Terminate("CRatioR4Reflection::CRatioR4Reflection", __FILE__, __LINE__,
					      "Cannot process this condition over entire domain.");
		}

		default:
			Terminate("CRatioR4Reflection::CRatioR4Reflection", __FILE__, __LINE__,
					      "Unknown process location.");
	}
}


CRatioR4Reflection::~CRatioR4Reflection
(
 void
)
 /*
  * Destructor for CRatioR4Reflection class, frees allocated memory.
  */
{

}


as3double CRatioR4Reflection::ComputeCoefficient
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSolver   *solver_container
)
 /*
	* Function which computes the reflection expression of R4. See header for definition.
	*/
{
	// Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

	// Abbreviation involving gamma.
	const as3double gm1  = GAMMA_MINUS_ONE;

	// Maximum reflection coefficient.
	as3double lmax = 0.0;

	// Extract node indices on this boundary.
	auto& NodeList   = element_container->GetIndexDOFsSol(iBoundary); 
	// Extract element indices on this boundary.
	auto& ElemList   = solver_container->GetBoundaryContainer(iBoundary)->GetElemIndexI();
  // Extract the current unit-normal for this boundary.
  auto& unitnorm   = geometry_container->GetUnitNormal(iBoundary);
	// Get the indices for the solution nodes on an element in the current boundary.
	auto& FaceIndexI = element_container->GetIndexDOFsSol(iBoundary);
	
	// Obtain the number of solution DOFs in 1D.
	const unsigned short nDOFsSol1D = element_container->GetnDOFsSol1D();

	// Get the lagrange interpolation function on SolDOFs, this is an identity matrix.
	auto *ell   = element_container->GetLagrangeSol1D();
	// Obtain the derivative over the face.
	auto* dells = element_container->GetDerLagrangeDOFsSol1DFace(iBoundary);
	// Get the transpose of the derivative of lagrange function on the solution DOFs in 1D.
	auto* dellT = element_container->GetDerLagrangeDOFsSol1DTranspose();
	
	// Extract unit-normals explicitly.
	const as3double nx = unitnorm[0];
	const as3double ny = unitnorm[1];

	// Total number of elements.
	unsigned long nElem = ElemList.size();

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel 
#endif
  {
		// Initialize dynamic array for computataion of the gradient and wave amplitudes.
		as3data1d<as3double> Lx(nVar, nullptr);
		as3data1d<as3double> Ly(nVar, nullptr);

		// Allocate dynamic memory per variable.
		for(int i=0; i<nVar; i++)
		{
			Lx[i] = new as3double[nDOFsSol1D]();
			Ly[i] = new as3double[nDOFsSol1D]();
			
			// Check if allocation failed.
			if( !Lx[i] || !Ly[i] )
				Terminate("CRatioR4Reflection::ComputeCoefficient", __FILE__, __LINE__,
						      "Could not allocate memory for Lx and Ly.");
		}

		// Loop over ther points and process the data.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), reduction(max:lmax)
#endif
		for(unsigned long i=0; i<nElem; i++)
		{
			// Explicitly extract element index.
			unsigned long iElem = ElemList[i];

			// Extract current solution data.
  		auto& sol = data_container[iElem]->GetDataDOFsSol();

			// Compute the derivative of the solution over the boundary.
			// Note, these are parametric gradients, however since we are only
			// interested in their ratio, then we can bypass the conversion into
			// Cartesian framework.
			TensorProductSolAndGradFace(iBoundary, nDOFsSol1D, nVar, nDOFsSol1D,
					                        FaceIndexI.data(), 
																	ell, dellT, dells,
																	sol.data(), nullptr,
																	Lx.data(), Ly.data());

			// Local nodal index counter for the Lx and Ly arrays.
			unsigned short idx = 0;

			// Loop over the points and process the data.
			for(unsigned short l : NodeList)
			{	
				// Extract primitive data.
  	  	const as3double rho   = sol[0][l];
  	  	const as3double ovrho = 1.0/rho;
  	  	const as3double u     = ovrho* sol[1][l];
  	  	const as3double v     = ovrho* sol[2][l];
  	  	const as3double p     = gm1*(  sol[3][l]
  	  	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

				// Compute the local speed of sound.
				const as3double a  = sqrt( GAMMA*p*ovrho );
				// Compute the normal velocity.
				const as3double un = nx*u + ny*v;
	    	// Kinetic energy.
  	  	const as3double ek = 0.5*(u*u + v*v);

				// Extract the derivative of the velocity and pressure in x-direction.
				const as3double dudx = ovrho*(    Lx[1][idx] 
						                 -          u*Lx[0][idx] );
				const as3double dvdx = ovrho*(    Lx[2][idx] 
						                 -          v*Lx[0][idx] );
				const as3double dpdx =   gm1*( ek*Lx[0][idx]
						                 -          u*Lx[1][idx]
														 -          v*Lx[2][idx]
														 +            Lx[3][idx] );

				// Extract the derivative of the velocity and pressure in y-direction.
				const as3double dudy = ovrho*(    Ly[1][idx] 
						                 -          u*Ly[0][idx] );
				const as3double dvdy = ovrho*(    Ly[2][idx] 
						                 -          v*Ly[0][idx] );
				const as3double dpdy =   gm1*( ek*Ly[0][idx]
						                 -          u*Ly[1][idx]
														 -          v*Ly[2][idx]
														 +            Ly[3][idx] );


				// Compute the normal velocity derivative in the normal direction.
				const as3double dundn = nx*dudx + ny*dvdy;

				// Abbreviation of repeating terms.
				const as3double ovradpdn = (ovrho/a)*( nx*dpdx + ny*dpdy );

				// Compute the acoustic wave amplitude: L(-). 
				// Note there is a "0.5" which is ommitted since we are only 
				// interested in the ratio of Lm/Lp.
				const as3double Lm = (un - a)*( ovradpdn - dundn ); 
				const as3double Lp = (un + a)*( ovradpdn + dundn ); 

				// Acoustic reflection coefficient.
				const as3double r  = fabs( Lm/Lp );
				
				// Determine the max acoustic reflection.
  	  	lmax = std::max(lmax, r);
			
				// Update local nodal coutner.
				idx++;
			}
		}

		// Free used memory.
		for(int i=0; i<Lx.size(); i++) if( Lx[i] ) delete [] Lx[i];
		for(int i=0; i<Ly.size(); i++) if( Ly[i] ) delete [] Ly[i];

	} // End of Open-MP parallel region.


	// Compute the reflection coefficient: R_4(t).
	const as3double R4 = lmax;

	// Return the computed reflection coefficient.
	return R4;
}



CRatioR5Reflection::CRatioR5Reflection
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CElement      *element_container,
 CInitial      *initial_container,
 CSolver       *solver_container,
 unsigned short iZone
)
	:
		CReflection
		(
		 config_container,
		 geometry_container,
		 element_container,
		 initial_container,
		 solver_container,
		 iZone
		)
 /*
  * Constructor, used to initialize CRatioR5Reflection in zone: iZone.
  */
{
	// Determine boundary ID needed for the current coefficient.
	switch( config_container->GetProcessLocation() )
	{
		case( PROCESS_LOCATION_XMIN ): { iBoundary = IDX_WEST;  break; }
		case( PROCESS_LOCATION_XMAX ): { iBoundary = IDX_EAST;  break; }
		case( PROCESS_LOCATION_YMIN ): { iBoundary = IDX_SOUTH; break; }
		case( PROCESS_LOCATION_YMAX ): { iBoundary = IDX_NORTH; break; }

		// ZONE_MAIN.
		case( PROCESS_LOCATION_DOMAIN ):
		{
			Terminate("CRatioR5Reflection::CRatioR5Reflection", __FILE__, __LINE__,
					      "Cannot process this condition over entire domain.");
		}

		default:
			Terminate("CRatioR5Reflection::CRatioR5Reflection", __FILE__, __LINE__,
					      "Unknown process location.");
	}
}


CRatioR5Reflection::~CRatioR5Reflection
(
 void
)
 /*
  * Destructor for CRatioR5Reflection class, frees allocated memory.
  */
{

}


as3double CRatioR5Reflection::ComputeCoefficient
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CElement  *element_container,
 CSolver   *solver_container
)
 /*
	* Function which computes the reflection expression of R5. See header for definition.
	*/
{
	// Extract data container of this zone.
  auto& data_container = solver_container->GetDataContainer();

	// Abbreviation involving gamma.
	const as3double gm1  = GAMMA_MINUS_ONE;

	// Maximum reflection coefficient.
	as3double lmax = 0.0;

	// Extract node indices on this boundary.
	auto& NodeList   = element_container->GetIndexDOFsSol(iBoundary); 
	// Extract element indices on this boundary.
	auto& ElemList   = solver_container->GetBoundaryContainer(iBoundary)->GetElemIndexI();
  // Extract the current unit-normal for this boundary.
  auto& unitnorm   = geometry_container->GetUnitNormal(iBoundary);
	// Get the indices for the solution nodes on an element in the current boundary.
	auto& FaceIndexI = element_container->GetIndexDOFsSol(iBoundary);
	
	// Obtain the number of solution DOFs in 1D.
	const unsigned short nDOFsSol1D = element_container->GetnDOFsSol1D();

	// Get the lagrange interpolation function on SolDOFs, this is an identity matrix.
	auto *ell   = element_container->GetLagrangeSol1D();
	// Obtain the derivative over the face.
	auto* dells = element_container->GetDerLagrangeDOFsSol1DFace(iBoundary);
	// Get the transpose of the derivative of lagrange function on the solution DOFs in 1D.
	auto* dellT = element_container->GetDerLagrangeDOFsSol1DTranspose();
	
	// Extract unit-normals explicitly.
	const as3double nx = unitnorm[0];
	const as3double ny = unitnorm[1];

	// Total number of elements.
	unsigned long nElem = ElemList.size();

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel 
#endif
  {
		// Initialize dynamic array for computataion of the gradient and wave amplitudes.
		as3data1d<as3double> Lx(nVar, nullptr);
		as3data1d<as3double> Ly(nVar, nullptr);

		// Allocate dynamic memory per variable.
		for(int i=0; i<nVar; i++)
		{
			Lx[i] = new as3double[nDOFsSol1D]();
			Ly[i] = new as3double[nDOFsSol1D]();
			
			// Check if allocation failed.
			if( !Lx[i] || !Ly[i] )
				Terminate("CRatioR4Reflection::ComputeCoefficient", __FILE__, __LINE__,
						      "Could not allocate memory for Lx and Ly.");
		}

		// Loop over ther points and process the data.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), reduction(max:lmax)
#endif
		for(unsigned long i=0; i<nElem; i++)
		{
			// Explicitly extract element index.
			unsigned long iElem = ElemList[i];

			// Extract current solution data.
  		auto& sol = data_container[iElem]->GetDataDOFsSol();

			// Compute the derivative of the solution over the boundary.
			TensorProductSolAndGradFace(iBoundary, nDOFsSol1D, nVar, nDOFsSol1D,
					                        FaceIndexI.data(), 
																	ell, dellT, dells,
																	sol.data(), nullptr,
																	Lx.data(), Ly.data());


    	// Extract current element size.
    	auto& ElemSize = geometry_container->GetGeometryZone(zoneID)->GetGeometryElem(iElem)->GetElemSize();
    	
			// Extract element sizes explicitly.
    	const as3double hx = ElemSize[0];
    	const as3double hy = ElemSize[1];

			// Transform the parametric gradients into Cartesian gradients.
			for(unsigned short i=0; i<nVar; i++)
			{
#pragma omp simd
				for(unsigned short l=0; l<nDOFsSol1D; l++)
				{
					Lx[i][l] *= 2.0/hx;
					Ly[i][l] *= 2.0/hy;
				}
			}


			// Local nodal index counter for the Lx and Ly arrays.
			unsigned short idx = 0;

			// Loop over the points and process the data.
			for(unsigned short l : NodeList)
			{	
				// Extract primitive data.
  	  	const as3double rho   = sol[0][l];
  	  	const as3double ovrho = 1.0/rho;
  	  	const as3double u     = ovrho* sol[1][l];
  	  	const as3double v     = ovrho* sol[2][l];
  	  	const as3double p     = gm1*(  sol[3][l]
  	  	                      - 0.5*(u*sol[1][l] + v*sol[2][l]) );

				// Compute the local speed of sound squared.
				const as3double a2 = GAMMA*p*ovrho;
				// Compute the normal velocity.
				const as3double un = nx*u + ny*v;
	    	// Kinetic energy.
  	  	const as3double ek = 0.5*(u*u + v*v);

				// Extract the derivative of the velocity and pressure in x-direction.
				const as3double drdx =            Lx[0][idx];
				const as3double dpdx =   gm1*( ek*Lx[0][idx]
						                 -          u*Lx[1][idx]
														 -          v*Lx[2][idx]
														 +            Lx[3][idx] );

				// Extract the derivative of the velocity and pressure in y-direction.
				const as3double drdy =            Ly[0][idx];
				const as3double dpdy =   gm1*( ek*Ly[0][idx]
						                 -          u*Ly[1][idx]
														 -          v*Ly[2][idx]
														 +            Ly[3][idx] );

				// Compute the normal derivative of the density and pressure.
				const as3double drdn = nx*drdx + ny*drdy;
				const as3double dpdn = nx*dpdx + ny*dpdy;

				// Compute the entropy derivative in the normal direction.
				const as3double dsdn = drdn - dpdn/a2;

				// Compute the entropy wave amplitude: L(s). 
				const as3double Ls = un*( dsdn ); 

				// Acoustic reflection coefficient.
				const as3double r  = fabs( Ls );
				
				// Determine the max acoustic reflection.
  	  	lmax = std::max(lmax, r);
			
				// Update local nodal coutner.
				idx++;
			}
		}

		// Free used memory.
		for(int i=0; i<Lx.size(); i++) if( Lx[i] ) delete [] Lx[i];
		for(int i=0; i<Ly.size(); i++) if( Ly[i] ) delete [] Ly[i];

	} // End of Open-MP parallel region.


	// Compute the reflection coefficient: R_5(t).
	const as3double R5 = lmax;

	// Return the computed reflection coefficient.
	return R5;
}


