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

	// Check if there need be any probes.
	if( config_container->GetProbeSpecified() )
	{
		// Extract probe locations.
		auto ProbeLocation = config_container->GetProbeLocation();
		// Extract probe variables.
		auto ProbeVariable = config_container->GetProbeVariable();
		
		// Consistency check.
		if( ProbeLocation.size()%2 != 0 )
			Terminate("CEEProcess::CEEProcess", __FILE__, __LINE__,
					      "Probes must be of a pair each, i.e. (x1, y1), (x2, y2).");

		// Another consistency check.
		if( ProbeLocation.size() == 0 )
			Terminate("CEEProcess::CEEProcess", __FILE__, __LINE__,
					      "No probe locations were specified.");

		// Total number of probes, division by nDim since these are pair of coordinates.
		unsigned short nProbe    = ProbeLocation.size()/nDim;
		// Toral number of variables to probe.
		unsigned short nProbeVar = ProbeVariable.size();
		// Use 2 items to process (pressure and v-velocity), per probe.
		unsigned short nItem     = nProbeVar*nProbe;

		// Initialize dimension of the probe.
		probe_container.resize(nProbe, nullptr);
		
		// Loop and initialize each probe.
		for(unsigned short iProbe=0; iProbe<nProbe; iProbe++)
		{
		  // Initialize current probe coordinates.
		  as3vector1d<as3double> probe = { ProbeLocation[iProbe*nDim  ], 
				                               ProbeLocation[iProbe*nDim+1] };
		  // Create current probe sensor.
		  probe_container[iProbe] = new CProbe(config_container,
		                                       geometry_container,
		                                       element_container[iZone],
		                                       probe, ProbeVariable, iZone);
		}
	}

	// Check if there need be any processing based on a reflection coefficient.
	if( config_container->GetTypeProcessData() != PROCESS_NOTHING )
	{
		switch( config_container->GetTypeProcessData() )
		{
			case( PROCESS_RATIO_P_U ):
			{
				reflection_container.resize(1, nullptr);
				reflection_container[0] = new CRatioR1Reflection(config_container,
						                                             geometry_container,
																											   element_container[iZone],
																											   initial_container[iZone],
																											   solver_container[iZone],
																											   iZone);
				break;
			}
		
			case( PROCESS_RATIO_P_M ):
			{
				reflection_container.resize(1, nullptr);
				reflection_container[0] = new CRatioR2Reflection(config_container,
						                                             geometry_container,
																											   element_container[iZone],
																											   initial_container[iZone],
																											   solver_container[iZone],
																											   iZone);
				break;
			}

			case( PROCESS_RATIO_WM_WP ):
			{
				reflection_container.resize(1, nullptr);
				reflection_container[0] = new CRatioR3Reflection(config_container,
						                                             geometry_container,
																											   element_container[iZone],
																											   initial_container[iZone],
																											   solver_container[iZone],
																											   iZone);
				break;
			}

			case( PROCESS_RATIO_LM_LP ):
			{
				reflection_container.resize(1, nullptr);
				reflection_container[0] = new CRatioR4Reflection(config_container,
						                                             geometry_container,
																											   element_container[iZone],
																											   initial_container[iZone],
																											   solver_container[iZone],
																											   iZone);
				break;
			}

			case( PROCESS_WAVE_ENTROPY ):
			{
				reflection_container.resize(1, nullptr);
				reflection_container[0] = new CRatioR5Reflection(config_container,
						                                             geometry_container,
																											   element_container[iZone],
																											   initial_container[iZone],
																											   solver_container[iZone],
																											   iZone);
				break;
			}


			default:
				Terminate("CEEProcess::CEEProcess", __FILE__, __LINE__,
						      "Unknown reflection coefficient for processing data.");
		}
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
  for(unsigned short i=0; i<probe_container.size(); i++)
    if( probe_container[i] ) delete probe_container[i];

	for(unsigned short i=0; i<reflection_container.size(); i++)
    if( reflection_container[i] ) delete reflection_container[i];
}


void CEEProcess::ProcessData
(
  CConfig   *config_container,
  CGeometry *geometry_container,
  CElement  *element_container,
  CSpatial  *spatial_container,
  CSolver   *solver_container,
  CInitial  *initial_container,
	COutput   *output_container,
  as3double  localTime
)
 /*
  * Function that selects what/how to process the data.
  */
{
	// Sample from probe location, if specified.
	if( config_container->GetProbeSpecified() )
	{
		// Total number of probe points.
		const unsigned short nProbe = probe_container.size(); 
		// Total number of probe variables.
		const unsigned short nSize  = config_container->GetProbeVariable().size();

		// Obtain the parent directory for this process.
		std::string dir = config_container->GetOutputProcessedDirectory();
		// Append the current base filename to the directory.
		dir += "probe_";

		// Loop over all the probe locations and sample the required data.
		for(unsigned short iProbe=0; iProbe<nProbe; iProbe++)
		{
			as3vector1d<as3double> data = probe_container[iProbe]->SampleProbeData(config_container,
					                                                                   geometry_container,
																							                               element_container,
																							                               solver_container,
																							                               initial_container);
			// Deduce the current probe filename.
			std::stringstream ss; ss << dir << iProbe;

			// Write and append data to file.
			output_container->WriteDataToFile(config_container, 
					                              geometry_container,
																				ss.str().c_str(),
																			  config_container->GetTypeFileFormatProcessedProbed(),
																				localTime, data);
		}
	}

	// If there is nothing else to process, then return.
	if( config_container->GetTypeProcessData() == PROCESS_NOTHING ) return;


	// Reserve memory for the data required.
	as3vector1d<as3double> data(reflection_container.size(), 0.0);

	// Loop over all the reflection functions and compute their coefficients.
	for(unsigned short i=0; i<reflection_container.size(); i++)
	{
		// Compute the respective metric.
		data[i] = reflection_container[i]->ComputeCoefficient(config_container,
				                                                  geometry_container,
																													element_container,
																													solver_container);
	}

	
	// Obtain the parent directory for this process.
	std::string dir = config_container->GetOutputProcessedDirectory();
	// Append the current base filename to the directory.
	dir += "profile";

	// Write and append data to file.
	output_container->WriteDataToFile(config_container, 
			                              geometry_container,
																		dir.c_str(),
																	  config_container->GetTypeFileFormatProcessedProbed(),
																		localTime, data);
}



