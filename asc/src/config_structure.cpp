#include "config_structure.hpp"



CConfig::CConfig
(
 const char *configFile
)
 /*
	* Constructor, reads the specified input configuration file.
	*/
{
  // Message stream.
	std::ostringstream message;

  // Check if file exists.
  std::ifstream inputFile(configFile);
  if( !inputFile.good() )
    Terminate("CConfig::CConfig", __FILE__, __LINE__,
              "File could not be opened!");


	if( !ReadGridOptions(configFile) ){
		message << "Failed to extract grid options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract solver options.
  if( !ReadSolverOptions(configFile) ){
		message << "Failed to read solver options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

 	// Extract input/output information.
	if( !ReadIOOptions(configFile) ){
		message << "Failed to read input/output options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

	// Extract boundary marker specification.
	if( !ReadBoundaryOptions(configFile) ){
		message << "Failed to read boundary marker options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract temporal information.
	if( !ReadTemporalOptions(configFile) ){
		message << "Failed to read temporal options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract flow information.
	if( !ReadFlowOptions(configFile) ){
		message << "Failed to read flow options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract initial condition information.
	if( !ReadICOptions(configFile) ){
		message << "Failed to read initial condition options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract data processing information.
  if( !ReadProcessingOptions(configFile) ){
    message << "Failed to read processing options from " << configFile;
    Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
  }

  // Extract modified boundary condition information.
  if( !ReadModifiedBCOptions(configFile) ){
    message << "Failed to read modified BC options from " << configFile;
    Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
  }

  // Close file.
  inputFile.close();
}


CConfig::~CConfig
(
 void
)
 /*
	* Destructor for CConfig, does nothing.
	*/
{

}

bool CConfig::ReadModifiedBCOptions
(
 const char *configFile
)
 /*
	* Function that reads the modified boundary condition specifications.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read which boundary condition to modify.
  AddVectorOption(paramFile, "BC_MODIFY_TYPE", NameTypeModifyBC, DefaultParam.NameTypeModifyBC, true);
  // Pad entries for NameTypeModifyBC.
  PadEntriesVectorData(NameTypeModifyBC, "BC_MODIFY_TYPE", nFace, 1);

  // Read temporal frequency of each modification per working variable.
  AddVectorOption(paramFile, "BC_MODIFY_FREQ", ModifyFreqBC, DefaultParam.ModifyFreqBC, true);
  // Pad entries for ModifyFreqBC.
  PadEntriesVectorData(ModifyFreqBC, "BC_MODIFY_FREQ", nVar, 1);

  // Read spatial width of each modification per working variable.
  AddVectorOption(paramFile, "BC_MODIFY_WIDTH", ModifyWidthBC, DefaultParam.ModifyWidthBC, true);
  // Pad entries for ModifyWidthBC.
  PadEntriesVectorData(ModifyWidthBC, "BC_MODIFY_WIDTH", nVar, 1);

  // Read spatial amplitude of each modification per working variable.
  AddVectorOption(paramFile, "BC_MODIFY_STRENGTH", ModifyStrengthBC, DefaultParam.ModifyStrengthBC, true);
  // Pad entries for ModifyStrengthBC.
  PadEntriesVectorData(ModifyStrengthBC, "BC_MODIFY_STRENGTH", nVar, 1);

  // Read pulse center(s).
  AddVectorOption(paramFile, "BC_MODIFY_CENTER", ModifyCenterBC, DefaultParam.ModifyCenterBC, true);

  // Read pulse center(s) shift per working variable.
  AddVectorOption(paramFile, "BC_MODIFY_SHIFT_CENTER", ModifyShiftCenterBC, DefaultParam.ModifyShiftCenterBC, true);
  // Pad entries for ModifyShiftCenterBC.
  PadEntriesVectorData(ModifyShiftCenterBC, "BC_MODIFY_SHIFT_CENTER", nVar, 1);

  // Assign TypeModifyBC map.
  MapTypeModifyBC();

  // Convert frequency to angular frequency.
  for(auto &omg : ModifyFreqBC)omg *= 2.0*PI_CONSTANT;

  // Make sure number of centers is a multiple of nDim.
  if(ModifyCenterBC.size()%2 != 0)
    Terminate("CConfig::ReadModifiedBCOptions", __FILE__, __LINE__,
              "ModifyCenterBC must be of size multiple of nDim,");


  // Close file.
  paramFile.close();

  // Return happily.
	return true;
}


bool CConfig::ReadTemporalOptions
(
 const char *configFile
)
 /*
	* Function that reads the temporal specifications.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

	// Read start time of the simulation.
	AddScalarOption(paramFile, "START_TIME", SimulationTime[0], true);
	// Read end time of the simulation.
	AddScalarOption(paramFile, "FINAL_TIME", SimulationTime[1], true);
	// Read time step selected.
	AddScalarOption(paramFile, "TIME_STEP", TimeStep, DefaultParam.TimeStep, true);
	// Read maximum time iteration.
	AddScalarOption(paramFile, "MAX_ITER", MaxIter, DefaultParam.MaxIter, true);
	// Read temporal scheme used in time marching.
	AddScalarOption(paramFile, "TIME_MARCHING", NameTemporalScheme, true);
  // Read input adaptive time step.
  AddBoolOption(paramFile, "ADAPT_TIME", AdaptTime, DefaultParam.NameAdaptTime, true);
  // Read Courant-Friedrichs-Lewy condition number.
  AddScalarOption(paramFile, "CFL_NUMBER", CFL, DefaultParam.CFL, true);
	// Assign TypeTemporalScheme map.
	MapTemporalScheme();

  // Close file.
  paramFile.close();

  // Return happily.
	return true;
}


bool CConfig::ReadGridOptions
(
 const char *configFile
)
 /*
	* Function that reads all grid information.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Number of zones expected.
	AddScalarOption(paramFile, "NUMBER_ZONE", nZone, DefaultParam.nZone, true);

  // Read zone regions.
	AddVectorOption(paramFile, "MARKER_ZONE", NameZoneMarker, DefaultParam.NameZoneMarker, true);
  // Pad entries for NameZoneMarker.
  PadEntriesVectorData(NameZoneMarker, "MARKER_ZONE", nZone);

  // Read type of DOFs.
	AddVectorOption(paramFile, "DOF_TYPE", NameNodalDOFs, DefaultParam.NameNodalDOFs, true);
  // Pad entries for NameNodalDOFs.
  PadEntriesVectorData(NameNodalDOFs, "DOF_TYPE", nZone, 1, 2);

  // Read solution polynomial orders in each zone.
	AddVectorOption(paramFile, "POLY_ORDER_SOL", nPolySolZone, true);
  // Pad entries for nPolySolZone.
  PadEntriesVectorData(nPolySolZone, "POLY_ORDER_SOL", nZone, 1, 2);

  // Read domain bounds in each zone.
	AddVectorOption(paramFile, "DOMAIN_BOUND", DomainBound, true);

  // Read number of elements in x-direction in each zone.
  AddVectorOption(paramFile, "NUMBER_XELEM", nxElemZone, true);
  // Pad entries for nxElemZone.
  PadEntriesVectorData(nxElemZone, "NUMBER_XELEM", nZone, 1, 2);

  // Read number of elements in y-direction in each zone.
  AddVectorOption(paramFile, "NUMBER_YELEM", nyElemZone, true);
  // Pad entries for nyElemZone.
  PadEntriesVectorData(nyElemZone, "NUMBER_YELEM", nZone, 1, 2);

  // Read type of riemann solver.
  AddVectorOption(paramFile, "RIEMANN_SOLVER", NameRiemannSolver, DefaultParam.NameRiemannSolver, true);
  // Pad entries for NameRiemannSolver.
  PadEntriesVectorData(NameRiemannSolver, "RIEMANN_SOLVER", nZone, 1, 2);

	// Read zone conformity, if specified.
  AddBoolOption(paramFile, "ZONE_CONFORMITY", ZoneConformity, DefaultParam.NameZoneConformity, true);
  // Read element ratio sizes.
  AddVectorOption(paramFile, "ELEMENT_RATIO", hElemRatioZone, DefaultParam.hElemRatioZone, true);

  // Read whether a uniform grid is used or not.
  AddBoolOption(paramFile, "UNIFORM_GRID_RESOLUTION", UniformGridResolution, DefaultParam.NameUniformGridResolution, true);

  // Read expansion ratios if a non-uniform grid is specified.
  if( !UniformGridResolution ){
    // Read block interface location.
    AddVectorOption(paramFile, "BLOCK_INTERFACE_LOCATION", BlockInterfaceLocation, true);
    // Read expansion ratio per blocks.
    AddVectorOption(paramFile, "BLOCK_EXPANSION_RATIO", DomainExpansionRatio, DefaultParam.DomainExpansionRatio, true);
    // Read number of element in x-direction per block.
    AddVectorOption(paramFile, "NUMBER_XELEM_BLOCK", nxBlockElem, true);
    // Read number of element in y-direction per block.
    AddVectorOption(paramFile, "NUMBER_YELEM_BLOCK", nyBlockElem, true);
  }


  // Assign TypeDOFs map.
  MapTypeDOFs();
  // Assign ZoneMarker.
  MapTypeZone();
  // Assign RiemannSolver.
  MapRiemannSolver();

  // Determine multizone strategy to use.
  DetermineMultizoneStrategy();
  // Check element ratio specified according to multizone strategy.
  CheckElementRatio();

  // Process zone conformity option so the solver overwrites input number of
  // elements in the different zones so the grid remains consistent with ZONE_MAIN.
  if( ZoneConformity ) ProcessZoneConformity();

  // Assign solver type. For now, use a EE only.
  TypeSolver.resize(nZone);
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    TypeSolver[iZone] = SOLVER_EE;

  // Close file.
  paramFile.close();

  // Return happily.
	return true;
}


bool CConfig::ReadFlowOptions
(
 const char *configFile
)
 /*
	* Function that reads the flow conditions.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read the free-stream Mach number, if specified.
  AddScalarOption(paramFile, "FREESTREAM_MACH", MachInf, DefaultParam.MachInf, true);
  // Read the flow angle, if specified.
  AddScalarOption(paramFile, "FLOW_ANGLE", FlowAngle, DefaultParam.FlowAngle, true);

  // Check whether a cross-flow is present of not.
  CrossFlow = ( fabs(sin(FlowAngle*PI_CONSTANT/180.0)) < 1e-10 ) ? false : true;

  // Close file.
  paramFile.close();

	// Return happily.
	return true;
}


bool CConfig::ReadICOptions
(
 const char *configFile
)
 /*
	* Function that reads the initial conditions.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

	// Read initial condition in regions.
	AddVectorOption(paramFile, "TYPE_IC", NameInitialCondition, true);
  // Pad entries for NameInitialCondition.
  PadEntriesVectorData(NameInitialCondition, "TYPE_IC", nZone, 1);
	// Assign TypeIC map.
	MapTypeIC();

  // Read the disturbance center, if specified.
  AddVectorOption(paramFile, "DISTURBANCE_CENTER", CenterX0, DefaultParam.CenterX0, true);

  // Read the disturbance peak percentage, if specified.
  AddScalarOption(paramFile, "DISTURBANCE_RATIO", DisturbanceRatio, DefaultParam.DisturbanceRatio, true);
  // Read the disturbance width, if specified.
  AddScalarOption(paramFile, "DISTURBANCE_WIDTH", DisturbanceWidth, DefaultParam.DisturbanceWidth, true);
  // Read angular frequency, if specified.
  AddScalarOption(paramFile, "FREQUENCY", Frequency, DefaultParam.Frequency, true);
	// Read frequency modification, if specified.
	AddBoolOption(paramFile, "CONSTANT_FREQUENCY", ConstantFrequency, DefaultParam.NameConstantFrequency, true);
	// Read the periodic source frequency parameters, if specified.
	AddVectorOption(paramFile, "PERIODIC_SOURCE_PARAM", SourceFrequencyParam, DefaultParam.SourceFrequencyParam, true);

	// Read periodic pulse, if specified.
  AddBoolOption(paramFile, "PERIODIC_PULSE", PeriodicPulse, DefaultParam.NamePeriodicPulse, true);
	// Read whether or not the periodic pulse is aligned with the flow.
	AddBoolOption(paramFile, "ALIGNED_PERIODIC_PULSE", AlignedPeriodicPulse, DefaultParam.NameAlignedPeriodicPulse, true);
	// Read the source-term frequency exponent, if specified.
	AddScalarOption(paramFile, "SOURCE_FREQUENCY_EXPONENT", SourceFrequencyExponent, DefaultParam.SourceFrequencyExponent, true);
	// Read whether or not the periodic pulse center varies in space and time.
	AddBoolOption(paramFile, "SOURCE_TERM_CENTER_FIXED", SourceTermCenterFixed, DefaultParam.NameSourceTermCenterFixed, true);
	// Read the periodic source center shift, if specified.
	AddVectorOption(paramFile, "PERIODIC_SOURCE_CENTER_SHIFT", SourceTermCenterShift, DefaultParam.SourceTermCenterShift, true);

  // Deduce the the angular frequency from the frequency.
  AngularFrequency = 2.0*PI_CONSTANT*Frequency;

	// Check source term frequency parameters specified has only 2 entries.
	if( SourceFrequencyParam.size() != 2 )
	  Terminate("CConfig::ReadICOptions", __FILE__, __LINE__,
	            "SourceFrequencyParam must have 2 inputs.");

	// Check source term center shift, must have only 2 entries.
	if( SourceTermCenterShift.size() != 2)
	  Terminate("CConfig::ReadICOptions", __FILE__, __LINE__,
	            "SourceTermCenterShift must have 2 inputs.");

	// Check disturbance center, must be a multiple of nDim.
	if( CenterX0.size()%2 != 0)
	  Terminate("CConfig::ReadICOptions", __FILE__, __LINE__,
	            "CenterX0 must be of size multiple of nDim,");


  // Close file.
  paramFile.close();

	// Return happily.
	return true;
}


bool CConfig::ReadIOOptions
(
 const char *configFile
)
 /*
	* Function that reads the output information.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read restart solution, if specified.
  AddBoolOption(paramFile, "RESTART_SOLUTION", RestartSolution, DefaultParam.NameRestartSolution, true);
  // Read output VTK file-writing frequency.
  AddScalarOption(paramFile, "WRITE_VTKDATA_FREQ", WriteVTKFreq, DefaultParam.WriteVTKFreq, true);
  // Read output restart data file-writing frequency.
  AddScalarOption(paramFile, "WRITE_RESTART_FREQ", WriteRestartFreq, DefaultParam.WriteRestartFreq, true);
  // Read output screen-monitoring frequency.
  AddScalarOption(paramFile, "OUTPUT_FREQ", OutputFreq, DefaultParam.OutputFreq, true);
	// Read output solution filename.
	AddScalarOption(paramFile, "OUTPUT_SOL_FILENAME", OutputSolFilename, true);
	// Read output solution filename.
	AddScalarOption(paramFile, "OUTPUT_VTK_FILENAME", OutputVTKFilename, true);
	// Read output VTK variables written.
	AddVectorOption(paramFile, "VTK_VARIABLE_WRITE",  NameOutputVTKVariable, DefaultParam.NameOutputVTKVariable, true);
	// Assign OutputVTKVariable.
	MapTypeOutputVTKVariable();

  // Read input solution restart filename.
	AddScalarOption(paramFile, "RESTART_FILENAME", RestartFilename, true);
  // Read output whether to write pml auxiliary data.
  AddBoolOption(paramFile, "OUTPUT_AUX_PML", WriteAuxiliaryDataPML, DefaultParam.NameWriteAuxiliaryDataPML, true);

	// Read whether to write a GNU-plot file.
	AddBoolOption(paramFile, "WRITE_GNUPLOT", WriteGNUplot, DefaultParam.NameWriteGNUplot, true);
	// Read GNU-plot filename.
	AddScalarOption(paramFile, "GNUPLOT_FILENAME", GNUplotFilename, DefaultParam.GNUplotFilename, true);

	// Read VTK output format file.
	AddScalarOption(paramFile, "VTK_FILE_FORMAT", NameFileFormatVTK, DefaultParam.NameFileFormatVTK, true);
	// Assign TypeFileFormatVTK map.
	MapTypeFileFormatVTK();
	// Read restart solution output format file.
	AddScalarOption(paramFile, "RESTART_SOLUTION_FILE_FORMAT", NameFileFormatSolution, DefaultParam.NameFileFormatSolution, true);
	// Assign TypeFileFormatSolution map.
	MapTypeFileFormatSolution();

  // Close file.
  paramFile.close();

  // Return happily.
	return true;
}


bool CConfig::ReadSolverOptions
(
 const char *configFile
)
 /*
	* Function that reads the solver specifications.
	*/
{
  // Read the integration rule information.
  ReadIntegrationRuleOptions(configFile);
  // Read the filtering information.
  ReadFilteringOptions(configFile);
  // Read the buffer layer information.
  ReadBufferLayerOptions(configFile);

	// Return happily.
	return true;
}


void CConfig::ReadIntegrationRuleOptions
(
 const char *configFile
)
 /*
	* Function that reads the integration rule specifications.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read integration factor used in over-integration.
	AddScalarOption(paramFile, "INTEGRATION_FACTOR", IntegrationFactor, DefaultParam.IntegrationFactor, true);

  // Close file.
  paramFile.close();
}


void CConfig::ReadFilteringOptions
(
 const char *configFile
)
 /*
	* Function that reads the filtering specifications.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read whether filtering is used.
	AddVectorOption(paramFile, "FILTER_SOLUTION", NameTypeFilterSolution, DefaultParam.NameTypeFilterSolution, true);
  // Pad entries for NameTypeFilterSolution.
  PadEntriesVectorData(NameTypeFilterSolution, "FILTER_SOLUTION", nZone, 1, 2);

  // Read filtering characteristic specified.
  AddVectorOption(paramFile, "FILTER_CHARACTERISTICS", FilterCharacteristics, DefaultParam.FilterCharacteristics, true);
  // Pad entries for FilterCharacteristics.
  PadEntriesVectorDefaultData(FilterCharacteristics, DefaultParam.FilterCharacteristics, "FILTER_CHARACTERISTICS");
  // Assign FilterSolution map.
  MapTypeFilterSolution();

  // Close file.
  paramFile.close();
}


void CConfig::ReadBufferLayerOptions
(
 const char *configFile
)
 /*
	* Function that reads the buffer-layer specifications.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read buffer layer type.
  AddVectorOption(paramFile, "TYPE_BUFFER_LAYER", NameTypeBufferLayer, DefaultParam.NameTypeBufferLayer, true);
  // Pad entries for NameTypeBufferLayer.
  PadEntriesVectorData(NameTypeBufferLayer, "TYPE_BUFFER_LAYER", nZone, 1, 2);
  // Assign BufferLayerType map.
  MapTypeBufferLayer();

  // If a PML buffer layer exists, set UsePML to true.
  UsePML = false; // default
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    if( TypeBufferLayer[iZone] == PML_LAYER ) UsePML = true;

  // Close file.
  paramFile.close();
}


bool CConfig::ReadProcessingOptions
(
 const char *configFile
)
 /*
	* Function that reads the processing specifications.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read processed output filename.
  AddScalarOption(paramFile, "PROCESSED_DATA_DIRECTORY", OutputProcessedDirectory, true);
	// Read type of processing data.
  AddScalarOption(paramFile, "PROCESS_DATA", NameTypeProcessData, DefaultParam.NameTypeProcessData, true);
	// Read location of processing data.
  AddScalarOption(paramFile, "PROCESS_LOCATION", NameProcessLocation, DefaultParam.NameProcessLocation, true);
	// Assign ProcessLocation.
	MapTypeProcessLocation();

	// Read processed/probed output format file.
	AddScalarOption(paramFile, "PROCESS_PROBE_FILE_FORMAT", NameFileFormatProcessedProbed, DefaultParam.NameFileFormatProcessedProbed, true);
	// Assign TypeFileFormatProcessedProbed map.
	MapTypeFileFormatProcessedProbed();

	// Read whether a probe is specified.
	AddBoolOption(paramFile, "PROBE_SPECIFIED", ProbeSpecified, DefaultParam.NameProbeSpecified, true);
	// Read probe locations.
  AddVectorOption(paramFile, "PROBE_LOCATION", ProbeLocation, DefaultParam.ProbeLocation, true);
  // Assign ProcessData.
  MapTypeProcessData();
	// Read probe variables.
  AddVectorOption(paramFile, "PROBE_VARIABLE", NameProbeVariable, DefaultParam.NameProbeVariable, true);
  // Assign ProbeVariable.
  MapTypeProbeVariable();

	// Read whether or not surface boundart data is sampled.
	AddBoolOption(paramFile, "SAMPLE_SURFACE_DATA", SampleSurfaceData, DefaultParam.NameSampleSurfaceData, true);
	// Read the surface boundaries that are sampled.
	AddVectorOption(paramFile, "MARKER_SAMPLE_BOUNDARY", NameMarkerSampleSurface, DefaultParam.NameMarkerSampleSurface, true);
	// Read sample surface data output filename.
	AddScalarOption(paramFile, "SAMPLE_SURFACE_DIRECTORY", OutputSampleSurfaceDirectory, DefaultParam.OutputSampleSurfaceDirectory, true);
	// Assign SampleDataBoundaryID.
	MapSampleDataBoundaryID();

  // Read whether or not zone data is written.
  AddBoolOption(paramFile, "SAMPLE_ZONE_DATA", SampleZoneData, DefaultParam.NameSampleZoneData, true);
  // Read zone data output filename.
  AddScalarOption(paramFile, "ZONE_DATA_FILENAME", OutputZoneDataFilename, DefaultParam.OutputZoneDataFilename, true);
  // Read zone data selected for output.
  AddScalarOption(paramFile, "MARKER_ZONE_DATA", NameMarkerZoneData, DefaultParam.NameMarkerZoneData, true);
  // Read output file-writing frequency for zone data.
  AddScalarOption(paramFile, "ZONE_DATA_WRITE_FREQ", WriteFreqZoneData, DefaultParam.WriteFreqZoneData, true);
  // Assign ZoneData.
  MapTypeZoneData();

	// Read zone sampling output format file.
	AddScalarOption(paramFile, "ZONE_FILE_FORMAT", NameFileFormatZone, DefaultParam.NameFileFormatZone, true);
	// Assign TypeFileFormatZone map.
	MapTypeFileFormatZone();


	// Pre-process sample boundary surface data subdirectories.
	if( SampleSurfaceData ){
	
		// Total size of the sampled boundaries.
		unsigned short nSampleData = SampleDataBoundaryID.size(); 

		for(int iSample=0; iSample<nSampleData; iSample++){

			// Create a lower-case version of the current boundary marker name.
			CreateLowerCase( NameMarkerSampleSurface[iSample] );
			std::string fn = NameMarkerSampleSurface[iSample];

			// Extract base path and file directory.
			fn = OutputSampleSurfaceDirectory + fn;
			std::cout << "\n";

#if __cplusplus > 201103L
		// This standard is higher than C++11.

		// Check if the current subdirectory exists, otherwise create one.
		if( !std::filesystem::is_directory( fn ) || !std::filesystem::exists( fn ) ){
			std::cout << "Directory: " << fn << ", does not exist, creating one." << std::endl;
			std::filesystem::create_directory( fn );	
		}

#else
		// The standard is equal or below C++11.
	
		// Check if the currect subdirectory exists, otherwise create one.
		struct stat info; bool DirExists = false;
		if( stat( fn.c_str(), &info ) != 0 ) 
			std::cout << "Cannot access: " << fn << std::endl;
		else if( info.st_mode & S_IFDIR )  
			DirExists = true;
		else
			std::cout << "Directory: " << fn << ", does not exist, creating one." << std::endl;
		
		// Create the subdirectory.
		if( !DirExists ){
			const int dir_err = mkdir(fn.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if( dir_err == -1 )
				std::cout << "Could not create the directory: " << fn << std::endl;
		}

#endif
		}
	}


  // Close file.
  paramFile.close();

	// Return happily.
	return true;
}


bool CConfig::ReadBoundaryOptions
(
 const char *configFile
)
 /*
	* Function that reads the boundary information.
	*/
{
  // Read boundary conditions specified.
  ReadBoundaryConditionOptions(configFile);
  // Read characteristic-boundary information.
  ReadCharacteristicBoundaryOptions(configFile);
  // Read sponge-damping information.
  ReadSpongeDampingOptions(configFile);
  // Read grid-stretching information.
  ReadGridStretchingOptions(configFile);
  // Read artificial-convection information.
  ReadArtificialConvectionOptions(configFile);
  // Read characteristic-layer information.
	ReadCharacteristicLayerOptions(configFile);
	// Process boundary conditions, based on the multizone strategy used.
  ProcessBoundaryConditions();

  // Return happily.
	return true;
}


void CConfig::ReadBoundaryConditionOptions
(
  const char *configFile
)
 /*
  * Function that reads the boundary condition information.
  */
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read requested boundary marker, if present.
  AddVectorOption(paramFile, "MARKER_BOUNDARY", NameBoundaryCondition, true);
  // Pad entries for NameBoundaryCondition.
  PadEntriesVectorData(NameBoundaryCondition, "MARKER_BOUNDARY", nFace);
  // Assign TypeBC map. Note, indices: (SOUTH, NORTH, WEST, EAST).
  MapTypeExternalBC();

  // Close file.
  paramFile.close();
}


void CConfig::ReadCharacteristicLayerOptions
(
 const char *configFile
)
 /*
	* Function that reads the characteristic layer information.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read the sponge-layer matching constant.
  AddVectorOption(paramFile, "CHARACTERISTIC_MATCHING_CONSTANT", CharacteristicConstant, DefaultParam.CharacteristicConstant, true);
  // Pad entries for CharacteristicConstant.
  PadEntriesVectorData(CharacteristicConstant, "CHARACTERISTIC_MATCHING_CONSTANT", nZone, 1, 2);

  // Read the characteristic-layer matching exponential.
  AddVectorOption(paramFile, "CHARACTERISTIC_MATCHING_EXPONENT", CharacteristicExponent, DefaultParam.CharacteristicExponent, true);
  // Pad entries for CharacteristicExponent.
  PadEntriesVectorData(CharacteristicExponent, "CHARACTERISTIC_MATCHING_EXPONENT", nZone, 1, 2);

  // Characteristic-matching indicator. Always off in the main zone.
  CharacteristicMatching.resize(nZone, false);
  for(unsigned short iZone=1; iZone<nZone; iZone++)
    if( fabs(CharacteristicConstant[iZone]) > 1.0e-5 )
      CharacteristicMatching[iZone] = true;

  // Close file.
  paramFile.close();
}


void CConfig::ReadCharacteristicBoundaryOptions
(
  const char *configFile
)
 /*
  * Function that reads the NSCBC information.
  */
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Check if an inlet NSCBC is prescribed.
  bool InletNSCBC = false;
  for(auto& BC : TypeExternalBC) 
		if ( (BC == BC_CBC_STATIC_INLET) || (BC == BC_CBC_TOTAL_INLET) ) InletNSCBC = true;

  // if an inlet NSCBC is found, extract the tuning parameters.
  if( InletNSCBC ){
    // Read the tuning coefficients for an inlet NSCBC.
    AddVectorOption(paramFile, "INLET_NSCBC_TUNING", ParamInletNSCBC, DefaultParam.ParamInletNSCBC, true);
    // Pad entries for ParamInletNSCBC.
    PadEntriesVectorDefaultData(ParamInletNSCBC, DefaultParam.ParamInletNSCBC, "INLET_NSCBC_TUNING");
		// Read whether or not to use a non-reflecting boundary inlet.
		AddBoolOption(paramFile, "INLET_NONREFLECTING", InletNRBC, DefaultParam.NameInletNRBC, true);
  }

  // Check if an outlet NSCBC is prescribed.
  bool OutletNSCBC = false;
  for(auto& BC : TypeExternalBC) if (BC == BC_CBC_OUTLET) OutletNSCBC = true;

  // If an outlet NSCBC is found, extract the tuning parameters.
  if( OutletNSCBC ){
    // Read the tuning coefficients for an outlet NSCBC.
    AddVectorOption(paramFile, "OUTLET_NSCBC_TUNING", ParamOutletNSCBC, DefaultParam.ParamOutletNSCBC, true);
    // Pad entries for ParamOutletNSCBC.
    PadEntriesVectorDefaultData(ParamOutletNSCBC, DefaultParam.ParamOutletNSCBC, "OUTLET_NSCBC_TUNING");
  }

  // Close file.
  paramFile.close();
}


void CConfig::ReadSpongeDampingOptions
(
  const char *configFile
)
 /*
  * Function that reads the sponge-damping information.
  */
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read the sponge-layer damping constant.
  AddVectorOption(paramFile, "SPONGE_DAMPING_CONSTANT", DampingConstant, DefaultParam.DampingConstant, true);
  // Pad entries for DampingConstant.
  PadEntriesVectorData(DampingConstant, "SPONGE_DAMPING_CONSTANT", nZone, 1, 2);

  // Read the sponge-damping exponential.
  AddVectorOption(paramFile, "SPONGE_DAMPING_EXPONENT", DampingExponent, DefaultParam.DampingExponent, true);
  // Pad entries for DampingExponent.
  PadEntriesVectorData(DampingExponent, "SPONGE_DAMPING_EXPONENT", nZone, 1, 2);

  // Close file.
  paramFile.close();
}


void CConfig::ReadGridStretchingOptions
(
  const char *configFile
)
 /*
  * Function that reads the grid-stretching information.
  */
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read the grid-stretching constant.
  AddVectorOption(paramFile, "GRID_STRETCHING_CONSTANT", GridStretchingConstant, DefaultParam.GridStretchingConstant, true);
  // Pad entries for GridStretchingConstant.
  PadEntriesVectorData(GridStretchingConstant, "GRID_STRETCHING_CONSTANT", nZone, 1, 2);

  // Read the grid-stretching exponential.
  AddVectorOption(paramFile, "GRID_STRETCHING_EXPONENT", GridStretchingExponent, DefaultParam.GridStretchingExponent, true);
  // Pad entries for GridStretchingExponent.
  PadEntriesVectorData(GridStretchingExponent, "GRID_STRETCHING_EXPONENT", nZone, 1, 2);

  // Grid-stretching indicator. Always off in the main zone.
  GridStretching.resize(nZone, false);
  for(unsigned short iZone=1; iZone<nZone; iZone++)
    if( fabs(GridStretchingConstant[iZone]) > 1.0e-5 )
      GridStretching[iZone] = true;

  // Close file.
  paramFile.close();
}


void CConfig::ReadArtificialConvectionOptions
(
  const char *configFile
)
 /*
  * Function that reads the artificial-convection information.
  */
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read the artificial-convection constant.
  AddVectorOption(paramFile, "ARTIFICIAL_CONVECTION_CONSTANT", ArtificialConvectionConstant, DefaultParam.ArtificialConvectionConstant);
  // Pad entries for ArtificialConvectionConstant.
  PadEntriesVectorData(ArtificialConvectionConstant, "ARTIFICIAL_CONVECTION_CONSTANT", nZone, 1, 2);

  // Read the artificial-convection exponential.
  AddVectorOption(paramFile, "ARTIFICIAL_CONVECTION_EXPONENT", ArtificialConvectionExponent, DefaultParam.ArtificialConvectionExponent);
  // Pad entries for ArtificialConvectionExponent.
  PadEntriesVectorData(ArtificialConvectionExponent, "ARTIFICIAL_CONVECTION_EXPONENT", nZone, 1, 2);

  // Artificial-convection indicator. Always off in the main zone.
  ArtificialConvection.resize(nZone, false);
  for(unsigned short iZone=1; iZone<nZone; iZone++)
    if( fabs(ArtificialConvectionConstant[iZone]) > 1.0e-5 )
      ArtificialConvection[iZone] = true;

  // Close file.
  paramFile.close();
}


void CConfig::MapTypeFileFormatProcessedProbed
(
 void
)
 /*
	* Function that maps NameFileFormatProcessedProbed to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["BINARY"] = BINARY_FORMAT;
	Mapper["ASCII"]  = ASCII_FORMAT;

  // Initialize to unknown.
	TypeFileFormatProcessedProbed = UNKNOWN_FORMAT;

  // Check if data abides by map convention.
	try {
		Mapper.at(NameFileFormatProcessedProbed);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapTypeFileFormatProcessedProbed", __FILE__, __LINE__,
							"Processed/probed file type does not follow associated map convention!");
	}
  // Assign data according to dedicated enum.
	TypeFileFormatProcessedProbed = Mapper.at(NameFileFormatProcessedProbed);
}


void CConfig::MapTypeFileFormatZone
(
 void
)
 /*
	* Function that maps NameFileFormatZone to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["BINARY"] = BINARY_FORMAT;
	Mapper["ASCII"]  = ASCII_FORMAT;

  // Initialize to unknown.
	TypeFileFormatZone = UNKNOWN_FORMAT;

  // Check if data abides by map convention.
	try {
		Mapper.at(NameFileFormatZone);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapTypeFileFormatZone", __FILE__, __LINE__,
							"Zone sampling file type does not follow associated map convention!");
	}
  // Assign data according to dedicated enum.
	TypeFileFormatZone = Mapper.at(NameFileFormatZone);
}


void CConfig::MapTypeFileFormatSolution
(
 void
)
 /*
	* Function that maps NameFileFormatSolution to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["BINARY"] = BINARY_FORMAT;
	Mapper["ASCII"]  = ASCII_FORMAT;

  // Initialize to unknown.
	TypeFileFormatSolution = UNKNOWN_FORMAT;

  // Check if data abides by map convention.
	try {
		Mapper.at(NameFileFormatSolution);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapTypeFileFormatSolution", __FILE__, __LINE__,
							"Restart solution file type does not follow associated map convention!");
	}
  // Assign data according to dedicated enum.
	TypeFileFormatSolution = Mapper.at(NameFileFormatSolution);
}


void CConfig::MapTypeFileFormatVTK
(
 void
)
 /*
	* Function that maps NameFileFormatVTK to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["BINARY"] = BINARY_FORMAT;
	Mapper["ASCII"]  = ASCII_FORMAT;

  // Initialize to unknown.
	TypeFileFormatVTK = UNKNOWN_FORMAT;

  // Check if data abides by map convention.
	try {
		Mapper.at(NameFileFormatVTK);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapTypeFileFormatVTK", __FILE__, __LINE__,
							"VTK file type does not follow associated map convention!");
	}
  // Assign data according to dedicated enum.
	TypeFileFormatVTK = Mapper.at(NameFileFormatVTK);
}


void CConfig::MapTypeOutputVTKVariable
(
 void
)
 /*
	* Function that maps NameOutputVTKVariable to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["PRESSURE"]    = VTK_VARIABLE_PRESSURE;
  Mapper["VELOCITY"]    = VTK_VARIABLE_VELOCITY;
  Mapper["VORTICITY"]   = VTK_VARIABLE_VORTICITY;
  Mapper["MACH"]        = VTK_VARIABLE_MACH;
	Mapper["TEMPERATURE"] = VTK_VARIABLE_TEMPERATURE;
	Mapper["ENTROPY"]     = VTK_VARIABLE_ENTROPY;

	// Always include the conservative variables.
	OutputVTKVariable.push_back( VTK_VARIABLE_DENSITY     );
	OutputVTKVariable.push_back( VTK_VARIABLE_MOMENTUM    );
	OutputVTKVariable.push_back( VTK_VARIABLE_TOTALENERGY );

  for(int iVar=0; iVar<NameOutputVTKVariable.size(); iVar++){

    // Check if data abides by map convention.
		try {
			Mapper.at(NameOutputVTKVariable[iVar]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapTypeOutputVTKVariable", __FILE__, __LINE__,
								"VTK variable does not follow associated map convention!");
		}
    // Assign data according to dedicated enum.
		OutputVTKVariable.push_back( Mapper.at(NameOutputVTKVariable[iVar]) );
  }
}


void CConfig::MapTypeExternalBC
(
 void
)
 /*
	* Function that maps NameBoundaryCondition to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["PERIODIC"]          = BC_INTERFACE;
  Mapper["CBC_OUTLET"]        = BC_CBC_OUTLET;
  Mapper["CBC_STATIC_INLET"]  = BC_CBC_STATIC_INLET;
	Mapper["CBC_TOTAL_INLET"]   = BC_CBC_TOTAL_INLET;
  Mapper["SYMMETRY"]          = BC_SYMMETRY;
  Mapper["STATIC_INLET"]      = BC_STATIC_INLET;
  Mapper["TOTAL_INLET"]       = BC_TOTAL_INLET;
  Mapper["STATIC_OUTLET"]     = BC_STATIC_OUTLET;
  Mapper["SUPERSONIC_OUTLET"] = BC_SUPERSONIC_OUTLET;
  Mapper["SUPERSONIC_INLET"]  = BC_SUPERSONIC_INLET;

  // Initialize actual mapped data.
  TypeExternalBC.resize(nFace);

  for(int iFace=0; iFace<nFace; iFace++){

    // Initialize to unknown.
		TypeExternalBC[iFace] = BC_UNKNOWN;

    // Check if data abides by map convention.
		try {
			Mapper.at(NameBoundaryCondition[iFace]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapTypeExternalBC", __FILE__, __LINE__,
								"BC data does not follow associated map convention!");
		}
    // Assign data according to dedicated enum.
		TypeExternalBC[iFace] = Mapper.at(NameBoundaryCondition[iFace]);
  }
}


void CConfig::MapTypeProbeVariable
(
 void
)
 /*
	* Function that maps NameProbeVariable to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["DENSITY"]     = PROBE_DENSITY;
	Mapper["XMOMENTUM"]   = PROBE_XMOMENTUM;
	Mapper["YMOMENTUM"]   = PROBE_YMOMENTUM;
	Mapper["TOTALENERGY"] = PROBE_TOTALENERGY;
	Mapper["XVELOCITY"]   = PROBE_XVELOCITY;
	Mapper["YVELOCITY"]   = PROBE_YVELOCITY;
	Mapper["PRESSURE"]    = PROBE_PRESSURE;

	// Extract number of input variables to probe.
	unsigned short nProbeVar = NameProbeVariable.size();

	// Initialize vector of probe variables.
	ProbeVariable.resize(nProbeVar, PROBE_NOTHING);

	for(int iVar=0; iVar<nProbeVar; iVar++){

		// Check if data abides by map convention.
	  try {
			Mapper.at(NameProbeVariable[iVar]);
		}
	  catch ( std::out_of_range& ) {
		  Terminate("CConfig::MapTypeProbeVariable", __FILE__, __LINE__,
		 		        "Probe variable does not follow associated map convention!");
		}
	  // Assign data according to dedicated enum.
	  ProbeVariable[iVar] = Mapper.at(NameProbeVariable[iVar]);
	}
}


void CConfig::MapSampleDataBoundaryID
(
 void
)
 /*
	* Function that maps NameMarkerSampleSurface to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["XMIN"] = IDX_WEST;
	Mapper["XMAX"] = IDX_EAST;
	Mapper["YMIN"] = IDX_SOUTH;
	Mapper["YMAX"] = IDX_NORTH;

	// Extract number of input sample data boundaries.
	unsigned short nSampleData = NameMarkerSampleSurface.size();

	// Initialize vector of boundary IDs.
	SampleDataBoundaryID.resize(nSampleData, NONE);

	for(int iSample=0; iSample<nSampleData; iSample++){

		// Check if data abides by map convention.
	  try {
			Mapper.at(NameMarkerSampleSurface[iSample]);
		}
	  catch ( std::out_of_range& ) {
		  Terminate("CConfig::MapSampleDataBoundaryID", __FILE__, __LINE__,
		 		        "Sample data boundary does not follow associated map convention!");
		}
	  // Assign data according to dedicated enum.
	  SampleDataBoundaryID[iSample] = Mapper.at(NameMarkerSampleSurface[iSample]);
	}
}


void CConfig::MapTypeProcessLocation
(
 void
)
 /*
	* Function that maps NameProcessLocation to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["XMIN"]   = PROCESS_LOCATION_XMIN;
	Mapper["XMAX"]   = PROCESS_LOCATION_XMAX;
	Mapper["YMIN"]   = PROCESS_LOCATION_YMIN;
	Mapper["YMAX"]   = PROCESS_LOCATION_YMAX;
	Mapper["DOMAIN"] = PROCESS_LOCATION_DOMAIN; 

	// Initialize process location.
	ProcessLocation = PROCESS_LOCATION_UNKNOWN;

	// Check if data abides by map convention.
	try {
		Mapper.at(NameProcessLocation);
	}
	catch ( std::out_of_range& ) {
	  Terminate("CConfig::MapTypeProcessLocation", __FILE__, __LINE__,
	 		        "Process location does not follow associated map convention!");
	}
	// Assign data according to dedicated enum.
	ProcessLocation = Mapper.at(NameProcessLocation);
}


void CConfig::MapTypeZoneData
(
 void
)
 /*
	* Function that maps NameMarkerZoneData to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
	Mapper["ZONE_MAIN"]     = ZONE_MAIN;
	Mapper["ZONE_WEST"]     = ZONE_WEST;
  Mapper["ZONE_EAST"]     = ZONE_EAST;
  Mapper["ZONE_SOUTH"]    = ZONE_SOUTH;
  Mapper["ZONE_NORTH"]    = ZONE_NORTH;
  Mapper["ZONE_CORNER_0"] = ZONE_CORNER_0;
  Mapper["ZONE_CORNER_1"] = ZONE_CORNER_1;
  Mapper["ZONE_CORNER_2"] = ZONE_CORNER_2;
  Mapper["ZONE_CORNER_3"] = ZONE_CORNER_3;

  // Initialize to unknown.
	TypeZoneData = ZONE_UNKNOWN;

  // Check if data abides by map convention.
	try {
		Mapper.at(NameMarkerZoneData);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapTypeZoneData", __FILE__, __LINE__,
							"Zone data type does not follow associated map convention!");
	}
  // Assign data according to dedicated enum.
	TypeZoneData = Mapper.at(NameMarkerZoneData);
}


void CConfig::MapTypeProcessData
(
 void
)
 /*
	* Function that maps NameTypeProcessData to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["RATIO_P_U"]              = PROCESS_RATIO_P_U;
	Mapper["RATIO_P_M"]              = PROCESS_RATIO_P_M;
	Mapper["RATIO_WM_WP"]            = PROCESS_RATIO_WM_WP; 
	Mapper["RATIO_LM_LP"]            = PROCESS_RATIO_LM_LP;
	Mapper["WAVE_AMPLITUDE_ENTROPY"] = PROCESS_WAVE_ENTROPY;
	Mapper["NOTHING"]                = PROCESS_NOTHING;

  // Initialize to unknown.
	TypeProcessData = PROCESS_NOTHING;

  // Check if data abides by map convention.
	try {
		Mapper.at(NameTypeProcessData);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapProcessData", __FILE__, __LINE__,
							"Processing data type does not follow associated map convention!");
	}
  // Assign data according to dedicated enum.
	TypeProcessData = Mapper.at(NameTypeProcessData);
}


void CConfig::MapTypeIC
(
 void
)
 /*
	* Function that maps NameInitialCondition to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["GAUSSIAN_PRESSURE"]      = IC_GAUSSIAN_PRESSURE;
  Mapper["ISENTROPIC_VORTEX"]      = IC_ISENTROPIC_VORTEX;
  Mapper["ENTROPY_WAVE"]           = IC_ENTROPY_WAVE;
  Mapper["VORTEX_ROLLUP"]          = IC_VORTEX_ROLLUP;
  Mapper["ACOUSTIC_PLANE_WAVE"]    = IC_ACOUSTIC_PLANE_WAVE;
	Mapper["GAUSSIAN_PRESSURE_1D_X"] = IC_GAUSSIAN_PRESSURE_1D_X;
	Mapper["GAUSSIAN_PRESSURE_1D_Y"] = IC_GAUSSIAN_PRESSURE_1D_Y;
	
	// Initialize actual mapped data.
	TypeIC.resize(nZone);

	// Iterate of data.
	for(int iZone=0; iZone<nZone; iZone++){

		// Initialize to unknown.
		TypeIC[iZone] = IC_UNKNOWN;

		// Check if data abides by map convention.
		try {
			Mapper.at(NameInitialCondition[iZone]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapTypeIC", __FILE__, __LINE__,
								"IC data does not follow associated map convention!");
		}
		// Assign data according to dedicated enum.
		TypeIC[iZone] = Mapper.at(NameInitialCondition[iZone]);
	}
}


void CConfig::MapTypeBufferLayer
(
 void
)
 /*
	* Function that maps NameTypeBufferLayer to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["NONE"]   = NO_LAYER;
  Mapper["PML"]    = PML_LAYER;
  Mapper["SPONGE"] = SPONGE_LAYER;

  // Initialize to no layer..
	TypeBufferLayer.resize(nZone, NO_LAYER);

  // Iterate of data. Note, main zone is skipped since it is always a physical domain.
	for(int iZone=1; iZone<nZone; iZone++){

    // Initialize to no layer.
    TypeBufferLayer[iZone] = NO_LAYER;

    // Check if data abides by map convention.
  	try {
  		Mapper.at(NameTypeBufferLayer[iZone]);
  	}
  	catch ( std::out_of_range& ) {
  		Terminate("CConfig::MapTypeBufferLayer", __FILE__, __LINE__,
  							"Buffer layer type does not follow associated map convention!");
  	}
    // Assign data according to dedicated enum.
  	TypeBufferLayer[iZone] = Mapper.at(NameTypeBufferLayer[iZone]);
  }
}


void CConfig::MapTemporalScheme
(
 void
)
 /*
	* Function that maps NameTemporalScheme to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["LSRK4"]  = TEMPORAL_SCHEME_LSRK4;
	Mapper["CRK4"]   = TEMPORAL_SCHEME_CRK4;
  Mapper["SSPRK3"] = TEMPORAL_SCHEME_SSPRK3;

	// Initialize to unknown.
	TypeTemporalScheme = TEMPORAL_SCHEME_UNKNOWN;

	// Check if data abides by map convention.
	try {
		Mapper.at(NameTemporalScheme);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapTemporalScheme", __FILE__, __LINE__,
							"Temporal data does not follow associated map convention!");
	}

	// Assign data according to dedicated enum.
	TypeTemporalScheme = Mapper.at(NameTemporalScheme);
}


void CConfig::MapTypeDOFs
(
 void
)
 /*
	* Function that maps NodalDOFs to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["EQD"] = TYPE_DOF_EQD;
	Mapper["LGL"] = TYPE_DOF_LGL;

	// Initialize actual mapped data.
	TypeDOFs.resize(nZone);

	// Iterate of data.
	for(int iZone=0; iZone<nZone; iZone++){

		// Initialize to unknown.
		TypeDOFs[iZone] = TYPE_DOF_UNKOWN;

		// Check if data abides by map convention.
		try {
			Mapper.at(NameNodalDOFs[iZone]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapTypeDOFs", __FILE__, __LINE__,
								"DOF data does not follow associated map convention!");
		}
		// Assign data according to dedicated enum.
		TypeDOFs[iZone] = Mapper.at(NameNodalDOFs[iZone]);
	}
}


void CConfig::MapTypeZone
(
 void
)
 /*
	* Function that maps TypeZone to its enum type.
	*/
{
	// Initialize actual mapped data.
	TypeZone.resize(nZone);

  // Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["ZONE_MAIN"]     = ZONE_MAIN;
	Mapper["ZONE_WEST"]     = ZONE_WEST;
  Mapper["ZONE_EAST"]     = ZONE_EAST;
  Mapper["ZONE_SOUTH"]    = ZONE_SOUTH;
  Mapper["ZONE_NORTH"]    = ZONE_NORTH;
  Mapper["ZONE_CORNER_0"] = ZONE_CORNER_0;
  Mapper["ZONE_CORNER_1"] = ZONE_CORNER_1;
  Mapper["ZONE_CORNER_2"] = ZONE_CORNER_2;
  Mapper["ZONE_CORNER_3"] = ZONE_CORNER_3;

  // Base zone name.
  std::string basename = "ZONE_";

	// Iterate of data.
	for(int iZone=0; iZone<nZone; iZone++){

    // Initialize to unknown.
		TypeZone[iZone] = ZONE_UNKNOWN;

		// Check if data abides by map convention.
		try {
			Mapper.at(NameZoneMarker[iZone]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapTypeZone", __FILE__, __LINE__,
								"Zone names do not follow associated map convention!");
		}
		// Assign data according to dedicated enum.
		TypeZone[iZone] = Mapper.at(NameZoneMarker[iZone]);
	}

	// Make sure input zone markers are unique.
	for(unsigned short jZone=0; jZone<nZone; jZone++)
		for(unsigned short iZone=jZone+1; iZone<nZone; iZone++)
			if( TypeZone[jZone] == TypeZone[iZone] )
				Terminate("CConfig::MapTypeZone", __FILE__, __LINE__,
									"Zone markers must be unique!");
}


void CConfig::MapRiemannSolver
(
  void
)
 /*
  * Function that maps the RiemannSolver to its enum type.
  */
{
  // Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["ROE"]       = RIEMANN_ROE;
  Mapper["RUSANOV"]   = RIEMANN_RUSANOV;
  Mapper["ROEISMAIL"] = RIEMANN_ROEISMAIL;

	// Initialize actual mapped data.
	RiemannSolver.resize(nZone);

	// Iterate of data.
	for(int iZone=0; iZone<nZone; iZone++){

		// Initialize to unknown.
		RiemannSolver[iZone] = RIEMANN_UNKNOWN;

		// Check if data abides by map convention.
		try {
			Mapper.at(NameRiemannSolver[iZone]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapRiemannSolver", __FILE__, __LINE__,
								"Riemann data does not follow associated map convention!");
		}
		// Assign data according to dedicated enum.
		RiemannSolver[iZone] = Mapper.at(NameRiemannSolver[iZone]);
	}
}


void CConfig::MapTypeFilterSolution
(
 void
)
 /*
	* Function that maps TypeFilterSolution to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["NONE"]        = NO_FILTER;
	Mapper["EXPONENTIAL"] = EXPONENTIAL_FILTER;

	// Initialize actual mapped data.
	TypeFilterSolution.resize(nZone);

	// Iterate of data.
	for(int iZone=0; iZone<nZone; iZone++){

		// Initialize to unknown.
		TypeFilterSolution[iZone] = NO_FILTER;

		// Check if data abides by map convention.
		try {
			Mapper.at(NameTypeFilterSolution[iZone]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapTypeFilterSolution", __FILE__, __LINE__,
								"Filtering data does not follow associated map convention!");
		}
		// Assign data according to dedicated enum.
		TypeFilterSolution[iZone] = Mapper.at(NameTypeFilterSolution[iZone]);
	}
}


void CConfig::MapTypeModifyBC
(
 void
)
 /*
	* Function that maps TypeModifyBC to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["NONE"]                = NO_BC_MODIFICATION;
	Mapper["GAUSSIAN_PRESSURE"]   = GAUSSIAN_PRESSURE_BC_MODIFICATION;
  Mapper["SINUSOIDAL_VELOCITY"] = SINUSOIDAL_VELOCITY_BC_MODIFICATION;

  // Initialize actual mapped data.
  TypeModifyBC.resize(nFace);

  // Iterate on data.
  for(unsigned short iFace=0; iFace<nFace; iFace++){

  	// Initialize to none.
  	TypeModifyBC[iFace] = NO_BC_MODIFICATION;

  	// Check if data abides by map convention.
  	try {
  		Mapper.at(NameTypeModifyBC[iFace]);
  	}
  	catch ( std::out_of_range& ) {
  		Terminate("CConfig::MapModifyTypeBC", __FILE__, __LINE__,
  							"TypeModifyBC data does not follow associated map convention!");
  	}

  	// Assign data according to dedicated enum.
  	TypeModifyBC[iFace] = Mapper.at(NameTypeModifyBC[iFace]);
  }
}


void CConfig::DetermineMultizoneStrategy
(
  void
)
 /*
  * Function that determines what strategy for a multizone simulation to adopt.
  */
{
  // Consistency check, input zone markers must be either 1, 2 or 9.
  if( (nZone != 1) && (nZone != 2) && (nZone != 3) && (nZone != 9) )
    Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
              "Multizone strategy demands nZone be: 1, 2, 3 or 9.");

  // To simplify things, make sure the first zone is always the main zone.
  if( TypeZone[0] != ZONE_MAIN )
    Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
              "First input zone must be ZONE_MAIN.");

  // Determine if this is strategy main.
  if( nZone == 1 ){
    // Assign multizone strategy.
    MultizoneStrategy = MULTIZONE_STRATEGY_MAIN;
    return;
  }

  // Determine if this is an full-zonal strategy.
  if( nZone == 9 ){
    // Assign multizone strategy.
    MultizoneStrategy = MULTIZONE_STRATEGY_ALL;
    // Make sure the element ratio input is correct.
    if( (hElemRatioZone[0] <= 0.0) || (hElemRatioZone.size() != 4) )
      Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
                "ELEMENT_RATIO must be positive and: (r(w), r(e), r(s), r(n).");
    return;
   }

   // Determine if this is either combination strategies. Note, the
   // first zone is always fixed as the main one so no need to check for that.
   if( nZone == 2 ){
     switch( TypeZone[1] ){
       case(ZONE_WEST):  MultizoneStrategy = MULTIZONE_STRATEGY_WEST;  break;
       case(ZONE_EAST):  MultizoneStrategy = MULTIZONE_STRATEGY_EAST;  break;
       case(ZONE_SOUTH): MultizoneStrategy = MULTIZONE_STRATEGY_SOUTH; break;
       case(ZONE_NORTH): MultizoneStrategy = MULTIZONE_STRATEGY_NORTH; break;

       default:
         Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
                   "Wrong second zone input, must be: WEST, EAST, SOUTH or NORTH.");
     }

     return;
   }

   // Determine if this is either combination strategies. Note, the first zone
   // is always fixed as the main one so no need to check for that.
   if( nZone == 3 ){

     // Make sure all zones are unique.
     if( TypeZone[1] == TypeZone[2] )
      Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
                "Second and third zone inputs must be unique.");

     // Check type of second zone.
     switch( TypeZone[1] ){

       // Check horizonal layers.
       case(ZONE_WEST): case(ZONE_EAST):
       {
         if( (TypeZone[2] == ZONE_WEST) || (TypeZone[2] == ZONE_EAST) )
           MultizoneStrategy = MULTIZONE_STRATEGY_HORIZONAL;

         break;
       }

       // Check vertical layers.
       case(ZONE_SOUTH): case(ZONE_NORTH):
       {
         if( (TypeZone[2] == ZONE_SOUTH) || (TypeZone[2] == ZONE_NORTH) )
           MultizoneStrategy = MULTIZONE_STRATEGY_VERTICAL;

         break;
       }

       default:
         Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
                   "Wrong second/third zone input combinations, must be vertical or horizontal.");
     }

     return;
   }
}


void CConfig::CheckElementRatio
(
  void
)
 /*
  * Function that checks for inconsistency errors in the element ratio
  * specified, w.r.t. the multizone strategy employed.
  */
{
  // In case this is a single-zone simulation, no need to check.
  if( MultizoneStrategy == MULTIZONE_STRATEGY_MAIN ) return;

  // Error flag.
  bool ErrorDetected = false;

  // If this is not a single-zone simulation, make sure all is consistent.
  if( hElemRatioZone.size() != nFace ) ErrorDetected = true;
  for(unsigned short i=0; i<nFace; i++)
    if( hElemRatioZone[i] <= 0.0)
      ErrorDetected = true;

  // Report error and exit, in case detected.
  if( ErrorDetected )
    Terminate("CConfig::CheckElementRatio", __FILE__, __LINE__,
              "ELEMENT_RATIO is inconsistent: must be positive and of size: 4");
}


void CConfig::ProcessZoneConformity
(
  void
)
 /*
  * Function that enforces the zone conformity if turned on.
  */
{
  // In case this is a single zone, do nothing and return.
  if( nZone == 1 ) return;

  // Determine what multizone strategy is employed.
  switch(MultizoneStrategy){

    // Combination of east or west zones.
    case(MULTIZONE_STRATEGY_EAST): case(MULTIZONE_STRATEGY_WEST):
    nyElemZone[1] = nyElemZone[0]; break;

    // Combination of south or north zones.
    case(MULTIZONE_STRATEGY_SOUTH): case(MULTIZONE_STRATEGY_NORTH):
    nxElemZone[1] = nxElemZone[0]; break;

    // All zones are utilized.
    case(MULTIZONE_STRATEGY_ALL):
    {
      // Consistency in naming check.
      assert( TypeZone[0] == ZONE_MAIN );

      for(unsigned short iZone=0; iZone<nZone; iZone++){
        switch(TypeZone[iZone]){
          case(ZONE_MAIN):  break;
          case(ZONE_WEST):  case(ZONE_EAST):  nyElemZone[iZone] = nyElemZone[ZONE_MAIN]; break;
          case(ZONE_SOUTH): case(ZONE_NORTH): nxElemZone[iZone] = nxElemZone[ZONE_MAIN]; break;
          case(ZONE_CORNER_0):
          {
            nxElemZone[iZone] = nxElemZone[ZONE_WEST];
            nyElemZone[iZone] = nyElemZone[ZONE_SOUTH];
            break;
          }
          case(ZONE_CORNER_1):
          {
            nxElemZone[iZone] = nxElemZone[ZONE_EAST];
            nyElemZone[iZone] = nyElemZone[ZONE_SOUTH];
            break;
          }
          case(ZONE_CORNER_2):
          {
            nxElemZone[iZone] = nxElemZone[ZONE_WEST];
            nyElemZone[iZone] = nyElemZone[ZONE_NORTH];
            break;
          }
          case(ZONE_CORNER_3):
          {
            nxElemZone[iZone] = nxElemZone[ZONE_EAST];
            nyElemZone[iZone] = nyElemZone[ZONE_NORTH];
            break;
          }
        }
      }

      break;
    }

    // Horizontal zones only.
    case(MULTIZONE_STRATEGY_HORIZONAL):
    {
      nyElemZone[1] = nyElemZone[0];
      nyElemZone[2] = nyElemZone[0];
      break;
    }

    // Vertical zones only.
    case(MULTIZONE_STRATEGY_VERTICAL):
    {
      nxElemZone[1] = nxElemZone[0];
      nxElemZone[2] = nxElemZone[0];
      break;
    }

  }
}


void CConfig::ProcessBoundaryConditions
(
  void
)
  /*
   * Function that processes the boundary conditions in each zone based on
   * the predetermined multizone strategy and their matching zone ID, in case
   * interface boundaries are prescribed.
   */
{
  // Reserve memory for the boundary conditions per zone.
  TypeBC.resize(nZone);
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    TypeBC[iZone].resize(nFace);

  // Reserve memory for the interface conditions per each boundary.
  InterfaceID.resize(nZone);
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    InterfaceID[iZone].resize(nFace);

  // Initialize everything as an interface BC and then modify accordingly, since
  // this is easier and more readable.
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    for(unsigned short iBoundary=0; iBoundary<nFace; iBoundary++)
      TypeBC[iZone][iBoundary] = BC_INTERFACE;

  // Determine what multizone strategy we are dealing with.
  switch( MultizoneStrategy ){

    // Single main zone.
    case(MULTIZONE_STRATEGY_MAIN):
    {
      // Make sure data is consistent.
      assert( (nZone == 1) && "Number of zones must be one!");

      // Assign boundaries.
      TypeBC[ZONE_MAIN][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[ZONE_MAIN][IDX_NORTH] = TypeExternalBC[IDX_NORTH];
      TypeBC[ZONE_MAIN][IDX_WEST]  = TypeExternalBC[IDX_WEST];
      TypeBC[ZONE_MAIN][IDX_EAST]  = TypeExternalBC[IDX_EAST];

      break;
    }

    // Combination of zones: main and west.
    case(MULTIZONE_STRATEGY_WEST):
    {
      // Make sure data is consistent.
      assert( (nZone == 2) && "Number of zones must be two!");

      // Assign boundaries: main zone.
      TypeBC[ZONE_MAIN][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[ZONE_MAIN][IDX_NORTH] = TypeExternalBC[IDX_NORTH];
      TypeBC[ZONE_MAIN][IDX_EAST]  = TypeExternalBC[IDX_EAST];

      // Assign boundaries: west zone.
      TypeBC[1][IDX_WEST]  = TypeExternalBC[IDX_WEST];
      TypeBC[1][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[1][IDX_NORTH] = TypeExternalBC[IDX_NORTH];

      break;
    }

    // Combination of zones: main and east.
    case(MULTIZONE_STRATEGY_EAST):
    {
      // Make sure data is consistent.
      assert( (nZone == 2) && "Number of zones must be two!");

      // Assign boundaries: main zone.
      TypeBC[ZONE_MAIN][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[ZONE_MAIN][IDX_NORTH] = TypeExternalBC[IDX_NORTH];
      TypeBC[ZONE_MAIN][IDX_WEST]  = TypeExternalBC[IDX_WEST];

      // Assign boundaries: east zone.
      TypeBC[1][IDX_EAST]  = TypeExternalBC[IDX_EAST];
      TypeBC[1][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[1][IDX_NORTH] = TypeExternalBC[IDX_NORTH];

      break;
    }

    // Combination of zones: main and south.
    case(MULTIZONE_STRATEGY_SOUTH):
    {
      // Make sure data is consistent.
      assert( (nZone == 2) && "Number of zones must be two!");

      // Assign boundaries: main zone.
      TypeBC[ZONE_MAIN][IDX_WEST]   = TypeExternalBC[IDX_WEST];
      TypeBC[ZONE_MAIN][IDX_NORTH]  = TypeExternalBC[IDX_NORTH];
      TypeBC[ZONE_MAIN][IDX_EAST]   = TypeExternalBC[IDX_EAST];

      // Assign boundaries: south zone.
      TypeBC[1][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[1][IDX_WEST]  = TypeExternalBC[IDX_WEST];
      TypeBC[1][IDX_EAST]  = TypeExternalBC[IDX_EAST];

      break;
    }

    // Combination of zones: main and north.
    case(MULTIZONE_STRATEGY_NORTH):
    {
      // Make sure data is consistent.
      assert( (nZone == 2) && "Number of zones must be two!");

      // Assign boundaries: main zone.
      TypeBC[ZONE_MAIN][IDX_WEST]   = TypeExternalBC[IDX_WEST];
      TypeBC[ZONE_MAIN][IDX_SOUTH]  = TypeExternalBC[IDX_SOUTH];
      TypeBC[ZONE_MAIN][IDX_EAST]   = TypeExternalBC[IDX_EAST];

      // Assign boundaries: north zone.
      TypeBC[1][IDX_NORTH] = TypeExternalBC[IDX_NORTH];
      TypeBC[1][IDX_WEST]  = TypeExternalBC[IDX_WEST];
      TypeBC[1][IDX_EAST]  = TypeExternalBC[IDX_EAST];

      break;
    }

    // All zones are used.
    case(MULTIZONE_STRATEGY_ALL):
    {
      // Make sure data is consistent.
      assert( (nZone == 9) && "Number of zones must be nine!");

      // Assign boundary: south zone.
      TypeBC[ZONE_SOUTH][IDX_SOUTH]    = TypeExternalBC[IDX_SOUTH];
      // Assign boundary: north zone.
      TypeBC[ZONE_NORTH][IDX_NORTH]    = TypeExternalBC[IDX_NORTH];
      // Assign boundary: east zone.
      TypeBC[ZONE_EAST][IDX_EAST]      = TypeExternalBC[IDX_EAST];
      // Assign boundary: west zone.
      TypeBC[ZONE_WEST][IDX_WEST]      = TypeExternalBC[IDX_WEST];

      // Assign boundaries: corner0 zone.
      TypeBC[ZONE_CORNER_0][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[ZONE_CORNER_0][IDX_WEST]  = TypeExternalBC[IDX_WEST];
      // Assign boundaries: corner1 zone.
      TypeBC[ZONE_CORNER_1][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[ZONE_CORNER_1][IDX_EAST]  = TypeExternalBC[IDX_EAST];
      // Assign boundaries: corner2 zone.
      TypeBC[ZONE_CORNER_2][IDX_WEST]  = TypeExternalBC[IDX_WEST];
      TypeBC[ZONE_CORNER_2][IDX_NORTH] = TypeExternalBC[IDX_NORTH];
      // Assign boundaries: corner3 zone.
      TypeBC[ZONE_CORNER_3][IDX_EAST]  = TypeExternalBC[IDX_EAST];
      TypeBC[ZONE_CORNER_3][IDX_NORTH] = TypeExternalBC[IDX_NORTH];

      break;
    }

    // Horizontal zones only.
    case(MULTIZONE_STRATEGY_HORIZONAL):
    {
      // Make sure data is consistent.
      assert( (nZone == 3) && "Number of zones must be three!");
      assert( (unsigned short)HORIZONTAL_ZONE_MAIN == ZONE_MAIN );

      // Assign boundaries: main zone.
      TypeBC[HORIZONTAL_ZONE_MAIN][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[HORIZONTAL_ZONE_MAIN][IDX_NORTH] = TypeExternalBC[IDX_NORTH];

      // Assign boundaries: west zone.
      TypeBC[HORIZONTAL_ZONE_WEST][IDX_WEST]  = TypeExternalBC[IDX_WEST];
      TypeBC[HORIZONTAL_ZONE_WEST][IDX_NORTH] = TypeExternalBC[IDX_NORTH];
      TypeBC[HORIZONTAL_ZONE_WEST][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];

      // Assign boundaries: east zone.
      TypeBC[HORIZONTAL_ZONE_EAST][IDX_EAST]  = TypeExternalBC[IDX_EAST];
      TypeBC[HORIZONTAL_ZONE_EAST][IDX_NORTH] = TypeExternalBC[IDX_NORTH];
      TypeBC[HORIZONTAL_ZONE_EAST][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];

      break;
    }

    // Vertical zones only.
    case(MULTIZONE_STRATEGY_VERTICAL):
    {
      // Make sure data is consistent.
      assert( (nZone == 3) && "Number of zones must be three!");
      assert( (unsigned short)VERTICAL_ZONE_MAIN == ZONE_MAIN );

      // Assign boundaries: main zone.
      TypeBC[VERTICAL_ZONE_MAIN][IDX_EAST] = TypeExternalBC[IDX_EAST];
      TypeBC[VERTICAL_ZONE_MAIN][IDX_WEST] = TypeExternalBC[IDX_WEST];

      // Assign boundaries: south zone.
      TypeBC[VERTICAL_ZONE_SOUTH][IDX_SOUTH] = TypeExternalBC[IDX_SOUTH];
      TypeBC[VERTICAL_ZONE_SOUTH][IDX_EAST]  = TypeExternalBC[IDX_EAST];
      TypeBC[VERTICAL_ZONE_SOUTH][IDX_WEST]  = TypeExternalBC[IDX_WEST];

      // Assign boundaries: north zone.
      TypeBC[VERTICAL_ZONE_NORTH][IDX_NORTH] = TypeExternalBC[IDX_NORTH];
      TypeBC[VERTICAL_ZONE_NORTH][IDX_EAST]  = TypeExternalBC[IDX_EAST];
      TypeBC[VERTICAL_ZONE_NORTH][IDX_WEST]  = TypeExternalBC[IDX_WEST];

      break;
    }


    default:
      Terminate("CConfig::ProcessBoundaryConditions", __FILE__, __LINE__,
                "Unknown/wrong multizone strategy detected.");
  }

  // Assign interface data information (jZone, jBoundary), if needed.
  for(unsigned short iZone=0; iZone<nZone; iZone++){
    for(unsigned short iBoundary=0; iBoundary<nFace; iBoundary++){
      if(TypeBC[iZone][iBoundary] == BC_INTERFACE){
        InterfaceID[iZone][iBoundary] = MatchInterface(iZone, iBoundary);
      }
    }
  }
}


as3vector1d<unsigned short> CConfig::MatchInterface
(
  unsigned short iZone,
  unsigned short iFace
)
 /*
  * Function that matches iZone iFace with its counterpart in jZone jFace.
  */
{
  // Matching data ID: (jZone, jFace).
  as3vector1d<unsigned short> jData(2);

  // Error flag.
  bool ErrorDetected = false;

  // Determine matching face.
  switch(iFace){
    case(IDX_SOUTH): jData[1] = IDX_NORTH; break;
    case(IDX_NORTH): jData[1] = IDX_SOUTH; break;
    case(IDX_WEST):  jData[1] = IDX_EAST;  break;
    case(IDX_EAST):  jData[1] = IDX_WEST;  break;

    default: ErrorDetected = true;
  }

  // Match zone.
  switch(MultizoneStrategy){

    // Main zone strategy.
    case(MULTIZONE_STRATEGY_MAIN): jData[0] = iZone; break;

    // West-buffer zone strategy.
    case(MULTIZONE_STRATEGY_WEST):
    {
      if(iZone == ZONE_MAIN){
        // This is the main zone.
        if( (iFace == IDX_SOUTH) || (iFace == IDX_NORTH) ) jData[0] = ZONE_MAIN;
        else                                               jData[0] = 1;
      }
      else {
        // This is the west zone.
        if( (iFace == IDX_SOUTH) || (iFace == IDX_NORTH) ) jData[0] = 1;
        else                                               jData[0] = ZONE_MAIN;
      }
      break;
    }

    // East-buffer zone strategy.
    case(MULTIZONE_STRATEGY_EAST):
    {
      if(iZone == ZONE_MAIN){
        // This is the main zone.
        if( (iFace == IDX_SOUTH) || (iFace == IDX_NORTH) ) jData[0] = ZONE_MAIN;
        else                                               jData[0] = 1;
      }
      else {
        // This is the east zone.
        if( (iFace == IDX_SOUTH) || (iFace == IDX_NORTH) ) jData[0] = 1;
        else                                               jData[0] = ZONE_MAIN;
      }
      break;
    }

    // South-buffer zone strategy.
    case(MULTIZONE_STRATEGY_SOUTH):
    {
      if(iZone == ZONE_MAIN){
        // This is the main zone.
        if( (iFace == IDX_EAST) || (iFace == IDX_WEST) ) jData[0] = ZONE_MAIN;
        else                                             jData[0] = 1;
      }
      else {
        // This is the south zone.
        if( (iFace == IDX_EAST) || (iFace == IDX_WEST) ) jData[0] = 1;
        else                                             jData[0] = ZONE_MAIN;
      }
      break;
    }

    // North-buffer zone strategy.
    case(MULTIZONE_STRATEGY_NORTH):
    {
      if(iZone == ZONE_MAIN){
        // This is the main zone.
        if( (iFace == IDX_EAST) || (iFace == IDX_WEST) ) jData[0] = ZONE_MAIN;
        else                                             jData[0] = 1;
      }
      else {
        // This is the south zone.
        if( (iFace == IDX_EAST) || (iFace == IDX_WEST) ) jData[0] = 1;
        else                                             jData[0] = ZONE_MAIN;
      }
      break;
    }

    // All zones buffer strategy.
    case(MULTIZONE_STRATEGY_ALL):
    {
      // Determine which zone we are pivoting from.
      switch(iZone){

        // Pivot from main zone.
        case(ZONE_MAIN):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_SOUTH; break;
            case(IDX_NORTH): jData[0] = ZONE_NORTH; break;
            case(IDX_WEST):  jData[0] = ZONE_WEST;  break;
            case(IDX_EAST):  jData[0] = ZONE_EAST;  break;
            default: ErrorDetected = true;
          }
          break;
        }

        // Pivot from west zone.
        case(ZONE_WEST):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_CORNER_0; break;
            case(IDX_NORTH): jData[0] = ZONE_CORNER_2; break;
            case(IDX_WEST):  jData[0] = ZONE_EAST;     break;
            case(IDX_EAST):  jData[0] = ZONE_MAIN;     break;
            default: ErrorDetected = true;
          }
          break;
        }

        // Pivot from east zone.
        case(ZONE_EAST):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_CORNER_1; break;
            case(IDX_NORTH): jData[0] = ZONE_CORNER_3; break;
            case(IDX_WEST):  jData[0] = ZONE_MAIN;     break;
            case(IDX_EAST):  jData[0] = ZONE_WEST;     break;
            default: ErrorDetected = true;
          }
          break;
        }

        // Pivot from south zone.
        case(ZONE_SOUTH):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_NORTH;    break;
            case(IDX_NORTH): jData[0] = ZONE_MAIN;     break;
            case(IDX_WEST):  jData[0] = ZONE_CORNER_0; break;
            case(IDX_EAST):  jData[0] = ZONE_CORNER_1; break;
            default: ErrorDetected = true;
          }
          break;
        }

        // Pivot from north zone.
        case(ZONE_NORTH):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_MAIN;     break;
            case(IDX_NORTH): jData[0] = ZONE_SOUTH;    break;
            case(IDX_WEST):  jData[0] = ZONE_CORNER_2; break;
            case(IDX_EAST):  jData[0] = ZONE_CORNER_3; break;
            default: ErrorDetected = true;
          }
          break;
        }

        // Pivot from corner0 zone.
        case(ZONE_CORNER_0):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_CORNER_2; break;
            case(IDX_NORTH): jData[0] = ZONE_WEST;     break;
            case(IDX_WEST):  jData[0] = ZONE_CORNER_1; break;
            case(IDX_EAST):  jData[0] = ZONE_SOUTH;    break;
            default: ErrorDetected = true;
          }
          break;
        }

        // Pivot from corner1 zone.
        case(ZONE_CORNER_1):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_CORNER_3; break;
            case(IDX_NORTH): jData[0] = ZONE_EAST;     break;
            case(IDX_WEST):  jData[0] = ZONE_SOUTH;    break;
            case(IDX_EAST):  jData[0] = ZONE_CORNER_0; break;
            default: ErrorDetected = true;
          }
          break;
        }

        // Pivot from corner2 zone.
        case(ZONE_CORNER_2):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_WEST;     break;
            case(IDX_NORTH): jData[0] = ZONE_CORNER_0; break;
            case(IDX_WEST):  jData[0] = ZONE_CORNER_3; break;
            case(IDX_EAST):  jData[0] = ZONE_NORTH;    break;
            default: ErrorDetected = true;
          }
          break;
        }

        // Pivot from corner3 zone.
        case(ZONE_CORNER_3):
        {
          // Assign matching zone depending on current iFace.
          switch(iFace){
            case(IDX_SOUTH): jData[0] = ZONE_EAST;     break;
            case(IDX_NORTH): jData[0] = ZONE_CORNER_1; break;
            case(IDX_WEST):  jData[0] = ZONE_NORTH;    break;
            case(IDX_EAST):  jData[0] = ZONE_CORNER_2; break;
            default: ErrorDetected = true;
          }
          break;
        }

        default: ErrorDetected = true;
      }

      break;
    }

    // Horizontal-buffer zone strategy.
    case(MULTIZONE_STRATEGY_HORIZONAL):
    {
      // Determine which zone this is.
      switch(iZone){

        // This is the main zone.
        case(HORIZONTAL_ZONE_MAIN):
        {
          // Determine which face this is.
          switch(iFace){
            case(IDX_SOUTH): case(IDX_NORTH): jData[0] = HORIZONTAL_ZONE_MAIN; break;
            case(IDX_EAST):                   jData[0] = HORIZONTAL_ZONE_EAST; break;
            case(IDX_WEST):                   jData[0] = HORIZONTAL_ZONE_WEST; break;
            default: ErrorDetected = true;
          }

          break;
        }

        // This is the west zone.
        case(HORIZONTAL_ZONE_WEST):
        {
          // Determine which face this is.
          switch(iFace){
            case(IDX_SOUTH): case(IDX_NORTH): jData[0] = HORIZONTAL_ZONE_WEST; break;
            case(IDX_EAST):                   jData[0] = HORIZONTAL_ZONE_MAIN; break;
            case(IDX_WEST):                   jData[0] = HORIZONTAL_ZONE_EAST; break;
            default: ErrorDetected = true;
          }

          break;
        }

        // This is the east zone.
        case(HORIZONTAL_ZONE_EAST):
        {
          // Determine which face this is.
          switch(iFace){
            case(IDX_SOUTH): case(IDX_NORTH): jData[0] = HORIZONTAL_ZONE_EAST; break;
            case(IDX_EAST):                   jData[0] = HORIZONTAL_ZONE_WEST; break;
            case(IDX_WEST):                   jData[0] = HORIZONTAL_ZONE_MAIN; break;
            default: ErrorDetected = true;
          }

          break;
        }

        default: ErrorDetected = true;
      }

      break;
    }

    // Vertical-buffer zone strategy.
    case(MULTIZONE_STRATEGY_VERTICAL):
    {
      // Determine which zone this is.
      switch(iZone){

        // This is the main zone.
        case(VERTICAL_ZONE_MAIN):
        {
          // Determine which face this is.
          switch(iFace){
            case(IDX_EAST): case(IDX_WEST): jData[0] = VERTICAL_ZONE_MAIN;  break;
            case(IDX_SOUTH):                jData[0] = VERTICAL_ZONE_NORTH; break;
            case(IDX_NORTH):                jData[0] = VERTICAL_ZONE_SOUTH; break;
            default: ErrorDetected = true;
          }

          break;
        }

        // This is the south zone.
        case(VERTICAL_ZONE_SOUTH):
        {
          // Determine which face this is.
          switch(iFace){
            case(IDX_EAST): case(IDX_WEST): jData[0] = VERTICAL_ZONE_SOUTH; break;
            case(IDX_SOUTH):                jData[0] = VERTICAL_ZONE_NORTH; break;
            case(IDX_NORTH):                jData[0] = VERTICAL_ZONE_MAIN;  break;
            default: ErrorDetected = true;
          }

          break;
        }

        // This is the north zone.
        case(VERTICAL_ZONE_NORTH):
        {
          // Determine which face this is.
          switch(iFace){
            case(IDX_EAST): case(IDX_WEST): jData[0] = VERTICAL_ZONE_NORTH; break;
            case(IDX_SOUTH):                jData[0] = VERTICAL_ZONE_MAIN;  break;
            case(IDX_NORTH):                jData[0] = VERTICAL_ZONE_SOUTH; break;
            default: ErrorDetected = true;
          }

          break;
        }

        default: ErrorDetected = true;
      }

      break;
    }

    default: ErrorDetected = true;
  }

  // Terminate immediately if an error is detected.
  if( ErrorDetected )
    Terminate("CConfig::MatchInterface", __FILE__, __LINE__,
              "Combination of iZone and iFace could not be matched!");


  // Return value of matching face.
  return jData;
}




