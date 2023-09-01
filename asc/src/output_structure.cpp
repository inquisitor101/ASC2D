#include "output_structure.hpp"




COutput::COutput
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 as3vector2d<unsigned long> &MapGlobalToLocal
)
 /*
	* Constructor, used to initiate COutput class.
	*/
{
	// Extract number of zones.
	nZone = config_container->GetnZone();

	// Initialize and pre-process the VTK object.
	Preprocess_VTK(config_container, geometry_container, MapGlobalToLocal);
	// Initialize and pre-process the restart file object.
	Preprocess_RestartSolution(config_container, geometry_container);

	// Reset binary data header files to false.
	WriteHeaderInfoDataFile = true;

	// Reset VTK file number to zero.
	FileNumberVTK = 0;
  // Reset zone data solution file number to zero.
  FileNumberZoneData = 0;
  // Reset data solution file number to zero.
  FileNumberDataSolution = 0;
}


COutput::~COutput
(
 void
)
 /*
	* Destructor for COutput class, frees allocated memory.
	*/
{
	if( vtk_container     != nullptr ) delete vtk_container;
	if( restart_container != nullptr ) delete restart_container;
}


void COutput::Preprocess_VTK
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 as3vector2d<unsigned long> &MapGlobalToLocal
)
 /*
	* Function which pre-processes the VTK object.
	*/
{
	// Determine the type of VTK file format specified.
	switch( config_container->GetTypeFileFormatVTK() ){
	
		// VTK file is in ASCII format.
		case( ASCII_FORMAT  ): 
		{ 
			vtk_container = new CASCIIFileVTK(config_container, geometry_container, MapGlobalToLocal);
			break; 
		}
		
		// VTK file is in binary format.
		case( BINARY_FORMAT ): 
		{ 
			vtk_container = new CBinaryFileVTK(config_container, geometry_container, MapGlobalToLocal);
			break; 
		}
	
		// Unknown VTK file format, report error and exit.
		default:
			Terminate("COutput::Preprocess_VTK", __FILE__, __LINE__, 
			 	       "Unknown VTK file format specified.");
	}
}


void COutput::Preprocess_RestartSolution
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function which pre-processes the restart solution object.
	*/
{
	// Determine the type of restart file format specified.
	switch( config_container->GetTypeFileFormatSolution() ){
	
		// Restart file is in ASCII format.
		case( ASCII_FORMAT  ): 
		{ 
			restart_container = new CASCIIRestart(config_container, geometry_container);
			break; 
		}
		
		// Restart file is in binary format.
		case( BINARY_FORMAT ): 
		{ 
			restart_container = new CBinaryRestart(config_container, geometry_container);
			break; 
		}
	
		// Unknown restart file format, report error and exit.
		default:
			Terminate("COutput::Preprocess_RestartSolution", __FILE__, __LINE__, 
			 	       "Unknown restart solution file format specified.");
	}
}


void COutput::WriteGNUplot
(
 CConfig                *config_container,
 const unsigned long     iIter,
 const as3double         SimTime,
 as3vector1d<as3double> &MonitoringData
)
 /*
	* Function that writes a GNU-plot output file.
	*/
{
  // Create GNU-plot file.
	std::ofstream GNU_File;

	// Deduce filename.
	std::stringstream ss;
	ss << config_container->GetGNUplotFilename() << ".dat";
	GNU_File.open(ss.str().c_str(), std::ios::app | std::ios::out);

	// Check if file can be open.
	if( !GNU_File.is_open() )
		Terminate("COutput::WriteGNUplot", __FILE__, __LINE__,
							"GNU-plot file could not be opened.");

	// Open file and write data.
	GNU_File.precision(10);
	GNU_File << std::scientific 
		        << SimTime           << ", "
		        << MonitoringData[0] << "\n";

  // Close file.
	GNU_File.close();
}


void COutput::WriteBoundaryDataToFile
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CElement      *element_container,
 CInitial      *initial_container,
 CSolver       *solver_container,
 unsigned short iSampleBoundary,
 as3double      SimTime,
 unsigned long  IterCount
)
 /*
	* Function which writes the data on a selected surface to file.
	*/
{
	// Abbreviation involving gamma.
	const as3double gm1 = GAMMA - 1.0;

	// Extract current zone geometry.
  auto* grid_zone          = geometry_container->GetGeometryZone(ZONE_MAIN);
	// Extract current boundary index.
	unsigned short iBoundary = config_container->GetSampleDataBoundaryID()[iSampleBoundary];

  // Get the indicial number for the solution on the face.
  auto& NodeList = element_container->GetIndexDOFsSol(iBoundary);
	// Extract current boundary element indices.
	auto& ElemList = solver_container->GetBoundaryContainer(iBoundary)->GetElemIndexI();  

	// Size of the total number of elements on this boundary.
	const unsigned long  nbElem = ElemList.size();
	// Size of the total number of nodes per element.
	const unsigned short nbNode = NodeList.size();

  // Extract the current unit-normal for this boundary.
  const as3double nx = geometry_container->GetUnitNormal(iBoundary)[0];
	const as3double ny = geometry_container->GetUnitNormal(iBoundary)[1];
	
	// Compute the tangential coordinate index, i.e. spatially-varying.
	const unsigned short sn = ( iBoundary == IDX_SOUTH || iBoundary == IDX_NORTH ) ? 0 : 1;   

	// Number of variables to write.
	// These are: xb, w(-), w(+).
	const unsigned short nVarWrite = 3;

	// // 
	// Assemble the buffer data to write.
	// //

	// Total number of data to write in buffer.
	const unsigned long nData = nbElem*nbNode*nVarWrite;

	// Initialize the buffer vector.
	as3vector1d<as3double> writeBuf(nData, 0.0);

	// Initialize index for the buffer data.
	unsigned long idx = 0;

	// Loop over all elements on the current boundary.
	for(unsigned long iElem : ElemList){
		
		// Extract current element solution   data.
		auto& data = solver_container->GetDataContainer(iElem)->GetDataDOFsSol();
		// Extract current element coordinate data.
    auto& coor = grid_zone->GetGeometryElem(iElem)->GetCoordSolDOFs();

		// Loop over all the surface nodes in each element.
		for(unsigned short l : NodeList){
			
			// Compute the primitive variables.
			const as3double rho   = data[0][l];
			const as3double ovrho = 1.0/rho;
      as3double u           = ovrho*  data[1][l];
      as3double v           = ovrho*  data[2][l];
      as3double p           = gm1*(   data[3][l]
                            - 0.5*( u*data[1][l] + v*data[2][l] ) );

			// Extract the relevant coordinates.
			as3vector1d<as3double> xy = { coor[0][l], coor[1][l] };

			// Initialize mean-state primitive vector.
			as3vector1d<as3double> Qinf(nVar, 0.0);

			// Compute the target-state in primitive form on the boundary. Note, this also
			// is the average state too, for the considered ICs used.
			initial_container->ComputeTargetStatePrimitivePerDOF(xy, Qinf);

			// Extract mean normal-velocity and pressure.
			const as3double uninf = nx*Qinf[1] 
				                    + ny*Qinf[2];
			const as3double pinf  =    Qinf[3];
			
			// Compute the local speed of sound.
			const as3double a  = sqrt( GAMMA*p*ovrho );
			// Compute the normal velocity.
			const as3double un = nx*u + ny*v; 

			// Compute the fluctuations of the normal velocity and pressure.
			const as3double dun = un - uninf;
			const as3double dp  = p  - pinf;

			// Assemble the incoming and outgoing characteristic.
			const as3double wp = dp + rho*a*dun;
			const as3double wm = dp - rho*a*dun;

			// Book-keep data values at this DOF.
			writeBuf[idx++] = coor[sn][l];
			writeBuf[idx++] = wm;
			writeBuf[idx++] = wp;
		}
	}


	// //
	// Write the data in binary format.
	// // // //

	// Deduce filename and number.
	std::stringstream ss;
	ss << config_container->GetOutputSampleSurfaceDirectory()
		 << config_container->GetNameMarkerSampleSurface()[iSampleBoundary]
		 << "/data"
     << "_" << IterCount << ".bin";

	// Open file.
	FILE *fh = std::fopen( ss.str().c_str(), "wb" ); 

	// Report error if file could not be opened.
	if( !fh )
		Terminate("COutput::WriteBoundaryDataToFile", __FILE__, __LINE__,
				      "Could not open file for writing data at boundary.");

	// Write the integer header.
	std::fwrite(&AS3_MagicNumber, 1, sizeof(int), fh);

	// Write the physical time in the simulation.
	std::fwrite(&SimTime, 1, sizeof(as3double), fh);

	// Write the total number of element on this boundary.
	std::fwrite(&nbElem, 1, sizeof(unsigned long), fh);

	// Write the total number of nodal indices on each element.
	std::fwrite(&nbNode, 1, sizeof(unsigned short), fh);

	// Write the total number of variables per DOFs.
	std::fwrite(&nVarWrite, 1, sizeof(unsigned short), fh);

	// Write the names of 3 global indices and the data.
	char varName[CGNS_STRING_SIZE];
  strncpy(varName, "xb", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

	strncpy(varName, "w(-)", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

	strncpy(varName, "w(+)", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

	// Write the actual processed data.
	std::fwrite(writeBuf.data(), writeBuf.size(), sizeof(as3double), fh);

	// Close the file.
	std::fclose(fh);
}


void COutput::WriteSolutionToFile
(
  CConfig    *config_container,
  CGeometry  *geometry_container,
  CElement  **element_container,
  CSolver   **solver_container,
  as3double   SimTime
)
 /*
  * Function that writes a restart solution file, depending on format specified.
  */
{
	// Write solution output.
	restart_container->WriteSolutionToFile(config_container,
			                                   geometry_container,
																				 element_container,
																				 solver_container,
																				 SimTime,
																				 FileNumberDataSolution);

	// Update file number.
	FileNumberDataSolution++;
}


void COutput::WriteDataToFile
(
  CConfig                      *config_container,
  CGeometry                    *geometry_container,
	const char                   *filename,
	const unsigned short          filetype,
  const as3double               time,
  const as3vector1d<as3double> &writebuf
)
 /*
  * Function that writes a data profile to a file.
  */
{
	if( filetype == ASCII_FORMAT )
	{
		// Specify filename.
		std::ofstream Data_File;
		std::stringstream ss;
		ss << filename << ".csv";
		
		// Open data file.
		Data_File.open(ss.str().c_str(), std::ios::app | std::ios::out);

		// Check if file could be opened.
		if( !Data_File.is_open() )
		Terminate("COutput::WriteDataToFile", __FILE__, __LINE__,
							"Data (processed/probed) file in ASCII format could not be opened.");

	  // Open file and write precision.
		Data_File.precision(10);
		
		// Write the temporal stamp.
		Data_File << std::scientific << time << ",";
		// Write the remaining variables.
		for(unsigned short i=0; i<writebuf.size(); i++)
			Data_File << std::scientific << writebuf[i] << ",";

		// Move to next line.
		Data_File << "\n";

		// Close file.
		Data_File.close();
	}
	else if( filetype == BINARY_FORMAT )
	{
		// Specify filename.
		std::stringstream ss;
		ss << filename << ".bin";

		// Open data file.
		FILE *fh = std::fopen( ss.str().c_str(), "ab" ); 

		// Report error if file could not be opened.
		if( !fh )
			Terminate("COutput::WriteBoundaryDataToFile", __FILE__, __LINE__,
					      "Data (processed/probed) file in binary format could not be opened.");

		// Write the integer header, only if this is the header.
		if( WriteHeaderInfoDataFile ) 
		{
			std::fwrite(&AS3_MagicNumber, 1, sizeof(int), fh);
			WriteHeaderInfoDataFile = false;
		}

		// Write the physical time in the simulation.
		std::fwrite(&time, 1, sizeof(as3double), fh);

		// Write the actual processed data.
		std::fwrite(writebuf.data(), writebuf.size(), sizeof(as3double), fh);

		// Close the file.
		std::fclose(fh);
	}
	else
	{
		// This is an unknown file type encountered, report error and exit.
		Terminate("COutput::WriteDataToFile", __FILE__, __LINE__, "Unknown data file type.");
	}
}


void COutput::WriteZoneDataToFile
(
  CConfig    *config_container,
  CGeometry  *geometry_container,
  CSolver   **solver_container,
  as3double   SimTime
)
 /*
  * Function that writes a zone data to a file.
  */
{
  // Extract required zone.
  unsigned short iZone = config_container->GetTypeZoneData();
  // Extract current zone geometry.
  auto* grid_zone      = geometry_container->GetGeometryZone(iZone);
  // Extract current zone solution.
  auto& data_zone      = solver_container[iZone]->GetDataContainer();
  // Extract total number of elements in zone.
  unsigned long nElem  = data_zone.size();
  // Extract total number of nodes per element in zone.
  unsigned short nNode = solver_container[iZone]->GetnDOFsSol2D();
  // Extract polynomial order of elements in zone.
  unsigned short nPoly = grid_zone->GetnPolySol();


	// Decide what type of output format used.
	// NOTE, in binary format, the coordinates are dismissed, along with nPoly and iZone.
	if( config_container->GetTypeFileFormatZone() == ASCII_FORMAT )
	{
  	// Report ASCII output.
  	std::cout << "----------------------------------------------"
								 "----------------------------------------------\n"
							<< "Writing zone data in ASCII format to file... ";

		// Open data file and write the header.
		std::ofstream Data_File;

  	// Create output file.
		std::stringstream ss;
		ss << config_container->GetOutputZoneDataFilename()
  	   << "_" << FileNumberZoneData << ".csv";
		Data_File.open(ss.str().c_str(), std::ios::out);

		// Check if file can be open.
		if( !Data_File.is_open() )
			Terminate("COutput::WriteZoneDataToFile", __FILE__, __LINE__,
					      "Could not open file for writing zone sampling in ASCII format.");

  	// Set precision of data.
  	Data_File.precision(15);
  	// Write time stamp and collective information.
  	Data_File << "# " << std::fixed
  	          << SimTime << ","
  	          << iZone   << ","
  	          << nElem   << ","
  	          << nNode   << ","
  	          << nPoly   << "\n";

  	// Write header information.
  	Data_File << "\"" << "X"          << "\"" << ","
  	          << "\"" << "Y"          << "\"" << ","
  	          << "\"" << "Density"    << "\"" << ","
  	          << "\"" << "Momentum:0" << "\"" << ","
  	          << "\"" << "Momentum:1" << "\"" << ","
  	          << "\"" << "Energy"     << "\"" << "\n";

  	// Write zone information to file.
  	for(unsigned long iElem=0; iElem<nElem; iElem++){

  	  // Extract current element geometry.
  	  auto* grid_element = grid_zone->GetGeometryElem(iElem);
  	  // Extract current element solution.
  	  auto* data_element = data_zone[iElem];

  	  // Obtain current coordinates of solution DOFs.
  	  auto& pos = grid_element->GetCoordSolDOFs();
  	  // Obtain current working variables of solution data.
  	  auto& var = data_element->GetDataDOFsSol();

  	  // Loop over every node and write working variables and coordinates.
  	  for(unsigned short l=0; l<nNode; l++){
  	    Data_File << std::scientific
  	              << pos[0][l] << ","
  	              << pos[1][l] << ","
  	              << var[0][l] << ","
  	              << var[1][l] << ","
  	              << var[2][l] << ","
  	              << var[3][l] << "\n";
  	  }
  	}

  	// Close file.
		Data_File.close();
	}
	else if( config_container->GetTypeFileFormatZone() == BINARY_FORMAT )
	{	
  	// Report binary output.
  	std::cout << "----------------------------------------------"
								 "----------------------------------------------\n"
							<< "Writing zone data in binary format to file... ";

		// Deduce filename and number.
		std::stringstream ss;
		ss << config_container->GetOutputZoneDataFilename()
  	   << "_" << FileNumberZoneData << ".bin";

		// Open file.
		FILE *fh = std::fopen( ss.str().c_str(), "wb" ); 

		// Report error if file could not be opened.
		if( !fh )
			Terminate("COutput::WriteZoneDataToFile", __FILE__, __LINE__,
					      "Could not open file for writing zone sampling in binary format.");

		// Write the integer header.
		std::fwrite(&AS3_MagicNumber, 1, sizeof(int), fh);

		// Write the current simulation time.
		std::fwrite(&SimTime, 1, sizeof(as3double), fh);

		// Write the total number of elements in this zone.
		std::fwrite(&nElem, 1, sizeof(unsigned long), fh);

		// Write the total number of nodes in each element.
		std::fwrite(&nNode, 1, sizeof(unsigned short), fh);

		// Loop over all the elements and write the data in this zone.
		for(unsigned long iElem=0; iElem<nElem; iElem++)
		{
      // Extract current solution data.
      auto& data = data_zone[iElem]->GetDataDOFsSol();

			// Write the solution in this element for each variable.
			for(unsigned short iVar=0; iVar<nVar; iVar++)
				std::fwrite(data[iVar], nNode, sizeof(as3double), fh);
		}

		// Close the file.
		std::fclose(fh);
	}
	else 
		Terminate("COutput::WriteZoneDataToFile", __FILE__, __LINE__,
				      "Unknown zone sampling file format specified.");

	// All is complete!
	std::cout << "Done." << std::endl;
  std::cout << "----------------------------------------------"
							 "----------------------------------------------\n" << std::endl;
  // Increment file number.
  FileNumberZoneData++;
}


void COutput::WriteFileVTK
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 CElement  **element_container,
 CSolver   **solver_container
)
 /*
	* Function that writes VTK data file, depending on the format.
	*/
{
	// Write file, depending on the format specified.
	vtk_container->WriteFileVTK(config_container, 
			                        geometry_container,
															element_container,
															solver_container, 
															FileNumberVTK);

	// Update file counter.
	FileNumberVTK++;
}





