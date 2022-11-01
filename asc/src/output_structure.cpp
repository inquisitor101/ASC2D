#include "output_structure.hpp"




COutput::COutput
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor, used to initiate COutput class.
	*/
{
	// Extract number of zones.
	nZone = config_container->GetnZone();

	// Extract output VTK filename.
	OutputVTKFilename = config_container->GetOutputVTKFilename();

	// Choose working data to write the output from. Note, for now
	// it is fixed as conservative working-variables.
	VariableDensity     = new CDensityConservative();
	VariableMomentum    = new CMomentumConservative();
	VariableEnergy      = new CEnergyConservative();
	VariablePressure    = new CPressureConservative();
  VariableTemperature = new CTemperatureConservative();
  VariableMachNumber  = new CMachNumberConservative();

	// Initialize local connectivity matrix.
	ConnLocal = nullptr;
	ConnLocal = new unsigned short*[nZone];

	// Populate connectivity matrix, per zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){

		// Extract polynomial order in zone.
		unsigned short nPoly    = geometry_container->GetGeometryZone(iZone)->GetnPolySol();
		// Deduce number of sub-elements, i.e. nPoly=1 elements within each element.
		unsigned short nSubElem = nPoly*nPoly;

		// Reserve memory for connectivity matrix in zone.
		ConnLocal[iZone] = new unsigned short[nSubElem*N_POINTS_QUADRILATERAL];

		// Number of sub-elements per dimension.
		unsigned short nxSubElem = sqrt(nSubElem);
		unsigned short nySubElem = nxSubElem;

		// Book-keeping indices.
		short idxLocal 						= -1;
		unsigned short ijSubElem 	=  0;

		// Assemble connectivity matrix.
		for(unsigned short jSubElem=0; jSubElem<nySubElem; jSubElem++){
			for(unsigned short iSubElem=0; iSubElem<nxSubElem; iSubElem++){
				ConnLocal[iZone][++idxLocal] = ijSubElem;
				ConnLocal[iZone][++idxLocal] = ijSubElem+1;
				ConnLocal[iZone][++idxLocal] = ijSubElem+1+nPoly+1;
				ConnLocal[iZone][++idxLocal] = ijSubElem+nPoly+1;
				ijSubElem++;
			}
			ijSubElem++;
		}
	}

	// Reset VTK file number to zero.
	FileNumberVTK = 0;
  // Reset zone data solution file number to zero.
  FileNumberZoneData = 0;
  // Resert processed data file number to zero.
  FileNumberProcessed = 0;
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
	if( VariableDensity     != nullptr ) delete VariableDensity;
	if( VariableMomentum    != nullptr ) delete VariableMomentum;
	if( VariableEnergy      != nullptr ) delete VariableEnergy;
	if( VariablePressure    != nullptr ) delete VariablePressure;
  if( VariableMachNumber  != nullptr ) delete VariableMachNumber;
  if( VariableTemperature != nullptr ) delete VariableTemperature;

	if(ConnLocal != nullptr){
		for(unsigned short iZone=0; iZone<nZone; iZone++)
			if(ConnLocal[iZone] != nullptr) delete [] ConnLocal[iZone];
		delete [] ConnLocal;
	}
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
  *
  */
{
  // Report output.
  std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution data to file... ";

  // Extract boundaries of main domain.
  auto& DomainBound = config_container->GetDomainBound();

  // Create data file.
	std::ofstream Data_File;

  // Deduce filename and number.
	std::stringstream ss;
	ss << config_container->GetOutputSolFilename()
     << "_" << FileNumberDataSolution << ".dat";
	Data_File.open(ss.str().c_str(), std::ios::out);

  // Open file and write header and info.
	Data_File.precision(10);
	Data_File << "# Data file for saving solution\n";
	Data_File << "\n\n";

	Data_File << "# physical time: \n";
	Data_File << "SimTime = " << SimTime << "\n";
	Data_File << "\n";

	Data_File << "# number of zones" << "\n";
	Data_File << "nZone = " << nZone << "\n";
	Data_File << "\n";

	Data_File << "# physical domain bounding coordinates\n";
	Data_File << std::showpos << "XMin = " << DomainBound[0] << "\n";
	Data_File << std::showpos << "XMax = " << DomainBound[1] << "\n";
	Data_File << std::showpos << "YMin = " << DomainBound[2] << "\n";
	Data_File << std::showpos << "YMax = " << DomainBound[3] << "\n";
	Data_File << "\n\n";

  // Loop over all zones and write data according to zone time.
  for(unsigned short iZone=0; iZone<nZone; iZone++){

    // Extract element information.
    const unsigned short TypeDOFs   = element_container[iZone]->GetTypeDOFs();
    const unsigned short nPolySol   = element_container[iZone]->GetnPolySol();
    const unsigned short nDOFsSol1D = element_container[iZone]->GetnDOFsSol1D();
    const unsigned short nDOFsSol2D = element_container[iZone]->GetnDOFsSol2D();
    const unsigned short nDOFsInt1D = element_container[iZone]->GetnDOFsInt1D();

    // Extract zone information.
    const unsigned long nxElem = config_container->GetnxElemZone(iZone);
    const unsigned long nyElem = config_container->GetnyElemZone(iZone);
    const unsigned long nElem  = solver_container[iZone]->GetnElem();

    // Extract type of solver.
    const unsigned short TypeSolver      = config_container->GetTypeSolver(iZone);
    // Extract type of buffer layer.
    const unsigned short TypeBufferLayer = config_container->GetTypeBufferLayer(iZone);
    // Extract type of zone.
    const unsigned short TypeZone        = config_container->GetTypeZone(iZone);
    // Extract type of IC.
    const unsigned short TypeIC          = config_container->GetTypeIC(iZone);

    // Extract data container in this zone.
    auto& data_container = solver_container[iZone]->GetDataContainer();

    // Consistency check.
    assert( nxElem*nyElem == nElem );

    // Skip line for readability.
    Data_File << "\n";

    // Current zone ID.
		Data_File << "# current zone\n";
		Data_File << "iZone = " << iZone << "\n";

    // Write solver type in this zone.
    Data_File << "# type of initial condition in this zone\n";
    Data_File << "TypeIC = " << TypeIC << "\n";

    // Write solver type in this zone.
    Data_File << "# type of solver in this zone\n";
    Data_File << "TypeSolver = " << TypeSolver << "\n";

    // Write buffer-layer type in this zone.
    Data_File << "# type of buffer layer in this zone\n";
    Data_File << "TypeBufferLayer = " << TypeBufferLayer << "\n";

    // Write current zone type.
    Data_File << "# type of current zone\n";
    Data_File << "TypeZone = " << TypeZone << "\n";

    // Write type of DOFs.
    Data_File << "# type of DOFs\n";
    Data_File << "TypeDOFs = " << TypeDOFs << "\n";

    // Write polynomial order.
    Data_File << "# polynomial order\n";
    Data_File << "nPolySol = " << nPolySol << "\n";

    // Write number of solution DOFs in 1D.
    Data_File << "# number of solution DOFs in 1D\n";
    Data_File << "nDOFsSol1D = " << nDOFsSol1D << "\n";

    // Write number of integration DOFs in 1D.
    Data_File << "# number of integration DOFs in 1D\n";
    Data_File << "nDOFsInt1D = " << nDOFsInt1D << "\n";

    // Write number of elements.
    Data_File << "# number of elements\n";
    Data_File << "nElem = " << nElem << "\n";

    // Write number of elements in x-direction.
    Data_File << "# number of elements in x-direction\n";
    Data_File << "nxElem = " << nxElem << "\n";

    // Write number of elements in y-direction.
    Data_File << "# number of elements in y-direction\n";
    Data_File << "nyElem = " << nyElem << "\n";

    // Write solution data in this zone.
    Data_File << "# solution data, per element, per variable, per solution node\n";

    // Loop over all elements and write the data.
    for(unsigned long iElem=0; iElem<nElem; iElem++){

      // Extract current solution data.
      auto& data = data_container[iElem]->GetDataDOFsSol();
      // Write data of current solution element.
      Data_File << "DataDOFsSol = " << iElem << "\n";

      // Write solution of current element.
      for(unsigned short iVar=0; iVar<nVar; iVar++){
        for(unsigned short iNode=0; iNode<nDOFsSol2D; iNode++){
          Data_File << std::scientific << data[iVar][iNode] << ", ";
        }
      }
      Data_File << "\n";
    }

    // If this is a PML zone, write auxiliary data too.
    if( TypeBufferLayer == PML_LAYER ){

      // Loop over all elements and write the data.
      for(unsigned long iElem=0; iElem<nElem; iElem++){

        // Extract current auxiliary data.
        auto** data = &data_container[iElem]->GetDataDOFsSol()[nVar];
        // Write data of current solution element.
        Data_File << "DataDOFsSolAux = " << iElem << "\n";

        // Write solution of current element.
        for(unsigned short iVar=0; iVar<nVar; iVar++){
          for(unsigned short iNode=0; iNode<nDOFsSol2D; iNode++){
            Data_File << std::scientific << data[iVar][iNode] << ", ";
          }
        }
        Data_File << "\n";
      }
    }

  }

  // Close file.
	Data_File.close();

	// Report progress.
	std::cout << "Done." << std::endl;
  std::cout << "----------------------------------------------"
							 "----------------------------------------------" << std::endl;

	// Update file number.
	FileNumberDataSolution++;
}


void COutput::WriteDataToFile
(
  CConfig                *config_container,
  CGeometry              *geometry_container,
  const char             *fileinfo,
  as3vector1d<as3double> &time,
  as3vector2d<as3double> &data
)
 /*
  * Function that writes a data profile to a file.
  */
{
  // Report output.
  std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing profile of processed data to file... ";

	// Open data file and write the header.
	std::ofstream Data_File;

	std::stringstream ss;
	ss << config_container->GetOutputProcessedFilename()
     << "_" << FileNumberProcessed << ".dat";
	Data_File.open(ss.str().c_str(), std::ios::out);

	Data_File.precision(10);
	Data_File << "# Data file for saving (post-)processed information: "
            << fileinfo;
	Data_File << "\n#\n";

  Data_File << "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #";
  Data_File << "\n";

  // Consistency check.
  for(unsigned short i=0; i<data.size(); i++) assert( time.size() == data[i].size() );

  // Write data to file.
  for(unsigned long i=0; i<time.size(); i++){
    Data_File << std::scientific << time[i] << ", ";
    for(unsigned short j=0; j<data.size(); j++){
      Data_File << std::scientific << data[j][i] << ", ";
    }
    Data_File << "\n";
  }

	// Close file.
	Data_File.close();

	// All is complete!
	std::cout << "Done." << std::endl;

  // Increment file number.
  FileNumberProcessed++;
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
  // Report output.
  std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing zone data to file... ";

	// Open data file and write the header.
	std::ofstream Data_File;

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

  // Create output file.
	std::stringstream ss;
	ss << config_container->GetOutputZoneDataFilename()
     << "_" << FileNumberZoneData << ".csv";
	Data_File.open(ss.str().c_str(), std::ios::out);


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
 CSolver   **solver_container
)
 /*
	* Function that writes VTK data file.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution in VTK format... " << std::endl;

	// Create output stream.
	std::ofstream Paraview_File;

	// Open string stream.
	std::stringstream ss;
	ss << OutputVTKFilename << "_" << FileNumberVTK << ".vtk";

	// Open ASCII file.
	Paraview_File.open(ss.str().c_str(), std::ios::out);

	// Check if file can be open.
	if( !Paraview_File.is_open() )
		Terminate("COutput::WriteFileVTK", __FILE__, __LINE__,
							"VTK file directory could not be located!");

	// Write header.
	Paraview_File.precision(6);
	Paraview_File << "# vtk DataFile Version 3.0\n";
	Paraview_File << "vtk output\n";
	Paraview_File << "ASCII\n";
	Paraview_File << "DATASET UNSTRUCTURED_GRID\n";

	// Initialize nPoly=1 element (index) size.
	unsigned short nNodeP1 = N_POINTS_QUADRILATERAL;

	// Extract total number of points assuming nPoly=1 elements.
	unsigned long nPoints = geometry_container->GetnPointSubElemP1();
	Paraview_File << "POINTS " << nPoints << " double\n";

	// Display progress.
	std::cout << "  writing grid data....... ";

	unsigned short idxElem   = 0;
	unsigned long  nDOFsGrid = 0;
	unsigned short nZone     = geometry_container->GetnZone();


	// Iterate on each zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){

		const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);

		unsigned short nPoly     = gridZone->GetnPolySol();
		unsigned short nSubElem  = nPoly*nPoly;
		unsigned long  nElem     = gridZone->GetnElem();

		for(unsigned long iElem=0; iElem<nElem; iElem++){
			const CGeometryElement *surfElem = gridZone->GetGeometryElem(iElem);
			as3data1d<as3double> elemNode = surfElem->GetCoordSolDOFs();

			idxElem = 0;
			for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
				for(unsigned short iNode=0; iNode<4; iNode++){
					for(unsigned short iDim=0; iDim<nDim; iDim++)
						Paraview_File << std::scientific << elemNode[iDim][ConnLocal[iZone][idxElem]] << "\t";
					// accout for z-coordinate
					Paraview_File << std::scientific << "0.0" << "\t";
					idxElem++;
					nDOFsGrid++;
				}
				Paraview_File << "\n";
			}
		}
	}


	// Register cell indices.
	unsigned long nSubElemP1Global 		 = nPoints/N_POINTS_QUADRILATERAL;
	unsigned long nGlobal_Elem_Storage = nSubElemP1Global*(N_POINTS_QUADRILATERAL+1);
	Paraview_File << "\nCELLS " << nSubElemP1Global << "\t" << nGlobal_Elem_Storage << "\n";

	unsigned long idxElemGlobal = 0;
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);
		unsigned long nElem = gridZone->GetnElem();

		unsigned short nPoly    = gridZone->GetnPolySol();
		unsigned short nSubElem = nPoly*nPoly;

		short idxLocal;
		for(unsigned long iElem=0; iElem<nElem; iElem++){
			idxLocal = 0;
			for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
				Paraview_File << N_POINTS_QUADRILATERAL << "\t";
				for(unsigned short iNode=0; iNode<nNodeP1; iNode++){
					Paraview_File << idxElemGlobal+idxLocal << "\t";
					idxLocal++;
				}
				Paraview_File << "\n";
			}
			idxElemGlobal += idxLocal;
		}
	}


	// Cell registration.
	Paraview_File << "\nCELL_TYPES " << nSubElemP1Global << "\n";
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);

		unsigned long  nElem 		= gridZone->GetnElem();
		unsigned short nPoly 	  = gridZone->GetnPolySol();
		unsigned short nSubElem = nPoly*nPoly;

		for(unsigned long iElem=0; iElem<nElem; iElem++){
			for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
				Paraview_File << QUADRILATERAL << "\t";
				Paraview_File << "\n";
			}
		}
	}
	Paraview_File << "\nPOINT_DATA " << nPoints << "\n";
	std::cout << "Done." << std::endl;


	// Register density data.
	std::cout << "  writing Density.........";
	Paraview_File << "\nSCALARS " << "Density" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableDensity);
	std::cout << " Done." << std::endl;

	// Register momentum data.
	std::cout << "  writing Momentum........";
	Paraview_File << "\nVECTORS " << "Momentum" << " double\n";
	WriteVector(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableMomentum);
	std::cout << " Done." << std::endl;

	// Register energy data.
	std::cout << "  writing Energy..........";
	Paraview_File << "\nSCALARS " << "Energy" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableEnergy);
	std::cout << " Done." << std::endl;

	// Register pressure data.
	std::cout << "  writing Pressure........";
	Paraview_File << "\nSCALARS " << "Pressure" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariablePressure);
	std::cout << " Done." << std::endl;

  // Register temperature data.
	std::cout << "  writing Temperature.....";
	Paraview_File << "\nSCALARS " << "Temperature" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableTemperature);
	std::cout << " Done." << std::endl;

  // Register Mach number data.
	std::cout << "  writing Mach............";
	Paraview_File << "\nSCALARS " << "Mach" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableMachNumber);
	std::cout << " Done." << std::endl;

  // If a PML zone exists, write its auxiliary variables.
  if( config_container->GetUsePML() ){

    // If user specifies auxiliary data included, write them.
    if( config_container->GetWriteAuxiliaryDataPML() ){

      // Register Q1 of density data.
      std::cout << "  writing Q1[Density].....";
      Paraview_File << "\nSCALARS " << "Q1[Density]" << " double 1\n";
      Paraview_File << "LOOKUP_TABLE default\n";
      WriteScalarAuxPML(config_container,
                        geometry_container,
                        solver_container,
                        Paraview_File,
                        CONTQ1_VAR);
      std::cout << " Done." << std::endl;

      // Register Q1 of momentum data.
      std::cout << "  writing Q1[Momentum]....";
      Paraview_File << "\nVECTORS " << "Q1[Momentum]" << " double\n";
      WriteVectorAuxPML(config_container,
                        geometry_container,
                        solver_container,
                        Paraview_File,
                        XMOMQ1_VAR);
      std::cout << " Done." << std::endl;

      // Register Q1 of energy data.
      std::cout << "  writing Q1[Energy]......";
      Paraview_File << "\nSCALARS " << "Q1[Energy]" << " double 1\n";
      Paraview_File << "LOOKUP_TABLE default\n";
      WriteScalarAuxPML(config_container,
                        geometry_container,
                        solver_container,
                        Paraview_File,
                        ENERQ1_VAR);
      std::cout << " Done." << std::endl;
    }
  }


	// Close file.
	Paraview_File.close();

	// Report: all is complete!
	std::cout << "Done." << std::endl;

	// Update file counter.
	FileNumberVTK++;
}


void COutput::WriteScalar
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File,
 CScalarVariable *Variable
)
 /*
	* Function that writes a scalar parameter.
	*/
{
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);

		unsigned long  nElem 		= gridZone->GetnElem();
		unsigned short nPoly 	  = gridZone->GetnPolySol();
		unsigned short nSubElem = nPoly*nPoly;

		for(unsigned long iElem=0; iElem<nElem; iElem++){

			// Extract element solution.
			const auto& variables = solver_container[iZone]->GetDataContainer(iElem)->GetDataDOFsSol();

			short idxLocal = 0;
			for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
				for(unsigned short iNode=0; iNode<N_POINTS_QUADRILATERAL; iNode++){
					// Compute scalar.
					auto scalar = Variable->GetValue(variables, ConnLocal[iZone][idxLocal]);
					// Write scalar to file.
					Paraview_File << std::scientific << scalar << "\t";
					idxLocal++;
				}
				Paraview_File << "\n";
			}
		}
	}
}


void COutput::WriteVector
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File,
 CVectorVariable *Variable
)
 /*
	* Function that writes a vector parameter.
	*/
{
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);

		unsigned long  nElem 		= gridZone->GetnElem();
		unsigned short nPoly 	  = gridZone->GetnPolySol();
		unsigned short nSubElem = nPoly*nPoly;

		for(unsigned long iElem=0; iElem<nElem; iElem++){

			// Extract element solution.
			const auto& variables = solver_container[iZone]->GetDataContainer(iElem)->GetDataDOFsSol();

			short idxLocal = 0;
			for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
				for(unsigned short iNode=0; iNode<N_POINTS_QUADRILATERAL; iNode++){
					// Compute vector.
					auto vector = Variable->GetValue(variables, ConnLocal[iZone][idxLocal]);
          // Write vector to file.
					Paraview_File << std::scientific << vector[0] << "\t"
																					 << vector[1] << "\t"
																					 << 0.0       << "\t";
					idxLocal++;
				}
				Paraview_File << "\n";
			}
		}
	}
}


void COutput::WriteScalarAuxPML
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File,
 unsigned short   iVar
)
 /*
	* Function that writes a PML auxiliary scalar parameter.
	*/
{
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);

		unsigned long  nElem 		= gridZone->GetnElem();
		unsigned short nPoly 	  = gridZone->GetnPolySol();
		unsigned short nSubElem = nPoly*nPoly;

    // If this is a non-PML zone, pad with zeros.
    if( config_container->GetTypeBufferLayer(iZone) != PML_LAYER ){
      for(unsigned long iElem=0; iElem<nElem; iElem++){
        short idxLocal = 0;
        for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
          for(unsigned short iNode=0; iNode<N_POINTS_QUADRILATERAL; iNode++){
            // Write scalar to file.
            Paraview_File << std::scientific << 0.0 << "\t";
            idxLocal++;
          }
          Paraview_File << "\n";
        }
      }
    }
    else{
      // This is a PML layer.
  		for(unsigned long iElem=0; iElem<nElem; iElem++){

  			// Extract element solution.
  			const auto& aux = solver_container[iZone]->GetDataContainer(iElem)->GetDataDOFsSol();

  			short idxLocal = 0;
  			for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
  				for(unsigned short iNode=0; iNode<N_POINTS_QUADRILATERAL; iNode++){
  					// Write scalar to file.
  					Paraview_File << std::scientific << aux[iVar][ConnLocal[iZone][idxLocal]] << "\t";
  					idxLocal++;
  				}
  				Paraview_File << "\n";
  			}
  		}
    }
	}
}


void COutput::WriteVectorAuxPML
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File,
 unsigned short   iVar
)
 /*
	* Function that writes a PML auxiliary vector parameter.
	*/
{
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);

		unsigned long  nElem 		= gridZone->GetnElem();
		unsigned short nPoly 	  = gridZone->GetnPolySol();
		unsigned short nSubElem = nPoly*nPoly;

    // If this is a non-PML zone, pad with zeros.
    if( config_container->GetTypeBufferLayer(iZone) != PML_LAYER ){
      for(unsigned long iElem=0; iElem<nElem; iElem++){
        short idxLocal = 0;
        for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
          for(unsigned short iNode=0; iNode<N_POINTS_QUADRILATERAL; iNode++){
            // Write scalar to file.
            Paraview_File << std::scientific << 0.0 << "\t"
                                             << 0.0 << "\t"
                                             << 0.0 << "\t";
            idxLocal++;
          }
          Paraview_File << "\n";
        }
      }
    }
    else {
      // This is a PML layer.
  		for(unsigned long iElem=0; iElem<nElem; iElem++){

  			// Extract element solution.
  			const auto& aux = solver_container[iZone]->GetDataContainer(iElem)->GetDataDOFsSol();

  			short idxLocal = 0;
  			for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
  				for(unsigned short iNode=0; iNode<N_POINTS_QUADRILATERAL; iNode++){
            unsigned short index = ConnLocal[iZone][idxLocal];
            // Write vector to file.
  					Paraview_File << std::scientific << aux[iVar  ][index]   << "\t"
  																					 << aux[iVar+1][index] << "\t"
  																					 << 0.0       << "\t";
  					idxLocal++;
  				}
  				Paraview_File << "\n";
  			}
  		}
    }
	}
}


