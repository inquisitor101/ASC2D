#include "restart_structure.hpp"




CRestart::CRestart
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor, used to initialize CRestart class.
	*/
{
	// Extract number of zones.
	nZone = config_container->GetnZone();
}


CRestart::~CRestart
(
 void
)
 /*
	* Destructor for CRestart class, frees allocated memory. 
	*/
{

}



CASCIIRestart::CASCIIRestart
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
	:
		CRestart
		(
		 config_container,
		 geometry_container
		)
 /*
	* Constructor, used to initialize CASCIIRestart class.
	*/
{

}


CASCIIRestart::~CASCIIRestart
(
 void
)
 /*
	* Destructor for CASCIIRestart class, frees allocated memory. 
	*/
{

}


void CASCIIRestart::WriteSolutionToFile
(
  CConfig             *config_container,
  CGeometry           *geometry_container,
  CElement           **element_container,
  CSolver            **solver_container,
  as3double            SimTime,
	const unsigned long  FileNumberDataSolution
)
 /*
	* Function that writes a restart solution in ASCII format.
	*/
{
  // Report output.
  std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution data to file in ASCII format... ";

  // Extract boundaries of main domain.
  auto& DomainBound = config_container->GetDomainBound();

  // Create data file.
	std::ofstream Data_File;

  // Deduce filename and number.
	std::stringstream ss;
	ss << config_container->GetOutputSolFilename()
     << "_" << FileNumberDataSolution << ".txt";
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

    // Write initial condition type in this zone.
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
}




CBinaryRestart::CBinaryRestart
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
	:
		CRestart
		(
		 config_container,
		 geometry_container
		)
 /*
	* Constructor, used to initialize CBinaryRestart class.
	*/
{

}


CBinaryRestart::~CBinaryRestart
(
 void
)
 /*
	* Destructor for CBinaryRestart class, frees allocated memory. 
	*/
{

}


void CBinaryRestart::WriteSolutionToFile
(
  CConfig             *config_container,
  CGeometry           *geometry_container,
  CElement           **element_container,
  CSolver            **solver_container,
  as3double            SimTime,
	const unsigned long  FileNumberDataSolution
)
 /*
	* Function that writes a restart solution in binary format.
	*/
{
  // Report output.
  std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution data to file in binary format... ";


  // Deduce filename and number.
	std::stringstream ss;
	ss << config_container->GetOutputSolFilename()
     << "_" << FileNumberDataSolution << ".bin";

	// Open file.
	FILE *fh = std::fopen( ss.str().c_str(), "wb" ); 

	// Report error if file could not be opened.
	if( !fh )
		Terminate("CBinaryRestart::WriteSolutionToFile", __FILE__, __LINE__,
				      "Could not open file for writing restart solution.");


  // Extract boundaries of main domain.
  auto& DomainBound = config_container->GetDomainBound();


	// Write the integer header.
	std::fwrite(&AS3_MagicNumber, 1, sizeof(int), fh);

	// Write the current simulation time.
	std::fwrite(&SimTime, 1, sizeof(as3double), fh);
	
	// Write the total number of zones.
	std::fwrite(&nZone, 1, sizeof(unsigned short), fh);

	// Write the physical domain bounding coordinates.
	std::fwrite(DomainBound.data(), 4, sizeof(as3double), fh);


	// Loop over all the zones and write the necessary information.
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

		// Write the current zone.
		std::fwrite(&iZone, 1, sizeof(unsigned short), fh);

		// Write the type of initial condition in this zone.
		std::fwrite(&TypeIC, 1, sizeof(unsigned short), fh);

		// Write the type of solver in this zone.
		std::fwrite(&TypeSolver, 1, sizeof(unsigned short), fh);

		// Write the type of buffer-layer in this zone.
		std::fwrite(&TypeBufferLayer, 1, sizeof(unsigned short), fh);

		// Write the type of zone in this zone.
		std::fwrite(&TypeZone, 1, sizeof(unsigned short), fh);

		// Write the type of DOFs in this zone.
		std::fwrite(&TypeDOFs, 1, sizeof(unsigned short), fh);

		// Write the polynomial order used in this zone.
		std::fwrite(&nPolySol, 1, sizeof(unsigned short), fh);

		// Write the number of solution DOFs in 1D in this zone.
		std::fwrite(&nDOFsSol1D, 1, sizeof(unsigned short), fh);

		// Write the number of integration points in 1D in this zone.
		std::fwrite(&nDOFsInt1D, 1, sizeof(unsigned short), fh);

		// Write the number of elements in this zone.
		std::fwrite(&nElem, 1, sizeof(unsigned long), fh);

		// Write the number of elements in the x-direction in this zone.
		std::fwrite(&nxElem, 1, sizeof(unsigned long), fh);

		// Write the number of elements in the y-direction in this zone.
		std::fwrite(&nyElem, 1, sizeof(unsigned long), fh);


		// Loop over all the elements and write the data in this zone.
		for(unsigned long iElem=0; iElem<nElem; iElem++){

      // Extract current solution data.
      auto& data = data_container[iElem]->GetDataDOFsSol();

			// Write the current element index.
			std::fwrite(&iElem, 1, sizeof(unsigned long), fh);

			// Write the solution in this element for each variable.
			for(unsigned short iVar=0; iVar<nVar; iVar++)
				std::fwrite(data[iVar], nDOFsSol2D, sizeof(as3double), fh);
		}

		// If this is a PML zone, write the associated auxiliary data as well.
		if( TypeBufferLayer == PML_LAYER ){

			// Loop over all the elements and write the auxiliary data in this zone.
			for(unsigned long iElem=0; iElem<nElem; iElem++){

        // Extract current auxiliary data.
        auto** data = &data_container[iElem]->GetDataDOFsSol()[nVar];

				// Write the current element index.
				std::fwrite(&iElem, 1, sizeof(unsigned long), fh);

				// Write the auxiliary solution in this element for each variable.
				for(unsigned short iVar=0; iVar<nVar; iVar++)
					std::fwrite(data[iVar], nDOFsSol2D, sizeof(as3double), fh);
			}
		}
	}


	// Close the file.
	std::fclose(fh);

	// Report progress.
	std::cout << "Done." << std::endl;
}



