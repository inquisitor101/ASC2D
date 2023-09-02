#include "vtk_structure.hpp"





CFileVTK::CFileVTK
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 as3vector2d<unsigned long> &MapGlobalToLocal
)
 /*
	* Constructor, used to initialize CFileVTK class.
	*/
{
	// Extract number of zones.
	nZone = config_container->GetnZone();

	// Extract output VTK filename.
	OutputVTKFilename = config_container->GetOutputVTKFilename();

	// Initialize the variables to write to zero.
	WriteDensity     = false;
	WriteMomentum    = false;
	WriteTotalEnergy = false;
	WritePressure    = false;
	WriteVelocity    = false;
	WriteVorticity   = false;
	WriteMach        = false;
	WriteTemperature = false;
	WriteEntropy     = false;

	// Extract the user-specified variables for writing.
	auto varwrite = config_container->GetOutputVTKVariable();
	
	// Loop over variables specified for writing.
	for(unsigned short i=0; i<varwrite.size(); i++)
	{
		if( varwrite[i] == VTK_VARIABLE_DENSITY     ) WriteDensity     = true;
		if( varwrite[i] == VTK_VARIABLE_MOMENTUM    ) WriteMomentum    = true;
		if( varwrite[i] == VTK_VARIABLE_TOTALENERGY ) WriteTotalEnergy = true;
		if( varwrite[i] == VTK_VARIABLE_PRESSURE    ) WritePressure    = true;
		if( varwrite[i] == VTK_VARIABLE_VELOCITY    ) WriteVelocity    = true;
		if( varwrite[i] == VTK_VARIABLE_VORTICITY   ) WriteVorticity   = true;
		if( varwrite[i] == VTK_VARIABLE_MACH        ) WriteMach        = true;
		if( varwrite[i] == VTK_VARIABLE_TEMPERATURE ) WriteTemperature = true;
		if( varwrite[i] == VTK_VARIABLE_ENTROPY     ) WriteEntropy     = true;
	}
}


CFileVTK::~CFileVTK
(
 void
)
 /*
	* Destructor for CFileVTK, frees allocated memory. 
	*/
{

}




CASCIIFileVTK::CASCIIFileVTK
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 as3vector2d<unsigned long> &MapGlobalToLocal
)
	:
		CFileVTK
		(
		 config_container,
		 geometry_container,
		 MapGlobalToLocal
		)
 /*
	* Constructor, used to initialize CASCIIFileVTK class.
	*/
{
	// Choose working data to write the output from. 
	if( WriteDensity     ) VariableDensity     = new CDensityConservative();
	if( WriteMomentum    ) VariableMomentum    = new CMomentumConservative();
	if( WriteTotalEnergy ) VariableEnergy      = new CEnergyConservative();
	if( WritePressure    ) VariablePressure    = new CPressureConservative();
	if( WriteTemperature ) VariableTemperature = new CTemperatureConservative();
	if( WriteMach        ) VariableMachNumber  = new CMachNumberConservative();
	if( WriteEntropy     ) VariableEntropy     = new CEntropyConservative();

	// Initialize local connectivity matrix.
	ConnLocal.resize(nZone);
	
	// Populate connectivity matrix, per zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){
	
		// Extract polynomial order in zone.
		unsigned short nPoly    = geometry_container->GetGeometryZone(iZone)->GetnPolySol();
		// Deduce number of sub-elements, i.e. nPoly=1 elements within each element.
		unsigned short nSubElem = nPoly*nPoly;
	
		// Reserve memory for connectivity matrix in zone.
		ConnLocal[iZone].resize( nSubElem*N_POINTS_QUADRILATERAL );
	
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
}


CASCIIFileVTK::~CASCIIFileVTK
(
 void
)
 /*
	* Destructor for CASCIIFileVTK, frees allocated memory. 
	*/
{
	if( VariableDensity     != nullptr ) delete VariableDensity;
	if( VariableMomentum    != nullptr ) delete VariableMomentum;
	if( VariableEnergy      != nullptr ) delete VariableEnergy;
	if( VariablePressure    != nullptr ) delete VariablePressure;
  if( VariableMachNumber  != nullptr ) delete VariableMachNumber;
  if( VariableTemperature != nullptr ) delete VariableTemperature;
	if( VariableEntropy     != nullptr ) delete VariableEntropy;
}


void CASCIIFileVTK::WriteFileVTK
(
 CConfig             *config_container,
 CGeometry           *geometry_container,
 CElement           **element_container,
 CSolver            **solver_container,
 const unsigned long  FileNumberVTK
)
 /*
	* Function that writes VTK data file in ASCII format.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution in VTK ASCII format... " << std::endl;

	// Create output stream.
	std::ofstream Paraview_File;

	// Open string stream.
	std::stringstream ss;
	ss << OutputVTKFilename << "_" << FileNumberVTK << ".vtk";

	// Open ASCII file.
	Paraview_File.open(ss.str().c_str(), std::ios::out);

	// Check if file can be open.
	if( !Paraview_File.is_open() )
		Terminate("CASCIIFileVTK::WriteFileVTK", __FILE__, __LINE__,
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
	if( WriteDensity )
	{
		std::cout << "  writing Density.........";
		Paraview_File << "\nSCALARS " << "Density" << " double 1\n";
		Paraview_File << "LOOKUP_TABLE default\n";
		WriteScalar(config_container,
								geometry_container,
								solver_container,
								Paraview_File,
								VariableDensity);
		std::cout << " Done." << std::endl;
	}

	// Register momentum data.
	if( WriteMomentum )
	{
		std::cout << "  writing Momentum........";
		Paraview_File << "\nVECTORS " << "Momentum" << " double\n";
		WriteVector(config_container,
								geometry_container,
								solver_container,
								Paraview_File,
								VariableMomentum);
		std::cout << " Done." << std::endl;
	}

	// Register energy data.
	if( WriteTotalEnergy )
	{
		std::cout << "  writing Energy..........";
		Paraview_File << "\nSCALARS " << "Energy" << " double 1\n";
		Paraview_File << "LOOKUP_TABLE default\n";
		WriteScalar(config_container,
								geometry_container,
								solver_container,
								Paraview_File,
								VariableEnergy);
		std::cout << " Done." << std::endl;
	}

	// Register pressure data.
	if( WritePressure )
	{
		std::cout << "  writing Pressure........";
		Paraview_File << "\nSCALARS " << "Pressure" << " double 1\n";
		Paraview_File << "LOOKUP_TABLE default\n";
		WriteScalar(config_container,
								geometry_container,
								solver_container,
								Paraview_File,
								VariablePressure);
		std::cout << " Done." << std::endl;
	}

  // Register temperature data.
	if( WriteTemperature )
	{
		std::cout << "  writing Temperature.....";
		Paraview_File << "\nSCALARS " << "Temperature" << " double 1\n";
		Paraview_File << "LOOKUP_TABLE default\n";
		WriteScalar(config_container,
								geometry_container,
								solver_container,
								Paraview_File,
								VariableTemperature);
		std::cout << " Done." << std::endl;
	}

  // Register Mach number data.
	if( WriteMach )
	{
		std::cout << "  writing Mach............";
		Paraview_File << "\nSCALARS " << "Mach" << " double 1\n";
		Paraview_File << "LOOKUP_TABLE default\n";
		WriteScalar(config_container,
								geometry_container,
								solver_container,
								Paraview_File,
								VariableMachNumber);
		std::cout << " Done." << std::endl;
	}

  // Register vorticity data.
	if( WriteVorticity )
	{
		std::cout << "  writing Vorticity.......";
		Paraview_File << "\nSCALARS " << "Vorticity" << " double 1\n";
		Paraview_File << "LOOKUP_TABLE default\n";
		WriteScalarVorticity(config_container,
								         geometry_container,
									       element_container,
								         solver_container,
								         Paraview_File);
		std::cout << " Done." << std::endl;
	}

	// Register entropy data.
	if( WriteEntropy )
	{
		std::cout << "  writing entropy.........";
		Paraview_File << "\nSCALARS " << "Entropy" << " double 1\n";
		Paraview_File << "LOOKUP_TABLE default\n";
		WriteScalar(config_container,
								geometry_container,
								solver_container,
								Paraview_File,
								VariableEntropy);
		std::cout << " Done." << std::endl;
	}

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
}


void CASCIIFileVTK::WriteScalar
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File,
 CScalarVariable *Variable
)
 /*
	* Function that writes a scalar parameter in ASCII format.
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


void CASCIIFileVTK::WriteVector
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File,
 CVectorVariable *Variable
)
 /*
	* Function that writes a vector parameter in ASCII format.
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


void CASCIIFileVTK::WriteScalarVorticity
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CElement       **element_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File
)
 /*
	* Function that writes a vorticity parameter in ASCII format.
	*/
{
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);

		unsigned long  nElem 		= gridZone->GetnElem();
		unsigned short nPoly 	  = gridZone->GetnPolySol();
		unsigned short nSubElem = nPoly*nPoly;
		unsigned short nNode    = gridZone->GetnDOFsSol2D();

		//Obtain number of solution DOFs in 1D.
		const unsigned short nDOFsSol1D = element_container[iZone]->GetnDOFsSol1D();
		// Obtain the 1D Lagrange differential matrix on the solution DOFs.
		auto* dell = element_container[iZone]->GetDerLagrangeDOFsSol1D();

		for(unsigned long iElem=0; iElem<nElem; iElem++){

			// Extract element solution.
			const auto& variables = solver_container[iZone]->GetDataContainer(iElem)->GetDataDOFsSol();
			
			// Extract current element size.
    	auto& ElemSize = gridZone->GetGeometryElem(iElem)->GetElemSize();

			// Jacobian in the x and y direction.
			const as3double drdx = 2.0/ElemSize[0];
			const as3double dsdy = 2.0/ElemSize[1];

			// Reserve memory for the required derivatives.
			as3vector1d<as3double> dudy(nNode, 0.0);
			as3vector1d<as3double> dvdx(nNode, 0.0);
			
			// Compute the Cartesian gradient of the velocity.
			for(unsigned short lj=0; lj<nDOFsSol1D; lj++)
			{
				// Starting r-index, based on level of s.
				unsigned short J0 = lj*nDOFsSol1D; 
				for(unsigned short li=0; li<nDOFsSol1D; li++)
				{
					// Index of the derivative in r.
					unsigned short I0 = li*nDOFsSol1D;
					// Current index in 2D.
					unsigned short II = J0+li;

					// Compute the parametric derivative dvdx in r-direction.
#pragma omp simd
					for(unsigned short kk=0; kk<nDOFsSol1D; kk++)
						dvdx[II] += ( variables[2][J0+kk]/variables[0][J0+kk] )*dell[I0+kk];

					// Compute the parametric derivative of dudy in s-direction.
#pragma omp simd
					for(unsigned short kk=0; kk<nDOFsSol1D; kk++)
						dudy[II] += ( variables[1][kk*nDOFsSol1D+li]/variables[0][kk*nDOFsSol1D+li] )*dell[J0+kk];
					
					// Convert the parametric derivatives into Cartesian ones.
					dvdx[II] *= drdx;
					dudy[II] *= dsdy;
				}
			}

			short idxLocal = 0;
			for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
				for(unsigned short iNode=0; iNode<N_POINTS_QUADRILATERAL; iNode++){
					// Compute vorticity.
					const as3double scalar = dvdx[ConnLocal[iZone][idxLocal]] 
						                     - dudy[ConnLocal[iZone][idxLocal]];
					// Write scalar to file.
					Paraview_File << std::scientific << scalar << "\t";
					idxLocal++;
				}
				Paraview_File << "\n";
			}
		}
	}
}


void CASCIIFileVTK::WriteScalarAuxPML
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File,
 unsigned short   iVar
)
 /*
	* Function that writes a PML auxiliary scalar parameter in ASCII format.
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


void CASCIIFileVTK::WriteVectorAuxPML
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		  **solver_container,
 std::ofstream   &Paraview_File,
 unsigned short   iVar
)
 /*
	* Function that writes a PML auxiliary vector parameter in ASCII format.
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



CBinaryFileVTK::CBinaryFileVTK
(
 CConfig                    *config_container,
 CGeometry                  *geometry_container,
 as3vector2d<unsigned long> &MapGlobalToLocal
)
	:
		CFileVTK
		(
		 config_container,
		 geometry_container,
		 MapGlobalToLocal
		)
 /*
	* Constructor, used to initialize CBinaryFileVTK class.
	*/
{
  // Check for big endian. If not, the bytes must be swapped when
  // writing a binary VTK file.
  union {int i; char c[4];} val;
  val.i = 0x76543210;
  if (val.c[0] == 0x10) BigEndian = false;
  else                  BigEndian = true;

	// Initialize the number of elements in each zone.
	nElemZone.resize(nZone, 0);
	// Initialize the number of nodes per element in each zone.
	nNodeZone.resize(nZone, 0);
	// Initialize the number of solution polynomial in each zone.
	nPolyZone.resize(nZone, 0);
	// Initialize the total DOFs per each zone.
	nDOFsTotZone.resize(nZone, 0);
	// Initialize the size of each sub-element of in each zone.
	nSubElemZone.resize(nZone, 0);
	// Obtain the actual information in each zone.
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		nElemZone[iZone] = geometry_container->GetGeometryZone(iZone)->GetnElem();
		nNodeZone[iZone] = geometry_container->GetGeometryZone(iZone)->GetnDOFsSol2D();
		nPolyZone[iZone] = geometry_container->GetGeometryZone(iZone)->GetnPolySol();
	
		// Note, ParaView expects the values in int.

		// Deduce the total number of DOFs in each zone. 
		nDOFsTotZone[iZone] = nElemZone[iZone]*nNodeZone[iZone];
		// Deduce the total number of sub-elements in each element of each zone.
		nSubElemZone[iZone] = nPolyZone[iZone]*nPolyZone[iZone];
	}

	// Initialize the number of total DOFs written in each previous zone.
	// This is needed in the connectivity array for multi-zone output.
	nDOFsTotWritten.resize(nZone, 0);
	// Loop over the previous zones and accumulate the data. 
	// Notice the index starts from iZone=1.
	for(unsigned short iZone=1; iZone<nZone; iZone++)
		nDOFsTotWritten[iZone] = nDOFsTotZone[iZone-1] + nDOFsTotWritten[iZone-1];

	// Retrieve the total number of elements in x and y directions, per zone.
	nxElemZone = config_container->GetnxElemZone();
	nyElemZone = config_container->GetnyElemZone();

	// Book-keep the address of the index mapper across all zones and elements.
	mMapGlobalToLocal = MapGlobalToLocal;

	// Initialize and define the variable names written.
	if( WriteDensity )
	{
		VariableNames.push_back("Density");
	}
	if( WriteMomentum )
	{
		VariableNames.push_back("Momentum_x");
		VariableNames.push_back("Momentum_y");
		VariableNames.push_back("Momentum_z");
	}
	if( WriteTotalEnergy )
	{
		VariableNames.push_back("Energy");
	}
	if( WritePressure )
	{
		VariableNames.push_back("Pressure");
	}
	if( WriteVelocity )
	{
		VariableNames.push_back("Velocity_x");
		VariableNames.push_back("Velocity_y");
		VariableNames.push_back("Velocity_z");
	}
	if( WriteTemperature )
	{
		VariableNames.push_back("Temperature");
	}
	if( WriteMach )
	{
		VariableNames.push_back("Mach");
	}
	if( WriteVorticity )
	{
		VariableNames.push_back("Vorticity");
	}
	if( WriteEntropy )
	{
		VariableNames.push_back("Entropy");
	}
}


CBinaryFileVTK::~CBinaryFileVTK
(
 void
)
 /*
	* Destructor for CBinaryFileVTK, frees allocated memory. 
	*/
{

}


void CBinaryFileVTK::WriteFileVTK
(
 CConfig             *config_container,
 CGeometry           *geometry_container,
 CElement           **element_container,
 CSolver            **solver_container,
 const unsigned long  FileNumberVTK
)
 /*
	* Function that writes a binary format of a VTK data file.
	*/
{
 	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution in binary VTK format... ";

	// Open string stream.
	std::ostringstream fn;
	fn << OutputVTKFilename << "_" << FileNumberVTK << ".vtk";


	// Define the maximum string length for the writing.
  const int MAX_STRING_LENGTH = 255;

	// Initialize the buffer data.
	as3vector2d<float> varbuf(nZone);
	as3vector2d<float> coorbuf(nZone);
	as3vector2d<int>   connbuf(nZone);
	as3vector2d<int>   elemtypebuf(nZone);
	

	// Assemble the data to write.
	for(unsigned short iZone=0; iZone<nZone; iZone++){
		
		// Extract the current zone information.
		const int nElem = nElemZone[iZone];
		const int nNode = nNodeZone[iZone];
		const int nPoly = nPolyZone[iZone];
		
		// Total DOFs and sub-elements in this zone.
		const int nDOFsTot = nDOFsTotZone[iZone];
		const int nSubElem = nSubElemZone[iZone];

		// Extract number of elements in x and y in this zone.
		const int nxElem = nxElemZone[iZone];
		const int nyElem = nyElemZone[iZone];

		// Deduce number of nodes for the DOFs in 1D.
		const int nDOFsSol1D = nPoly + 1;

		// Obtain the current element's geometry.
		const CGeometryZone *grid = geometry_container->GetGeometryZone(iZone);

		// Initialize and specify the coordinates. Note, ParaView expects them in 3D.
		coorbuf[iZone].resize( 3*nDOFsTot );
		// Initialize the global connectivity matrix.
		// Note, the addition of one is for the number of points per element.
		connbuf[iZone].resize( nElem*nSubElem*(N_POINTS_QUADRILATERAL+1) );
		// Initialize the element types. For this code, just specify a quadrilateral.
		elemtypebuf[iZone].resize( nElem*nSubElem, QUADRILATERAL );

		// Initialize the buffer data.
  	for(int ijElem=0; ijElem<nElem; ++ijElem)
  	{
			// Deduce the local element indices in x and y.
			const int jElem = ijElem/nxElem;
			const int iElem = ijElem - jElem*nxElem;


			// Deduce the current element coordinate array location.
			float *coorelem = coorbuf[iZone].data() + ijElem*3*nNode;  
			// Deduce the current element connectivity array location.
			int *connelem   = connbuf[iZone].data() + ijElem*(N_POINTS_QUADRILATERAL+1)*nSubElem;


			// Loop over every DOF on this element and store the coordinate.
			int idx = 0;
			for(int l=0; l<nNode; l++){
				coorelem[idx++] = (float) grid->GetGeometryElem(ijElem)->GetCoordSolDOFs()[0][l];
				coorelem[idx++] = (float) grid->GetGeometryElem(ijElem)->GetCoordSolDOFs()[1][l];
				coorelem[idx++] = (float) 0.0;
			}

			// Deduce the start of the current element's connectivity array.
			const int offset = ijElem*nNode + nDOFsTotWritten[iZone];
			
			// Connectivity array.
			idx = 0;
			for(int j=0; j<nPoly; j++){
				const int joff = offset + j*nDOFsSol1D;
				for(int i=0; i<nPoly; i++){

					// Starting index of each nPoly=1 sub-element.
					const int n0 = joff + i; 
					// Number of points in this element.
					connelem[idx++] = N_POINTS_QUADRILATERAL;
					// Actual connectivity array.
					connelem[idx++] = n0;
					connelem[idx++] = n0 + 1;
					connelem[idx++] = n0 + 1 + nDOFsSol1D;
					connelem[idx++] = n0 + nDOFsSol1D;
				}
			}
		}


		// Allocate the size of the written data variables.
		varbuf[iZone].resize(nDOFsTot*VariableNames.size(), 0.0);
	}

	// Compute and store the required data for visualizatoin. 
	// Note, this is done over all elements in all zones simulataneoously.
	DetermineVisualizationData(geometry_container, 
			                       element_container, 
														 solver_container, 
														 VariableNames, varbuf); 
	
	// Check if there need be any swapping, since ParaView expects data in big endian format.
	if( !BigEndian ){
		for(unsigned short iZone=0; iZone<nZone; iZone++){
			SwapBytes( coorbuf[iZone].data(),     sizeof(float),   coorbuf[iZone].size() );
			SwapBytes( varbuf[iZone].data(),      sizeof(float),    varbuf[iZone].size() );
			SwapBytes( connbuf[iZone].data(),     sizeof(int),     connbuf[iZone].size() );
			SwapBytes( elemtypebuf[iZone].data(), sizeof(int), elemtypebuf[iZone].size() );
		}
	}


	// Determine the total number of DOFs in all zones combined.
	int nPoints = 0;
	for(int iZone=0; iZone<nZone; iZone++) 
		nPoints += nDOFsTotZone[iZone];

	// Determine total number of nPoly=1 elements in all the zones combined.
	int nSubElemTot = 0;
	for(int iZone=0; iZone<nZone; iZone++) 
		nSubElemTot += nElemZone[iZone]*nSubElemZone[iZone]; 

	// Open the visualization file for binary writing.
	FILE *fh = std::fopen( fn.str().c_str(), "wb" );
	if( !fh ) 
		Terminate("CBinaryFileVTK::WriteFileVTK", __FILE__, __LINE__,
				      "Visualization file in binary format could not be openned.");


	// Write the header of the visualization file.
  char str_buf[MAX_STRING_LENGTH];
  std::strcpy(str_buf, "# vtk DataFile Version 3.0\n");
	std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "vtk output\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "BINARY\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);


  // Write the coordinates.
  std::sprintf(str_buf, "POINTS %i float\n", nPoints);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  for(int iZone=0; iZone<nZone; iZone++) 
		std::fwrite(coorbuf[iZone].data(), sizeof(float), coorbuf[iZone].size(), fh);


  // Write the connectivity data.
  std::sprintf(str_buf, "\nCELLS %i %i\n", nSubElemTot, (N_POINTS_QUADRILATERAL+1)*nSubElemTot);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  for(int iZone=0; iZone<nZone; iZone++)
		std::fwrite(connbuf[iZone].data(), sizeof(int), connbuf[iZone].size(), fh);


  // Write the element type data.
  std::sprintf(str_buf, "\nCELL_TYPES %i\n", nSubElemTot);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  for(int iZone=0; iZone<nZone; iZone++)
		std::fwrite(elemtypebuf[iZone].data(), sizeof(int), elemtypebuf[iZone].size(), fh);


  // Write the ASCII line for the point data.
  std::sprintf(str_buf, "\nPOINT_DATA %i\n", nPoints);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);


	// Loop over each variables and write output to file.
	for(int iVar=0; iVar<VariableNames.size(); iVar++)
	{

		// Copy variable name, since it will get modified.
		std::string varname = VariableNames[iVar];

		// Check whether this is a scalar or vector variable.
    // Note that the y- and z-components of a vector are skipped.
    bool writevar = true, isvector = false;
    size_t found = varname.find("_x");
    if(found != std::string::npos) isvector = true;

    found = varname.find("_y");
    if(found != std::string::npos) writevar = false;

    found = varname.find("_z");
    if(found != std::string::npos) writevar = false;

    // Check for a vector field.
    if( isvector )
    {
      // Vector variable. Remove the trailing "_x".
      varname.erase(varname.end()-2, varname.end());

      // Write the ASCII line with the information.
      std::sprintf(str_buf, "\nVECTORS %s float\n", varname.c_str());
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      // Write the vector data.
			for(int iZone=0; iZone<nZone; iZone++)
      	std::fwrite(varbuf[iZone].data() + iVar*nDOFsTotZone[iZone], 
						sizeof(float), 3*nDOFsTotZone[iZone], fh);
    }
    else if( writevar )
    {
      // Scalar variable. Write the ASCII lines with the information.
      std::sprintf(str_buf, "\nSCALARS %s float 1\n", varname.c_str());
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      std::sprintf(str_buf, "LOOKUP_TABLE default\n");
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      // Write the scalar data.
			for(int iZone=0; iZone<nZone; iZone++)
      	std::fwrite(varbuf[iZone].data() + iVar*nDOFsTotZone[iZone], 
						sizeof(float), nDOFsTotZone[iZone], fh);
    }
  }

  // Close the file again.
  std::fclose(fh);

	// Report: all is complete!
	std::cout << "Done." << std::endl;
}

	
void CBinaryFileVTK::DetermineVisualizationData
(
 CGeometry                 *geometry_container,
 CElement                 **element_container,
 CSolver                  **solver_container,
 as3vector1d<std::string>  &varnames,
 as3vector2d<float>        &varbuf
)
 /*
	* Function that computes the required data for visualization in binary format.
	*/
{
	// Abbreviation involving gamma.
	const as3double gm1 = GAMMA_MINUS_ONE;

	// Deduce the total number of elements across all zones.
	const unsigned long nElemTotal = mMapGlobalToLocal.size();

	// For parallelization efficiency reasons, the computation of the variables are
	// done by joining the zones and elements. This yields an efficient distribution
	// of work over the threads.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(unsigned long iii=0; iii<nElemTotal; iii++){

		// Extract local zone number.
		const int iZone    = (int) mMapGlobalToLocal[iii][0];
		// Extract local zone element index.
		const int ijElem   = (int) mMapGlobalToLocal[iii][1];

		// Extract the current zone information.
		const int nNode    = (int) nNodeZone[iZone];
		// Total DOFs and sub-elements in this zone.
		const int nDOFsTot = (int) nDOFsTotZone[iZone];
		// Extract number of elements in x-direction in this zone.
		const int nxElem   = (int) nxElemZone[iZone];

		// Deduce the local element indices in x and y.
		const int jElem = ijElem/nxElem;
		const int iElem = ijElem - jElem*nxElem;
	
		// Deduce the start of the current element's connectivity array.
		const int offset = ijElem*nNode;
	
		// Extract current element's data solution.
		auto& sol = solver_container[iZone]->GetDataContainer(ijElem)->GetDataDOFsSol();
	
		// Allocate data for the primitive variables on the solution DOFs for this element.
		// These contain: 1/rho, u, v, p. 
		// Note, the inverse of rho is computed for convenience.
		as3vector2d<as3double> primvar( nVar, as3vector1d<as3double>(nNode, 0.0) );
	
		// Loop over the element DOFs and compute the primitive variables needed.
		for(int l=0; l<nNode; l++){
			const as3double ovrho =   1.0/sol[0][l];                       
			const as3double u     = ovrho*sol[1][l];                     
			const as3double v     = ovrho*sol[2][l];                     
			const as3double p     = gm1*( sol[3][l] 
					                  -   0.5*sol[0][l]*(u*u + v*v) );
	
			// Store the values.
			primvar[0][l] = ovrho;
			primvar[1][l] = u;
			primvar[2][l] = v;
			primvar[3][l] = p;
		}
	
	
		// Determine which variables to compute.
		for(int iVar=0; iVar<varnames.size(); iVar++)
		{
	
			if( varnames[iVar] == "Density" )
			{
				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
			
				// Compute and store the density in each DOF.
				for(int l=0; l<nNode; l++){
					buf[l] = (float) sol[0][l]; 
				}
	
			}
			else if( varnames[iVar] == "Momentum_x" )
			{
				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode*3;
	
				// Compute and store the momentum variables.
				for(int l=0; l<nNode; l++){
					const int l3 = 3*l;
					buf[l3  ] = (float) sol[1][l];
					buf[l3+1] = (float) sol[2][l];
					buf[l3+2] = (float) 0.0;
				}
	
				// Increment the number of variables by 2, to offset the other 2 components
				// of this vector, which are already computed and stored here.
				iVar += 2;
			}
			else if( varnames[iVar] == "Energy" )
			{
				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
	
				// Compute and store the energy variable.
				for(int l=0; l<nNode; l++){
					buf[l] = (float) sol[3][l];
				}
			}
			else if( varnames[iVar] == "Pressure" )
			{
				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
				
				// Compute and store the pressure variable.
				for(int l=0; l<nNode; l++){
					buf[l] = (float) primvar[3][l];
				}
			}
			else if( varnames[iVar] == "Velocity_x" )
			{
				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode*3;
	
				// Compute and store the momentum variables.
				for(int l=0; l<nNode; l++){
					const int l3 = 3*l;
					const as3double ovrho = 1.0/sol[0][l];
					buf[l3  ] = (float) ovrho*sol[1][l];
					buf[l3+1] = (float) ovrho*sol[2][l];
					buf[l3+2] = (float) 0.0;
				}
	
				// Increment the number of variables by 2, to offset the other 2 components
				// of this vector, which are already computed and stored here.
				iVar += 2;
			}
			else if( varnames[iVar] == "Temperature" )
			{
				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
	
				// Abbreviation for 1/R.
				const as3double ovR = 1.0/GAS_CONSTANT;
	
				// Compute and store the temperature variable.
				for(int l=0; l<nNode; l++){
					buf[l] = (float) ovR*primvar[0][l]*primvar[3][l];
				}
			}
			else if( varnames[iVar] == "Mach" )
			{
				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
	
				// Compute and store the Mach number variable.
				for(int l=0; l<nNode; l++){
					const as3double a2    = GAMMA*primvar[0][l]*primvar[3][l];
					const as3double umag2 = primvar[1][l]*primvar[1][l] 
						                    + primvar[2][l]*primvar[2][l]; 
					buf[l] = (float) sqrt( umag2/a2 ); 
				}
			}
			//
			else if( varnames[iVar] == "Entropy" )
			{
				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
	
				// Compute and store the Mach number variable.
				for(int l=0; l<nNode; l++){
					buf[l] = (float) logf( primvar[3][l]*pow( primvar[0][l], GAMMA ) ); 
				}
			}		
			//
			else if( varnames[iVar] == "Vorticity" )
			{
				//Obtain number of solution DOFs in 1D.
				const unsigned short nDOFsSol1D = element_container[iZone]->GetnDOFsSol1D();
				// Obtain the 1D Lagrange differential matrix on the solution DOFs.
				auto* dell = element_container[iZone]->GetDerLagrangeDOFsSol1D();

				// Extract current element size.
    		auto& ElemSize = geometry_container->GetGeometryZone(iZone)->GetGeometryElem(iElem)->GetElemSize();

				// Jacobian in the x and y direction.
				const as3double drdx = 2.0/ElemSize[0];
				const as3double dsdy = 2.0/ElemSize[1];

				// Reserve memory for the required derivatives.
				as3vector1d<as3double> dudy(nNode, 0.0);
				as3vector1d<as3double> dvdx(nNode, 0.0);

				// Compute the Cartesian gradient of the velocity.
				for(unsigned short lj=0; lj<nDOFsSol1D; lj++)
				{
					// Starting r-index, based on level of s.
					unsigned short J0 = lj*nDOFsSol1D; 
					for(unsigned short li=0; li<nDOFsSol1D; li++)
					{
						// Index of the derivative in r.
						unsigned short I0 = li*nDOFsSol1D;
						// Current index in 2D.
						unsigned short II = J0+li;

						// Compute the parametric derivative dvdx in r-direction.
						for(unsigned short kk=0; kk<nDOFsSol1D; kk++)
							dvdx[II] += primvar[2][J0+kk]*dell[I0+kk];

						// Compute the parametric derivative of dudy in s-direction.
						for(unsigned short kk=0; kk<nDOFsSol1D; kk++)
							dudy[II] += primvar[1][kk*nDOFsSol1D+li]*dell[J0+kk];
						
						// Convert the parametric derivatives into Cartesian ones.
						dvdx[II] *= drdx;
						dudy[II] *= dsdy;
					}
				}

				// Set the pointer to the correct location.
				float *buf = varbuf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;

				// Compute and store the Vorticity variable.
				for(int l=0; l<nNode; l++){
					buf[l] = (float) ( dvdx[l] - dudy[l] ); 
				}
			}
			else
			{
#ifdef HAVE_OPENMP
				if(omp_get_thread_num() == 0)
#endif
				Terminate("CBinaryFileVTK::DetermineVisualizationData", __FILE__, __LINE__,
						      "Unknown visualization variable name encountered.");
			}
		}
	} // End of OpenMP loop.
}


