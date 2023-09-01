#include "import_structure.hpp"





CImport::CImport
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor, used to initialize CImport class. 
	*/
{
	// Extract number of zones.
	nZone = config_container->GetnZone();
}


CImport::~CImport
(
 void
)
 /*
	* Destructor for CImport, frees allocated memory.
	*/
{

}




CASCIIImport::CASCIIImport
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
	:
		CImport
		(
		 config_container,
		 geometry_container
		)
 /*
	* Constructor, used to initialize CASCIIImport class. 
	*/
{

}


CASCIIImport::~CASCIIImport
(
 void
)
 /*
	* Destructor for CASCIIImport, frees allocated memory.
	*/
{

}


void CASCIIImport::ReadSolutionRestartFile
(
  CConfig    *config_container,
  CGeometry  *geometry_container,
  CElement  **element_container,
  CSolver   **solver_container,
  as3double  &SimTime
)
 /*
  * Function that reads the solution from a restart file in ASCII format.
  */
{
  // Report output.
  std::cout << "\n  reading ASCII data from restart file... ";

  // Open restart file.
  std::ifstream Restart_File( config_container->GetRestartFilename() );

  // Check if file exists.
  if(!Restart_File.is_open())
		Terminate("CASCIIImport::ReadSolutionRestartFile", __FILE__, __LINE__,
							"ASCII restart file could not be opened!");

  // Create temporary input data.
  as3double input_XMin, input_XMax, input_YMin, input_YMax;
  unsigned short input_nZone;

  // Read header information.
  AddScalarOption(Restart_File, "SimTime", SimTime);
  AddScalarOption(Restart_File, "nZone",   input_nZone);
  AddScalarOption(Restart_File, "XMin",    input_XMin);
  AddScalarOption(Restart_File, "XMax",    input_XMax);
  AddScalarOption(Restart_File, "YMin",    input_YMin);
  AddScalarOption(Restart_File, "YMax",    input_YMax);

  // Extract main domain bound.
  auto& DomainBound = config_container->GetDomainBound();

  // Consistency check.
  assert( input_nZone == nZone );
  assert( input_XMin  == DomainBound[0] );
  assert( input_XMax  == DomainBound[1] );
  assert( input_YMin  == DomainBound[2] );
  assert( input_YMax  == DomainBound[3] );

  // Read data in every zone.
  for(unsigned short iZone=0; iZone<nZone; iZone++){

    // Create temporary input data.
    unsigned long  input_nElem, input_nxElem, input_nyElem, input_iElem;
    unsigned short input_iZone, input_TypeSolver, input_TypeZone,
                   input_TypeDOFs, input_nPolySol, input_nDOFsSol1D,
                   input_nDOFsInt1D, input_TypeIC, input_TypeBufferLayer;

    // Temporary storage of the solution being read in.
    as3data1d<as3double> storage(nVar, nullptr);

    // Extract data container in this zone.
    auto& data_container = solver_container[iZone]->GetDataContainer();

    // Boolean to indicate whether or not an interpolation is needed. Default, false.
    bool Interpolate = false;

    // Read general information.
    AddScalarOption(Restart_File, "iZone",            input_iZone);
    AddScalarOption(Restart_File, "TypeIC",           input_TypeIC);
    AddScalarOption(Restart_File, "TypeSolver",       input_TypeSolver);
    AddScalarOption(Restart_File, "TypeBufferLayer",  input_TypeBufferLayer);
    AddScalarOption(Restart_File, "TypeZone",         input_TypeZone);
    AddScalarOption(Restart_File, "TypeDOFs",         input_TypeDOFs);
    AddScalarOption(Restart_File, "nPolySol",         input_nPolySol);
    AddScalarOption(Restart_File, "nDOFsSol1D",       input_nDOFsSol1D);
    AddScalarOption(Restart_File, "nDOFsInt1D",       input_nDOFsInt1D);
    AddScalarOption(Restart_File, "nElem",            input_nElem);
    AddScalarOption(Restart_File, "nxElem",           input_nxElem);
    AddScalarOption(Restart_File, "nyElem",           input_nyElem);

    // Extract current data from config file.
    unsigned short TypeIC          = config_container->GetTypeIC(iZone);
    unsigned short TypeSolver      = config_container->GetTypeSolver(iZone);
    unsigned short TypeBufferLayer = config_container->GetTypeBufferLayer(iZone);
    unsigned short TypeZone        = config_container->GetTypeZone(iZone);
    unsigned short TypeDOFs        = config_container->GetTypeDOFs(iZone);
    unsigned long  nxElem          = config_container->GetnxElemZone(iZone);
    unsigned long  nyElem          = config_container->GetnyElemZone(iZone);
    unsigned long  nElem           = solver_container[iZone]->GetnElem();
    unsigned short nPolySol        = element_container[iZone]->GetnPolySol();
    unsigned short nDOFsSol1D      = element_container[iZone]->GetnDOFsSol1D();

    // Consistency check.
    assert( iZone           == input_iZone );
    assert( TypeIC          == input_TypeIC );
    assert( TypeSolver      == input_TypeSolver );
    assert( TypeBufferLayer == input_TypeBufferLayer );
    assert( TypeZone        == input_TypeZone );
    assert( nElem           == input_nElem );
    assert( nxElem          == input_nxElem );
    assert( nyElem          == input_nyElem );


   	// Interpolating polynomial that converts from the restart solution to the
    // current specified/expected solution polynomial.
    as3vector1d<as3double> lagrange1DTranspose;

    // Check if there need be no interpolation.
    if( (nPolySol != input_nPolySol) || (TypeDOFs != input_TypeDOFs) ){

      // Interpolation is required.
      Interpolate = true;

      // Extract basis of current data.
      auto& rDOFsSol1D = element_container[iZone]->GetrDOFsSol1D();

      // Reserve memory for the basis used in the restart file.
      as3vector1d<as3double> rDOFsSolInput1D(input_nDOFsSol1D);

      // Create basis for the solution from restart file.
      element_container[iZone]->LocationDOFs1D(input_TypeDOFs, rDOFsSolInput1D);

      // Reserve (transposed) interpolation polynomial between the restart file and current data.
      lagrange1DTranspose.resize(input_nDOFsSol1D*nDOFsSol1D);

      // Compute lagrange polynomial that interpolates from the nDOFsSol1D of the
      // restart solution to this current zone.
      element_container[iZone]->LagrangeBasisFunctions(rDOFsSolInput1D, rDOFsSol1D,
                                                       nullptr, nullptr,
                                                       lagrange1DTranspose.data(),
                                                       nullptr, nullptr, nullptr);
    }

    // Deduce the number of solution DOFs in 2D.
    const unsigned short input_nDOFsSol2D = input_nDOFsSol1D*input_nDOFsSol1D;

    // Reserve memory for the intermediary solution storage data. Note, this does not
    // matter whether this is an interpolation procedure or not, since this is the input data.
    for(unsigned short iVar=0; iVar<nVar; iVar++)
      storage[iVar] = new as3double[input_nDOFsSol2D]();

    // Loop over all elements and read the data.
    for(unsigned long iElem=0; iElem<nElem; iElem++){

      // Extract current data.
      auto& data = data_container[iElem]->GetDataDOFsSol();
      // Find the starting point for the data.
      AddScalarOption(Restart_File, "DataDOFsSol", input_iElem);
      // Consistency check.
      assert( iElem == input_iElem );

      // Temporary storage.
      std::string line;

      // Update line read.
      getline(Restart_File, line);
      std::stringstream ss(line);

      // Read and copy the data in this element.
      for(unsigned short iVar=0; iVar<nVar; iVar++){
        for(unsigned short iNode=0; iNode<input_nDOFsSol2D; iNode++){

          // Read new value and check the buffer is good.
          std::string tmp; getline(ss, tmp, ','); if( !ss.good() ) break;

          // Convert the string value read into a double.
          std::stringstream convertor(tmp); convertor >> storage[iVar][iNode];

        }
      }

      // Check whether an interpolation is needed or not. If not, just copy the data.
      if( Interpolate )
        // Interpolate solution.
        TensorProductSolAndGradVolume(nDOFsSol1D, nVar, input_nDOFsSol1D,
                                      lagrange1DTranspose.data(),
                                      nullptr, storage.data(),
                                      data.data(),
                                      nullptr, nullptr);
      else
        // Copy solution.
        for(unsigned short i=0; i<nVar; i++)
          for(unsigned short l=0; l<input_nDOFsSol2D; l++)
            data[i][l] = storage[i][l];

    }

    // If this is a PML zone, read the auxiliary data too.
    if( TypeBufferLayer == PML_LAYER ){

      // Loop over all elements and read the data.
      for(unsigned long iElem=0; iElem<nElem; iElem++){

        // Extract current data.
        auto** data = &data_container[iElem]->GetDataDOFsSol()[nVar];
        // Find the starting point for the data.
        AddScalarOption(Restart_File, "DataDOFsSolAux", input_iElem);
        // Consistency check.
        assert( iElem == input_iElem );

        // Temporary storage.
        std::string line;

        // Update line read.
        getline(Restart_File, line);
        std::stringstream ss(line);

        // Read and copy the data in this element.
        for(unsigned short iVar=0; iVar<nVar; iVar++){
          for(unsigned short iNode=0; iNode<input_nDOFsSol2D; iNode++){

            // Read new value and check the buffer is good.
            std::string tmp; getline(ss, tmp, ','); if( !ss.good() ) break;

            // Convert the string value read into a double.
            std::stringstream convertor(tmp); convertor >> storage[iVar][iNode];

          }
        }

        // Check whether an interpolation is needed or not. If not, just copy the data.
        if( Interpolate )
          // Interpolate solution.
          TensorProductSolAndGradVolume(nDOFsSol1D, nVar, input_nDOFsSol1D,
                                        lagrange1DTranspose.data(),
                                        nullptr, storage.data(),
                                        data,
                                        nullptr, nullptr);
        else
          // Copy solution.
          for(unsigned short i=0; i<nVar; i++)
            for(unsigned short l=0; l<input_nDOFsSol2D; l++)
              data[i][l] = storage[i][l];

      }
    }

    // Peek to next character (used for determining end of file assertion).
		Restart_File.peek();

    // delete the temporary storage used.
    for(unsigned short i=0; i<storage.size(); i++)
      if( storage[i] ) delete [] storage[i];
  }

  // Make sure end of file is reached.
	assert( Restart_File.eof() == true );

  // Close file.
  Restart_File.close();

  // Report progress.
	std::cout << "Done." << std::endl;
}




CBinaryImport::CBinaryImport
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
	:
		CImport
		(
		 config_container,
		 geometry_container
		)
 /*
	* Constructor, used to initialize CBinaryImport class. 
	*/
{

}


CBinaryImport::~CBinaryImport
(
 void
)
 /*
	* Destructor for CBinaryImport, frees allocated memory.
	*/
{

}


void CBinaryImport::ReadSolutionRestartFile
(
  CConfig    *config_container,
  CGeometry  *geometry_container,
  CElement  **element_container,
  CSolver   **solver_container,
  as3double  &SimTime
)
 /*
  * Function that reads the solution from a restart file in binary format.
  */
{
  // Report output.
  std::cout << "\n  reading binary data from restart file... ";

	// Open the file for binary reading.
	FILE *fh = std::fopen(config_container->GetRestartFilename(), "rb");

	// Check if the file can be opened. 
	if( !fh ){
		// Log message.
		std::ostringstream message;
		message << "Could not open file: \n" 
			      << "'" 
						<< config_container->GetRestartFilename() 
						<< "'";
		// Report error and exit.
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__, message.str());
	}

	// Initialize header information.
	int header[1];

	// Open and read the first integer.
	if( std::fread(&header, sizeof(int), 1, fh) != 1 ) 
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
				       "File header could not be read.");

	// Check if byte swapping must be applied.
	bool byteswap = false;
	if( header[0] != AS3_MagicNumber ){
		byteswap = true;
		SwapBytes(header, sizeof(int), 1);
	}

	// Check if this indeed is an AS3 binary file.
	if( header[0] != AS3_MagicNumber ) 
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
				      "Imported file is not an AS3 file.");


  // Create temporary input data.
  as3double input_XMin, input_XMax, input_YMin, input_YMax;
  unsigned short input_nZone;

  // Extract main domain bound.
  auto& DomainBound = config_container->GetDomainBound();


	// Read the physical simulation time in the current file.
	if( std::fread(&SimTime, sizeof(as3double), 1, fh) != 1 ) 
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__, 
				      "Could not read SimTime.");
	// Swap bytes, if necessary.
	if(byteswap) SwapBytes(&SimTime, sizeof(as3double), 1); 

	// Read the total number of zones.
	if( std::fread(&input_nZone, sizeof(unsigned short), 1, fh) != 1 )
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
				      "Could not read nZone.");
	// Swap bytes, if necessary.
	if(byteswap) SwapBytes(&input_nZone, sizeof(unsigned short), 1);

	// Read the physical domain bounding coordinate: xmin.
	if( std::fread(&input_XMin, sizeof(as3double), 1, fh) != 1 )
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
				      "Could not read xmin.");
	// Swap bytes, if necessary.
	if(byteswap) SwapBytes(&input_XMin, sizeof(as3double), 1);

	// Read the physical domain bounding coordinate: xmax.
	if( std::fread(&input_XMax, sizeof(as3double), 1, fh) != 1 )
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
				      "Could not read xmax.");
	// Swap bytes, if necessary.
	if(byteswap) SwapBytes(&input_XMax, sizeof(as3double), 1);

	// Read the physical domain bounding coordinate: ymin.
	if( std::fread(&input_YMin, sizeof(as3double), 1, fh) != 1 )
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
				      "Could not read ymin.");
	// Swap bytes, if necessary.
	if(byteswap) SwapBytes(&input_YMin, sizeof(as3double), 1);

	// Read the physical domain bounding coordinate: ymax.
	if( std::fread(&input_YMax, sizeof(as3double), 1, fh) != 1 )
		Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
				      "Could not read ymax.");
	// Swap bytes, if necessary.
	if(byteswap) SwapBytes(&input_YMax, sizeof(as3double), 1);


  // Consistency check.
  assert( input_nZone == nZone );
  assert( fabs( input_XMin - DomainBound[0] ) <= 1.0e-5 );
  assert( fabs( input_XMax - DomainBound[1] ) <= 1.0e-5 );
  assert( fabs( input_YMin - DomainBound[2] ) <= 1.0e-5 );
  assert( fabs( input_YMax - DomainBound[3] ) <= 1.0e-5 );


  // Read data in every zone.
  for(unsigned short iZone=0; iZone<nZone; iZone++){

    // Create temporary input data.
    unsigned long  input_nElem, input_nxElem, input_nyElem, input_iElem;
    unsigned short input_iZone, input_TypeSolver, input_TypeZone,
                   input_TypeDOFs, input_nPolySol, input_nDOFsSol1D,
                   input_nDOFsInt1D, input_TypeIC, input_TypeBufferLayer;

    // Temporary storage of the solution being read in.
    as3data1d<as3double> storage(nVar, nullptr);
 
		// Extract data container in this zone.
    auto& data_container = solver_container[iZone]->GetDataContainer();

    // Boolean to indicate whether or not an interpolation is needed. Default, false.
    bool Interpolate = false;


		// Read the current zone index.
		if( std::fread(&input_iZone, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read iZone.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_iZone, sizeof(unsigned short), 1);

		// Read the current type of initial condition specified.
		if( std::fread(&input_TypeIC, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read TypeIC.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_TypeIC, sizeof(unsigned short), 1);

		// Read the current type of solver specified.
		if( std::fread(&input_TypeSolver, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read TypeSolver.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_TypeSolver, sizeof(unsigned short), 1);

		// Read the current type of buffer layer specified.
		if( std::fread(&input_TypeBufferLayer, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read TypeBufferLayer.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_TypeBufferLayer, sizeof(unsigned short), 1);

		// Read the current type of zone specified.
		if( std::fread(&input_TypeZone, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read TypeZone.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_TypeZone, sizeof(unsigned short), 1);

		// Read the current type of DOFs specified.
		if( std::fread(&input_TypeDOFs, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read TypeDOFs.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_TypeDOFs, sizeof(unsigned short), 1);

		// Read the current polynomial order specified.
		if( std::fread(&input_nPolySol, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read nPolySol.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_nPolySol, sizeof(unsigned short), 1);

		// Read the current solution DOFs in 1D specified.
		if( std::fread(&input_nDOFsSol1D, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read nDOFsSol1D.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_nDOFsSol1D, sizeof(unsigned short), 1);

		// Read the current integration nodes in 1D specified.
		if( std::fread(&input_nDOFsInt1D, sizeof(unsigned short), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read nDOFsInt1D.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_nDOFsInt1D, sizeof(unsigned short), 1);

		// Read the current total number of elements specified.
		if( std::fread(&input_nElem, sizeof(unsigned long), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read nElem.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_nElem, sizeof(unsigned long), 1);

		// Read the current number of elements in x-direction specified.
		if( std::fread(&input_nxElem, sizeof(unsigned long), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read nxElem.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_nxElem, sizeof(unsigned long), 1);

		// Read the current number of elements in y-direction specified.
		if( std::fread(&input_nyElem, sizeof(unsigned long), 1, fh) != 1 )
			Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
					      "Could not read nyElem.");
		// Swap bytes, if necessary.
		if(byteswap) SwapBytes(&input_nyElem, sizeof(unsigned long), 1);

    // Extract current data from config file.
    unsigned short TypeIC          = config_container->GetTypeIC(iZone);
    unsigned short TypeSolver      = config_container->GetTypeSolver(iZone);
    unsigned short TypeBufferLayer = config_container->GetTypeBufferLayer(iZone);
    unsigned short TypeZone        = config_container->GetTypeZone(iZone);
    unsigned short TypeDOFs        = config_container->GetTypeDOFs(iZone);
    unsigned long  nxElem          = config_container->GetnxElemZone(iZone);
    unsigned long  nyElem          = config_container->GetnyElemZone(iZone);
    unsigned long  nElem           = solver_container[iZone]->GetnElem();
    unsigned short nPolySol        = element_container[iZone]->GetnPolySol();
    unsigned short nDOFsSol1D      = element_container[iZone]->GetnDOFsSol1D();

    // Consistency check.
    assert( iZone           == input_iZone );
    assert( TypeIC          == input_TypeIC );
    assert( TypeSolver      == input_TypeSolver );
    assert( TypeBufferLayer == input_TypeBufferLayer );
    assert( TypeZone        == input_TypeZone );
    assert( nElem           == input_nElem );
    assert( nxElem          == input_nxElem );
    assert( nyElem          == input_nyElem );
	

   	// Interpolating polynomial that converts from the restart solution to the
    // current specified/expected solution polynomial.
    as3vector1d<as3double> lagrange1DTranspose;

    // Check if there need be interpolation.
    if( (nPolySol != input_nPolySol) || (TypeDOFs != input_TypeDOFs) ){

      // Interpolation is required.
      Interpolate = true;

      // Extract basis of current data.
      auto& rDOFsSol1D = element_container[iZone]->GetrDOFsSol1D();

      // Reserve memory for the basis used in the restart file.
      as3vector1d<as3double> rDOFsSolInput1D(input_nDOFsSol1D);

      // Create basis for the solution from restart file.
      element_container[iZone]->LocationDOFs1D(input_TypeDOFs, rDOFsSolInput1D);

      // Reserve (transposed) interpolation polynomial between the restart file and current data.
      lagrange1DTranspose.resize(input_nDOFsSol1D*nDOFsSol1D);

      // Compute lagrange polynomial that interpolates from the nDOFsSol1D of the
      // restart solution to this current zone.
      element_container[iZone]->LagrangeBasisFunctions(rDOFsSolInput1D, rDOFsSol1D,
                                                       nullptr, nullptr,
                                                       lagrange1DTranspose.data(),
                                                       nullptr, nullptr, nullptr);
    }

    // Deduce the number of solution DOFs in 2D.
    const unsigned short input_nDOFsSol2D = input_nDOFsSol1D*input_nDOFsSol1D;

    // Reserve memory for the intermediary solution storage data. Note, this does not
    // matter whether this is an interpolation procedure or not, since this is the input data.
	  for(unsigned short iVar=0; iVar<nVar; iVar++)
      storage[iVar] = new as3double[input_nDOFsSol2D]();

		// Loop over all the elements and read the data in this zone.
		for(unsigned long iElem=0; iElem<nElem; iElem++){

      // Extract current solution data.
      auto& data = data_container[iElem]->GetDataDOFsSol();

			// Read the current element index.
			if( std::fread(&input_iElem, sizeof(unsigned long), 1, fh) != 1 )
				Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
						      "Could not read iElem.");
			// Swap bytes, if necessary.
			if(byteswap) SwapBytes(&input_iElem, sizeof(unsigned long), 1);

			// Consistency check.
      assert( iElem == input_iElem );

			// Read the current element solution.
			for(unsigned short iVar=0; iVar<nVar; iVar++){
				if( std::fread(storage[iVar], sizeof(as3double), input_nDOFsSol2D, fh) != input_nDOFsSol2D )
					Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
							      "Could not read physical solution.");
				// Swap bytes, if necessary.
				if(byteswap) SwapBytes(storage[iVar], sizeof(as3double), input_nDOFsSol2D);
			}

			
      // Check whether an interpolation is needed or not. If not, just copy the data.
      if( Interpolate )
			{
        // Interpolate solution.
        TensorProductSolAndGradVolume(nDOFsSol1D, nVar, input_nDOFsSol1D,
                                      lagrange1DTranspose.data(),
                                      nullptr, storage.data(),
                                      data.data(),
                                      nullptr, nullptr);
			}
      else 
			{
        // Copy solution.
        for(unsigned short i=0; i<nVar; i++)
          for(unsigned short l=0; l<input_nDOFsSol2D; l++)
            data[i][l] = storage[i][l];
			}
		}


		// If this is a PML zone, then read also the auxiliary data in it as well.
		if( TypeBufferLayer == PML_LAYER ){

      // Loop over all elements and read the data.
      for(unsigned long iElem=0; iElem<nElem; iElem++){

        // Extract current auxiliary data.
        auto** data = &data_container[iElem]->GetDataDOFsSol()[nVar];

				// Read the current element index.
				if( std::fread(&input_iElem, sizeof(unsigned long), 1, fh) != 1 )
					Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
							      "Could not read iElem.");
				// Swap bytes, if necessary.
				if(byteswap) SwapBytes(&input_iElem, sizeof(unsigned long), 1);

				// Consistency check.
      	assert( iElem == input_iElem );

				// Read the current element auxiliary solution.
				for(unsigned short iVar=0; iVar<nVar; iVar++){
					if( std::fread(storage[iVar], sizeof(as3double), input_nDOFsSol2D, fh) != input_nDOFsSol2D )
						Terminate("CBinaryImport::ReadSolutionRestartFile", __FILE__, __LINE__,
								      "Could not read auxiliary solution.");
					// Swap bytes, if necessary.
					if(byteswap) SwapBytes(storage[iVar], sizeof(as3double), input_nDOFsSol2D);
				}

			
      	// Check whether an interpolation is needed or not. If not, just copy the data.
      	if( Interpolate )
				{
      	  // Interpolate solution.
      	  TensorProductSolAndGradVolume(nDOFsSol1D, nVar, input_nDOFsSol1D,
      	                                lagrange1DTranspose.data(),
      	                                nullptr, storage.data(),
      	                                data,
      	                                nullptr, nullptr);
				}
      	else 
				{
      	  // Copy solution.
      	  for(unsigned short i=0; i<nVar; i++)
      	    for(unsigned short l=0; l<input_nDOFsSol2D; l++)
      	      data[i][l] = storage[i][l];
				}
			}
    }

    // delete the temporary storage used.
    for(unsigned short i=0; i<storage.size(); i++)
      if( storage[i] ) delete [] storage[i];
	}

	// Close the file.
	std::fclose(fh);

  // Report progress.
	std::cout << "Done." << std::endl;
}



