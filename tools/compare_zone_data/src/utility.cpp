#include "utility.hpp"


// Define the global expected variables.
unsigned long  nElemExpected;
unsigned short nNodeExpected;


as3vector2d<unsigned long> PreprocessDataDistribution
(
 const char         *dirA,
 const char         *dirB,
 const char         *filename,
 const unsigned long nData
)
 /*
	* Function which pre-processes the data and decides how to distribute the 
	* computations via different chunk sizes.
	*/
{
	// Report progress.
	std::cout << "---------------------------------------------------" << std::endl;

#ifdef HAVE_OPENMP
  // Get max number of threads specified.
  unsigned short nThreads = omp_get_max_threads();
  std::cout << "\nThis is a parallel implementation using: "
            << nThreads << " threads." << std::endl;
#else
	std::cout << "\nThis is a serial implementation." << std::endl;
#endif

	// Report progress.
	std::cout << "---------------------------------------------------" << std::endl;
	std::cout << "Pre-processing information...";

	// Check if all files exist in both directories.
	CheckFilesExist(dirA, filename, nData);
	CheckFilesExist(dirB, filename, nData);

	// Open the first file from the first directory.
	unsigned long nInfoA = DeduceHeaderInformation(dirA, filename, 0);
	// Open the first file from the second directory.
	unsigned long nInfoB = DeduceHeaderInformation(dirB, filename, 0);

	// Make sure the header information is the same in both directories.
	if( nInfoA != nInfoB ) ERROR("Inconsistent header information between directories.");

	// Make sure that the expected values match the read values -- consistency check.
	if( nInfoA != (nElemExpected*nNodeExpected) ) ERROR("Inconsistency in expected values.");

	// Compute total number of as3double-size memory per file.
	unsigned long nMemoryFile = (unsigned long) nInfoA*nVar;

	// Memory size of 1 GB in bytes.
	const unsigned long mem_1gb_byte = 1073741824;

	// Current data type.
	as3double datasize;

	// Memory size of 1 GB in as3double.
	unsigned long mem_1gb_as3 = mem_1gb_byte/sizeof(datasize);

	// Approximate number of files that can fit 1GB of memory.
	// Note, the division by two is because there are two 
	// directories extracted per time step.
	unsigned long nMaxFiles = ( mem_1gb_as3/nMemoryFile )/2;

	// Reduce the maximum allowed number of files by a factor of two (default), 
	// such that the processing memory is also factored in.
	// Note, this quantity can be adjusted. 
	nMaxFiles /= 2;

	// Initialize the Chunks vector to 1.
	as3vector2d<unsigned long> Chunks; 

	// If all the files fit into memory, then the chunk size is one.
	if( nMaxFiles >= nData )
	{
		// Initialize the chunk.
		Chunks.resize(1); Chunks[0] = {0, nData};
	
		// Report progress.
		std::cout << " Done." << std::endl;
		std::cout << "---------------------------------------------------" << std::endl;
	
		// Return the information.
		return Chunks;
	}

	// Not all data can fit into the memory requirements at once. 
	// Decide on an equal distribution.
	unsigned short nChunk = std::max( nData/nMaxFiles, (unsigned long) 2 );

	// Initialize the chunks.
	Chunks.resize(nChunk);

	// Estimate the equal chunk sizes.
	unsigned long nFileChunk = nData/nChunk;

	for(int i=0; i<nChunk; i++)
	{
		// Deduce the starting and ending indices.
		unsigned long I0 = i*nFileChunk;
		unsigned long I1 = I0 + nFileChunk;

		// Book-keep starting and ending indices on this chunk.
		Chunks[i] = { I0, I1};
	}

	// Correct for the remaining entry, just in case it is not an integer 
	// multiple of the chunk sizes.
	if( Chunks[nChunk-1][1] < nData ) Chunks[nChunk-1][1] = nData;

	// Report progress.
	std::cout << "Done." << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;

	// Return the number of chunks needed.
	return Chunks;
}


void CheckFilesExist
(
 const char   *directory,
 const char   *filename,
 unsigned long nFile
)
 /*
	* Function, which checks if all the input files exist.
	*/
{
	// Pre-assemble basename of directory and file.
	std::string basename;
	basename += directory;

	// If there is no slash in the directory as specified by the user, add one.
	if( basename.back() != '/' ) basename += "/";
	// Form rest of basename.
	basename += filename;
	basename += "_";

	// Log message.
	std::ostringstream message;

	// Loop over all files and check if they can be opened.
	for(unsigned long iFile=0; iFile<nFile; iFile++)
	{
		// Deduce the current filename.	
		std::stringstream fn; fn << basename << iFile << ".bin";

		// Check if file exists.
		std::ifstream file(fn.str());
		if( !file.good() )
		{
			message << "Could not open file: \n" 
				      << "'" << fn.str() << "'";
			// Report error and exit.
			ERROR(message.str());
		}
	}

	// Check if there are still some files available, not specified.
	std::stringstream fn; fn << basename << nFile << ".bin";
	// Check if file exists.
	std::ifstream file(fn.str());
	if( file.good() )
	{
		message << "There exists additional files, such as: \n" 
			      << "'" << fn.str() << "'";
		// Report error and exit.
		ERROR(message.str());
	}
}


unsigned long DeduceHeaderInformation
(
 const char   *directory,
 const char   *filename,
 unsigned long iFile
)
 /*
	* Function, which computes from a reference file the necessary information.
	* Moreover, it also returns the total number of information: nElem*nNode.
	*/
{
	// Initialize remaining header information.
	as3double time; 

	// Pre-assemble basename of directory and file.
	std::stringstream fn;
	fn << directory;

	// If there is no slash in the directory as specified by the user, add one.
	if( fn.str().back() != '/' ) fn << "/";
	// Form rest of filename.
	fn << filename << "_" << iFile << ".bin";

	// Open the file for binary reading.
	FILE *file = std::fopen(fn.str().c_str(), "rb");

	// Check if the file can be opened. 
	if( !file )
	{
		// Log message.
		std::ostringstream message;
		message << "Could not open file: \n" 
			      << "'" << fn.str() << "'";
		// Report error and exit.
		ERROR(message.str());
	}

	// Initialize header information.
	int header[1];
	unsigned long  nElem;
	unsigned short nNode;

	// Open and read the first integer.
	if( std::fread(&header, sizeof(int), 1, file) != 1 ) ERROR("File header could not be read.");

	// Check if byte swapping must be applied.
	bool byteswap = false;
	if( header[0] != AS3_MagicNumber )
	{
		byteswap = true;
		SwapBytes(header, sizeof(int), 1);
	}

	// Check if this indeed is an AS3 binary file.
	if( header[0] != AS3_MagicNumber ) ERROR("Imported file is not an AS3 file.");

	// Read the physical simulation time in the current file.
	if( std::fread(&time, sizeof(as3double), 1, file) != 1 ) ERROR("Could not read time.");
	if(byteswap) SwapBytes(&time, sizeof(as3double), 1); 

	// Read the total number of elements in this file.
	if( std::fread(&nElem, sizeof(unsigned long), 1, file) != 1 ) ERROR("Could not read nElem.");
	if(byteswap) SwapBytes(&nElem, sizeof(unsigned long), 1);

	// Read the total number of nodes in each element in this file.
	if( std::fread(&nNode, sizeof(unsigned short), 1, file) != 1 ) ERROR("Could not read nNode.");
	if(byteswap) SwapBytes(&nNode, sizeof(unsigned short), 1);

	// Close file.
	std::fclose(file);

	// Book-keep the expected values.
	nElemExpected = nElem;
	nNodeExpected = nNode;

	// Compute the total number of DOFs contained in all this file.
	unsigned long nInfo = (unsigned long) nElem*nNode;

	// Return total number of information needed per variable.
	return nInfo;
}




void Terminate
(
  const char        *functionName,
  const char        *fileName,
  const int          lineNumber,
  const std::string &errorMessageIn
)
 /*
  * Function which prints an error message and exits.
  */
{
  // Insert a comment character and a space on places where a newline
  // character occurs in the string, such that the error message
  // looks nicer.
  std::string  errorMessage  = errorMessageIn;
  std::string  tmpString     = errorMessage;
  std::string::size_type off = 1;

  for(;;)
  {
    std::string::size_type loc = tmpString.find("\n");
    if(loc == std::string::npos) break;
    errorMessage.insert(loc+off, "# ");
    off += loc+3;
    tmpString.erase(0,loc+1);
  }

  // Header of the error message.
	std::cout << std::endl;
	std::cout << "#" << std::endl
            << "#============================= !!! Error !!! "
            << "==============================" << std::endl
            << "#" << std::endl;

  // Write the function name, file name and line number.
  std::cout << "#* Run-time error in " << functionName << std::endl;
  std::cout << "#* File: " << fileName
            << ", Line: "  << lineNumber << std::endl
            << "#" << std::endl;

  // Write the error message and the terminating line.
  std::cout << "# " << errorMessage << std::endl
            << "#" << std::endl;
  std::cout << "#============================================"
            << "==============================" << std::endl
            << "#" << std::endl << std::flush;

  // And exit.
  exit(1);
}


void SwapBytes
(
 void   *buffer,
 size_t  nBytes,
 int     nItems
)
 /*
	* Function, which swaps bytes.
	*/
{
	// Store half the number of bytes in kk and cast the buffer
	// to a character buffer.
	char *buf       = (char *) buffer;
	const size_t kk = nBytes/2;
	
	// Loop over the number of items in the buffer.
	for(int j=0; j<nItems; ++j)
	{
	  // Initialize ii and jj, which are used to store the
	  // indices of the bytes to be swapped.
	  size_t ii = j*nBytes;
	  size_t jj = ii + nBytes - 1;
	
	  // Swap the bytes.
	  for(size_t i=0; i<kk; ++i, ++ii, --jj)
	  {
	    const char tmp = buf[jj];
	    buf[jj] = buf[ii];
	    buf[ii] = tmp;
	  }
	}
}



