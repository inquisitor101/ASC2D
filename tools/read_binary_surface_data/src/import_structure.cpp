#include "import_structure.hpp"



CImport::CImport
(
 const char   *directory,
 const char   *filename,
 unsigned long nfiles
)
	:
		nFile(nfiles)
 /*
	* Constructor, that initializes an import container.
	*/
{
	// Report user specification.
	std::cout << "directory : " << directory << "\n"
		        << "filename  : " << filename  << "\n"
						<< "# of files: " << nFile     << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
	
	
	// Check if all files exist.
	CheckFilesExist(directory, filename, nFile);

	// Deduce the header information (nElem, nNode, nVarRead, SimTime... etc).
	// Assign reference data file at initial step: 0.
	DeduceHeaderInformation(directory, filename, 0);

	// Initialize the total number of physical simulation time needed.
	SimTime.resize(nFile, 0.0);

	// Initialize the importad data.
	InputData.resize(nFile);
	for(unsigned long iFile=0; iFile<nFile; iFile++)
	{
		InputData[iFile].resize(nElem);
		for(unsigned long iElem=0; iElem<nElem; iElem++)
		{
			InputData[iFile][iElem].resize(nNode);
			for(unsigned short iNode=0; iNode<nNode; iNode++)
				InputData[iFile][iElem][iNode].resize(nVarRead, 0.0);
		}
	}

	// Import all the data in the directory.
	ImportDataFromDirectory(directory, filename);
}


CImport::~CImport
(
 void
)
 /*
	* Destructor for CImport, frees dynamic memory.
	*/
{

}


void CImport::CheckFilesExist
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


void CImport::DeduceHeaderInformation
(
 const char   *directory,
 const char   *filename,
 unsigned long iFile
)
 /*
	* Function, which computes from a reference file the necessary information.
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

	// Read the total number of variables to read per DOF in this file.
	if( std::fread(&nVarRead, sizeof(unsigned short), 1, file) != 1 ) ERROR("Could not read nVarRead.");
	if(byteswap) SwapBytes(&nVarRead, sizeof(unsigned short), 1);

	// Read the names of the written variables.
	for(unsigned short i=0; i<nVarRead; i++)
	{
		char varname[CGNS_STRING_SIZE];
		if( std::fread(varname, sizeof(char), CGNS_STRING_SIZE, file) != CGNS_STRING_SIZE ) 
			ERROR("Could not read variable name.");
	}

	// Close file.
	std::fclose(file);
}


void CImport::ImportDataFromDirectory
(
 const char *directory,
 const char *filename
)
 /*
	* Function, which imports all the data from the entire directory.
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

	// Report progress.
	std::cout << "Importing files.... ";

	// Loop over all the files.
	for(unsigned long iFile=0; iFile<nFile; iFile++)
	{
		// Assemble current file name.
		std::stringstream fn;
		fn << basename << iFile << ".bin";

		// Import data in this file.
		ImportDataFromFile(fn.str().c_str(), iFile); 
	}

	// Report progress.
	std::cout << "Done." << std::endl;
}


void CImport::ImportDataFromFile
(
 const char   *fn,
 unsigned long iFile
)
 /*
	* Function, which extract data from a single file.
	*/
{
	// Initialize remaining header information.
	as3double time; 
	unsigned long  nelem;
	unsigned short nnode;
	unsigned short nvarread;

	// Open the file for binary reading.
	FILE *file = std::fopen(fn, "rb");

	// Check if the file can be opened. 
	if( !file )
	{
		// Log message.
		std::ostringstream message;
		message << "Could not open file: \n" 
			      << "'" << fn << "'";
		// Report error and exit.
		ERROR(message.str());
	}

	// Initialize header information.
	int header[1];

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
	if( std::fread(&nelem, sizeof(unsigned long), 1, file) != 1 ) ERROR("Could not read nElem.");
	if(byteswap) SwapBytes(&nelem, sizeof(unsigned long), 1);

	// Read the total number of nodes in each element in this file.
	if( std::fread(&nnode, sizeof(unsigned short), 1, file) != 1 ) ERROR("Could not read nNode.");
	if(byteswap) SwapBytes(&nnode, sizeof(unsigned short), 1);

	// Read the total number of variables to read per DOF in this file.
	if( std::fread(&nvarread, sizeof(unsigned short), 1, file) != 1 ) ERROR("Could not read nVarRead.");
	if(byteswap) SwapBytes(&nvarread, sizeof(unsigned short), 1);

	// Consistency check.
	if( nelem    != nElem    ) ERROR("nElem is not consistent.");
	if( nnode    != nNode    ) ERROR("nNode is not consistent.");
	if( nvarread != nVarRead ) ERROR("nVarRead is not consistent.");

	// Save the current time.
	SimTime[iFile] = time;

	// Read the names of the written variables.
	for(unsigned short i=0; i<nVarRead; i++)
	{
		char varname[CGNS_STRING_SIZE];
		if( std::fread(varname, sizeof(char), CGNS_STRING_SIZE, file) != CGNS_STRING_SIZE ) 
			ERROR("Could not read variable name.");
	}

	// Read data on each elementi.
	for(unsigned long iElem=0; iElem<nElem; iElem++)
	{
		// Determine the number of items to read per element.
		const unsigned int sizeread = nNode*nVarRead;

		// Initialize the buffer that reads data for each element.
		as3vector1d<as3double> readbuf(sizeread, 0.0);

		// Read the data belonging to this element.
		if( std::fread(readbuf.data(), sizeof(as3double), sizeread, file) != sizeread ) 
			ERROR("Data could not be read.");

		// Check if there need be any byte swapping done.
		if(byteswap) SwapBytes(readbuf.data(), sizeof(as3double), sizeread);

		// Reset local nodal index.
		unsigned short ii = 0;

		// Loop over all the nodes and assign the data to this element.
		for(unsigned short iNode=0; iNode<nNode; iNode++)
		{
			for(unsigned short iVarRead=0; iVarRead<nVarRead; iVarRead++, ii++)
				InputData[iFile][iElem][iNode][iVarRead] = readbuf[ii];
		}
	}

	// Close file.
	std::fclose(file);
}






