#include "import_structure.hpp"



CImport::CImport
(
 const char         *directory,
 const char         *filename,
 const unsigned long I0,
 const unsigned long I1,
 const unsigned long nFileTotal
)
 /*
	* Constructor, that initializes an import container.
	*/
{
	// Deduce the total number of files.
	nFile = I1 - I0;

	// Deduce percentage complete, thus far.
	const as3double perc0 = 100.0*( (as3double) I0/nFileTotal );
	const as3double perc1 = 100.0*( (as3double) I1/nFileTotal );

	// Report user specification.
	std::cout << "\n"
		        << "*) directory : " << directory << "\n"
		        << "   filename  : " << filename  << "\n"
						<< "   from index: " << I0        << "\n"
						<< "   to   index: " << I1-1      << "\n"
						<< "   # of files: " << nFile     << "\n";
	std::cout << std::fixed 
		        << std::setprecision(2) 
						<< std::setw(2)
						<< "   percentage: " << perc0 << "% - " 
						                     << perc1 << "%\n" 
																 << std::endl;

	// Assign total number of elements and nodes.
	nElem = nElemExpected;
	nNode = nNodeExpected;

	// Reserve memory for required temporal data.
	SimTime.resize(nFile, 0.0);

	// Reserve memory for required solution data.
	ImportedData.resize(nFile);

	// Allocate the memory on each DOF.
	for(unsigned long iTime=0; iTime<ImportedData.size(); iTime++)
	{
		ImportedData[iTime].resize(nElem);
		for(unsigned long iElem=0; iElem<ImportedData[iTime].size(); iElem++)
		{
			ImportedData[iTime][iElem].resize(nVar);
			for(unsigned short iVar=0; iVar<ImportedData[iTime][iElem].size(); iVar++)
			{
				// Allocate actual memory for each variable.
				ImportedData[iTime][iElem][iVar].resize(nNode, 0.0);
			}
		}
	}

	// Import data into this class.
	ImportDataFromDirectory(directory, filename, I0, I1);

	// Report progress.
	std::cout << "---------------------------------------------------" << std::endl;
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


void CImport::ImportDataFromDirectory
(
 const char         *directory,
 const char         *filename,
 const unsigned long I0,
 const unsigned long I1
)
 /*
	* Function, which imports all the data from the entire directory 
	* from index I0 (inclusive) to I1 (exclusive).
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
	std::cout << "   Importing files.... ";

	// Loop over all the files.
	for(unsigned long iFile=0; iFile<nFile; iFile++)
	{
		// Assemble current file name.
		std::stringstream fn;
		fn << basename << I0+iFile << ".bin";

		// Import data in this file.
		ImportDataFromFile(fn.str().c_str(), iFile); 
	}

	// Consistency check.
	if( I0+nFile != I1 ) ERROR("Inconsistency in total number of files read.");

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

	// Consistency check.
	if( nElem != nelem ) ERROR("nElem is not consistent.");
	if( nNode != nnode ) ERROR("nNode is not consistent.");

	// Save the current time.
	SimTime[iFile] = time;

	// Determine the number of items to read per entire zone.
	const unsigned long sizeread = (unsigned long) nElem*nNode*nVar;

	// Initialize the buffer that reads the entire data.
	as3vector1d<as3double> readbuf(sizeread, 0.0);

	// Read the actual data into the buffer.
	if( std::fread(readbuf.data(), sizeof(as3double), sizeread, file) != sizeread ) 
		ERROR("Data could not be read.");

	// Check if there need be any byte swapping done.
	if(byteswap) SwapBytes(readbuf.data(), sizeof(as3double), sizeread);

	// Copy the data on each element, according to format.
	for(unsigned long iElem=0; iElem<nElem; iElem++)
	{
		// Reset local nodal index.
		unsigned short ii = 0;

		// Loop over all the variables.
		for(unsigned short iVar=0; iVar<nVar; iVar++)
		{
			// Loop over all the nodes and extract the data.
			for(unsigned short iNode=0; iNode<nNode; iNode++, ii++)
				ImportedData[iFile][iElem][iVar][iNode] = readbuf[ii];
		}
	}

	// Close file.
	std::fclose(file);
}


