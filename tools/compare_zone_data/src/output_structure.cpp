#include "output_structure.hpp"



COutput::COutput
(
 const char     *directory,
 const char     *filename,
 const CProcess *process_container
)
 /*
	* Constructor, that initializes an output container.
	*/
{
	// Pre-assemble basename of directory and file.
	std::string basename;
	basename += directory;
	
	// If there is no slash in the directory as specified by the user, add one.
	if( basename.back() != '/' ) basename += "/";
	// Add the current zone processing directory extension.
	basename += "../zone_proc/";

	// Depending on the C++ standard used, create the direcory if it does not exist.
#if __cplusplus > 201103L
	// This standard is higher than C++11.

	// Check if the current subdirectory exists, otherwise create one.
	if( !std::filesystem::is_directory( basename ) || !std::filesystem::exists( basename ) )
	{
		std::cout << "Directory: " << basename << ", does not exist, creating one." << std::endl;
		std::filesystem::create_directory( basename );	
	}
#else
	// The standard is equal or below C++11.
	
	// Check if the currect subdirectory exists, otherwise create one.
	struct stat info; bool DirExists = false;
	if( stat( basename.c_str(), &info ) != 0 ) 
		std::cout << "Cannot access: " << basename << std::endl;
	else if( info.st_mode & S_IFDIR )  
		DirExists = true;
	else
		std::cout << "Directory does not exist: '" << basename << "'" << std::endl;
	
	// Create the subdirectory.
	if( !DirExists )
	{
		std::cout << "Creating directory: '" << basename << "' ... ";
		const int dir_err = mkdir(basename.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if( dir_err == -1 ) std::cout << "Could not create the directory." << std::endl;
		else                std::cout << "Done." << std::endl;
	}
#endif

	// For now, initialize the output files to the 3 norms.
	OutputFileName.resize(3, basename);
	OutputFileName[0] += "L1_norm.bin";
	OutputFileName[1] += "L2_norm.bin";
	OutputFileName[2] += "Linf_norm.bin";

	// Report progress.
	std::cout << "*) Writing files: \n";
	for(unsigned short i=0; i<OutputFileName.size(); i++)
		std::cout << "   file #" << i+1 << ": " << OutputFileName[i] << std::endl;
}


COutput::~COutput
(
 void
)
 /*
	* Destructor for COutput, frees dynamic memory.
	*/
{

}


void COutput::WriteProcessedDataBinary
(
 const CProcess *process_container
)
 /*
	* Function that writes the processed output in binary format.
	*/
{
	// Report progress.
	std::cout << "\n   Writing data in binary format...";

	// Extract processed data and simulation time.
	auto& data = process_container->GetData();
	auto& time = process_container->GetTime();

	// Consistency check.
	if( OutputFileName.size() != data[0].size() ) ERROR("Inconsistent number of files and metrics.");
	if( data.size()           != time.size()    ) ERROR("Inconsistent size of data and time steps.");

	// Deduce total number of time steps.
	unsigned long  nTime     = data.size();
	unsigned short nMetric   = data[0].size();
	unsigned short nVarWrite = data[0][0].size();

	// Total number of data to write in every file.
	// Note, the additional variable corresponds to time.
	unsigned long nData = (unsigned long) nTime*(1+nVarWrite); 


	// Loop over every file and write it separately.
	for(unsigned short iFile=0; iFile<nMetric; iFile++)
	{
		// Reserve memory for the write buffer.
		as3vector1d<as3double> writebuf(nData, 0.0);

		// Data item counter.
		unsigned long ii = 0;

		// Loop over every time step.
		for(unsigned long iTime=0; iTime<data.size(); iTime++)
		{
			writebuf[ii++] = time[iTime];
			for(unsigned short iVar=0; iVar<nVarWrite; iVar++, ii++)
				writebuf[ii] = data[iTime][iFile][iVar];
		}

		// Open file.
		FILE *fh = std::fopen( OutputFileName[iFile].c_str(), "wb" ); 

		// Report error if file could not be opened.
		if( !fh ) ERROR("Could not open file for writing data.");

		// Write the integer header.
		std::fwrite(&AS3_MagicNumber, 1, sizeof(int), fh);

		// Write the actual processed data.
		std::fwrite(writebuf.data(), writebuf.size(), sizeof(as3double), fh);

		// Close the file.
		std::fclose(fh);
	}


	// Report progress.
	std::cout << " Done." << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;
}



