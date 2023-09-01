#include "process_structure.hpp"




CProcess::CProcess
(
 const char    *directory,
 const CImport *import_container
)
 /*
	* Constructor, that initializes a process container. 
	*/
{
	// Extract the imported data.
	auto& data = import_container->GetInputData();

	// Extract current data information.
	nTime = data.size();
	nElem = data[0].size();
	nNode = data[0][0].size();
	nItem = data[0][0][0].size();

	// Deduce the total number of spatial probes, i.e. DOFs on boundary.
	nProbe = nElem*nNode;

	// Write the rearranged data in all directory.
	WriteTimeDomainData(directory, import_container);
}


CProcess::~CProcess
(
 void
)
 /*
	* Destructor for CProcess, frees dynamic memory.
	*/
{

}


void CProcess::WriteTimeDomainData
(
 const char    *directory,
 const CImport *import_container
)
 /*
	* Function, which converts the input data into time-domain format for 
	* all the directory.
	*/
{
	// Pre-assemble basename of directory and file.
	std::string basename;
	basename += directory;
	
	// If there is no slash in the directory as specified by the user, add one.
	if( basename.back() != '/' ) basename += "/";
	
	// Remove last character in basename, i.e. "/".
	basename.pop_back();
	// Form rest of basename.
	basename += "_proc/";


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


	// Append the current filename to the basename.
	basename += "xy_";

	// Report progress.
	std::cout << "Processing files... ";

	// Extract input data and simulation time.
	auto& data = import_container->GetInputData();
	auto& time = import_container->GetSimTime();

	//// Initialize probe counter.
	//unsigned long iProbe = 0;

	//// Loop over all the elements and nodes, which constitute the probes.
	//for(unsigned long iElem=0; iElem<nElem; iElem++)
	//{
	//	for(unsigned short iNode=0; iNode<nNode; iNode++, iProbe++)
	//	{
	//		// Create the current filename.
	//		std::stringstream fn;
	//		fn << basename << iProbe << ".csv";

  //		// Create data file.
	//		std::ofstream Data_File;

	//		// Deduce filename and number.
	//		Data_File.open(fn.str().c_str(), std::ios::out);

	//		// Open file and write header and info.
	//		Data_File.precision(10);
	//		Data_File << "# Probe location: " << data[0][iElem][iNode][0] << "\n";
	//		Data_File << "#\n";

	//		// Write output format.
	//		Data_File << "#      time      ,       w(-)      ,       w(+)\n";


	//		// Loop over all the temporal samples.
	//		for(unsigned long iTime=0; iTime<nTime; iTime++)
	//		{
	//			// Write data, according to time-domain format.
	//			Data_File << std::scientific << std::showpos 
	//				        << time[iTime]                  << ","
	//								<< data[iTime][iElem][iNode][1] << ","
	//								<< data[iTime][iElem][iNode][2] << "\n";
	//		}


	//		// Close file.
	//		Data_File.close();
	//	}
	//}

	//// Report progress.
	//std::cout << "Done." << std::endl;


	// //
	// Write in binary format.
	// //

	// Determine total number of items per file(probe).
	unsigned long nData = nTime*nItem;

	// Initialize probe counter.
	unsigned long iProbe = 0;

	// Loop over all the elements and nodes, which constitute the probes.
	for(unsigned long iElem=0; iElem<nElem; iElem++)
	{
		for(unsigned short iNode=0; iNode<nNode; iNode++, iProbe++)
		{
			// Initialize the buffer vector.
			as3vector1d<as3double> writeBuf(nData, 0.0);

			// Reset the data counter.
			unsigned long iData = 0;

			// Loop over all the temporal samples.
			for(unsigned long iTime=0; iTime<nTime; iTime++)
			{
				writeBuf[iData++] = time[iTime];
				writeBuf[iData++] = data[iTime][iElem][iNode][1];
				writeBuf[iData++] = data[iTime][iElem][iNode][2];
			}

			// Create the current filename.
			std::stringstream fn;
			fn << basename << iProbe << ".bin";

			// Open file.
			FILE *fh = std::fopen( fn.str().c_str(), "wb" ); 

			// Report error if file could not be opened.
			if( !fh ) ERROR("Could not open file for writing data.");

			// Write the integer header.
			std::fwrite(&AS3_MagicNumber, 1, sizeof(int), fh);

			// Write the actual processed data.
			std::fwrite(writeBuf.data(), writeBuf.size(), sizeof(as3double), fh);

			// Close the file.
			std::fclose(fh);
		}
	}

	// Report progress.
	std::cout << "Done." << std::endl;
}







