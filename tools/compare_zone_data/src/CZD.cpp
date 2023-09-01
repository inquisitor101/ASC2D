#include "CZD.hpp"


int main(int argc, char **argv)
{
	// Check if the correct usage is used.
	std::ostringstream message;
	if( argc != 5 )
	{
		// Create instruction on how to use the program.
		message << "Executable: "     << argv[0] << "\n\n"
			      << "Usage: \n"
			      << "<executable>"              << "  "     // argv[0]
						<< "<input simulation:A dir>"  << "  "     // argv[1]
						<< "<input simulation:B dir>"  << "  "     // argv[2]
						<< "<input target filenames>"  << "  "     // argv[3]
						<< "<# of files to process>"   << "\n"     // argv[4]
						<< "... Note: sim:B is the reference.";
		
		// Report error and exit.
		ERROR(message.str());
	}

	// Convert total number of files into a numeric value.
	std::istringstream iss(argv[4]); unsigned long nFileTotal; iss >> nFileTotal;
	// Check if conversion failed.
	if( iss.fail() ) ERROR( "Could not convert nFile into unsigned long.");


	// Pre-process number of chunks required.
	as3vector2d<unsigned long> Chunks = PreprocessDataDistribution(argv[1], argv[2], argv[3], nFileTotal);	

	// Initialize process container to null.
	CProcess *process_container = nullptr;

	// Create process object.
	process_container = new CProcess( nFileTotal );
	

	// Import data based on specified chunks, in order not to run out of memory.
	for(unsigned short iChunk=0; iChunk<Chunks.size(); iChunk++)
	{
		// Initialize required containers to null.
		CImport *import_A_container = nullptr; 
		CImport *import_B_container = nullptr;

		// Starting index of the files in this chunk.
		unsigned long I0 = Chunks[iChunk][0];
		// Ending   index of the files in this chunk.
		unsigned long I1 = Chunks[iChunk][1];	

		// Create import object for dir: A, according to user input.
		import_A_container = new CImport(argv[1], argv[3], I0, I1, nFileTotal); 
		// Create import object for dir: B, according to user input.
		import_B_container = new CImport(argv[2], argv[3], I0, I1, nFileTotal);

		// Create process object.
		process_container->ComputeMetrics(import_A_container,
				                              import_B_container,
																	    I0, I1);
		
		// Delete objects.
		if( import_A_container != nullptr ) delete import_A_container;
		if( import_B_container != nullptr ) delete import_B_container;	
	}

	// Create output object.
	COutput *output_container = new COutput( argv[1], argv[3], process_container );

	// Write the output of the processed data.
	output_container->WriteProcessedDataBinary( process_container );


	// Delete objects.
	if( process_container != nullptr ) delete process_container;
	if( output_container  != nullptr ) delete output_container;

	// All good, exit.
	return 0;
}
