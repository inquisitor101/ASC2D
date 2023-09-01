#include "ASD.hpp"


int main(int argc, char **argv)
{
	// Check if the correct usage is used.
	std::ostringstream message;
	if( argc != 4 )
	{
		// Create instruction on how to use the program.
		message << "Executable: "     << argv[0] << "\n\n"
			      << "Usage: \n"
			      << "<executable>"     << "  "
						<< "<input dir>"      << "  "
						<< "<input filename>" << "  "
						<< "<# of files>";
		
		// Report error and exit.
		ERROR(message.str());
	}

	// Convert total number of files into a numeric value.
	std::istringstream iss(argv[3]); unsigned long nFile; iss >> nFile;
	// Check if conversion failed.
	if( iss.fail() ) ERROR( "Could not convert nFile into unsigned long.");


	// Initialize required containers to null.
	CImport  *import_container  = nullptr;
	CProcess *process_container = nullptr;

	// Create new import container.
	import_container  = new CImport(argv[1], argv[2], nFile);

	// Create new post process container.
	process_container = new CProcess(argv[1], import_container);

	// Delete objects.
	if( import_container  != nullptr ) delete import_container;
	if( process_container != nullptr ) delete process_container;

	// All good, exit.
	return 0;
}
