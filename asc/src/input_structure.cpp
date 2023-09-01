#include "input_structure.hpp"




CInput::CInput
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor, used to initialize CInput class.
	*/
{
	// Extract number of zones.
	nZone = config_container->GetnZone();

	// Initialize and pre-process the import container, only if a restart is specified.
	if( config_container->GetRestartSolution() )
		Preprocess_ImportSolution(config_container, geometry_container);
}


CInput::~CInput
(
 void
)
 /*
	* Deconstructor for CInput class.
	*/
{
	if( import_container != nullptr ) delete import_container;
}


void CInput::Preprocess_ImportSolution
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function which pre-processes the import solution object.
	*/
{
	// Extract the restart filename in string format.
	std::string fn = config_container->GetRestartFilename();

	// Deduce the extension specified. Do not worry, if there is no extension, 
	// it already will be reported in the below check.
	std::string ext = fn.substr(fn.find(".") + 1);

	// Check what the input restart file format is and assign the necessary import object.
	if      ( ext == "bin" ) import_container = new CBinaryImport(config_container, geometry_container);
	else if ( ext == "txt" ) import_container = new CASCIIImport( config_container, geometry_container);
	else 
		Terminate("CInput::Preprocess_ImportSolution", __FILE__, __LINE__,
				      "Unknown file extension for input restart solution.");
}



void CInput::ReadSolutionRestartFile
(
  CConfig    *config_container,
  CGeometry  *geometry_container,
  CElement  **element_container,
  CSolver   **solver_container,
  as3double  &SimTime
)
 /*
  * Function that reads the solution from a restart file, depending on its input format.
  */
{
	// Read and initialize the current solution from the specified file.
	import_container->ReadSolutionRestartFile(config_container,
			                                      geometry_container,
																						element_container,
																						solver_container,
																						SimTime);
}


