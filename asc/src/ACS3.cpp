#include "ACS3.hpp"



int main(int argc, char **argv){

  // Check if correct usage is executed.
  std::ostringstream message;
  if (argc != 2){
    message << "Usage: \n" << argv[0]
            << "\n <input config file>";
    Terminate("main", __FILE__, __LINE__, message.str());
  }


  // Create configuration file.
  CConfig *config = nullptr; config = new CConfig(argv[1]);

  // Create a pointer to the main driver.
  CDriver *driver = nullptr; driver = new CDriver(config);


	// Start solver.
	driver->StartSolver();

	// Free memory.
  if( driver != nullptr ) delete driver;
	// if( config != nullptr ) delete config;

	// Exit happily.
	return 0;
}

