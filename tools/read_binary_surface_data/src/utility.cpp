#include "utility.hpp"



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





