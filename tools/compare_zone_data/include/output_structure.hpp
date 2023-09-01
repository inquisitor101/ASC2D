#pragma once 

#include "utility.hpp"
#include "process_structure.hpp"

#if __cplusplus > 201103L
	#include <filesystem>
#else
	#include <sys/types.h>
	#include <sys/stat.h>
#endif



class COutput 
{
	public:
		// Default constructor.
		COutput(const char     *directory,
				    const char     *filename, 
				    const CProcess *process_container);
		
		// Destructor.
		~COutput(void);
	
		// Function that writes the processed output in binary format.
		void WriteProcessedDataBinary(const CProcess *process_container);

	protected:

	private:
		// File names to write.
		// For now, these are fixed to 3 corresponding to:
		// L1 norm, L2 norm and Linf norm.
		as3vector1d<std::string> OutputFileName;
};


