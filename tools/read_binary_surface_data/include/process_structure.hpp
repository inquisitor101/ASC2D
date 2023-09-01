#pragma once 

#include "utility.hpp"
#include "import_structure.hpp"

#if __cplusplus > 201103L
	#include <filesystem>
#else
	#include <sys/types.h>
	#include <sys/stat.h>
#endif



class CProcess
{
	public:
		// Default constructor.
		CProcess(const char    *directory,
				     const CImport *import_container);
		
		// Destructor.
		~CProcess(void);

		// Function, which writes all the data of the directory in time-domain.
		// i.e. each output [iSpace] file has: [iTime][iData].
		void WriteTimeDomainData(const char    *directory,
				                     const CImport *import_container);

	protected:


	private:
		// Total number of time steps.
		unsigned long  nTime;
		// Total number of spatial probes.
		unsigned long  nProbe;
		// Total number of elements in each file.
		unsigned long  nElem;
		// Total number of nodes in each element.
		unsigned short nNode;
		// Total number of item in each node
		unsigned short nItem;


};
