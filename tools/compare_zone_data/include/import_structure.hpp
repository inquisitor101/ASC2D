#pragma once 

#include "utility.hpp"
#include <iomanip>



class CImport 
{
	public:
		// Default constructor.
		CImport(const char         *directory,
				    const char         *filename,
						const unsigned long I0,
						const unsigned long I1,
						const unsigned long nFileTotal);
		
		// Destructor.
		~CImport(void);

		// Function, which reads and extracts all the data from the files.
		void ImportDataFromDirectory(const char          *directory,
				                         const char          *filename,
																 const unsigned long  I0,
																 const unsigned long  I1);

		// Getter: returns ImportData.
		const as3vector4d<as3double> &GetImportedData(void) const {return ImportedData;}
		// Getter: returns SimTime.
		const as3vector1d<as3double> &GetSimTime(void)      const {return SimTime;}

	protected:

	private:
		// Total number of files to import.
		unsigned long  nFile;
		// Total number of elements in each file.
		unsigned long  nElem;
		// Total number of nodes in each element of each file.
		unsigned short nNode;

		// Simulation time used in every file.
		// [iTime].
		as3vector1d<as3double> SimTime;

		// Imported data. 
		// [iTime][iElem][iVar][iNode].
		as3vector4d<as3double> ImportedData;

		// Function, which imports the data from a single file.
		void ImportDataFromFile(const char   *fn,
				                    unsigned long iFile);

};
