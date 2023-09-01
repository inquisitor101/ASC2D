#pragma once 

#include "utility.hpp"
#include <fstream>



class CImport 
{
	public:
		// Default constructor.
		CImport(const char   *directory,
				    const char   *filename,
						unsigned long nfiles);
		
		 // Default destructor, frees any allocated memory.
		~CImport(void);

		// Function, which reads and extracts all the data from the files.
		void ImportDataFromDirectory(const char *directory,
				                         const char *filename);

		// Getter: returns InputData.
		const as3vector4d<as3double> &GetInputData(void)  const {return InputData;}
		// Getter: returns SimTime.
		const as3vector1d<as3double> &GetSimTime(void)    const {return SimTime;}

	protected:


	private:
		// Total number of files to import.
		unsigned long  nFile;
		// Total number of elements in each file.
		unsigned long  nElem;
		// Total number of nodes in each element.
		unsigned short nNode;
		// Total number of variables in each node
		unsigned short nVarRead;

		// Simulation time used in every file.
		// [iTime].
		as3vector1d<as3double> SimTime;

		// Imported data. 
		// [iTime][iElem][iNode][iVar].
		as3vector4d<as3double> InputData;

		// Function, which checks if all the specified files exist.
		void CheckFilesExist(const char   *directory,
				                 const char   *filename,
												 unsigned long nFile);

		// Function, which reads from a reference input file the 
		// necessary header information.
		void DeduceHeaderInformation(const char    *directory,
				                         const char    *filename,
																 unsigned long  iFile);

		// Function, which imports the data from a single file.
		void ImportDataFromFile(const char   *fn,
				                    unsigned long iFile);

};
