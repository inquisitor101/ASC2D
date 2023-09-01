#pragma once

/*!
 * @file vtk_structure.hpp
 * @brief The file containing all the VTK output functionalities.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "variable_structure.hpp"


/*!
 * @brief An interface class used for initializing a generic VTK class.
 */
class CFileVTK {

	public:
		/*!
		 * @brief Default constructor of CFileVTK, which initializes a generic VTK class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] MapGlobalToLocal reference to the indices of the elements.
		 */	
		CFileVTK(CConfig                    *config_container,
				     CGeometry                  *geometry_container,
						 as3vector2d<unsigned long> &MapGlobalToLocal);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~CFileVTK(void);

		/*!
		 * @brief Pure virtual function that writes VTK data file. 
		 * Must be implemented by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] FileNumberVTK file number used in writing sequence. 
		 */
		virtual void WriteFileVTK(CConfig             *config_container,
											        CGeometry           *geometry_container,
															CElement           **element_container,
											        CSolver            **solver_container,
															const unsigned long  FileNumberVTK) = 0;
	protected:
		unsigned short nZone;              ///< Total number of zones.
		const char    *OutputVTKFilename;  ///< Output VTK filename.
		bool           WriteDensity;       ///< Option for writing a density variable.
		bool           WriteMomentum;      ///< Option for writing a momentum variable.
		bool           WriteTotalEnergy;   ///< Option for writing a total energy variable.
		bool           WritePressure;      ///< Option for writing a pressure variable.
		bool           WriteVelocity;      ///< Option for writing a velocity variable.
		bool           WriteVorticity;     ///< Option for writing a vorticity variable.
		bool           WriteMach;          ///< Option for writing a Mach number variable.
		bool           WriteTemperature;   ///< Option for writing a temperature variable.
		bool           WriteEntropy;       ///< Option for writing a specific entropy variable.

	private:

};


/*!
 * @brief A class used for initializing an ASCII type VTK class.
 */
class CASCIIFileVTK : public CFileVTK {

	public:
		/*!
		 * @brief Default constructor of CASCIIFileVTK, which initializes an ASCII type VTK class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] MapGlobalToLocal reference to the indices of the elements.
		 */	
		CASCIIFileVTK(CConfig                    *config_container,
				          CGeometry                  *geometry_container,
									as3vector2d<unsigned long> &MapGlobalToLocal);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CASCIIFileVTK(void);

		/*!
		 * @brief Function that writes VTK data file in ASCII format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] FileNumberVTK file number used in writing sequence. 
		 */
		void WriteFileVTK(CConfig             *config_container,
											CGeometry           *geometry_container,
											CElement           **element_container,
											CSolver            **solver_container,
											const unsigned long  FileNumberVTK);
	protected:

	private:
		as3vector2d<unsigned short> ConnLocal;           ///< Local indicial connectivity matrix.
		CScalarVariable *VariableDensity     = nullptr;  ///< Type of density scalar variable used.
		CScalarVariable *VariableEnergy      = nullptr;  ///< Type of energy scalar variable used.
		CScalarVariable *VariablePressure    = nullptr;  ///< Type of pressure scalar variable used.
    CScalarVariable *VariableTemperature = nullptr;  ///< Type of temperature scalar variable used.
    CScalarVariable *VariableMachNumber  = nullptr;  ///< Type of Mach number scalar variable used.
		CVectorVariable *VariableMomentum    = nullptr;  ///< Type of momentum vector variable used.
		CScalarVariable *VariableEntropy     = nullptr;  ///< Type of entropy scalar variable used.


		/*!
		 * @brief Function that writes a scalar parameter in ASCII format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] Paraview_File reference to input Paraview file.
		 * @param[in] Variable pointer to type of scalar variable to write.
		 */
		void WriteScalar(CConfig         *config_container,
										 CGeometry       *geometry_container,
										 CSolver     	  **solver_container,
										 std::ofstream   &Paraview_File,
										 CScalarVariable *Variable);

		/*!
		 * @brief Function that writes a vector parameter in ASCII format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] Paraview_File reference to input Paraview file.
		 * @param[in] Variable pointer to type of vector variable to write.
		 */
		void WriteVector(CConfig         *config_container,
										 CGeometry       *geometry_container,
										 CSolver     	  **solver_container,
										 std::ofstream   &Paraview_File,
										 CVectorVariable *Variable);
 
		/*!
		 * @brief Function that writes a vorticity (scalar) parameter in ASCII format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] Paraview_File reference to input Paraview file.
		 */
		void WriteScalarVorticity(CConfig         *config_container,
										          CGeometry       *geometry_container,
												      CElement       **element_container,
										          CSolver     	 **solver_container,
										          std::ofstream   &Paraview_File);

		/*!
		 * @brief Function that writes a PML auxiliary (scalar) parameter in ASCII format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] Paraview_File reference to input Paraview file.
		 */
		void WriteScalarAuxPML(CConfig         *config_container,
                           CGeometry       *geometry_container,
                           CSolver     	  **solver_container,
                           std::ofstream   &Paraview_File,
                           unsigned short   iVar);

 		/*!
		 * @brief Function that writes a PML auxiliary (vector) parameter in ASCII format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] Paraview_File reference to input Paraview file.
		 */
		void WriteVectorAuxPML(CConfig         *config_container,
                           CGeometry       *geometry_container,
                           CSolver     	  **solver_container,
                           std::ofstream   &Paraview_File,
                           unsigned short   iVar);
};


/*!
 * @brief A class used for initializing a binary type VTK class.
 */
class CBinaryFileVTK : public CFileVTK {

	public:
		/*!
		 * @brief Default constructor of CBinaryFileVTK, which initializes a binary type VTK class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] MapGlobalToLocal reference to the indices of the elements.
		 */	
		CBinaryFileVTK(CConfig                    *config_container,
				           CGeometry                  *geometry_container,
									 as3vector2d<unsigned long> &MapGlobalToLocal);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CBinaryFileVTK(void);

		/*!
		 * @brief Function that writes VTK data file in binary format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] FileNumberVTK file number used in writing sequence. 
		 */
		void WriteFileVTK(CConfig             *config_container,
											CGeometry           *geometry_container,
											CElement           **element_container,
											CSolver            **solver_container,
											const unsigned long  FileNumberVTK);
	protected:

	private:
    as3vector2d<unsigned long>  mMapGlobalToLocal;  ///< Map data from [iZone][iElemZone] to [iElem]. 
																										///< This is used for parallelization efficiency.
																										///< Dimension: [iElem][iData], where:[iData] is: [0]: iZone, [1]: iElemZone.

		bool BigEndian;                                 ///< Flag that determines whether this machine uses big endian or not.
		as3vector1d<unsigned long>  nElemZone;          ///< The number of elements in each zone.
		as3vector1d<unsigned short> nNodeZone;          ///< The number of nodes per element in each zone.
		as3vector1d<unsigned short> nPolyZone;          ///< The number of solution polynomial in each zone.
		as3vector1d<int>            nDOFsTotZone;       ///< The total DOFs per each zone.
		as3vector1d<int>            nSubElemZone;       ///< The size of each sub-element of in each zone.
		as3vector1d<int>            nDOFsTotWritten;    ///< Total number of DOFs written in previous zones.
		as3vector1d<unsigned long>  nxElemZone;         ///< Total number of elements in x-direction in each zone.
		as3vector1d<unsigned long>  nyElemZone;         ///< Total number of elements in y-direction in each zone.
		as3vector1d<std::string>    VariableNames;      ///< Variable names for writing.

		/*!
		 * @brief Function that computes and stores the data required for writing a binary VTK file.
		 *
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] varnames reference to the variable string names.
		 * @param[in] varbuf reference to the variables buffer data.
		 */
	  void DetermineVisualizationData(CGeometry                 *geometry_container,
				                            CElement                 **element_container,
				                            CSolver                  **solver_container,
				                            as3vector1d<std::string>  &varnames,
																		as3vector2d<float>        &varbuf);
};

