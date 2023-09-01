#pragma once

/*!
 * @file output_structure.hpp
 * @brief The file containing all the output functionalities.
 */

#include "option_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "vtk_structure.hpp"
#include "restart_structure.hpp"


/*!
 * @brief A class used for writing the output information.
 */
class COutput {

	public:
		/*!
		 * @brief Default constructor of COutput, which initializes an output class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] MapGlobalToLocal reference to the indices of the elements.
		 */
		COutput(CConfig                    *config_container,
						CGeometry                  *geometry_container,
						as3vector2d<unsigned long> &MapGlobalToLocal);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~COutput(void);


		/*!
		 * @brief Function that writes a GNU-plot file in ASCII format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] iIter current iteration number.
		 * @param[in] SimTime current physical time.
		 * @param[in] MonitoringData reference to the data monitored.
		 */
		void WriteGNUplot(CConfig                *config_container,
				              const unsigned long     iIter,
				              const as3double         SimTime,
				              as3vector1d<as3double> &MonitoringData);

		/*!
		 * @brief Function that writes VTK data file, depending on the format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 */
		void WriteFileVTK(CConfig    *config_container,
											CGeometry  *geometry_container,
											CElement  **element_container,
											CSolver   **solver_container);

    /*!
		 * @brief Function that writes temporal data to a file.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] filename pointer to the output file name.
		 * @param[in] filetype file format type.
		 * @param[in] time written physical time.
		 * @param[in] writebuf reference to buffer that stores the written data.
		 */
    void WriteDataToFile(CConfig                      *config_container,
                         CGeometry                    *geometry_container,
												 const char                   *filename,
												 const unsigned short          filetype,
                         const as3double               time,
                         const as3vector1d<as3double> &writebuf);

    /*!
		 * @brief Function that writes the data to a file.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] SimTime current physical simulation time.
		 */
    void WriteSolutionToFile(CConfig    *config_container,
                             CGeometry  *geometry_container,
                             CElement  **element_container,
                             CSolver   **solver_container,
                             as3double   SimTime);

    /*!
		 * @brief Function that writes zone data to a file.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] SimTime current physical simulation time.
		 */
    void WriteZoneDataToFile(CConfig    *config_container,
                             CGeometry  *geometry_container,
                             CSolver   **solver_container,
                             as3double   SimTime);

		/*!
		 * @brief Function that writes the boundary data to a file.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial condition container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] iSampleBoundary input boundary ID. 
		 * @param[in] SimTime current physical simulation time.
		 * @param[in] IterCount current (temporal-)iteration number.
		 */
		void WriteBoundaryDataToFile(CConfig       *config_container,
				                         CGeometry     *geometry_container,
																 CElement      *element_container,
																 CInitial      *initial_container,
																 CSolver       *solver_container,
																 unsigned short iSampleBoundary,
																 as3double      SimTime,
																 unsigned long  IterCount);

	protected:

	private:
		unsigned short nZone;                       ///< Total number of zones.
		bool           WriteHeaderInfoDataFile;     ///< Option to write a header in binary data files.
		unsigned long  FileNumberVTK;               ///< VTK file number used in writing sequence.
    unsigned long  FileNumberDataSolution;      ///< Written file number for solution data files.
    unsigned long  FileNumberZoneData;          ///< Written file number for zone solution data files.
		CFileVTK      *vtk_container     = nullptr; ///< Object for VTK file format containers.
		CRestart      *restart_container = nullptr; ///< Object for restart solution containers.

		/*!
		 * @brief Function which pre-processes the VTK object.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] MapGlobalToLocal reference to the indices of the elements.
		 */
		void Preprocess_VTK(CConfig                    *config_container,
				                CGeometry                  *geometry_container,
												as3vector2d<unsigned long> &MapGlobalToLocal);

		/*!
		 * @brief Function which pre-processes the restart solution object.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */
		void Preprocess_RestartSolution(CConfig   *config_container,
				                            CGeometry *geometry_container);
};


