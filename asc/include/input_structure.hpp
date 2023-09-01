#pragma once

/*!
 * @file input_structure.hpp
 * @brief The file containing functionalities for reading input data files.
 */

#include "option_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "element_structure.hpp"
#include "import_structure.hpp"


/*!
 * @brief A class used for reading input files.
 */
class CInput {

	public:
		/*!
		 * @brief Default constructor of CInput, which initializes an input class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */
		CInput(CConfig   *config_container,
					 CGeometry *geometry_container);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CInput(void);

    /*!
		 * @brief Function that reads the solution from a restart file.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[out] SimTime reference to the imported simulation time.
		 */
    void ReadSolutionRestartFile(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CElement  **element_container,
                                 CSolver   **solver_container,
                                 as3double  &SimTime);
	protected:

	private:
		unsigned short nZone;                        ///< Total number of zones.
		CImport       *import_container = nullptr;   ///< Object for importing files.

		/*!
		 * @brief Function that pre-processes the import file.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */
		void Preprocess_ImportSolution(CConfig   *config_container,
				                           CGeometry *geometry_container);
};



