#pragma once 

/*!
 * @file import_structure.hpp
 * @brief The file containing all the import functionalities.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"


/*!
 * @brief An interface class used for importing a solution restart file.
 */
class CImport {

	public:
		/*!
		 * @brief Default constructor of CImport, which initializes a generic import class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */
		CImport(CConfig   *config_container,
				    CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~CImport(void);

    /*!
		 * @brief Pure virtual function that reads the solution from a restart file.
		 * Must be implemented by a derived class. 
		 *
		 * @param[in] config_container pointer to input configuration/dictionary container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[out] SimTime reference to the imported simulation time.
		 */
    virtual void ReadSolutionRestartFile(CConfig    *config_container,
                                         CGeometry  *geometry_container,
                                         CElement  **element_container,
                                         CSolver   **solver_container,
                                         as3double  &SimTime) = 0;
	protected:
		unsigned short nZone;  ///< Total number of zones.

 	private:

};


/*!
 * @brief A class used for importing a solution restart file in ASCII format.
 */
class CASCIIImport : public CImport {

	public:
		/*!
		 * @brief Default constructor of CASCIIImport, which initializes an ASCII import class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */
		CASCIIImport(CConfig   *config_container,
				         CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CASCIIImport(void);

    /*!
		 * @brief Function that reads the solution from a restart file in ASCII format.
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

};


/*!
 * @brief A class used for importing a solution restart file in binary format.
 */
class CBinaryImport : public CImport {

	public:
		/*!
		 * @brief Default constructor of CBinaryImport, which initializes a binary import class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */	
		CBinaryImport(CConfig   *config_container,
				          CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CBinaryImport(void);

    /*!
		 * @brief Function that reads the solution from a restart file in binary format.
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

};

