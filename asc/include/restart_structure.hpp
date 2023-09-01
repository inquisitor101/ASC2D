#pragma once

/*!
 * @file restart_structure.hpp
 * @brief The file containing all the restart files writing functionalities.
 */

#include "option_structure.hpp"
#include "element_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"


/*!
 * @brief An interface class used for initializing a generic restart class.
 */
class CRestart {

	public:
		/*!
		 * @brief Default constructor of CRestart, which initializes a generic restart class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */
		CRestart(CConfig   *config_container,
				     CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		virtual ~CRestart(void);
						 
    /*!
		 * @brief Pure virtual function that writes the data to a file.
		 * Must be implemented by a derived class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] SimTime current physical (simulation) time.
		 * @param[in] FileNumberDataSolution file number index.
		 */
    virtual void WriteSolutionToFile(CConfig             *config_container,
                                     CGeometry           *geometry_container,
                                     CElement           **element_container,
                                     CSolver            **solver_container,
                                     as3double            SimTime,
																		 const unsigned long  FileNumberDataSolution) = 0;
	protected:
		unsigned short nZone;  ///< Total number of zones.

	private:

};


/*!
 * @brief A class used for initializing an ASCII restart class.
 */
class CASCIIRestart : public CRestart {

	public:
		/*!
		 * @brief Default constructor of CASCIIRestart, which initializes an ASCII restart class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */	
		CASCIIRestart(CConfig   *config_container,
				          CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */		
		~CASCIIRestart(void);
						
    /*!
		 * @brief Function that writes the data to a file in ASCII format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] SimTime current physical (simulation) time.
		 * @param[in] FileNumberDataSolution file number index.
		 */
    void WriteSolutionToFile(CConfig             *config_container,
                             CGeometry           *geometry_container,
                             CElement           **element_container,
                             CSolver            **solver_container,
                             as3double            SimTime,
														 const unsigned long  FileNumberDataSolution);
	protected:

	private:

};


/*!
 * @brief A class used for initializing a binary restart class.
 */
class CBinaryRestart : public CRestart {

	public:
		/*!
		 * @brief Default constructor of CBinaryRestart, which initializes a binary restart class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 */		
		CBinaryRestart(CConfig  *config_container,
				          CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */		
		~CBinaryRestart(void);
						 
    /*!
		 * @brief Function that writes the data to a file in binary format.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] solver_container pointer to input solver container.
		 * @param[in] SimTime current physical (simulation) time.
		 * @param[in] FileNumberDataSolution file number index.
		 */
    void WriteSolutionToFile(CConfig             *config_container,
                             CGeometry           *geometry_container,
                             CElement           **element_container,
                             CSolver            **solver_container,
                             as3double            SimTime,
														 const unsigned long  FileNumberDataSolution);
	protected:

	private:

};
