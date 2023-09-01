#pragma once 

/*!
 * @file reflection_structure.hpp
 * @brief The file containing all the reflection coefficient definitions.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "initial_structure.hpp"
#include "solver_structure.hpp"


/*!
 * @brief An interface class used for initializing a generic reflection coefficient.
 */
class CReflection {

	public:
		/*!
		 * @brief Default constructor of CReflection, which initializes a generic reflection coefficient.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] iZone input zone ID.
		 */	
		CReflection(CConfig       *config_container,
		            CGeometry     *geometry_container,
				        CElement      *element_container,
				        CInitial      *initial_container,
				        CSolver       *solver_container,
				        unsigned short iZone);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		virtual ~CReflection(void);

		/*!
		 * @brief Pure virtual function which computes the reflection coefficient. 
		 * Must be overwritten by a derived class.
		 *
 		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] solver_container pointer to current zone solver container.
		 *
		 * @return reflection coefficient
		 */
		virtual as3double ComputeCoefficient(CConfig   *config_container,
				                                 CGeometry *geometry_container,
				                                 CElement  *element_container,
															           CSolver   *solver_container) = 0;

	protected:
		unsigned short zoneID;  ///< Current zone ID.

	private:

};


/*!
 * @brief A class used for initializing the reflection coefficient: R1.
 */
class CRatioR1Reflection : public CReflection {

	public:
		/*!
		 * @brief Default constructor of CRatioR1Reflection, which initializes the reflection coefficient: R1.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] iZone input zone ID.
		 */	
		CRatioR1Reflection(CConfig       *config_container,
				               CGeometry     *geometry_container,
								       CElement      *element_container,
								       CInitial      *initial_container,
								       CSolver       *solver_container,
								       unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CRatioR1Reflection(void) override;

		/*!
		 * @brief Function which computes the reflection coefficient: R1.
		 * R_1(t) = max( p(t)/pInf ) / max( u(t)/uInf ).
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] solver_container pointer to current zone solver container.
		 *
		 * @return reflection coefficient: R1
		 */
		as3double ComputeCoefficient(CConfig   *config_container,
		                             CGeometry *geometry_container,
		                             CElement  *element_container,
											           CSolver   *solver_container);
	
	private:
		as3double pInf;  ///< Average pressure.
		as3double uInf;  ///< Average u-velocity.
};


/*!
 * @brief A class used for initializing the reflection coefficient: R2.
 */
class CRatioR2Reflection : public CReflection {

	public:
		/*!
		 * @brief Default constructor of CRatioR2Reflection, which initializes the reflection coefficient: R2.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] iZone input zone ID.
		 */	
		CRatioR2Reflection(CConfig       *config_container,
				               CGeometry     *geometry_container,
								       CElement      *element_container,
								       CInitial      *initial_container,
								       CSolver       *solver_container,
								       unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CRatioR2Reflection(void) override;

		/*!
		 * @brief Function which computes the reflection coefficient: R2.
		 * R_2(t) = max( p(t)/pInf ) / max( M(t)/Minf ).
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] solver_container pointer to current zone solver container.
		 *
		 * @return reflection coefficient: R2
		 */
		as3double ComputeCoefficient(CConfig   *config_container,
		                             CGeometry *geometry_container,
		                             CElement  *element_container,
											           CSolver   *solver_container);
	
	private:
		as3double pInf;  ///< Average pressure.
		as3double Minf;  ///< Average Mach number.
};


/*!
 * @brief A class used for initializing the reflection coefficient: R3.
 */
class CRatioR3Reflection : public CReflection {

	public:
		/*!
		 * @brief Default constructor of CRatioR3Reflection, which initializes the reflection coefficient: R3.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] iZone input zone ID.
		 */	
		CRatioR3Reflection(CConfig       *config_container,
				               CGeometry     *geometry_container,
								       CElement      *element_container,
								       CInitial      *initial_container,
								       CSolver       *solver_container,
								       unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CRatioR3Reflection(void) override;

		/*!
		 * @brief Function which computes the reflection coefficient: R3. 
		 * R_3(t) = max( w(-) / w(+) ).
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] solver_container pointer to current zone solver container.
		 *
		 * @return reflection coefficient: R3
		 */
		as3double ComputeCoefficient(CConfig   *config_container,
		                             CGeometry *geometry_container,
		                             CElement  *element_container,
											           CSolver   *solver_container);
	
	private:
		as3double      pInf;      ///< Average pressure.
		as3double      uInf;      ///< Average u-velocity.
		as3double      vInf;      ///< Average v-velocity.
		unsigned short iBoundary; ///< Current boundary ID used for processing.
};


/*!
 * @brief A class used for initializing the reflection coefficient: R4.
 */
class CRatioR4Reflection : public CReflection {

	public:
		/*!
		 * @brief Default constructor of CRatioR4Reflection, which initializes the reflection coefficient: R4.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] iZone input zone ID.
		 */
		CRatioR4Reflection(CConfig       *config_container,
				               CGeometry     *geometry_container,
								       CElement      *element_container,
								       CInitial      *initial_container,
								       CSolver       *solver_container,
								       unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CRatioR4Reflection(void) override;

		/*!
		 * @brief Function which computes the reflection coefficient: R4.
		 * R_4(t) = max( L(-) / L(+) ).
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] solver_container pointer to current zone solver container.
		 *
		 * @return reflection coefficient: R4
		 */
		as3double ComputeCoefficient(CConfig   *config_container,
		                             CGeometry *geometry_container,
		                             CElement  *element_container,
											           CSolver   *solver_container);
	
	private:
		unsigned short iBoundary;  ///< Current boundary ID used for processing.
};


/*!
 * @brief A class used for initializing the reflection coefficient: R5.
 */
class CRatioR5Reflection : public CReflection {

	public:
		/*!
		 * @brief Default constructor of CRatioR5Reflection, which initializes the reflection coefficient: R5.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] initial_container pointer to current zone initial conditions container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] iZone input zone ID.
		 */
		CRatioR5Reflection(CConfig       *config_container,
				               CGeometry     *geometry_container,
								       CElement      *element_container,
								       CInitial      *initial_container,
								       CSolver       *solver_container,
								       unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CRatioR5Reflection(void) override;

		/*!
		 * @brief Function which computes the reflection coefficient: R5.
		 * R_5(t) = max( L(s) ).
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] solver_container pointer to current zone solver container.
		 *
		 * @return reflection coefficient: R5
		 */
		as3double ComputeCoefficient(CConfig   *config_container,
		                             CGeometry *geometry_container,
		                             CElement  *element_container,
											           CSolver   *solver_container);
	
	private:
		unsigned short iBoundary;  ///< Current boundary ID used for processing. 
		
};
