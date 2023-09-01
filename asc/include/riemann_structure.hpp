#pragma once

/*!
 * @file riemann_structure.hpp
 * @brief The file containing all the Riemann solver implementation.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"


/*!
 * @brief An interface class used for initializing a generic Riemann class.
 */
class CRiemann {

  public:
 		/*!
		 * @brief Default constructor of CRiemann, which initializes a generic Riemann class.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 */
		CRiemann(CConfig *config_container);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */      
		virtual ~CRiemann(void);

		/*!
		 * @brief Function that determines the unique state on a boundary via upwinding.
		 *
		 * @param[in] UnitNormal reference to the unit-vector on this boundary.
		 * @param[in] weights reference to the integration weights on this boundary.
		 * @param[in] hElem reference to the element size at this boundary.
		 * @param[in] Velocity relevant convective velocity (scalar) needed for the upwinding. 
		 * @param[in] VarI pointer to the solution at integration points (1D) from the current element. 
		 * @param[in] VarJ pointer to the solution at integration points (1D) from the neighboring element.
		 * @param[out] VarF pointer to the unique flux computed using a purely upwinding method.
		 */
		void ComputeVariableStateUpwinding(const as3vector1d<as3double>  &UnitNormal,
		                                   const as3vector1d<as3double>  &weights,
		                                   const as3vector1d<as3double>  &hElem,
		                                   as3double                     Velocity,
		                                   as3double                   **VarI,
		                                   as3double                   **VarJ,
		                                   as3double                   **VarF);

    /*!
		 * @brief Pure virtual function that determines the unique state of the flux at a face.
		 * Note, must be overriden by a derived class.
		 *
		 * @param[in] UnitNormal reference to the unit-vector on this boundary.
		 * @param[in] weights reference to the integration weights on this boundary.
		 * @param[in] hElem reference to the element size at this boundary.
		 * @param[in] VarI pointer to the solution at integration points (1D) from the current element. 
		 * @param[in] VarJ pointer to the solution at integration points (1D) from the neighboring element.
		 * @param[out] VarF pointer to the unique flux computed using a Riemann solver.
		 */
    virtual void ComputeFluxState(const as3vector1d<as3double>  &UnitNormal,
                                  const as3vector1d<as3double>  &weights,
                                  const as3vector1d<as3double>  &hElem,
                                  as3double                    **VarI,
                                  as3double                    **VarJ,
                                  as3double                    **Flux) = 0;

  protected:
    as3double gm1;  ///< Abbreviation: gamma minus one.

  private:

};


/*!
 * @brief A class used for initializing Roe's Riemann class.
 */
class CRoeRiemann : public CRiemann {

  public:
 		/*!
		 * @brief Default constructor of CRoeRiemann, which initializes a Roe Riemann solver.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 */
		CRoeRiemann(CConfig *config_container);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */    
		~CRoeRiemann(void) final;

  protected:
    /*!
		 * @brief Function that determines the unique state of the flux at a face using Roe's Riemann solver.
		 *
		 * @param[in] UnitNormal reference to the unit-vector on this boundary.
		 * @param[in] weights reference to the integration weights on this boundary.
		 * @param[in] hElem reference to the element size at this boundary.
		 * @param[in] VarI pointer to the solution at integration points (1D) from the current element. 
		 * @param[in] VarJ pointer to the solution at integration points (1D) from the neighboring element.
		 * @param[out] VarF pointer to the unique flux computed.
		 */
    void ComputeFluxState(const as3vector1d<as3double>  &UnitNormal,
                          const as3vector1d<as3double>  &weights,
                          const as3vector1d<as3double>  &hElem,
                          as3double                    **VarI,
                          as3double                    **VarJ,
                          as3double                    **Flux) final;
  private:
    as3double Delta;  ///< Entropy fix constant.
};


/*!
 * @brief A class used for initializing Rusanov's Riemann class.
 */
class CRusanovRiemann : public CRiemann {

  public:
 		/*!
		 * @brief Default constructor of CRusanovRiemann, which initializes a Rusanov Riemann solver.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 */
    CRusanovRiemann(CConfig *config_container);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */ 
		~CRusanovRiemann(void) final;

  protected:
    /*!
		 * @brief Function that determines the unique state of the flux at a face using Rusanov's Riemann solver.
		 *
		 * @param[in] UnitNormal reference to the unit-vector on this boundary.
		 * @param[in] weights reference to the integration weights on this boundary.
		 * @param[in] hElem reference to the element size at this boundary.
		 * @param[in] VarI pointer to the solution at integration points (1D) from the current element. 
		 * @param[in] VarJ pointer to the solution at integration points (1D) from the neighboring element.
		 * @param[out] VarF pointer to the unique flux computed.
		 */   
		void ComputeFluxState(const as3vector1d<as3double>  &UnitNormal,
                          const as3vector1d<as3double>  &weights,
                          const as3vector1d<as3double>  &hElem,
                          as3double                    **VarI,
                          as3double                    **VarJ,
                          as3double                    **Flux) final;
  private:

};


/*!
 * @brief A class used for initializing Ismail-Roe's Riemann class.
 */
class CRoeIsmailRiemann : public CRiemann {

  public:
 		/*!
		 * @brief Default constructor of CRoeIsmailRiemann, which initializes an Ismail-Roe Riemann solver.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 */   
		CRoeIsmailRiemann(CConfig *config_container);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */   
		~CRoeIsmailRiemann(void) final;

  protected:
    /*!
		 * @brief Function that determines the unique state of the flux at a face using Ismail-Roe's Riemann solver.
		 *
		 * @param[in] UnitNormal reference to the unit-vector on this boundary.
		 * @param[in] weights reference to the integration weights on this boundary.
		 * @param[in] hElem reference to the element size at this boundary.
		 * @param[in] VarI pointer to the solution at integration points (1D) from the current element. 
		 * @param[in] VarJ pointer to the solution at integration points (1D) from the neighboring element.
		 * @param[out] VarF pointer to the unique flux computed.
		 */   
		void ComputeFluxState(const as3vector1d<as3double>  &UnitNormal,
                          const as3vector1d<as3double>  &weights,
                          const as3vector1d<as3double>  &hElem,
                          as3double                    **VarI,
                          as3double                    **VarJ,
                          as3double                    **Flux) final;
  private:
    // Abbreviations involving gamma.
    as3double ovgm1;    ///< Abbreviation: 1/(gamma-1).
    as3double gp1Ovg;   ///< Abbreviation: (gamma+1)/gamma.
    as3double gm1Ovg;   ///< Abbreviation: (gamma-1)/gamma.
    as3double beta;     ///< A scaling coefficient for the acoustic eigenvalue. 
    as3double alphaMax; ///< A scaling coefficient for the acoustic eigenvalue
};




