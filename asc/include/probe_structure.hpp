#pragma once

/*!
 * @file probe_structure.hpp
 * @brief The file responsible for defining probes.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "initial_structure.hpp"


/*!
 * @brief A class used for initializing all probe sensors.
 */
class CProbe {

  public:
		/*!
		 * @brief Default constructor of CProbe, which initializes the probe class.
		 *
		 * @param[in] config input configuration/dictionary container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] probe probe coordinates.
		 * @param[in] vars type of variables to sample via probe.
		 * @param[in] iZone input zone ID.
		 */
    CProbe(CConfig                    *config_container,
           CGeometry                  *geometry_container,
           CElement                   *element_container,
           as3vector1d<as3double>      probe,
					 as3vector1d<unsigned short> vars,
           unsigned short              iZone);
    
 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */   
		~CProbe(void);

    /*!
		 * @brief Getter function which returns the value of ProbeLocation.
		 *
		 * @return ProbeLocation
		 */
    const as3vector1d<as3double> &GetProbeLocation(void) const {return ProbeLocation;}
    /*!
		 * @brief Getter function which returns the value of Interpolation.
		 *
		 * @return Interpolation
		 */
    const as3vector1d<as3double> &GetInterpolation(void) const {return Interpolation;}

    /*!
		 * @brief Function that samples the current data on this probe location.
		 *
		 * @param[in] config input configuration/dictionary container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to current zone standard element container.
		 * @param[in] solver_container pointer to current zone solver container.
		 * @param[in] initial_container pointer to current zone initial condition container.
		 *
		 * @return the interpolated variables at the probe location
		 */
    as3vector1d<as3double> SampleProbeData(CConfig   *config_container,
                                           CGeometry *geometry_container,
                                           CElement  *element_container,
                                           CSolver   *solver_container,
                                           CInitial  *initial_container);


  private:
    unsigned long          iElemProbe;      ///< Global element index of the probe ownership.
    as3vector1d<as3double> ProbeLocation;   ///< Probe coordinates.
    as3vector1d<as3double> Interpolation;   ///< Lagrange interpolation over this element.

		as3vector1d<as3double (*)(as3vector1d<as3double> &)> variables; ///< Vector of function pointers.
																																		///< This is used to obtain different variables. 

    /*!
		 * @brief Function that computes the Lagrange interpolation entry in 1D for a given basis.
		 *
		 * @param[in] xPoint evaluation point.
		 * @param[in] rBasis reference to the polynomial basis.
		 * @param[in] iDegree degree of the Lagrange polynomial.
		 *
		 * @return the evaluation of a Lagrange polynomial at the evaluated point. 
		 */
    as3double EvaluateEntryLagrangePolynomial(const as3double               xPoint,
                                              const as3vector1d<as3double> &rBasis,
                                              const unsigned short          iDegree);

		/*!
		 * @brief Function which extracts the density from a conservative vector.
		 *
		 * @param[in] U reference to a conservative variables vector.
		 *
		 * @return density
		 */
		static inline as3double GetValueDensity(as3vector1d<as3double> &U)
		{
			return U[0];
		}
		/*!
		 * @brief Function which extracts the x-momentum from a conservative vector.
		 *
		 * @param[in] U reference to a conservative variables vector.
		 *
		 * @return x-momentum
		 */
		static inline as3double GetValueXMomentum(as3vector1d<as3double> &U)
		{
			return U[1];
		}
		/*!
		 * @brief Function which extracts the y-momentum from a conservative vector.
		 *
		 * @param[in] U reference to a conservative variables vector.
		 *
		 * @return y-momentum
		 */
		static inline as3double GetValueYMomentum(as3vector1d<as3double> &U)
		{
			return U[2];
		}
		/*!
		 * @brief Function which extracts the total energy from a conservative vector.
		 *
		 * @param[in] U reference to a conservative variables vector.
		 *
		 * @return total energy
		 */
		static inline as3double GetValueTotalEnergy(as3vector1d<as3double> &U)
		{
			return U[3];
		}
		/*!
		 * @brief Function which extracts the u-velocity from a conservative vector.
		 *
		 * @param[in] U reference to a conservative variables vector.
		 *
		 * @return u-velocity
		 */
		static inline as3double GetValueXVelocity(as3vector1d<as3double> &U)
		{
			return U[1]/U[0];
		}
		/*!
		 * @brief Function which extracts the v-velocity from a conservative vector.
		 *
		 * @param[in] U reference to a conservative variables vector.
		 *
		 * @return v-velocity
		 */	
		static inline as3double GetValueYVelocity(as3vector1d<as3double> &U)
		{
			return U[2]/U[0];
		}
		/*!
		 * @brief Function which extracts the pressure from a conservative vector.
		 *
		 * @param[in] U reference to a conservative variables vector.
		 *
		 * @return pressure
		 */	
		static inline as3double GetValuePressure(as3vector1d<as3double> &U)
		{
			return GAMMA_MINUS_ONE*( U[3] - 0.5*( U[1]*U[1] + U[2]*U[2] )/U[0] );
		}
};

