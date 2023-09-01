#pragma once

/*!
 * @file initial_structure.hpp
 * @brief The file containing all the initial condition specifications.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"


/*!
 * @brief An interface class used for initializing a solution.
 */
class CInitial {

	public:
		/*!
		 * @brief Default constructor of CInitial, which initializes a generic initial condition.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */
		CInitial(CConfig    	 *config_container,
						 CGeometry  	 *geometry_container,
						 CElement   	**element_container,
						 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~CInitial(void);

		/*!
		 * @brief Pure virtual function that sets the initial condition.
		 * Must be overriden by one of the derived classes.
		 *
		 * @param[in] grid_nodes reference to the input coordinates.
		 * @param[in] data_nodes reference to the input (conservative) solution data.
		 * @param[in] nNode number of points, size of grid_nodes and data_nodes.
		 * @param[in] time input physical time.
		 */
		virtual void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                                     as3data1d<as3double>       &data_nodes,
                                     unsigned short              nNode,
                                     as3double                   time) = 0;

    /*!
		 * @brief Pure virtual function that computes the target state in primitive form.
		 * Must be overriden by one of the derived classes.
		 *
		 * @param[in] Coords reference to the input coordinates.
		 * @param[out] TargetState reference to the target (primitive) solution.
		 * @param[in] nData input number of points, size of TargetState.
		 */
    virtual void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                             as3data1d<as3double>       &TargetState,
                                             unsigned short              nData) = 0;

		/*!
		 * @brief Pure virtual function that computes the target state in primitive form on a single DOFs. 
		 * Must be overriden by one of the derived classes.
		 *
		 * @param[in] xy reference to the input coordinates.
		 * @param[out] Qinf reference to the computed target state in primitive form on a single point.
		 */
		virtual void ComputeTargetStatePrimitivePerDOF(const as3vector1d<as3double> &xy,
				                                           as3vector1d<as3double>       &Qinf) = 0;

    /*!
		 * @brief Getter function which returns the value of betaPML.
		 *
		 * @return betaPML
		 */
    as3double GetBetaPML(void)                       const {return betaPML;}
    /*! 
		 * @brief Getter function which returns the value of kappax.
		 *
		 * @return 0.0
		 */
    virtual as3double GetKappax(void)                const {return 0.0;}
    /*!
		 * @brief Getter function which returns the value of kappat.
		 *
		 * @return 0.0
		 */
    virtual as3double GetKappat(void)                const {return 0.0;}
    /*!
		 * @brief Pure virtual getter function which returns the value of aInf. 
		 * Must be implemented in a derived class.
		 */
    virtual as3double GetReferenceSpeedOfSound(void) const = 0;
    /*!
		 * @brief Pure virtual getter function which returns the value of Tinf. 
		 * Must be implemented in a derived class.
		 */
    virtual as3double GetTinf(void)                  const = 0;
    /*!
		 * @brief Pure virtual getter function which returns the value of pInf. 
		 * Must be implemented in a derived class.
		 */
    virtual as3double GetPinf(void)                  const = 0;
    /*!
		 * @brief Pure virtual getter function which returns the value of A0. 
		 * Must be implemented in a derived class.
		 */
    virtual as3double GetA0(void)                    const = 0;
    /*!
		 * @brief Pure virtual getter function which returns the value of Mach. 
		 * Must be implemented in a derived class.
		 */
    virtual as3double GetMach(void)                  const = 0;
    /*!
		 * @brief Pure virtual getter function which returns the value of t0. 
		 * Must be implemented in a derived class.
		 */
    virtual as3double Gett0(void)                    const = 0;
    /*!
		 * @brief Pure virtual getter function which returns true/false. 
		 * Must be implemented in a derived class.
		 */
    virtual bool GetDimensionalProblem(void)         const = 0;
		/*!
		 * @brief Pure virtual getter function which returns the value of uInf. 
		 * Must be implemented in a derived class.
		 */
    virtual as3double GetUinf(void)                  const = 0;
		/*!
		 * @brief Pure virtual getter function which returns the value of vInf. 
		 * Must be implemented in a derived class.
		 */
    virtual as3double GetVinf(void)                  const = 0;

	protected:
		unsigned short zoneID;  ///< Current zone ID.
    as3double      betaPML; ///< Dispersion-relaxation correction parameter. Used in PML formulation.

	private:

};


/*!
 * @brief A class used for initializing a solution based on a Gaussian pressure pulse.
 */
class CGaussianInitial : public CInitial {

	public:
		/*!
		 * @brief Default constructor of CGaussianInitial, which initializes a Gaussian pressure pulse IC.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */	
		CGaussianInitial(CConfig    	 *config_container,
						 				 CGeometry  	 *geometry_container,
						 				 CElement   	**element_container,
						 				 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CGaussianInitial(void) final;

		/*!
		 * @brief Function that sets a Gaussian initial condition.
		 *
		 * @param[in] grid_nodes reference to the input coordinates.
		 * @param[in] data_nodes reference to the input (conservative) solution data.
		 * @param[in] nNode number of points, size of grid_nodes and data_nodes.
		 * @param[in] time input physical time.
		 */	
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    /*!
		 * @brief Function that computes the target state in primitive form.
		 *
		 * @param[in] Coords reference to the input coordinates.
		 * @param[out] TargetState reference to the target (primitive) solution.
		 * @param[in] nData input number of points, size of TargetState.
		 */
		void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

		/*!
		 * @brief Function that computes the target state in primitive form on a single DOFs.
		 *
		 * @param[in] xy reference to the input coordinates.
		 * @param[out] Qinf reference to the computed target state in primitive form on a single point.
		 */
		void ComputeTargetStatePrimitivePerDOF(const as3vector1d<as3double> &xy,
		                                       as3vector1d<as3double>       &Qinf);

    /*!
		 * @brief Getter function which returns the value of aInf.
		 *
		 * @return aInf
		 */
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    /*!
		 * @brief Getter function which returns the value of Tinf.
		 *
		 * @return Tinf
		 */
    as3double GetTinf(void)                  const final {return Tinf;}
    /*!
		 * @brief Getter function which returns the value of pInf.
		 *
		 * @return pInf
		 */
    as3double GetPinf(void)                  const final {return pInf;}
    /*!
		 * @brief Getter function which returns the value of A0.
		 *
		 * @return A0
		 */
    as3double GetA0(void)                    const final {return A0;}
    /*!
		 * @brief Getter function which returns the value of Mach.
		 *
		 * @return Mach
		 */
    as3double GetMach(void)                  const final {return Mach;}
    /*!
		 * @brief Getter function which returns the value of t0. 
		 * Note, this is not needed in this IC.
		 *
		 * @return 0.0
		 */
    as3double Gett0(void)                    const final {return 0.0;}
    /*!
		 * @bried Getter function which returns the value true.
		 *
		 * @return true
		 */
    bool GetDimensionalProblem(void)         const final {return true;}
		/*!
		 * @brief Getter function which returns the value of uInf.
		 *
		 * @return uInf
		 */
    as3double GetUinf(void)                  const final {return uInf;}
		/*!
		 * @brief Getter function which returns the value of vInf.
		 *
		 * @return vInf
		 */
    as3double GetVinf(void)                  const final {return vInf;}

	protected:

	private:
		as3double sigmax;       ///< x-dimension constant of the pulse.
		as3double sigmay;       ///< y-dimension constant of the pulse.
		as3double x0;           ///< Pulse x-coordinate center location.
		as3double y0;           ///< Pulse y-coordinate center location.
		as3double A0;           ///< Pulse strength.
		as3double b;            ///< Pulse width.      
		as3double kappa;        ///< Pulse attenuation.
    as3double Mach;         ///< Background Mach number.
    as3double theta;        ///< Background flow direction [degrees].
    as3double VelInf;       ///< Free-stream velocity magnitude.
		as3double uInf;         ///< Free-stream u-velocity.
		as3double vInf;         ///< Free-stream v-velocity.
    as3double pInf;         ///< Free-stream pressure. 
    as3double Tinf;         ///< Free-stream temperature.
		as3double rhoInf;       ///< Free-stream density.
    as3double aInf;         ///< Free-stream speed of sound.
		bool      OneWayPulse;  ///< Option whether this is a one-way 1D pulse.
};


/*!
 * @brief A class used for initializing a solution based on an isentropic vortex.
 */
class CIsentropicVortexInitial : public CInitial {

	public:
		/*!
		 * @brief Default constructor of CIsentropicVortexInitial, which initializes an isentropic vortex IC.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */	
		CIsentropicVortexInitial(CConfig    	 *config_container,
        						 				 CGeometry  	 *geometry_container,
        						 				 CElement   	**element_container,
        						 				 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CIsentropicVortexInitial(void) final;

		/*!
		 * @brief Function that sets an isentropic vortex initial condition.
		 *
		 * @param[in] grid_nodes reference to the input coordinates.
		 * @param[in] data_nodes reference to the input (conservative) solution data.
		 * @param[in] nNode number of points, size of grid_nodes and data_nodes.
		 * @param[in] time input physical time.
		 */
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    /*!
		 * @brief Function that computes the target state in primitive form.
		 *
		 * @param[in] Coords reference to the input coordinates.
		 * @param[out] TargetState reference to the target (primitive) solution.
		 * @param[in] nData input number of points, size of TargetState.
		 */   
		void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

		/*!
		 * @brief Function that computes the target state in primitive form on a single DOFs.
		 *
		 * @param[in] xy reference to the input coordinates.
		 * @param[out] Qinf reference to the computed target state in primitive form on a single point.
		 */	
		void ComputeTargetStatePrimitivePerDOF(const as3vector1d<as3double> &xy,
		                                       as3vector1d<as3double>       &Qinf);

    /*!
		 * @brief Getter function which returns the velue of aInf.
		 *
		 * @return aInf
		 */
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    /*!
		 * @brief Getter function which returns the value of Tinf.
		 *
		 * @return Tinf
		 */
    as3double GetTinf(void)                  const final {return Tinf;}
    /*!
		 * @brief Getter function which returns the value of pInf.
		 *
		 * @return pInf
		 */
    as3double GetPinf(void)                  const final {return pInf;}
    /*!
		 * @brief Getter function which returns the value of A0.
		 *
		 * @return A0
		 */
    as3double GetA0(void)                    const final {return A0;}
    /*!
		 * @brief Getter function which returns the value of Mach.
		 *
		 * @return Mach
		 */
    as3double GetMach(void)                  const final {return Mach;}
    /*!
		 * @brief Getter function which returns the value of t0. 
		 * Note, this is not needed in this IC.
		 *
		 * @return 0.0
		 */
    as3double Gett0(void)                    const final {return 0.0;}
    /*!
		 * @brief Getter function which returns the value of DimensionalProblem.
		 *
		 * @return true
		 */
    bool GetDimensionalProblem(void)         const final {return true;}
		/*!
		 * @brief Getter function which returns the value of uInf.
		 *
		 * @return uInf
		 */
    as3double GetUinf(void)                  const final {return uInf;}
		/*!
		 * @brief Getter function which returns the value of vInf.
		 *
		 * @return vInf
		 */
    as3double GetVinf(void)                  const final {return vInf;}

	protected:

	private:
		as3double x0;     ///< Vortex center x-coordinate.
		as3double y0;     ///< Vortex center y-coordinate.
		as3double A0;     ///< Vortex strength.
		as3double Rv;     ///< Vortex radius.
    as3double Mach;   ///< Background Mach number.
    as3double theta;  ///< Background flow direction [degrees].
    as3double VelInf; ///< Free-stream velocity magnitude.
		as3double uInf;   ///< Free-stream u-velocity.
		as3double vInf;   ///< Free-stream v-velocity.
    as3double pInf;   ///< Free-stream pressure.
    as3double Tinf;   ///< Free-stream temperature.
		as3double rhoInf; ///< Free-stream density.
    as3double aInf;   ///< Free-stream speed of sound.
};


/*!
 * @brief A class used for initializing a solution based on an entropy wave.
 */
class CEntropyWave : public CInitial {

	public:
		/*!
		 * @brief Default constructor of CEntropyWave, which initializes an entropy wave IC.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */	
		CEntropyWave(CConfig    	 *config_container,
				 				 CGeometry  	 *geometry_container,
				 				 CElement   	**element_container,
				 				 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CEntropyWave(void) final;

		/*!
		 * @brief Function that sets an entropy wave initial condition.
		 *
		 * @param[in] grid_nodes reference to the input coordinates.
		 * @param[in] data_nodes reference to the input (conservative) solution data.
		 * @param[in] nNode number of points, size of grid_nodes and data_nodes.
		 * @param[in] time input physical time.
		 */
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    /*!
		 * @brief Function that computes the target state in primitive form.
		 *
		 * @param[in] Coords reference to the input coordinates.
		 * @param[out] TargetState reference to the target (primitive) solution.
		 * @param[in] nData input number of points, size of TargetState.
		 */      
		void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

		/*!
		 * @brief Function that computes the target state in primitive form on a single DOFs.
		 *
		 * @param[in] xy reference to the input coordinates.
		 * @param[out] Qinf reference to the computed target state in primitive form on a single point.
		 */		
		void ComputeTargetStatePrimitivePerDOF(const as3vector1d<as3double> &xy,
		                                       as3vector1d<as3double>       &Qinf);

    /*!
		 * @brief Getter function which returns the value of aInf.
		 *
		 * @return aInf
		 */
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    /*!
		 * @brief Getter function which returns the value of Tinf.
		 *
		 * @return Tinf
		 */
    as3double GetTinf(void)                  const final {return Tinf;}
    /*!
		 * @brief Getter function which returns the value of pInf.
		 *
		 * @return pInf
		 */
    as3double GetPinf(void)                  const final {return pInf;}
    /*!
		 * @brief Getter function which returns the value of A0.
		 *
		 * @return A0
		 */
    as3double GetA0(void)                    const final {return A0;}
    /*!
		 * @brief Getter function which returns the value of Mach.
		 *
		 * @return Mach
		 */
    as3double GetMach(void)                  const final {return Mach;}
    /*!
		 * @brief Getter function which returns the value of t0. 
		 * Note, this is not needed in this IC.
		 *
		 * @return 0.0
		 */
    as3double Gett0(void)                    const final {return 0.0;}
    /*!
		 * @brief Getter function which returns the value of DimensionalProblem.
		 *
		 * @return true
		 */
    bool GetDimensionalProblem(void)         const final {return true;}
		/*!
		 * @brief Getter function which returns the value of uInf.
		 *
		 * @return uInf
		 */
    as3double GetUinf(void)                  const final {return uInf;}
		/*!
		 * @brief Getter function which returns the value of vInf.
		 *
		 * @return vInf
		 */
    as3double GetVinf(void)                  const final {return vInf;}

	protected:

	private:
    as3double x0;     ///< Entropy wave center x-coordinate.
		as3double y0;     ///< Entropy wave center y-coordinate.
		as3double A0;     ///< Gaussian strength.
		as3double b;      ///< Gaussian width.
		as3double kappa;  ///< Gaussian attenuation.
    as3double Mach;   ///< Background Mach number.
    as3double theta;  ///< Background flow direction [degrees].
    as3double VelInf; ///< Free-stream velocity magnitude.
		as3double uInf;   ///< Free-stream u-velocity.
		as3double vInf;   ///< Free-stream v-velocity.
    as3double pInf;   ///< Free-stream pressure.
    as3double Tinf;   ///< Free-stream temperature.
		as3double rhoInf; ///< Free-stream density.
    as3double aInf;   ///< Free-stream speed of sound.
};


/*!
 * @brief A class used for initializing a solution based on a shear-flow vortex roll-up.
 */
class CVortexRollup : public CInitial {

	public:
		/*!
		 * @brief Default constructor of CVortexRollup, which initializes a vortex roll-up IC.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */	
		CVortexRollup(CConfig    	  *config_container,
				 				  CGeometry  	  *geometry_container,
				 				  CElement     **element_container,
				 				  unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CVortexRollup(void) final;

		/*!
		 * @brief Function that sets a vortex roll-up initial condition.	
		 *
		 * @param[in] grid_nodes reference to the input coordinates.
		 * @param[in] data_nodes reference to the input (conservative) solution data.
		 * @param[in] nNode number of points, size of grid_nodes and data_nodes.
		 * @param[in] time input physical time.
		 */
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    /*!
		 * @brief Function that computes the target state in primitive form.
		 *
		 * @param[in] Coords reference to the input coordinates.
		 * @param[out] TargetState reference to the target (primitive) solution.
		 * @param[in] nData input number of points, size of TargetState.
		 */      
		void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

		/*!
		 * @brief Function that computes the target state in primitive form on a single DOFs.
		 *
		 * @param[in] xy reference to the input coordinates.
		 * @param[out] Qinf reference to the computed target state in primitive form on a single point.
		 */			
		void ComputeTargetStatePrimitivePerDOF(const as3vector1d<as3double> &xy,
		                                       as3vector1d<as3double>       &Qinf);

    /*!
		 * @brief Getter function which returns the value of aInf. 
		 * Note, this is based on the (1) state.
		 *
		 * @return aInf
		 */
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    /*!
		 * @brief Getter function which returns the value of Tinf. 
		 * Note, this is based on the (1) state.
		 *
		 * @return Tinf
		 */
    as3double GetTinf(void)                  const final {return Tinf;}
    /*!
		 * @brief Getter function which returns the value of pInf. 
		 * Note, this is based on the (1) state.
		 *
		 * @return pInf
		 */
    as3double GetPinf(void)                  const final {return pInf;}
    /*!
		 * @brief Getter function which returns the value of A0. 
		 * Note, it does not exist in this IC.
		 *
		 * @return A0
		 */
    as3double GetA0(void)                    const final {return A0;}
    /*!
		 * @brief Getter function which returns the value of Mach. 
		 * Note, this is based on the (1) state.
		 *
		 * @return Mach
		 */
    as3double GetMach(void)                  const final {return Mach;}
    /*!
		 * @brief Getter function which returns the value of t0. 
		 * Note, this is not needed in this IC.
		 *
		 * @return 0.0
		 */
    as3double Gett0(void)                    const final {return 0.0;}
    /*!
		 * @brief Getter function which returns the value of DimensionalProblem.
		 *
		 * @return false
		 */
    bool GetDimensionalProblem(void)         const final {return false;}
		/*!
		 * @brief Getter function which returns the value of uInf.
		 *
		 * @return uInf
		 */
    as3double GetUinf(void)                  const final {return uInf;}
		/*!
		 * @brief Getter function which returns the value of vInf.
		 *
		 * @return vInf
		 */
    as3double GetVinf(void)                  const final {return vInf;}

	protected:

	private:
    // NOTE, this is a non-dimensional problem, as defined in the paper by Hu.
    // If there need be dimensionalization, then the dispersion-relation needs
    // to be taken care of and not use the 1/4 ratio defined by Hu.

    as3double delta;  ///< Flow width around shear layer center..
		as3double A0;     ///< Gaussian strength.
    as3double T1;     ///< Temperature on top.
    as3double T2;     ///< Temperature on bottom.
    as3double U1;     ///< (u-)Velocity on top.
    as3double U2;     ///< (u-)Velocity on bottom.
    as3double Mach;   ///< Background Mach number (not used in this IC).
    as3double theta;  ///< Free-stream flow direction [degrees] (not used in this IC).
    as3double VelInf; ///< Free-stream velocity magnitude.
    as3double Uavg;   ///< Free-stream average of u-velocity.
    as3double Ujmp;   ///< Free-stream jump of u-velocity.
		as3double uInf;   ///< Free-stream based on (1)-state of u-velocity.
		as3double vInf;   ///< Free-stream based on (1)-state of v-velocity.
    as3double pInf;   ///< Free-stream based on (1)-state of pressure.
    as3double Tinf;   ///< Free-stream based on (1)-state of temperature.
		as3double rhoInf; ///< Free-stream based on (1)-state of density.
    as3double aInf;   ///< Free-stream based on (1)-state of speed of sound.

		// Some abbreviations.
		as3double tovd;   ///< Abbreviation: 2/delta.
		as3double ovujmp; ///< Abbreviation: 1/Ujmp.
		as3double ovgm1;  ///< Abbreviation: 1/(gamma-1).
		as3double gm1ov2; ///< Abbreviation: (gamma-1)/2.
};


/*!
 * @brief A class used for initializing a solution based on an acoustic plane.
 */
class CAcousticPlane : public CInitial {

	public:
		/*!
		 * @brief Default constructor of CAcousticPlane, which initializes an acoustic pulse IC.
		 *
		 * @param[in] config_container pointer to input configuration container.
		 * @param[in] geometry_container pointer to input geometry container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone input zone ID.
		 */	
		CAcousticPlane(CConfig    	 *config_container,
					 				 CGeometry  	 *geometry_container,
					 				 CElement     **element_container,
					 				 unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CAcousticPlane(void) final;

		/*!
		 * @brief Function that sets an acoustic Gaussian plane pulse initial condition.
		 *
		 * @param[in] grid_nodes reference to the input coordinates.
		 * @param[in] data_nodes reference to the input (conservative) solution data.
		 * @param[in] nNode number of points, size of grid_nodes and data_nodes.
		 * @param[in] time input physical time.
		 */
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    /*!
		 * @brief Function that computes the target state in primitive form.
		 *
		 * @param[in] Coords reference to the input coordinates.
		 * @param[out] TargetState reference to the target (primitive) solution.
		 * @param[in] nData input number of points, size of TargetState.
		 */      
		void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

		/*!
		 * @brief Function that computes the target state in primitive form on a single DOFs.
		 *
		 * @param[in] xy reference to the input coordinates.
		 * @param[out] Qinf reference to the computed target state in primitive form on a single point.
		 */			
		void ComputeTargetStatePrimitivePerDOF(const as3vector1d<as3double> &xy,
		                                       as3vector1d<as3double>       &Qinf);

    /*!
		 * @brief Getter function which returns the value of aInf.
		 *
		 * @return aInf
		 */
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    /*!
		 * @brief Getter function which returns the value of Tinf.
		 *
		 * @return Tinf
		 */
    as3double GetTinf(void)                  const final {return Tinf;}
    /*!
		 * @brief Getter function which returns the value of pInf.
		 *
		 * @return pInf
		 */
    as3double GetPinf(void)                  const final {return pInf;}
    /*!
		 * @brief Getter function which returns the value of A0.
		 *
		 * @return A0
		 */
    as3double GetA0(void)                    const final {return A0;}
    /*!
		 * @brief Getter function which returns the value of Mach.
		 *
		 * @return Mach
		 */
    as3double GetMach(void)                  const final {return Mach;}
    /*!
		 * @brief Getter function which returns the value of t0.
		 *
		 * @return t0
		 */
    as3double Gett0(void)                    const final {return t0;}
    /*!
		 * @brief Getter function which returns the value of kappax.
		 *
		 * @return kappax
		 */
    as3double GetKappax(void)                const final {return kappax;}
    /*!
		 * @brief Getter function which returns the value of kappat.
		 *
		 * @return kappat
		 */
    as3double GetKappat(void)                const final {return kappat;}
    /*!
		 * @brief Getter function which returns the value of DimensionalProblem.
		 *
		 * @return true
		 */
    bool GetDimensionalProblem(void)         const final {return true;}
		/*!
		 * @brief Getter function which returns the value of uInf.
		 *
		 * @return uInf
		 */
    as3double GetUinf(void)                  const final {return uInf;}
		/*!
		 * @brief Getter function which returns the value of vInf.
		 *
		 * @return vInf
		 */
    as3double GetVinf(void)                  const final {return vInf;}

	protected:

	private:
		as3double x0;     ///< Pulse center x-coordinate.
		as3double y0;     ///< Pulse center y-coordinate.
    as3double t0;     ///< Time lag parameter.
    as3double pa;     ///< Pressure disturbance coefficient.
    as3double st;     ///< Gaussian temporal width.
    as3double sx;     ///< Gaussian spatial width.
		as3double A0;     ///< Gaussian amplitude coefficient.
		as3double kappax; ///< Pulse spatial attenuation.
		as3double kappat; ///< Pulse temporal attenuation.
    as3double Mach;   ///< Background Mach number.
    as3double theta;  ///< Background flow direction [degrees].
    as3double VelInf; ///< Free-stream velocity magnitude.
		as3double uInf;   ///< Free-stream u-velocity.
		as3double vInf;   ///< Free-stream v-velocity.
    as3double pInf;   ///< Free-stream pressure.
    as3double Tinf;   ///< Free-stream temperature.
		as3double rhoInf; ///< Free-stream density.
    as3double aInf;   ///< Free-stream speed of sound.
};
