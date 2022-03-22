#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"


class CInitial {

	public:
		// Constructor.
		CInitial(CConfig    	 *config_container,
						 CGeometry  	 *geometry_container,
						 CElement   	**element_container,
						 unsigned short iZone);

		// Destructor.
		virtual ~CInitial(void);

		// Pure virtual function that sets the initial condition.
		// Must be overriden by one of the derived classes.
		virtual void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                                     as3data1d<as3double>       &data_nodes,
                                     unsigned short              nNode,
                                     as3double                   time) = 0;

    // Pure virtual function that computes the target state in primitive form.
    // Must be overriden by one of the derived classes.
    virtual void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                             as3data1d<as3double>       &TargetState,
                                             unsigned short              nData) = 0;

    // Getter: returns betaPML.
    as3double GetBetaPML(void)                       const {return betaPML;}
    // Getter: returns kappax.
    virtual as3double GetKappax(void)                const {return 0.0;}
    // Getter: returns kappat.
    virtual as3double GetKappat(void)                const {return 0.0;}
    // Pure virtual getter: returns aInf. Must be implemented in a derived class.
    virtual as3double GetReferenceSpeedOfSound(void) const = 0;
    // Pure virtual getter: returns V0. Must be implemented in a derived class.
    virtual as3double GetVelocityTransverse(void)    const = 0;
    // Pure virtual getter: returns U0. Must be implemented in a derived class.
    virtual as3double GetVelocityNormal(void)        const = 0;
    // Pure virtual getter: returns Tinf. Must be implemented in a derived class.
    virtual as3double GetTinf(void)                  const = 0;
    // Pure virtual getter: returns pInf. Must be implemented in a derived class.
    virtual as3double GetPinf(void)                  const = 0;
    // Pure virtual getter: returns A0. Must be implemented in a derived class.
    virtual as3double GetA0(void)                    const = 0;
    // Pure virtual getter: returns Mach. Must be implemented in a derived class.
    virtual as3double GetMach(void)                  const = 0;
    // Pure virtual getter: returns t0. Must be implemented in a derived class.
    virtual as3double Gett0(void)                    const = 0;
    // Pure virtual getter: returns true/false. Must be implemented in a derived class.
    virtual bool GetDimensionalProblem(void)         const = 0;

	protected:
		// Zone ID.
		unsigned short zoneID;

    // Dispersion-relaxation correction parameter. Used in PML formulation.
    as3double betaPML;

	private:

};


class CGaussianInitial : public CInitial {

	public:
		// Constructor.
		CGaussianInitial(CConfig    	 *config_container,
						 				 CGeometry  	 *geometry_container,
						 				 CElement   	**element_container,
						 				 unsigned short iZone);

		// Destructor.
		~CGaussianInitial(void) final;

		// Function that sets a Gaussian initial condition.
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    // Function that computes the target state in primitive form.
    void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

    // Getter: returns aInf.
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    // Getter: returns V0.
    as3double GetVelocityTransverse(void)    const final {return vInf;}
    // Getter: returns U0.
    as3double GetVelocityNormal(void)        const final {return uInf;}
    // Getter: returns Tinf.
    as3double GetTinf(void)                  const final {return Tinf;}
    // Getter: returns pInf.
    as3double GetPinf(void)                  const final {return pInf;}
    // Getter: returns A0.
    as3double GetA0(void)                    const final {return A0;}
    // Getter: returns Mach.
    as3double GetMach(void)                  const final {return Mach;}
    // Getter: returns t0. Not needed in this IC.
    as3double Gett0(void)                    const final {return 0.0;}
    // Getter: returns true/false.
    bool GetDimensionalProblem(void)         const final {return true;}

	protected:

	private:
		// Pulse center.
		as3double x0;
		as3double y0;
		// Pulse strength.
		as3double A0;
		// Pulse width.
		as3double b;
		// Pulse attenuation.
		as3double kappa;
    // Background Mach number.
    as3double Mach;
    // Flow direction [degrees].
    as3double theta;
    // Velocity magnitude free-stream.
    as3double VelInf;
		// Free-stream values.
		as3double uInf;
		as3double vInf;
    as3double pInf;
    as3double Tinf;
		as3double rhoInf;
    as3double aInf;
};


class CIsentropicVortexInitial : public CInitial {

	public:
		// Constructor.
		CIsentropicVortexInitial(CConfig    	 *config_container,
        						 				 CGeometry  	 *geometry_container,
        						 				 CElement   	**element_container,
        						 				 unsigned short iZone);

		// Destructor.
		~CIsentropicVortexInitial(void) final;

		// Function that sets an isentropic vortex initial condition.
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    // Function that computes the target state in primitive form.
    void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

    // Getter: returns aInf.
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    // Getter: returns V0.
    as3double GetVelocityTransverse(void)    const final {return vInf;}
    // Getter: returns U0.
    as3double GetVelocityNormal(void)        const final {return uInf;}
    // Getter: returns Tinf.
    as3double GetTinf(void)                  const final {return Tinf;}
    // Getter: returns pInf.
    as3double GetPinf(void)                  const final {return pInf;}
    // Getter: returns A0.
    as3double GetA0(void)                    const final {return A0;}
    // Getter: returns Mach.
    as3double GetMach(void)                  const final {return Mach;}
    // Getter: returns t0. Not needed in this IC.
    as3double Gett0(void)                    const final {return 0.0;}
    // Getter: returns DimensionalProblem.
    bool GetDimensionalProblem(void)         const final {return true;}

	protected:

	private:
		// Vortex center.
		as3double x0;
		as3double y0;
		// Vortex strength.
		as3double A0;
		// Vortex radius.
		as3double Rv;
    // Background Mach number.
    as3double Mach;
    // Flow direction [degrees].
    as3double theta;
    // Velocity magnitude free-stream.
    as3double VelInf;
		// Free-stream values.
		as3double uInf;
		as3double vInf;
    as3double pInf;
    as3double Tinf;
		as3double rhoInf;
    as3double aInf;
};


class CEntropyWave : public CInitial {

	public:
		// Constructor.
		CEntropyWave(CConfig    	 *config_container,
				 				 CGeometry  	 *geometry_container,
				 				 CElement   	**element_container,
				 				 unsigned short iZone);

		// Destructor.
		~CEntropyWave(void) final;

		// Function that sets an entropy wave initial condition.
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    // Function that computes the target state in primitive form.
    void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

    // Getter: returns aInf.
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    // Getter: returns V0.
    as3double GetVelocityTransverse(void)    const final {return vInf;}
    // Getter: returns U0.
    as3double GetVelocityNormal(void)        const final {return uInf;}
    // Getter: returns Tinf.
    as3double GetTinf(void)                  const final {return Tinf;}
    // Getter: returns pInf.
    as3double GetPinf(void)                  const final {return pInf;}
    // Getter: returns A0.
    as3double GetA0(void)                    const final {return A0;}
    // Getter: returns Mach.
    as3double GetMach(void)                  const final {return Mach;}
    // Getter: returns t0. Not needed in this IC.
    as3double Gett0(void)                    const final {return 0.0;}
    // Getter: returns DimensionalProblem.
    bool GetDimensionalProblem(void)         const final {return true;}

	protected:

	private:
		// Entropy wave center.
    as3double x0;
		as3double y0;
		// Gaussian strength.
		as3double A0;
		// Gaussian width.
		as3double b;
		// Gaussian attenuation.
		as3double kappa;
    // Background Mach number.
    as3double Mach;
    // Flow direction [degrees].
    as3double theta;
    // Velocity magnitude free-stream.
    as3double VelInf;
		// Free-stream values.
		as3double uInf;
		as3double vInf;
    as3double pInf;
    as3double Tinf;
		as3double rhoInf;
    as3double aInf;
};


class CVortexRollup : public CInitial {

	public:
		// Constructor.
		CVortexRollup(CConfig    	  *config_container,
				 				  CGeometry  	  *geometry_container,
				 				  CElement     **element_container,
				 				  unsigned short iZone);

		// Destructor.
		~CVortexRollup(void) final;

		// Function that sets a vortex roll-up initial condition.
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    // Function that computes the target state in primitive form.
    void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

    // Getter: returns aInf. Note, based on the (1) state.
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    // Getter: returns V0. Note, based on the (1) state.
    as3double GetVelocityTransverse(void)    const final {return vInf;}
    // Getter: returns U0.
    as3double GetVelocityNormal(void)        const final {return uInf;}
    // Getter: returns Tinf. Note, based on the (1) state.
    as3double GetTinf(void)                  const final {return Tinf;}
    // Getter: returns pInf. Note, based on the (1) state.
    as3double GetPinf(void)                  const final {return pInf;}
    // Getter: returns A0. Note, it does not exist in this IC.
    as3double GetA0(void)                    const final {return A0;}
    // Getter: returns Mach. Note, based on the (1) state.
    as3double GetMach(void)                  const final {return Mach;}
    // Getter: returns t0. Not needed in this IC.
    as3double Gett0(void)                    const final {return 0.0;}
    // Getter: returns DimensionalProblem.
    bool GetDimensionalProblem(void)         const final {return false;}

	protected:

	private:
    // NOTE, this is a non-dimensional problem, as defined in the paper by Hu.
    // If there need be dimensionalization, then the dispersion-relation needs
    // to be taken care of and not use the 1/4 ratio defined by Hu.

    // Duct normalizing dimension.
    as3double delta;
		// Gaussian strength.
		as3double A0;
    // Temperature on top.
    as3double T1;
    // Temperature in bottom.
    as3double T2;
    // Velocity on top.
    as3double U1;
    // Velocity in bottom.
    as3double U2;
    // Background Mach number.
    as3double Mach;
    // Flow direction [degrees].
    as3double theta;
    // Velocity magnitude free-stream.
    as3double VelInf;
    // Average of u-velocity.
    as3double Uavg;
    // Jump of u-velocity.
    as3double Ujmp;
    // Free-stream values. These are based off the U1, T1 state.
		as3double uInf;
		as3double vInf;
    as3double pInf;
    as3double Tinf;
		as3double rhoInf;
    as3double aInf;
};


class CAcousticPlane : public CInitial {

	public:
		// Constructor.
		CAcousticPlane(CConfig    	 *config_container,
					 				 CGeometry  	 *geometry_container,
					 				 CElement     **element_container,
					 				 unsigned short iZone);

		// Destructor.
		~CAcousticPlane(void) final;

		// Function that sets an acoustic Gaussian plane pulse initial condition.
		void SetInitialCondition(const as3data1d<as3double> &grid_nodes,
                             as3data1d<as3double>       &data_nodes,
                             unsigned short              nNode,
                             as3double                   time) final;

    // Function that computes the target state in primitive form.
    void ComputeTargetStatePrimitive(const as3data1d<as3double> &Coords,
                                     as3data1d<as3double>       &TargetState,
                                     unsigned short              nData) final;

    // Getter: returns aInf.
    as3double GetReferenceSpeedOfSound(void) const final {return aInf;}
    // Getter: returns V0.
    as3double GetVelocityTransverse(void)    const final {return vInf;}
    // Getter: returns U0.
    as3double GetVelocityNormal(void)        const final {return uInf;}
    // Getter: returns Tinf.
    as3double GetTinf(void)                  const final {return Tinf;}
    // Getter: returns pInf.
    as3double GetPinf(void)                  const final {return pInf;}
    // Getter: returns A0.
    as3double GetA0(void)                    const final {return A0;}
    // Getter: returns Mach.
    as3double GetMach(void)                  const final {return Mach;}
    // Getter: returns t0.
    as3double Gett0(void)                    const final {return t0;}
    // Getter: returns kappax.
    as3double GetKappax(void)                const final {return kappax;}
    // Getter: returns kappat.
    as3double GetKappat(void)                const final {return kappat;}
    // Getter: returns true/false.
    bool GetDimensionalProblem(void)         const final  {return true;}

	protected:

	private:
		// Pulse center.
		as3double x0;
		as3double y0;
    // Time lag.
    as3double t0;
    // Pressure disturbance coefficient.
    as3double pa;
    // Gaussian temporal width.
    as3double st;
		// Gaussian spatial width.
    as3double sx;
    // Gaussian amplitude coefficient.
		as3double A0;
		// Pulse spatial attenuation.
		as3double kappax;
    // Pulse temporal attenuation.
		as3double kappat;
    // Background Mach number.
    as3double Mach;
    // Flow direction [degrees].
    as3double theta;
    // Velocity magnitude free-stream.
    as3double VelInf;
		// Free-stream values.
		as3double uInf;
		as3double vInf;
    as3double pInf;
    as3double Tinf;
		as3double rhoInf;
    as3double aInf;
};