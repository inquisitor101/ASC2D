#pragma once

#include "option_structure.hpp"
#include <map>



class CConfig {

	public:
		// Constructor.
		CConfig(const char *configFile);

		// Destructor.
		~CConfig(void);

    // Getter: returns CFL.
    as3double GetCFL(void)                                                                            const {return CFL;}
    // Getter: returns AdaptTime.
    bool GetAdaptTime(void)                                                                           const {return AdaptTime;}
    // Getter: returns OutputFreq.
    unsigned long GetOutputFreq(void)                                                                 const {return OutputFreq;}
    // Getter: returns WriteFreq.
    unsigned long GetWriteFreq(void)                                                                  const {return WriteFreq;}
    // Getter: returns OutputSolFilename.
		const char *GetOutputSolFilename(void)                                                            const {return OutputSolFilename.c_str();}
		// Getter: returns OutputVTKFilename.
		const char *GetOutputVTKFilename(void)                                                            const {return OutputVTKFilename.c_str();}
    // Getter: returns OutputProcessedFilename.
    const char *GetOutputProcessedFilename(void)                                                      const {return OutputProcessedFilename.c_str();}
    // Getter: returns OutputZoneDataFilename.
    const char *GetOutputZoneDataFilename(void)                                                       const {return OutputZoneDataFilename.c_str();}
    // Getter: returns RestartFilename.
    const char *GetRestartFilename(void)                                                              const {return RestartFilename.c_str();}
    // Getter: returns nZone.
		unsigned short GetnZone(void)                                                                     const {return nZone;}
    // Getter: returns MultizoneStrategy.
    unsigned short GetMultizoneStrategy(void)                                                         const {return MultizoneStrategy;}
    // Getter: returns nPolySolZone, per input zone.
		unsigned short GetnPolySolZone(unsigned short iZone)                                              const {return nPolySolZone[iZone];}
    // Getter: returns nPolySolZone.
		const as3vector1d<unsigned short> &GetnPolySolZone(void)                                          const {return nPolySolZone;}
    // Getter: returns BlockInterfaceLocation.
    const as3vector1d<as3double> &GetBlockInterfaceLocation(void)                                     const {return BlockInterfaceLocation;}
    // Getter: returns DomainExpansionRatio.
    const as3vector1d<as3double> &GetDomainExpansionRatio(void)                                       const {return DomainExpansionRatio;}
    // Getter: returns nxBlockElem.
    const as3vector1d<unsigned long> &GetnxBlockElem(void)                                            const {return nxBlockElem;}
    // Getter: returns nyBlockElem.
    const as3vector1d<unsigned long> &GetnyBlockElem(void)                                            const {return nyBlockElem;}
    // Getter: returns UniformGridResolution.
    bool GetUniformGridResolution(void)                                                               const {return UniformGridResolution;}
    // Getter: returns ZoneMarker, per input zone.
    unsigned short GetTypeZone(unsigned short iZone)                                                  const {return TypeZone[iZone];}
    // Getter: returns ZoneMarker.
		const as3vector1d<unsigned short> &GetTypeZone(void)                                              const {return TypeZone;}
    // Getter: returns TypeDOFs, per input zone.
		unsigned short GetTypeDOFs(unsigned short iZone)                                                  const {return TypeDOFs[iZone];}
    // Getter: returns TypeDOFs.
		const as3vector1d<unsigned short> &GetTypeDOFs(void)                                              const {return TypeDOFs;}
    // Getter: returns TypeSolver, per input zone.
    unsigned short GetTypeSolver(unsigned short iZone)                                                const {return TypeSolver[iZone];}
    // Getter: returns TypeSolver.
    const as3vector1d<unsigned short> &GetTypeSolver(void)                                            const {return TypeSolver;}
    // Getter: returns TypeBufferLayer.
    const as3vector1d<unsigned short> &GetTypeBufferLayer(void)                                       const {return TypeBufferLayer;}
    // Getter: returns TypeBufferLayer, per input zone.
    unsigned short GetTypeBufferLayer(unsigned short iZone)                                           const {return TypeBufferLayer[iZone];}
    // Getter: returns UsePML.
    bool GetUsePML(void)                                                                              const {return UsePML;}
    // Getter: returns WriteAuxiliaryDataPML.
    bool GetWriteAuxiliaryDataPML(void)                                                               const {return WriteAuxiliaryDataPML;}
    // Getter: returns integrationFactor.
		unsigned short GetIntegrationFactor(void)                                                         const {return IntegrationFactor;}
		// Getter: returns SimulationTime[0].
		as3double GetSimulationStartTime(void)                                                            const {return SimulationTime[0];}
		// Getter: returns SimulationTime[1].
		as3double GetSimulationFinalTime(void)                                                            const {return SimulationTime[1];}
		// Getter: returns TypeTemporalScheme.
		unsigned short GetTypeTemporalScheme(void)                                                        const {return TypeTemporalScheme;}
		// Getter: returns TypeIC.
		const as3vector1d<unsigned short> &GetTypeIC(void)                                                const {return TypeIC;}
    // Getter: returns TypeIC, per input zone.
    unsigned short GetTypeIC(unsigned short iZone)                                                    const {return TypeIC[iZone];}
    // Getter: returns FlowAngle.
    as3double GetFlowAngle(void)                                                                      const {return FlowAngle;}
    // Getter: returns MachInf.
    as3double GetMachInf(void)                                                                        const {return MachInf;}
    // Getter: returns CrossFlow.
    bool GetCrossFlow(void)                                                                           const {return CrossFlow;}
    // Getter: returns CenterX0.
    const as3vector1d<as3double> &GetCenterX0(void)                                                   const {return CenterX0;}
    // Getter: returns DisturbanceRatio.
    as3double GetDisturbanceRatio(void)                                                               const {return DisturbanceRatio;}
    // Getter: returns DisturbanceWidth.
    as3double GetDisturbanceWidth(void)                                                               const {return DisturbanceWidth;}
    // Getter: returns AngularFrequency.
    as3double GetAngularFrequency(void)                                                               const {return AngularFrequency;}
    // Getter: returns ConstantFrequency.
    bool GetConstantFrequency(void)                                                                   const {return ConstantFrequency;}
    // Getter: returns TypeModifyBC, per input boundary.
    unsigned short GetTypeModifyBC(unsigned short iBoundary)                                          const {return TypeModifyBC[iBoundary];}
    // Getter: returns ModifyFreqBC.
    const as3vector1d<as3double> GetModifyFreqBC(void)                                                const {return ModifyFreqBC;}
    // Getter: returns ModifyStrengthBC.
    const as3vector1d<as3double> GetModifyStrengthBC(void)                                            const {return ModifyStrengthBC;}
    // Getter: returns ModifyWidthBC.
    const as3vector1d<as3double> GetModifyWidthBC(void)                                               const {return ModifyWidthBC;}
    // Getter: returns ModifyCenterBC.
    const as3vector1d<as3double> GetModifyCenterBC(void)                                              const {return ModifyCenterBC;}
    // Getter: returns ModifyShiftCenterBC.
    const as3vector1d<as3double> GetModifyShiftCenterBC(void)                                         const {return ModifyShiftCenterBC;}
    // Getter: returns RiemannSolver.
    const as3vector1d<unsigned short> &GetRiemannSolver(void)                                         const {return RiemannSolver;}
    // Getter: returns RiemannSolver, per input zone.
    unsigned short GetRiemannSolver(unsigned short iZone)                                             const {return RiemannSolver[iZone];}
    // Getter: returns AlignedPeriodicPulse.
    bool GetAlignedPeriodicPulse(void)                                                                const {return AlignedPeriodicPulse;}
		// Getter: returns SourceFrequencyParam.
		const as3vector1d<as3double> &GetSourceFrequencyParam(void)                                       const {return SourceFrequencyParam;}
		// Getter: returns SourceTermCenterShift.
		const as3vector1d<as3double> &GetSourceTermCenterShift(void)                                      const {return SourceTermCenterShift;}
		// Getter: returns SourceFrequencyExponent.
		as3double GetSourceFrequencyExponent(void)                                                        const {return SourceFrequencyExponent;}
		// Getter: returns SourceTermCenterFixed.
		bool GetSourceTermCenterFixed(void)                                                               const {return SourceTermCenterFixed;}
		// Getter: returns RestartSolution.
		bool GetRestartSolution(void)                                                                     const {return RestartSolution;}
		// Getter: returns TimeStep.
		as3double GetTimeStep(void)                                                                       const {return TimeStep;}
		// Getter: returns MaxIter.
		unsigned long GetMaxIter(void)                                                                    const {return MaxIter;}
    // Getter: returns DomainBound.
    const as3vector1d<as3double> &GetDomainBound(void)                                                const {return DomainBound;}
    // Getter: returns nxElemZone.
    const as3vector1d<unsigned long> &GetnxElemZone(void)                                             const {return nxElemZone;}
    // Getter: returns nxElemZone, per input zone.
    unsigned long GetnxElemZone(unsigned short iZone)                                                 const {return nxElemZone[iZone];}
    // Getter: returns nyElemZone.
    const as3vector1d<unsigned long> &GetnyElemZone(void)                                             const {return nyElemZone;}
    // Getter: returns nyElemZone, per input zone.
    unsigned long GetnyElemZone(unsigned short iZone)                                                 const {return nyElemZone[iZone];}
    // Getter: returns hRatioZone.
    const as3vector1d<as3double> &GethElemRatioZone(void)                                             const {return hElemRatioZone;}
    // Getter: returns ParamOutletNSCBC.
    const as3vector1d<as3double> &GetParamOutletNSCBC(void)                                           const {return ParamOutletNSCBC;}
    // Getter: returns ParamInletNSCBC.
    const as3vector1d<as3double> &GetParamInletNSCBC(void)                                            const {return ParamInletNSCBC;}
    // Getter: returns DampingConstant, per input zone.
    as3double GetDampingConstant(unsigned short iZone)                                                const {return DampingConstant[iZone];}
    // Getter: returns DampingExponent, per input zone.
    as3double GetDampingExponent(unsigned short iZone)                                                const {return DampingExponent[iZone];}
    // Getter: returns TypeProcessData.
    unsigned short GetTypeProcessData(void)                                                           const {return TypeProcessData;}
    // Getter: returns ProbeLocation.
    const as3vector1d<as3double> &GetProbeLocation(void)                                              const {return ProbeLocation;}
    // Getter: returns SampleZoneData.
    bool GetSampleZoneData(void)                                                                      const {return SampleZoneData;}
    // Getter: returns TypeZoneData.
    unsigned short GetTypeZoneData(void)                                                              const {return TypeZoneData;}
    // Getter: returns WriteFreqZoneData.
    unsigned long GetWriteFreqZoneData(void)                                                          const {return WriteFreqZoneData;}
    // Getter: returns TypeFilterSolution, per input zone.
    unsigned short GetTypeFilterSolution(unsigned short iZone)                                        const {return TypeFilterSolution[iZone];}
    // Getter: returns FilterCharacteristics.
    const as3vector1d<unsigned short> &GetFilterCharacteristics(void)                                 const {return FilterCharacteristics;}
    // Getter: returns GridStretching, per input zone.
    bool GetGridStretching(unsigned short iZone)                                                      const {return GridStretching[iZone];}
    // Getter: returns GridStretchingConstant, per input zone.
    as3double GetGridStretchingConstant(unsigned short iZone)                                         const {return GridStretchingConstant[iZone];}
    // Getter: returns GridStretchingExponent, per input zone.
    as3double GetGridStretchingExponent(unsigned short iZone)                                         const {return GridStretchingExponent[iZone];}
    // Getter: returns ArtificialConvection, per input zone.
    bool GetArtificialConvection(unsigned short iZone)                                                const {return ArtificialConvection[iZone];}
    // Getter: returns ArtificialConvectionConstant, per input zone.
    as3double GetArtificialConvectionConstant(unsigned short iZone)                                   const {return ArtificialConvectionConstant[iZone];}
    // Getter: returns ArtificialConvectionExponent, per input zone.
    as3double GetArtificialConvectionExponent(unsigned short iZone)                                   const {return ArtificialConvectionExponent[iZone];}
    // Getter: returns PeriodicPulse.
    bool GetPeriodicPulse(void)                                                                       const {return PeriodicPulse;}
    // Getter: returns TypeBC.
    const as3vector2d<unsigned short> &GetTypeBC(void)                                                const {return TypeBC;}
    // Getter: returns TypeBC, per input zone.
    const as3vector1d<unsigned short> &GetTypeBC(unsigned short iZone)                                const {return TypeBC[iZone];}
    // Getter: returns TypeBC, per input zone, per input boundary.
    unsigned short GetTypeBC(unsigned short iZone, unsigned short iBoundary)                          const {return TypeBC[iZone][iBoundary];}
    // Getter: returns InterfaceID.
    const as3vector3d<unsigned short> &GetInterfaceID(void)                                           const {return InterfaceID;}
    // Getter: returns InterfaceID, per input zone.
    const as3vector2d<unsigned short> &GetInterfaceID(unsigned short iZone)                           const {return InterfaceID[iZone];}
    // Getter: returns InterfaceID, per input zone, per input boundary.
    const as3vector1d<unsigned short> &GetInterfaceID(unsigned short iZone, unsigned short iBoundary) const {return InterfaceID[iZone][iBoundary];}

    private:
      // Default data parameters.
  		struct {
  			// Number of zones.
  			unsigned short nZone = 1;
  			// Integration factor.
  			unsigned short IntegrationFactor = 3;
  			// Polynomial order.
  			as3vector1d<unsigned short> nPoly = { 1 };
        // Type of DOFs.
        as3vector1d<std::string> NameNodalDOFs;
  			// Zone marker.
  			as3vector1d<std::string> NameZoneMarker = { "ZONE_MAIN" };
  			// Markers.
  			as3vector1d<std::string> NameMarker = {};
  			// Restart solution or not.
  			std::string NameRestartSolution = "false";
  			// Time step.
  			as3double TimeStep = -1.0;
  			// Maximum iterations in time.
  			unsigned long MaxIter = 100000;
        // Ratio of element sizes in (WEST, EAST, SOUTH, NORTH).
        as3vector1d<as3double> hElemRatioZone = { 1.0, 1.0, 1.0, 1.0 };
        // Zone conformity w.r.t. ZONE_MAIN.
        std::string NameZoneConformity = "true";
        // Type of riemann solver.
        as3vector1d<std::string> NameRiemannSolver = { "ROE" };
        // Output writing frequency.
        unsigned long WriteFreq  = 10;
        // Output monitoring frequency.
        unsigned long OutputFreq = 100;
        // CFL number.
        as3double CFL = 1.0;
        // Adaptive time step.
        std::string NameAdaptTime = "false";
        // Relaxation parameters in the case of an outlet NSCBC.
        // Dimension: [0]: sigma, [1]: beta_l, [2]: beta_t, [3]: len, [4]: eta.
        as3vector1d<as3double> ParamOutletNSCBC = { 0.1, -1.0, -1.0, 1.0, 0.8 };
        // Relaxation parameters in the case of an inlet NSCBC.
        // Dimension: [0]: sigma, [1]: beta_t, [2]: len, [3]: eta.
        as3vector1d<as3double> ParamInletNSCBC = { 1.0, -1.0, 1.0, 0.8 };
        // Sponge-layer damping exponential coefficient.
        as3vector1d<as3double> DampingExponent = { 2.0 };
        // Sponge-layer damping constant.
        as3vector1d<as3double> DampingConstant = { 2.0 };
        // Free-stream Mach number.
        as3double MachInf = 0.5;
        // Flow angle in degrees (w.r.t. x-direction).
        as3double FlowAngle = 0.0;
        // Process data.
        std::string NameTypeProcessData = "PROCESS_NOTHING";
        // Disturbance center.
        as3vector1d<as3double> CenterX0 = { 0.0, 0.0 };
        // Ratio of disturbance w.r.t. background flow.
        as3double DisturbanceRatio = 0.2;
        // Disturbance width.
        as3double DisturbanceWidth = 0.25;
        // frequency.
        as3double Frequency = 1.0;
        // Type of filter applied on the solution.
        as3vector1d<std::string> NameTypeFilterSolution = { "NONE" };
        // Filtering parameters.
        // Dimension: [0]: freq, [1]: Nc, [2]: s, [3]: alpha.
        as3vector1d<unsigned short> FilterCharacteristics = { 5, 4, 16, 36 };
        // Periodic pulse.
        std::string NamePeriodicPulse = "false";
        // Type of buffer layer.
        as3vector1d<std::string> NameTypeBufferLayer = { "NONE" };
        // Grid-stretching constant.
        as3vector1d<as3double> GridStretchingConstant = { 1.0 };
        // Grid-stretching exponent.
        as3vector1d<as3double> GridStretchingExponent = { 1.0 };
        // Grid-stretching used.
        bool GridStretching = false;
        // Artificial-convection constant (max Mach number at boundary).
        as3vector1d<as3double> ArtificialConvectionConstant = { 1.15 };
        // Artificial-convection exponent.
        as3vector1d<as3double> ArtificialConvectionExponent = { 2.0 };
        // Artificial-convection used.
        bool ArtificialConvection = false;
        // Probe locations as 1D array of (x1,y1, x2,y2, ... etc.).
        as3vector1d<as3double> ProbeLocation = { };
        // Whether or not to sample and write zone data.
        std::string NameSampleZoneData = "false";
        // Zone data output filename.
        std::string OutputZoneDataFilename = "zonedata";
        // Zone data selected for writing.
        std::string NameMarkerZoneData = "ZONE_MAIN";
        // Output writing frequency of zone data.
        unsigned long WriteFreqZoneData = 100;
        // Face of BCs to modify.
        as3vector1d<std::string> NameTypeModifyBC = { "NONE" };
        // Angular frequency of modified BC of each working variable.
        as3vector1d<as3double> ModifyFreqBC = { 1.0 };
        // Pulse spatial width of modified BC of each working variable.
        as3vector1d<as3double> ModifyWidthBC = { 0.25 };
        // Pulse spatial amplitude w.r.t. background flow of modified BC of each working variable.
        as3vector1d<as3double> ModifyStrengthBC = { 0.2 };
        // Pulse center(s) for BC modification.
        as3vector1d<as3double> ModifyCenterBC = { -1.0, 0.0 };
        // Pulse center(s) shoft for BC modification of each working variable.
        as3vector1d<as3double> ModifyShiftCenterBC = { 0.25 };
        // Write PML auxiliary data.
        std::string NameWriteAuxiliaryDataPML = "false";
        // Periodic flow direction aligned with flow.
        std::string NameAlignedPeriodicPulse = "false";
        // Constant frequency in source term.
        std::string NameConstantFrequency = "true";
				// Periodic source center fixed or not.
				std::string NameSourceTermCenterFixed = "true";
				// Periodic source varying spatial center shift: (dx, dy).
				as3vector1d<as3double> SourceTermCenterShift = { 1.0, 1.0 };
				// Periodic source frequency exponent alpha: sin^alpha.
				as3double SourceFrequencyExponent = 1.0;
				// Source frequency parameters f = f0 * A^(t/Tau).
				// Dimension: [0]: A, [1]: Tau.
				as3vector1d<as3double> SourceFrequencyParam = { 1.0, 0.1 };
        // Domain expansion ratios: ( rxb1, rxb2, ryb1, ryb2 ).
        as3vector1d<as3double> DomainExpansionRatio = { 1.0, 1.0, 1.0, 1.0 };
        // Uniform equidistant grid resolution or not.
        std::string NameUniformGridResolution = "true";

  		} DefaultParam;

      // Output writing frequency.
      unsigned long WriteFreq;
      // Output monitoring frequency.
      unsigned long OutputFreq;
      // Output solution filename.
  		std::string OutputSolFilename;
  		// Output visualization filename.
  		std::string OutputVTKFilename;
      // Restart solution filename.
      std::string RestartFilename;

      // Number of zones.
  		unsigned short nZone;
  		// Polynomial order of solution per zone.
  		as3vector1d<unsigned short> nPolySolZone;
      // Multizone strategy used in the simulation.
      unsigned short MultizoneStrategy;

      // Domain bounds (physical domain).
      // Convention: (WEST, EAST, SOUTH, NORTH).
      as3vector1d<as3double> DomainBound;
      // Number of elements in x-direction in each zone.
      as3vector1d<unsigned long> nxElemZone;
      // Number of elements in y-direction in each zone.
      as3vector1d<unsigned long> nyElemZone;
      // Element ratio per extension of zone: (WEST, EAST, SOUTH, NORTH).
      as3vector1d<as3double> hElemRatioZone;
      // Use zone conformity between adjacent elements on zone.
      bool ZoneConformity;

      // Block interface location.
      as3vector1d<as3double> BlockInterfaceLocation;
      // Domain expansion ratios: ( rxb1, rxb2, ryb1, ryb2 ).
      as3vector1d<as3double> DomainExpansionRatio;
      // Number of sub-elements in x-direction in main zone.
      as3vector1d<unsigned long> nxBlockElem;
      // Number of sub-elements in y-direction in main zone.
      as3vector1d<unsigned long> nyBlockElem;
      // Uniform equidistant grid resolution or not.
      bool UniformGridResolution;

      // Zone markers, string name.
      as3vector1d<std::string> NameZoneMarker;
      // Zone type, enum values.
      as3vector1d<unsigned short> TypeZone;

      // Nodal type of DOFs, string name.
      as3vector1d<std::string> NameNodalDOFs;
      // Nodal type of DOFs, enum values.
      as3vector1d<unsigned short> TypeDOFs;
      // Type of solver per zone.
      as3vector1d<unsigned short> TypeSolver;

      // Type of temporal discretization, string name.
  		std::string NameTemporalScheme;
  		// Type of temporal discretization.
  		unsigned short TypeTemporalScheme;
  		// Simulation start and end times.
  		// index: [0]: start, [1]: end.
  		as3double SimulationTime[2];
  		// Time step input.
  		as3double TimeStep;
  		// Maximum temporal iterations.
  		unsigned long MaxIter;
      // CFL number.
      as3double CFL;
      // Adaptive time step.
      bool AdaptTime;

      // Integration factor (multiplied by polynomial of solution).
  		unsigned short IntegrationFactor;
      // Type of filtering used, string name.
      as3vector1d<std::string> NameTypeFilterSolution;
      // Type of filtering used.
      as3vector1d<unsigned short> TypeFilterSolution;
      // Filtering parameters.
      // Dimension: [0]: freq, [1]: Nc, [2]: s, [3]: alpha.
      as3vector1d<unsigned short> FilterCharacteristics;

      // Type of initial conditions, string name.
  		as3vector1d<std::string> NameInitialCondition;
  		// Type of initial conditions, enum values.
  		as3vector1d<unsigned short> TypeIC;
  		// Decide whether or not to restart a simulation.
  		bool RestartSolution;

      // Type of data processing, string name.
      std::string NameTypeProcessData;
      // Type of data processing, enum values.
      unsigned short TypeProcessData;
      // Output processed data filename.
      std::string OutputProcessedFilename;
      // Probe locations, array of coordinates: (x1,y1, x2,y2, ... etc.).
      as3vector1d<as3double> ProbeLocation;
      // Zone data output filename.
      std::string OutputZoneDataFilename;
      // Zone data selected for writing.
      std::string NameMarkerZoneData;
      // Type of zone selected for writing, enum values.
      unsigned short TypeZoneData;
      // Whether a PML zone exists in all zones or not.
      bool UsePML;
      // Write PML auxiliary data.
      bool WriteAuxiliaryDataPML;
      // Output writing frequency of zone data.
      unsigned long WriteFreqZoneData;
      // Whether or not to sample and write zone data.
      bool SampleZoneData;

      // Whether or not to use grid-stretching.
      as3vector1d<bool> GridStretching;
      // Grid-stretching constant.
      as3vector1d<as3double> GridStretchingConstant;
      // Grid-stretching exponent.
      as3vector1d<as3double> GridStretchingExponent;

      // Free-stream Mach number.
      as3double MachInf;
      // Flow angle in degrees (w.r.t. x-direction).
      as3double FlowAngle;
      // Presence of a cross-flow.
      bool CrossFlow;
      // Disturbance center.
      as3vector1d<as3double> CenterX0;
      // Ratio of disturbance w.r.t. background flow.
      as3double DisturbanceRatio;
      // Disturbance width.
      as3double DisturbanceWidth;
      // Frequency.
      as3double Frequency;
      // Angular frequency.
      as3double AngularFrequency;
      // Presence of periodic pulse.
      bool PeriodicPulse;
      // Use a constant frequency in source term.
      bool ConstantFrequency;
      // Periodic flow direction aligned with flow.
      bool AlignedPeriodicPulse;
			// Periodic source varying frequency parameters: f = f0 * A^(t/Tau).
			// Dimension: [0]: A, [1]: Tau.
			as3vector1d<as3double> SourceFrequencyParam;
			// Periodic source varying spatial center shift: (dx, dy).
			as3vector1d<as3double> SourceTermCenterShift;
			// Periodic source frequency exponent alpha: sin^alpha.
			as3double SourceFrequencyExponent;
			// Periodic source center fixed or not.
			bool SourceTermCenterFixed;

      // Name of type of BCs to modify per face.
      as3vector1d<std::string> NameTypeModifyBC;
      // Type of modified BC per face.
      as3vector1d<unsigned short> TypeModifyBC;
      // Angular frequency of modified BC of each working variable.
      as3vector1d<as3double> ModifyFreqBC;
      // Pulse spatial width of modified BC of each working variable.
      as3vector1d<as3double> ModifyWidthBC;
      // Pulse spatial amplitude w.r.t. background flow of modified BC of each working variable.
      as3vector1d<as3double> ModifyStrengthBC;
      // Pulse center(s) for BC modification.
      as3vector1d<as3double> ModifyCenterBC;
      // Pulse center(s) shoft for BC modification of each working variable.
      as3vector1d<as3double> ModifyShiftCenterBC;

      // Type of buffer layer, string name.
      as3vector1d<std::string> NameTypeBufferLayer;
      // Type of buffer layer, enum values.
      as3vector1d<unsigned short> TypeBufferLayer;

      // Type of riemann solver selected, string name.
      as3vector1d<std::string> NameRiemannSolver;
      // Type of riemann solver, enum values.
      as3vector1d<unsigned short> RiemannSolver;

      // Type of boundary conditions, string name:
      // ... indices: (SOUTH, NORTH, WEST, EAST).
      as3vector1d<std::string> NameBoundaryCondition;
      // Type of external boundary conditions over all zones:
      // ... indices: (SOUTH, NORTH, WEST, EAST).
      as3vector1d<unsigned short> TypeExternalBC;
      // Type of boundary conditions per each zone.
      // Dimension: [iZone][iBoundary], where iBoundary
      // has indices: (SOUTH, NORTH, WEST, EAST).
      as3vector2d<unsigned short> TypeBC;
      // Interface ID for periodic/zonal interface boundaries.
      // Dimension: [iZone][iBoundary][iData], with indices:
      // [iData] = [0]: jZone, [1]: jFace.
      as3vector3d<unsigned short> InterfaceID;

      // Relaxation parameters in the case of an outlet NSCBC.
      // Dimension: [0]: sigma, [1]: beta_l, [2]: beta_t, [3]: len, [4]: eta.
      as3vector1d<as3double> ParamOutletNSCBC;
      // Relaxation parameters in the case of an inlet NSCBC.
      // Dimension: [0]: sigma, [1]: beta_t, [2]: len, [3]: eta.
      as3vector1d<as3double> ParamInletNSCBC;

      // Sponge-layer damping exponential coefficient.
      as3vector1d<as3double> DampingExponent;
      // Sponge-layer damping constant.
      as3vector1d<as3double> DampingConstant;

      // Artificial-convection constant (max Mach number at boundary).
      as3vector1d<as3double> ArtificialConvectionConstant;
      // Artificial-convection exponent.
      as3vector1d<as3double> ArtificialConvectionExponent;
      // Artificial-convection used.
      as3vector1d<bool> ArtificialConvection;

      // Function that maps TypeZoneData from string to enum.
      void MapTypeZoneData(void);
      // Function that maps TypeDOFs from string to enum.
      void MapTypeDOFs(void);
      // Function that maps TypeZone from string to enum.
      void MapTypeZone(void);
      // Function that maps RiemannSolver from string to enum.
      void MapRiemannSolver(void);
      // Function that maps TemporalScheme from string to enum.
      void MapTemporalScheme(void);
      // Function that maps TypeIC from string to enum.
      void MapTypeIC(void);
      // Function that maps TypeExternalBC from string to enum.
      void MapTypeExternalBC(void);
      // Function that maps TypeProcessData from string to enum.
      void MapTypeProcessData(void);
      // Function that maps TypeFilterSolution from string to enum.
      void MapTypeFilterSolution(void);
      // Function that maps TypeBufferLayer from string to enum.
      void MapTypeBufferLayer(void);
      // Function that maps TypeModifyBC from string to enum.
      void MapTypeModifyBC(void);

      // Function that determines the multizone strategy adopted.
      void DetermineMultizoneStrategy(void);
      // Function that checks the element ratios specified according
      // to specified multizone strategy.
      void CheckElementRatio(void);
      // Function that processes the zone conformity in the specified elements.
      void ProcessZoneConformity(void);
      // Function that processes the boundary conditions per each zone.
      void ProcessBoundaryConditions(void);
      // Function that matches each iZone iBoundary with its
      // counterpart in jZone jBoundary.
      as3vector1d<unsigned short> MatchInterface(unsigned short iZone,
                                                 unsigned short iFace);


      // Function that reads the Integration rule data.
      void ReadIntegrationRuleOptions(const char *configFile);
      // Function that reads the filtering data.
      void ReadFilteringOptions(const char *configFile);
      // Function that reads the buffer-layer data.
      void ReadBufferLayerOptions(const char *configFile);
      // Function that reads the boundary conditions data.
      void ReadBoundaryConditionOptions(const char *configFile);
      // Function that reads the sponge-damping data.
      void ReadSpongeDampingOptions(const char *configFile);
      // Function that reads the grid-stretching data.
      void ReadGridStretchingOptions(const char *configFile);
      // Function that reads the artificial-convection data.
      void ReadArtificialConvectionOptions(const char *configFile);
      // Function that reads the NSCBC data.
      void ReadCharacteristicBoundaryOptions(const char *configFile);

      // Function that reads grid options.
      bool ReadGridOptions(const char *configFile);
      // Function that reads input/output options.
  		bool ReadIOOptions(const char *configFile);
  		// Function that reads solver options.
  		bool ReadSolverOptions(const char *configFile);
  		// Function that reads boundary options.
  		bool ReadBoundaryOptions(const char *configFile);
  		// Function that reads temporal options.
  		bool ReadTemporalOptions(const char *configFile);
  		// Function that reads initial condition options.
  		bool ReadICOptions(const char *configFile);
      // Function that reads flow characteristics options.
      bool ReadFlowOptions(const char *configFile);
      // Function that reads processing information.
      bool ReadProcessingOptions(const char *configFile);
      // Function that reads modified boundary condition information.
      bool ReadModifiedBCOptions(const char *configFile);
};





