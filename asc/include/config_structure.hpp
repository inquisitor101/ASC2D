#pragma once

/*!
 * @file config_structure.hpp
 * @brief The file containing all the user input options.
 */

#include "option_structure.hpp"
#include <map>

#if __cplusplus > 201103L
	#include <filesystem>
#else
	#include <sys/types.h>
	#include <sys/stat.h>
#endif

/*!
 * @brief A class used for storing and processing the user-specified input options. 
 */
class CConfig {

	public:
		/*!
		 * @brief Default constructor of CConfig, which reads the input dictionary.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
		CConfig(const char *configFile);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CConfig(void);


		/*!
		 * @brief Getter function which returns the value of InletNRBC.
		 *
		 * @return InletNRBC
		 */
		bool GetInletNRBC(void)                                                                           const {return InletNRBC;}
		/*!
		 * @brief Getter function which returns the value of OutputVTKVariable.
		 *
		 * @return OutputVTKVariable
		 */
		as3vector1d<unsigned short> GetOutputVTKVariable(void)                                            const {return OutputVTKVariable;}
    /*!
		 * @brief Getter function which returns the value of ProcessLocation.
		 *
		 * @return ProcessLocation
		 */	
		unsigned short GetProcessLocation(void)                                                           const {return ProcessLocation;}
    /*!
		 * @brief Getter function which returns the value of TypeFileFormatProcessedProbed.
		 *
		 * @return TypeFileFormatProcessedProbed
		 */	
		unsigned short GetTypeFileFormatProcessedProbed(void)                                             const {return TypeFileFormatProcessedProbed;}
    /*!
		 * @brief Getter function which returns the value of GNUplotFilename.
		 *
		 * @return GNUplotFilename
		 */	
		const char *GetGNUplotFilename(void)                                                              const {return GNUplotFilename.c_str();}
    /*!
		 * @brief Getter function which returns the value of WriteGNUplot.
		 *
		 * @return WriteGNUplot
		 */	
		bool GetWriteGNUplot(void)                                                                        const {return WriteGNUplot;}
	  /*!
		 * @brief Getter function which returns the value of TypeFileFormatZone.
		 *
		 * @return TypeFileFormatZone
		 */		
		unsigned short GetTypeFileFormatZone(void)                                                        const {return TypeFileFormatZone;}
	  /*!
		 * @brief Getter function which returns the value of TypeFileFormatSolution.
		 *
		 * @return TypeFileFormatSolution
		 */	
		unsigned short GetTypeFileFormatSolution(void)                                                    const {return TypeFileFormatSolution;}
	  /*!
		 * @brief Getter function which returns the value of TypeFileFormatVTK.
		 *
		 * @return TypeFileFormatVTK
		 */	
		unsigned short GetTypeFileFormatVTK(void)                                                         const {return TypeFileFormatVTK;}
 	  /*!
		 * @brief Getter function which returns the value of CFL.
		 *
		 * @return CFL
		 */	
    as3double GetCFL(void)                                                                            const {return CFL;}
 	  /*!
		 * @brief Getter function which returns the value of AdaptTime.
		 *
		 * @return AdaptTime
		 */	
    bool GetAdaptTime(void)                                                                           const {return AdaptTime;}
 	  /*!
		 * @brief Getter function which returns the value of OutputFreq.
		 *
		 * @return OutputFreq
		 */	
    unsigned long GetOutputFreq(void)                                                                 const {return OutputFreq;}
 	  /*!
		 * @brief Getter function which returns the value of WriteVTKFreq.
		 *
		 * @return WriteVTKFreq
		 */	
    unsigned long GetWriteVTKFreq(void)                                                               const {return WriteVTKFreq;}
 	  /*!
		 * @brief Getter function which returns the value of WriteRestartFreq.
		 *
		 * @return WriteRestartFreq
		 */
    unsigned long GetWriteRestartFreq(void)                                                           const {return WriteRestartFreq;}
 	  /*!
		 * @brief Getter function which returns the value of OutputSolFilename.
		 *
		 * @return OutputSolFilename
		 */   
		const char *GetOutputSolFilename(void)                                                            const {return OutputSolFilename.c_str();}
 	  /*!
		 * @brief Getter function which returns the value of OutputVTKFilename.
		 *
		 * @return OutputVTKFilename
		 */ 
		const char *GetOutputVTKFilename(void)                                                            const {return OutputVTKFilename.c_str();}
 	  /*!
		 * @brief Getter function which returns the value of OutputProcessedDirectory.
		 *
		 * @return OutputProcessedDirectory
		 */ 
    const char *GetOutputProcessedDirectory(void)                                                     const {return OutputProcessedDirectory.c_str();}
 	  /*!
		 * @brief Getter function which returns the value of OutputZoneDataFilename.
		 *
		 * @return OutputZoneDataFilename
		 */ 
    const char *GetOutputZoneDataFilename(void)                                                       const {return OutputZoneDataFilename.c_str();}
 	  /*!
		 * @brief Getter function which returns the value of RestartFilename.
		 *
		 * @return RestartFilename
		 */ 
    const char *GetRestartFilename(void)                                                              const {return RestartFilename.c_str();}
 	  /*!
		 * @brief Getter function which returns the value of OutputSampleSurfaceDirectory.
		 *
		 * @return OutputSampleSurfaceDirectory
		 */ 
		const char *GetOutputSampleSurfaceDirectory(void)                                                 const {return OutputSampleSurfaceDirectory.c_str();}
 	  /*!
		 * @brief Getter function which returns the value of nZone.
		 *
		 * @return nZone
		 */ 
		unsigned short GetnZone(void)                                                                     const {return nZone;}
 	  /*!
		 * @brief Getter function which returns the value of MultizoneStrategy.
		 *
		 * @return MultizoneStrategy
		 */ 
    unsigned short GetMultizoneStrategy(void)                                                         const {return MultizoneStrategy;}
 	  /*!
		 * @brief Getter function which returns the value of nPolySolZone, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return nPolySolZone[iZone]
		 */ 
		unsigned short GetnPolySolZone(unsigned short iZone)                                              const {return nPolySolZone[iZone];}
 	  /*!
		 * @brief Getter function which returns the value of nPolySolZone.
		 *
		 * @return nPolySolZone
		 */ 
		const as3vector1d<unsigned short> &GetnPolySolZone(void)                                          const {return nPolySolZone;}
 	  /*!
		 * @brief Getter function which returns the value of BlockInterfaceLocation.
		 *
		 * @return BlockInterfaceLocation
		 */    
    const as3vector1d<as3double> &GetBlockInterfaceLocation(void)                                     const {return BlockInterfaceLocation;}
 	  /*!
		 * @brief Getter function which returns the value of DomainExpansionRatio.
		 *
		 * @return DomainExpansionRatio
		 */ 
    const as3vector1d<as3double> &GetDomainExpansionRatio(void)                                       const {return DomainExpansionRatio;}
 	  /*!
		 * @brief Getter function which returns the value of nxBlockElem.
		 *
		 * @return nxBlockElem
		 */ 
    const as3vector1d<unsigned long> &GetnxBlockElem(void)                                            const {return nxBlockElem;}
 	  /*!
		 * @brief Getter function which returns the value of nyBlockElem.
		 *
		 * @return nyBlockElem
		 */ 
    const as3vector1d<unsigned long> &GetnyBlockElem(void)                                            const {return nyBlockElem;}
 	  /*!
		 * @brief Getter function which returns the value of UniformGridResolution.
		 *
		 * @return UniformGridResolution
		 */
    bool GetUniformGridResolution(void)                                                               const {return UniformGridResolution;}
 	  /*!
		 * @brief Getter function which returns the value of SampleSurfaceData.
		 *
		 * @return SampleSurfaceData
		 */
		bool GetSampleSurfaceData(void)                                                                   const {return SampleSurfaceData;}
 	  /*!
		 * @brief Getter function which returns the value of SampleDataBoundaryID.
		 *
		 * @return SampleDataBoundaryID
		 */
		as3vector1d<unsigned short> GetSampleDataBoundaryID(void)                                         const {return SampleDataBoundaryID;}
 	  /*!
		 * @brief Getter function which returns the value of NameMarkerSampleSurface.
		 *
		 * @return NameMarkerSampleSurface
		 */
		as3vector1d<std::string> GetNameMarkerSampleSurface(void)                                         const {return NameMarkerSampleSurface;}
 	  /*!
		 * @brief Getter function which returns the value of TypeZone, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return TypeZone[iZone]
		 */
    unsigned short GetTypeZone(unsigned short iZone)                                                  const {return TypeZone[iZone];}
 	  /*!
		 * @brief Getter function which returns the value of TypeZone.
		 *
		 * @return TypeZone
		 */
		const as3vector1d<unsigned short> &GetTypeZone(void)                                              const {return TypeZone;}
 	  /*!
		 * @brief Getter function which returns the value of TypeDOFs, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return TypeDOFs[iZone]
		 */   
		unsigned short GetTypeDOFs(unsigned short iZone)                                                  const {return TypeDOFs[iZone];}
 	  /*!
		 * @brief Getter function which returns the value of TypeDOFs.
		 *
		 * @return TypeDOFs
		 */ 
		const as3vector1d<unsigned short> &GetTypeDOFs(void)                                              const {return TypeDOFs;}
 	  /*!
		 * @brief Getter function which returns the value of TypeSolver, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return TypeSolver[iZone]
		 */
    unsigned short GetTypeSolver(unsigned short iZone)                                                const {return TypeSolver[iZone];}
 	  /*!
		 * @brief Getter function which returns the value of TypeSolver.
		 *
		 * @return TypeSolver
		 */
    const as3vector1d<unsigned short> &GetTypeSolver(void)                                            const {return TypeSolver;}
 	  /*!
		 * @brief Getter function which returns the value of TypeBufferLayer, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return TypeBufferLayer[iZone]
		 */ 
    unsigned short GetTypeBufferLayer(unsigned short iZone)                                           const {return TypeBufferLayer[iZone];}
		/*!
		 * @brief Getter function which returns the value of TypeBufferLayer.
		 *
		 * @return TypeBufferLayer
		 */   
    const as3vector1d<unsigned short> &GetTypeBufferLayer(void)                                       const {return TypeBufferLayer;}
  	/*!
		 * @brief Getter function which returns the value of UsePML.
		 *
		 * @return UsePML
		 */   
    bool GetUsePML(void)                                                                              const {return UsePML;}
  	/*!
		 * @brief Getter function which returns the value of WriteAuxiliaryDataPML.
		 *
		 * @return WriteAuxiliaryDataPML
		 */
    bool GetWriteAuxiliaryDataPML(void)                                                               const {return WriteAuxiliaryDataPML;}
  	/*!
		 * @brief Getter function which returns the value of IntegrationFactor.
		 *
		 * @return IntegrationFactor
		 */
		unsigned short GetIntegrationFactor(void)                                                         const {return IntegrationFactor;}
  	/*!
		 * @brief Getter function which returns the value of SimulationTime[0].
		 *
		 * @return SimulationTime[0]
		 */
		as3double GetSimulationStartTime(void)                                                            const {return SimulationTime[0];}
	 	/*!
		 * @brief Getter function which returns the value of SimulationTime[1].
		 *
		 * @return SimulationTime[1]
		 */	
		as3double GetSimulationFinalTime(void)                                                            const {return SimulationTime[1];}
	 	/*!
		 * @brief Getter function which returns the value of TypeTemporalScheme.
		 *
		 * @return TypeTemporalScheme
		 */	
		unsigned short GetTypeTemporalScheme(void)                                                        const {return TypeTemporalScheme;}
 	 	/*!
		 * @brief Getter function which returns the value of TypeIC, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return TypeIC[iZone]
		 */
    unsigned short GetTypeIC(unsigned short iZone)                                                    const {return TypeIC[iZone];}
		/*!
		 * @brief Getter function which returns the value of TypeIC.
		 *
		 * @return TypeIC
		 */
		const as3vector1d<unsigned short> &GetTypeIC(void)                                                const {return TypeIC;}
		/*!
		 * @brief Getter function which returns the value of FlowAngle.
		 *
		 * @return FlowAngle
		 */
    as3double GetFlowAngle(void)                                                                      const {return FlowAngle;}
  	/*!
		 * @brief Getter function which returns the value of MachInf.
		 *
		 * @return MachInf
		 */  
    as3double GetMachInf(void)                                                                        const {return MachInf;}
  	/*!
		 * @brief Getter function which returns the value of CrossFlow.
		 *
		 * @return CrossFlow
		 */
    bool GetCrossFlow(void)                                                                           const {return CrossFlow;}
  	/*!
		 * @brief Getter function which returns the value of CenterX0.
		 *
		 * @return CenterX0
		 */
    const as3vector1d<as3double> &GetCenterX0(void)                                                   const {return CenterX0;}
  	/*!
		 * @brief Getter function which returns the value of DisturbanceRatio.
		 *
		 * @return DisturbanceRatio
		 */
    as3double GetDisturbanceRatio(void)                                                               const {return DisturbanceRatio;}
  	/*!
		 * @brief Getter function which returns the value of DisturbanceWidth.
		 *
		 * @return DisturbanceWidth
		 */
    as3double GetDisturbanceWidth(void)                                                               const {return DisturbanceWidth;}
  	/*!
		 * @brief Getter function which returns the value of AngularFrequency.
		 *
		 * @return AngularFrequency
		 */
    as3double GetAngularFrequency(void)                                                               const {return AngularFrequency;}
  	/*!
		 * @brief Getter function which returns the value of ConstantFrequency.
		 *
		 * @return ConstantFrequency
		 */
    bool GetConstantFrequency(void)                                                                   const {return ConstantFrequency;}
  	/*!
		 * @brief Getter function which returns the value of TypeModifyBC, per input boundary.
		 *
		 * @param[in] iBoundary boundary ID.
		 *
		 * @return TypeModifyBC[iBoundary]
		 */
    unsigned short GetTypeModifyBC(unsigned short iBoundary)                                          const {return TypeModifyBC[iBoundary];}
  	/*!
		 * @brief Getter function which returns the value of ModifyFreqBC.
		 *
		 * @return ModifyFreqBC
		 */
    const as3vector1d<as3double> GetModifyFreqBC(void)                                                const {return ModifyFreqBC;}
   	/*!
		 * @brief Getter function which returns the value of ModifyStrengthBC.
		 *
		 * @return ModifyStrengthBC
		 */   
    const as3vector1d<as3double> GetModifyStrengthBC(void)                                            const {return ModifyStrengthBC;}
   	/*!
		 * @brief Getter function which returns the value of ModifyWidthBC.
		 *
		 * @return ModifyWidthBC
		 */
    const as3vector1d<as3double> GetModifyWidthBC(void)                                               const {return ModifyWidthBC;}
   	/*!
		 * @brief Getter function which returns the value of ModifyCenterBC.
		 *
		 * @return ModifyCenterBC
		 */
    const as3vector1d<as3double> GetModifyCenterBC(void)                                              const {return ModifyCenterBC;}
    /*!
		 * @brief Getter function which returns the value of ModifyShiftCenterBC.
		 *
		 * @return ModifyShiftCenterBC
		 */   
    const as3vector1d<as3double> GetModifyShiftCenterBC(void)                                         const {return ModifyShiftCenterBC;}
    /*!
		 * @brief Getter function which returns the value of RiemannSolver, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return RiemannSolver[iZone]
		 */    
    unsigned short GetRiemannSolver(unsigned short iZone)                                             const {return RiemannSolver[iZone];}
		/*!
		 * @brief Getter function which returns the value of RiemannSolver.
		 *
		 * @return RiemannSolver
		 */   
    const as3vector1d<unsigned short> &GetRiemannSolver(void)                                         const {return RiemannSolver;}
		/*!
		 * @brief Getter function which returns the value of AlignedPeriodicPulse.
		 *
		 * @return AlignedPeriodicPulse
		 */
    bool GetAlignedPeriodicPulse(void)                                                                const {return AlignedPeriodicPulse;}
		/*!
		 * @brief Getter function which returns the value of SourceFrequencyParam.
		 *
		 * @return SourceFrequencyParam
		 */
		const as3vector1d<as3double> &GetSourceFrequencyParam(void)                                       const {return SourceFrequencyParam;}
		/*!
		 * @brief Getter function which returns the value of SourceFrequencyExponent.
		 *
		 * @return SourceFrequencyExponent
		 */
		as3double GetSourceFrequencyExponent(void)                                                        const {return SourceFrequencyExponent;}
		/*!
		 * @brief Getter function which returns the value of SourceTermCenterShift.
		 *
		 * @return SourceTermCenterShift
		 */
		const as3vector1d<as3double> &GetSourceTermCenterShift(void)                                      const {return SourceTermCenterShift;}	
		/*!
		 * @brief Getter function which returns the value of SourceTermCenterFixed.
		 *
		 * @return SourceTermCenterFixed
		 */
		bool GetSourceTermCenterFixed(void)                                                               const {return SourceTermCenterFixed;}
		/*!
		 * @brief Getter function which returns the value of RestartSolution.
		 *
		 * @return RestartSolution
		 */	
		bool GetRestartSolution(void)                                                                     const {return RestartSolution;}
		/*!
		 * @brief Getter function which returns the value of TimeStep.
		 *
		 * @return TimeStep
		 */
		as3double GetTimeStep(void)                                                                       const {return TimeStep;}
		/*!
		 * @brief Getter function which returns the value of MaxIter.
		 *
		 * @return MaxIter
		 */
		unsigned long GetMaxIter(void)                                                                    const {return MaxIter;}
 		/*!
		 * @brief Getter function which returns the value of DomainBound.
		 *
		 * @return DomainBound
		 */
    const as3vector1d<as3double> &GetDomainBound(void)                                                const {return DomainBound;}
 		/*!
		 * @brief Getter function which returns the value of nxElemZone, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return nxElemZone[iZone]
		 */   
    unsigned long GetnxElemZone(unsigned short iZone)                                                 const {return nxElemZone[iZone];}
		/*!
		 * @brief Getter function which returns the value of nxElemZone.
		 *
		 * @return nxElemZone
		 */
    const as3vector1d<unsigned long> &GetnxElemZone(void)                                             const {return nxElemZone;}
 		/*!
		 * @brief Getter function which returns the value of nyElemZone, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return nyElemZone[iZone]
		 */    
    unsigned long GetnyElemZone(unsigned short iZone)                                                 const {return nyElemZone[iZone];}
		/*!
		 * @brief Getter function which returns the value of nyElemZone.
		 *
		 * @return nyElemZone
		 */
    const as3vector1d<unsigned long> &GetnyElemZone(void)                                             const {return nyElemZone;}
		/*!
		 * @brief Getter function which returns the value of hElemRatioZone.
		 *
		 * @return hElemRatioZone
		 */
    const as3vector1d<as3double> &GethElemRatioZone(void)                                             const {return hElemRatioZone;}
  	/*!
		 * @brief Getter function which returns the value of ParamOutletNSCBC.
		 *
		 * @return ParamOutletNSCBC
		 */  
    const as3vector1d<as3double> &GetParamOutletNSCBC(void)                                           const {return ParamOutletNSCBC;}
  	/*!
		 * @brief Getter function which returns the value of ParamInletNSCBC.
		 *
		 * @return ParamInletNSCBC
		 */
    const as3vector1d<as3double> &GetParamInletNSCBC(void)                                            const {return ParamInletNSCBC;}
   	/*!
		 * @brief Getter function which returns the value of DampingConstant, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return DampingConstant[iZone] 
		 */   
    as3double GetDampingConstant(unsigned short iZone)                                                const {return DampingConstant[iZone];}
   	/*!
		 * @brief Getter function which returns the value of DampingExponent, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return DampingExponent[iZone]
		 */   
    as3double GetDampingExponent(unsigned short iZone)                                                const {return DampingExponent[iZone];}
   	/*!
		 * @brief Getter function which returns the value of CharacteristicConstant, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return CharacteristicConstant[iZone]
		 */   
    as3double GetCharacteristicConstant(unsigned short iZone)                                         const {return CharacteristicConstant[iZone];}
   	/*!
		 * @brief Getter function which returns the value of CharacteristicExponent, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return CharacteristicExponent[iZone]
		 */
    as3double GetCharacteristicExponent(unsigned short iZone)                                         const {return CharacteristicExponent[iZone];}
   	/*!
		 * @brief Getter function which returns the value of CharacteristicMatching, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return CharacteristicMatching[iZone]
		 */   
    bool GetCharacteristicMatching(unsigned short iZone)                                              const {return CharacteristicMatching[iZone];}
   	/*!
		 * @brief Getter function which returns the value of TypeProcessData.
		 *
		 * @return TypeProcessData
		 */    
    unsigned short GetTypeProcessData(void)                                                           const {return TypeProcessData;}
   	/*!
		 * @brief Getter function which returns the value of ProbeSpecified.
		 *
		 * @return ProbeSpecified
		 */ 
		bool GetProbeSpecified(void)                                                                      const {return ProbeSpecified;}
   	/*!
		 * @brief Getter function which returns the value of ProbeLocation.
		 *
		 * @return ProbeLocation
		 */
    const as3vector1d<as3double> &GetProbeLocation(void)                                              const {return ProbeLocation;}
   	/*!
		 * @brief Getter function which returns the value of ProbeVariable.
		 *
		 * @return ProbeVariable
		 */
		const as3vector1d<unsigned short> &GetProbeVariable(void)                                         const {return ProbeVariable;}
   	/*!
		 * @brief Getter function which returns the value of SampleZoneData.
		 *
		 * @return SampleZoneData
		 */
    bool GetSampleZoneData(void)                                                                      const {return SampleZoneData;}
   	/*!
		 * @brief Getter function which returns the value of TypeZoneData.
		 *
		 * @return TypeZoneData
		 */   
    unsigned short GetTypeZoneData(void)                                                              const {return TypeZoneData;}
   	/*!
		 * @brief Getter function which returns the value of WriteFreqZoneData.
		 *
		 * @return WriteFreqZoneData
		 */   
    unsigned long GetWriteFreqZoneData(void)                                                          const {return WriteFreqZoneData;}
   	/*!
		 * @brief Getter function which returns the value of TypeFilterSolution, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return TypeFilterSolution[iZone]
		 */
    unsigned short GetTypeFilterSolution(unsigned short iZone)                                        const {return TypeFilterSolution[iZone];}
   	/*!
		 * @brief Getter function which returns the value of FilterCharacteristics.
		 *
		 * @return FilterCharacteristics
		 */
    const as3vector1d<unsigned short> &GetFilterCharacteristics(void)                                 const {return FilterCharacteristics;}
   	/*!
		 * @brief Getter function which returns the value of GridStretching, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return GridStretching[iZone]
		 */   
    bool GetGridStretching(unsigned short iZone)                                                      const {return GridStretching[iZone];}
   	/*!
		 * @brief Getter function which returns the value of GridStretchingConstant, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return GridStretchingConstant[iZone]
		 */   
    as3double GetGridStretchingConstant(unsigned short iZone)                                         const {return GridStretchingConstant[iZone];}
   	/*!
		 * @brief Getter function which returns the value of GridStretchingExponent, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return GridStretchingExponent[iZone]
		 */   
    as3double GetGridStretchingExponent(unsigned short iZone)                                         const {return GridStretchingExponent[iZone];}
   	/*!
		 * @brief Getter function which returns the value of ArtificialConvection, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return ArtificialConvection[iZone]
		 */    
    bool GetArtificialConvection(unsigned short iZone)                                                const {return ArtificialConvection[iZone];}
   	/*!
		 * @brief Getter function which returns the value of ArtificialConvectionConstant, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return ArtificialConvectionConstant[iZone]
		 */
    as3double GetArtificialConvectionConstant(unsigned short iZone)                                   const {return ArtificialConvectionConstant[iZone];}
   	/*!
		 * @brief Getter function which returns the value of ArtificialConvectionExponent, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return ArtificialConvectionExponent[iZone]
		 */
    as3double GetArtificialConvectionExponent(unsigned short iZone)                                   const {return ArtificialConvectionExponent[iZone];}
   	/*!
		 * @brief Getter function which returns the value of PeriodicPulse.
		 *
		 * @return PeriodicPulse
		 */   
    bool GetPeriodicPulse(void)                                                                       const {return PeriodicPulse;}
   	/*!
		 * @brief Getter function which returns the value of TypeBC, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return TypeBC[iZone]
		 */ 
    const as3vector1d<unsigned short> &GetTypeBC(unsigned short iZone)                                const {return TypeBC[iZone];}
   	/*!
		 * @brief Getter function which returns the value of TypeBC, per input zone, per input boundary.
		 *
		 * @param[in] iZone zone ID.
		 * @param[in] iBoundary boundary ID.
		 *
		 * @return TypeBC[iZone][iBoundary]
		 */    
    unsigned short GetTypeBC(unsigned short iZone, unsigned short iBoundary)                          const {return TypeBC[iZone][iBoundary];}
		/*!
		 * @brief Getter function which returns the value of TypeBC.
		 *
		 * @return TypeBC
		 */    
    const as3vector2d<unsigned short> &GetTypeBC(void)                                                const {return TypeBC;}
   	/*!
		 * @brief Getter function which returns the value of InterfaceID, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return InterfaceID[iZone]
		 */    
    const as3vector2d<unsigned short> &GetInterfaceID(unsigned short iZone)                           const {return InterfaceID[iZone];}
   	/*!
		 * @brief Getter function which returns the value of InterfaceID, per input zone, per input boundary.
		 *
		 * @param[in] iZone zone ID.
		 * @param[in] iBoundary boundary ID.
		 *
		 * @return InterfaceID[iZone][iBoundary]
		 */    
    const as3vector1d<unsigned short> &GetInterfaceID(unsigned short iZone, unsigned short iBoundary) const {return InterfaceID[iZone][iBoundary];}
		/*!
		 * @brief Getter function which returns the value of InterfaceID.
		 *
		 * @return InterfaceID
		 */ 
    const as3vector3d<unsigned short> &GetInterfaceID(void)                                           const {return InterfaceID;}
 
	private:
		/*!
		 * @brief Struct that sets the default values for a dictionary file.
		 */    
  	struct {
  		
			unsigned short nZone = 1;                                                ///< Default number of zones.
  		unsigned short IntegrationFactor = 3;                                    ///< Default (over-)integration factor.
  		as3vector1d<unsigned short> nPoly = { 1 };                               ///< Default polynomial order.
      as3vector1d<std::string> NameNodalDOFs = { "LGL" };                      ///< Default Type of DOFs.
  		as3vector1d<std::string> NameZoneMarker = { "ZONE_MAIN" };               ///< Default zone marker.
  		as3vector1d<std::string> NameMarker = {};                                ///< Default marker names.
  		std::string NameRestartSolution = "false";                               ///< Default restart solution option.
  		as3double TimeStep = -1.0;                                               ///< Default time step.
  		unsigned long MaxIter = 100000;                                          ///< Default maximum number of temporal iterations.
      as3vector1d<as3double> hElemRatioZone = { 1.0, 1.0, 1.0, 1.0 };          ///< Default element size ratio in (West, East, South, North).
      std::string NameZoneConformity = "true";                                 ///< Default zone conformity w.r.t. ZONE_MAIN.
      as3vector1d<std::string> NameRiemannSolver = { "ROE" };                  ///< Default type of Riemann solver.
      unsigned long WriteVTKFreq     = 100000;                                 ///< Default output VTK file writing frequency.
      unsigned long WriteRestartFreq = 100000;                                 ///< Default output solution file writing frequency.
      unsigned long OutputFreq       = 100;                                    ///< Default output (runtime-)monitoring frequency.
      as3double CFL = 1.0;                                                     ///< Default CFD number.
      std::string NameAdaptTime = "false";                                     ///< Default adaptive time option.
      as3vector1d<as3double> ParamOutletNSCBC = { 0.1, -1.0, -1.0, 1.0, 0.8 }; ///< Default outlet NSCBC relaxation parameters.
																																							 ///< Dimension: [0]: sigma, [1]: beta_l, [2]: beta_t, [3]: length, [4]: eta.
      as3vector1d<as3double> ParamInletNSCBC = { 1.0, 1.0 };                   ///< Default inlet NSCBC relaxation parameters.
																																							 ///< Dimension: [0]: sigma, [1]: length.
      as3vector1d<as3double> DampingExponent = { 2.0 };                        ///< Default sponge-layer damping exponential coefficient.
      as3vector1d<as3double> DampingConstant = { 2.0 };                        ///< Default sponge-layer damping constant coefficient.
      as3double MachInf = 0.5;                                                 ///< Default free-stream Mach number.
      as3double FlowAngle = 0.0;                                               ///< Default free-stream flow angle w.r.t x-direction (in degrees).
      std::string NameTypeProcessData = "NOTHING";                             ///< Default option to process data.
      as3vector1d<as3double> CenterX0 = { 0.0, 0.0 };                          ///< Default disturbance center coordinates.
      as3double DisturbanceRatio = 0.2;                                        ///< Default disturbance ratio w.r.t. background flow.
      as3double DisturbanceWidth = 0.25;                                       ///< Default disturbance width.
      as3double Frequency = 1.0;                                               ///< Default disturbance frequency.
      as3vector1d<std::string> NameTypeFilterSolution = { "NONE" };            ///< Default type of solution filter.
      as3vector1d<unsigned short> FilterCharacteristics = { 5, 4, 16, 36 };    ///< Default filtering parameters.
																																							 ///< Dimension: [0]: frequency, [1]: Nc, [2]: s, [3]: alpha.
																																							 ///< For info, see Hesthaven DG book or Nordstrom paper on filtering in JSC 2021.
      std::string NamePeriodicPulse = "false";                                 ///< Default option for periodic pulse.
      as3vector1d<std::string> NameTypeBufferLayer = { "NONE" };               ///< Default type of buffer layer.
      as3vector1d<as3double> GridStretchingConstant = { 1.0 };                 ///< Default grid-stretching constant(s).
      as3vector1d<as3double> GridStretchingExponent = { 1.0 };                 ///< Default grid-stretching exponential(s).
      bool GridStretching = false;                                             ///< Default option for grid-stretching.
      as3vector1d<as3double> ArtificialConvectionConstant = { 1.15 };          ///< Default artificial convection constant(s) (max Mach at boundary).
      as3vector1d<as3double> ArtificialConvectionExponent = { 2.0 };           ///< Default artificial convection exponential(s).
      bool ArtificialConvection = false;                                       ///< Default option for artificial convection.
      as3vector1d<as3double> ProbeLocation = { };                              ///< Default location for probes.
																																							 ///< Format: (x1,y1, x2,y2, ..., etc).
      std::string NameSampleZoneData = "false";                                ///< Default option for samping zone data.
      std::string OutputZoneDataFilename = "zone/sample";                      ///< Default zone samping filename.
      std::string NameMarkerZoneData = "ZONE_MAIN";                            ///< Default zone for samping zone data.
      unsigned long WriteFreqZoneData = 100;                                   ///< Default zone samping frequency (in terms of iterations).
      as3vector1d<std::string> NameTypeModifyBC = { "NONE" };                  ///< Default faces of boundaries to modify.
      as3vector1d<as3double> ModifyFreqBC = { 1.0 };                           ///< Default boundary condition frequency modification on each working variable [Hz].
      as3vector1d<as3double> ModifyWidthBC = { 0.25 };                         ///< Default boundary condition modification perturbation width.
      as3vector1d<as3double> ModifyStrengthBC = { 0.2 };                       ///< Default boundary condition modification perturbation amplitude strength.
																																							 ///< Note, these are done w.r.t. background flow of each working variable.
      as3vector1d<as3double> ModifyCenterBC = { -1.0, 0.0 };                   ///< Default boundary condition modification perturbation center(s).
      as3vector1d<as3double> ModifyShiftCenterBC = { 0.25 };                   ///< Default boundary condition modification perturbation center(s) spatial shift.
      std::string NameWriteAuxiliaryDataPML = "false";                         ///< Default option for writing PML auxiliary data.
      std::string NameAlignedPeriodicPulse = "false";                          ///< Default option for an aligned periodic pulse test-case.
      std::string NameConstantFrequency = "true";                              ///< Default option for a constant frequency in a temporal perturbation source term.
			std::string NameSourceTermCenterFixed = "true";                          ///< Default option for a spatially fixed source term in a temporal perturbation.
			as3vector1d<as3double> SourceTermCenterShift = { 1.0, 1.0 };             ///< Default periodic source-term varying spatial center shift values: (dx, dy).
			as3double SourceFrequencyExponent = 1.0;                                 ///< Default periodic source-term frequency amplitude alpha, such that: sin^alpha.
			as3vector1d<as3double> SourceFrequencyParam = { 1.0, 0.1 };              ///< Default periodic source-term frequency paramters, based on: f = f0 * A^(t/tau).
																																							 ///< Dimension: [0]: A, [1]: tau.
      as3vector1d<as3double> DomainExpansionRatio = { 1.0, 1.0, 1.0, 1.0 };    ///< Default domain expansion ratios.
																																							 ///< Format: (rxb1, rxb2, ryb1, ryb2).
      std::string NameUniformGridResolution = "true";                          ///< Default option for equidistant grid resolution.
			as3vector1d<std::string> NameMarkerSampleSurface = { "NONE" };           ///< Default boundary markers for sampling boundary data.
    	std::string OutputSampleSurfaceDirectory = "surf/";                      ///< Default surface samping data directory name.
			std::string NameSampleSurfaceData = "false";                             ///< Default option for surface/boundary data sampling.
			std::string NameFileFormatVTK = "BINARY";                                ///< Default option for VTK output file format.
			std::string NameFileFormatSolution = "BINARY";                           ///< Default option for restart solution file format.
			std::string NameFileFormatZone = "BINARY";                               ///< Default option for zone samping file format.
			std::string NameWriteGNUplot = "false";                                  ///< Default option for GNU-plot output (process-time) file.
			std::string GNUplotFilename = "gnuplot/log";                             ///< Default GNU-plot output filename.
      as3vector1d<as3double> CharacteristicExponent = { 2.0 };                 ///< Default characteristic matching layer exponential coefficient(s).
      as3vector1d<as3double> CharacteristicConstant = { 0.0 };                 ///< Default characteristic matching layer constant coefficient(s).
			std::string NameProbeSpecified = "false";                                ///< Default option for probe data sampling.
			as3vector1d<std::string> NameProbeVariable = { "PRESSURE" };             ///< Default probe data variables sampled.
			std::string NameFileFormatProcessedProbed = "BINARY";                    ///< Default format of processed/probed data files.
			std::string NameProcessLocation = "DOMAIN";                              ///< Default location of data processed.
			as3vector1d<std::string> NameOutputVTKVariable = { };                    ///< Default VTK variables written.
			std::string NameInletNRBC = "true";                                      ///< Default option for a non-reflecting inlet boundary condition (NRBC).

  	} DefaultParam;

    unsigned long  WriteVTKFreq;                       ///< Output VTK file writing frequency. 
    unsigned long  WriteRestartFreq;                   ///< Output solution file writing frequency. 
    unsigned long  OutputFreq;                         ///< Output (runtime-)monitoring frequency.
  	std::string    OutputSolFilename;                  ///< Output solution filename.
  	std::string    OutputVTKFilename;                  ///< Output VTK filename.
    std::string    RestartFilename;                    ///< Input restart solution filename.
		std::string    NameFileFormatVTK;                  ///< VTK output file format string name.
		unsigned short TypeFileFormatVTK;                  ///< VTK output file format.

		as3vector1d<std::string>    NameOutputVTKVariable; ///< Name of VTK variables written in output file.
		as3vector1d<unsigned short> OutputVTKVariable;     ///< VTK variables written in output file.
		
		std::string    NameFileFormatSolution;             ///< Solution output file format string name.
		unsigned short TypeFileFormatSolution;             ///< Solution output file format.
		std::string    NameFileFormatProcessedProbed;      ///< Processed/probed output file format string name.
		unsigned short TypeFileFormatProcessedProbed;      ///< Processed/probed output file format.
		
		bool        WriteGNUplot;                          ///< Option for writing GNU-plot file.     
		std::string GNUplotFilename;                       ///< GNU-plot output filename.

  	unsigned short              nZone;                 ///< Number of zones.
  	as3vector1d<unsigned short> nPolySolZone;          ///< Solution polynomial order per each zone.
    unsigned short              MultizoneStrategy;     ///< Multizone strategy used.

    as3vector1d<as3double>     DomainBound;            ///< (physical) Domain bound, format: (West, East, South, North).
    as3vector1d<unsigned long> nxElemZone;             ///< Number of element in the x-direction in each zone.
    as3vector1d<unsigned long> nyElemZone;             ///< Number of element in the y-direction in each zone.
    as3vector1d<as3double>     hElemRatioZone;         ///< Element geometric expansion ratio per each zone: (West, East, South, Noth).
    bool                       ZoneConformity;         ///< Option to specify a conforming zone.

    as3vector1d<as3double>     BlockInterfaceLocation; ///< Block interface location.
    as3vector1d<as3double>     DomainExpansionRatio;   ///< Domain expantion ratios: (rxb1, rxb2, ryb1. ryb2).
    as3vector1d<unsigned long> nxBlockElem;            ///< Number of elements in x-direction in the main zone.
    as3vector1d<unsigned long> nyBlockElem;            ///< Number of elements in y-direction in the main zone.
    bool                       UniformGridResolution;  ///< Option to use a uniform/equidistant grid resolution.

    as3vector1d<std::string>    NameZoneMarker;        ///< Name of the zone marker. 
    as3vector1d<unsigned short> TypeZone;              ///< Type of each zone.

    as3vector1d<std::string>    NameNodalDOFs;         ///< Type of the DOFs in each zone, string name.
    as3vector1d<unsigned short> TypeDOFs;              ///< Type of DOFs in each zone.
    as3vector1d<unsigned short> TypeSolver;            ///< Type of solver in each zone.

  	std::string    NameTemporalScheme;                 ///< Name of temporal discretization.
  	unsigned short TypeTemporalScheme;                 ///< Type of temporal discretization.
  	as3double      SimulationTime[2];                  ///< Simulation physical starting [0] and ending [1] times.
  	as3double      TimeStep;                           ///< Physical time step.
  	unsigned long  MaxIter;                            ///< Maximum number of temporal iterations.
    as3double      CFL;                                ///< CFL number.
    bool           AdaptTime;                          ///< Option for using an adaptive time step.

  	unsigned short              IntegrationFactor;      ///< (over-)Integration factor.
    as3vector1d<std::string>    NameTypeFilterSolution; ///< Type of solution filtering, string name.
    as3vector1d<unsigned short> TypeFilterSolution;     ///< Type of solution filtering.
    as3vector1d<unsigned short> FilterCharacteristics;  ///< Filter parameters.
																											  ///< Dimension: [0]: frequency, [1]: Nc, [2]: s, [3]: alpha.
																											  ///< For info, see Hesthaven DG book or Nordstrom paper on filtering in JSC 2021.
 
		as3vector1d<std::string>    NameMarkerSampleSurface;      ///< Boundary marker for sampling data, string name.
    std::string                 OutputSampleSurfaceDirectory; ///< Output samping surface data directory.
		bool                        SampleSurfaceData;            ///< Option for samping data from a surface/boundary.
		as3vector1d<unsigned short> SampleDataBoundaryID;         ///< Boundary ID for the sampled surface data.

 		
		as3vector1d<std::string>    NameInitialCondition;         ///< Type of initial condition per each zone, string name.
  	as3vector1d<unsigned short> TypeIC;                       ///< Type of initial condition per each zone.
  	bool                        RestartSolution;              ///< Option for running a restart solution.

    std::string                 NameTypeProcessData;          ///< Type of data processing, string name.
    unsigned short              TypeProcessData;              ///< Type of data processing.
		std::string                 NameProcessLocation;          ///< Location of the data processed, string name.
		unsigned short              ProcessLocation;              ///< Location of the data processed.
    std::string                 OutputProcessedDirectory;     ///< Processed data output directory.
		bool                        ProbeSpecified;               ///< Option to use a probe.
    as3vector1d<as3double>      ProbeLocation;                ///< Probe(s) location. Format: (x1,y1, x2,y2, ...etc).
		as3vector1d<std::string>    NameProbeVariable;            ///< Name of the probe variables.
		as3vector1d<unsigned short> ProbeVariable;                ///< Probe variables.
    
    std::string                 OutputZoneDataFilename;       ///< Zone samping data output filename.
    std::string                 NameMarkerZoneData;           ///< Zone marker for sampling zone data.
    unsigned short              TypeZoneData;                 ///< Type of zone selected for samping zone data.
    bool                        UsePML;                       ///< Option to use a PML.
    bool                        WriteAuxiliaryDataPML;        ///< Write auxiliary data in a PML layer.
    unsigned long               WriteFreqZoneData;            ///< Writing frequency of output zone sampling.
		std::string                 NameFileFormatZone;           ///< Type of zone sampling file format, string name.
		unsigned short              TypeFileFormatZone;           ///< Type of zone sampling file format.
    bool                        SampleZoneData;               ///< Option for samping zone data.

    as3vector1d<bool>      GridStretching;                    ///< Option to use grid-stretching in each zone.
    as3vector1d<as3double> GridStretchingConstant;            ///< Grid-stretching constant per each zone.
    as3vector1d<as3double> GridStretchingExponent;            ///< Grid-stretching exponent per each zone.

    as3double              MachInf;                           ///< Free-stream Mach number.
    as3double              FlowAngle;                         ///< Free-stream flow angle w.r.t. x-direction (in degrees).
    bool                   CrossFlow;                         ///< Option for the existence of cross(oblique)-flow.
    
    as3vector1d<as3double> CenterX0;                          ///< Disturbance source center.
    as3double              DisturbanceRatio;                  ///< Disturbance strength ratio w.r.t. background flow.
    as3double              DisturbanceWidth;                  ///< Disturbance width.
    as3double              Frequency;                         ///< Disturbance frequency [Hz].
    as3double              AngularFrequency;                  ///< Disturbance angular frequency.
    bool                   PeriodicPulse;                     ///< Option for using a periodic pulse disturbance.
    bool                   ConstantFrequency;                 ///< Option for having a constant frequency in the disturbance.
    bool                   AlignedPeriodicPulse;              ///< Option for having a periodic flow direction aligned with the flow.
		as3vector1d<as3double> SourceFrequencyParam;              ///< Periodic source-term frequency paramters, based on: f = f0 * A^(t/tau).
																															///< Dimension: [0]: A, [1]: tau. 
		as3vector1d<as3double> SourceTermCenterShift;             ///< Periodic source-term varying spatial center shift values: (dx, dy).
		as3double              SourceFrequencyExponent;           ///< Periodic source-term frequency amplitude alpha, such that: sin^alpha.
		bool                   SourceTermCenterFixed;             ///< Option for a spatially fixed source term in a temporal perturbation.

    as3vector1d<std::string>    NameTypeModifyBC;             ///< Type of boundary conditions to modify per each face, string name.
    as3vector1d<unsigned short> TypeModifyBC;                 ///< Type of boundary conditions to modify per each face.
    as3vector1d<as3double>      ModifyFreqBC;                 ///< Angular frequency in the modified boundary condition term.
    as3vector1d<as3double>      ModifyWidthBC;                ///< Pulse width in the modified boundary condition of each working variable.
    as3vector1d<as3double>      ModifyStrengthBC;             ///< Pulse amplitude ratio w.r.t. background flow of the modified boundary condition of each working variable.
    as3vector1d<as3double>      ModifyCenterBC;               ///< Pulse center location(s) in the modified boundary conditions.
    as3vector1d<as3double>      ModifyShiftCenterBC;          ///< Pulse center(s) shift in the modified boundary conditions.

    as3vector1d<std::string>    NameTypeBufferLayer;          ///< Type of buffer layer, string name.
    as3vector1d<unsigned short> TypeBufferLayer;              ///< Type of buffer layer.

    as3vector1d<std::string>    NameRiemannSolver;            ///< Type of Riemann solver per each zone, string name.
    as3vector1d<unsigned short> RiemannSolver;                ///< Type of Riemann solver per each zone.

    as3vector1d<std::string>    NameBoundaryCondition;        ///< Type of boundary conditions specified, string name.
																															///< Format: (SOUTH, NORTH, WEST, EAST).
    as3vector1d<unsigned short> TypeExternalBC;               ///< Type of external or terminating boundary conditions over all zones.
																															///< Format: (SOUTH, NORTH, WEST, EAST).
    as3vector2d<unsigned short> TypeBC;                       ///< Type of boundary conditions per each zone per each face.
																															///< Dimension: [iZone][iBoundary], where iBoundary has indices: (SOUTH, NORTH, WEST, EAST).
    as3vector3d<unsigned short> InterfaceID;                  ///< Interface ID for periodic/zonal interface boundaries.
																															///< Dimension: [iZone][iBoundary][iData], with indices: [iData] = [0]: jZone, [1]: jFace.

    as3vector1d<as3double> ParamOutletNSCBC;                  ///< Relaxation parameters in the case of an outlet NSCBC.
																															///< Dimension: [0]: sigma, [1]: beta_l, [2]: beta_t, [3]: len, [4]: eta.
    as3vector1d<as3double> ParamInletNSCBC;                   ///< Relaxation parameters in the case of an inlet NSCBC.
																															///< Dimension: [0]: sigma, [1]: len. 
		bool                   InletNRBC;                         ///< Option for using a non-reflective NSCBC inlet (NRBC).

    as3vector1d<as3double> DampingExponent;                   ///< Sponge-layer damping exponential coefficient(s).
    as3vector1d<as3double> DampingConstant;                   ///< Sponge-layer damping constant coefficient(s).


    as3vector1d<bool>      CharacteristicMatching;            ///< Option for using a characteristic matching layer.
    as3vector1d<as3double> CharacteristicExponent;            ///< Characteristic matching layer exponential coefficient(s).
    as3vector1d<as3double> CharacteristicConstant;            ///< Characteristic matching layer constant coefficient(s).

    as3vector1d<as3double> ArtificialConvectionConstant;      ///< Artificial-convection constant coefficient(s) (max Mach at boundary).
    as3vector1d<as3double> ArtificialConvectionExponent;      ///< Artificial-convection exponential coefficient(s).
    as3vector1d<bool>      ArtificialConvection;              ///< Option for using an artificial-convection (supersonic) layer.

		/*!
		 * @brief Function that maps TypeZoneData from string to enum.
		 */
    void MapTypeZoneData(void);
    /*!
		 * @brief Function that maps TypeDOFs from string to enum.
		 */
    void MapTypeDOFs(void);
    /*!
		 * @brief Function that maps NameMarkerZoneData from string to enum.
		 */
    void MapTypeZone(void);
    /*!
		 * @brief Function that maps RiemannSolver from string to enum.
		 */
    void MapRiemannSolver(void);
    /*!
		 * @brief Function that maps TemporalScheme from string to enum.
		 */
    void MapTemporalScheme(void);
    /*!
		 * @brief Function that maps TypeIC from string to enum.
		 */
    void MapTypeIC(void);
    /*!
		 * @brief Function that maps TypeExternalBC from string to enum.
		 */
    void MapTypeExternalBC(void);
    /*!
		 * @brief Function that maps TypeProcessData from string to enum.
		 */
    void MapTypeProcessData(void);
    /*!
		 * @brief Function that maps TypeFilterSolution from string to enum.
		 */
    void MapTypeFilterSolution(void);
    /*!
		 * @brief Function that maps TypeBufferLayer from string to enum.
		 */
    void MapTypeBufferLayer(void);
    /*!
		 * @brief Function that maps TypeModifyBC from string to enum.
		 */
    void MapTypeModifyBC(void);
		/*!
		 * @brief Function that maps NameMarkerSampleSurface from string to enum.
		 */
		void MapSampleDataBoundaryID(void);
		/*!
		 * @brief Function that maps NameFileFormatVTK from string to enum.
		 */
		void MapTypeFileFormatVTK(void);
		/*!
		 * @brief Function that maps NameFileFormatSolution from string to enum.
		 */
		void MapTypeFileFormatSolution(void);
		/*!
		 * @brief Function that maps NameFileFormatZone from string to enum.
		 */
		void MapTypeFileFormatZone(void);
		/*!
		 * @brief Function that maps ProbeVariable from string to enum.
		 */
		void MapTypeProbeVariable(void);
		// Function that maps NameFileFormatProcessedProbed from string to enum.
		void MapTypeFileFormatProcessedProbed(void);
		// Function that maps NameProcessLocation from string to enum.
		void MapTypeProcessLocation(void);
		// Function that maps NameOutputVTKVariable from string to enum.
		void MapTypeOutputVTKVariable(void);

    /*!
		 * @brief Function that determines the multizone strategy adopted.
		 */
    void DetermineMultizoneStrategy(void);
    /*!
		 * @brief Function that checks the element ratios specified according to specified multizone strategy.
		 */
    void CheckElementRatio(void);
    /*!
		 * @brief Function that processes the zone conformity in the specified elements.
		 */
    void ProcessZoneConformity(void);
    /*!
		 * @brief Function that processes the boundary conditions per each zone.
		 */
    void ProcessBoundaryConditions(void);
    /*!
		 * @brief Function that matches each iZone iBoundary with its counterpart in jZone jBoundary.
		 *
		 * @param[in] iZone zone ID.
		 * @param[in] iFace boundary/face ID.
		 *
		 * @return the corresponding matching face which interconnects the zones and boundaries.
		 */
    as3vector1d<unsigned short> MatchInterface(unsigned short iZone,
                                               unsigned short iFace);

    /*!
		 * @brief Function that reads the Integration rule data.
		 *
		 * @param[in] configFile input configuration/dictionary file. 
		 */
    void ReadIntegrationRuleOptions(const char *configFile);
    /*!
		 * @brief Function that reads the filtering data.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
    void ReadFilteringOptions(const char *configFile);
    /*!
		 * @brief Function that reads the buffer-layer data.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
    void ReadBufferLayerOptions(const char *configFile);
    /*!
		 * @brief Function that reads the boundary conditions data.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
    void ReadBoundaryConditionOptions(const char *configFile);
    /*!
		 * @brief Function that reads the sponge-damping data.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
    void ReadSpongeDampingOptions(const char *configFile);
    /*!
		 * @brief Function that reads the grid-stretching data.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
    void ReadGridStretchingOptions(const char *configFile);
    /*! @brief Function that reads the artificial-convection data.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
    void ReadArtificialConvectionOptions(const char *configFile);
    /*!
		 * @brief Function that reads the NSCBC data.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
    void ReadCharacteristicBoundaryOptions(const char *configFile);
		/*!
		 * @brief Function that reads the characteristic layer data.
		 *
		 * @param[in] configFile input configuration/dictionary file.
		 */
		void ReadCharacteristicLayerOptions(const char *configFile);

    /*!
		 * @brief Function that reads grid options.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
    bool ReadGridOptions(const char *configFile);
    /*!
		 * @brief Function that reads input/output options.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
  	bool ReadIOOptions(const char *configFile);
  	/*!
		 * @brief Function that reads solver options.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
  	bool ReadSolverOptions(const char *configFile);
  	/*!
		 * @brief Function that reads boundary options.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
  	bool ReadBoundaryOptions(const char *configFile);
  	/*!
		 * @brief Function that reads temporal options.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
  	bool ReadTemporalOptions(const char *configFile);
  	/*!
		 * @brief Function that reads initial condition options.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
  	bool ReadICOptions(const char *configFile);
    /*!
		 * @brief Function that reads flow characteristics options.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
    bool ReadFlowOptions(const char *configFile);
    /*!
		 * @brief Function that reads processing information.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
    bool ReadProcessingOptions(const char *configFile);
    /*!
		 * @brief Function that reads modified boundary condition information.
		 *
		 * @return bool whether operation failed/succeeded.
		 */
    bool ReadModifiedBCOptions(const char *configFile);
};





