#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "spatial_structure.hpp"
#include "output_structure.hpp"
#include "initial_structure.hpp"


// Forward declaration to avoid compiler problems.
class CSensorData;


class CProcess {

  public:
    // Constructor.
    CProcess(CConfig       *config_container,
             CGeometry     *geometry_container,
             COutput       *output_container,
             CElement     **element_container,
             CInitial     **initial_container,
             CSolver      **solver_container,
             CSpatial     **spatial_container,
             unsigned short iZone);

    // Destructor.
    virtual ~CProcess(void);

    // Function that filters the solution in this zone.
    void FilterSolution(CConfig   *config_container,
                        CGeometry *geometry_container,
                        CElement  *element_container,
                        CSolver   *solver_container,
                        CSpatial  *spatial_container);

    // Pure virtual function that selects what/how to process the data.
    // Must be overriden by a derived class.
    virtual void ProcessData(CConfig   *config_container,
                             CGeometry *geometry_container,
                             CElement  *element_container,
                             CSpatial  *spatial_container,
                             CSolver   *solver_container,
                             CInitial  *initial_container,
                             as3double  localTime) = 0;

    // Pure virtual function that decides what data to write as output.
    // Must be overriden by a derived class.
    virtual void WriteProcessedData(CConfig   *config_container,
                                    CGeometry *geometry_container,
                                    COutput   *output_container) = 0;

  protected:
    // Zone ID.
    unsigned short zoneID;

  private:

};


class CEEProcess : public CProcess {

  public:
    // Constructor.
    CEEProcess(CConfig       *config_container,
               CGeometry     *geometry_container,
               COutput       *output_container,
               CElement     **element_container,
               CInitial     **initial_container,
               CSolver      **solver_container,
               CSpatial     **spatial_container,
               unsigned short iZone);

    // Destructor.
    ~CEEProcess(void) override;

    // Function that selects what/how to process the data.
    void ProcessData(CConfig   *config_container,
                     CGeometry *geometry_container,
                     CElement  *element_container,
                     CSpatial  *spatial_container,
                     CSolver   *solver_container,
                     CInitial  *initial_container,
                     as3double  localTime) override;

    // Function that decides what data to write as output.
    void WriteProcessedData(CConfig   *config_container,
                            CGeometry *geometry_container,
                            COutput   *output_container);

  protected:
    // Function that computes the reflection coeffcient as a function of (P).
    // Expression: R = max( p'(t) )/pinf.
    void ComputeReflectionCoefficientFuncP(CConfig   *config_container,
                                           CGeometry *geometry_container,
                                           CElement  *element_container,
                                           CSpatial  *spatial_container,
                                           CSolver   *solver_container,
                                           CInitial  *initial_container,
                                           as3double  localTime);

    // Function that computes the reflection coeffcient as a function of (P, T).
    // Expression: R = ( max( p'(t) )/pinf ) / ( max( T'(t=0) )/ Tinf ).
    void ComputeReflectionCoefficientFuncPT(CConfig   *config_container,
                                            CGeometry *geometry_container,
                                            CElement  *element_container,
                                            CSpatial  *spatial_container,
                                            CSolver   *solver_container,
                                            CInitial  *initial_container,
                                            as3double  localTime);

    // Function that computes the reflection coeffcient as a function of (P, M).
    // Expression: R = ( max( p'(t) )/pinf ) / ( max( M'(t) )/ Minf ).
    void ComputeReflectionCoefficientFuncPM(CConfig   *config_container,
                                            CGeometry *geometry_container,
                                            CElement  *element_container,
                                            CSpatial  *spatial_container,
                                            CSolver   *solver_container,
                                            CInitial  *initial_container,
                                            as3double  localTime);

    // Function that samples data from probe points of a domain.
    void SampleProbePoints(CConfig   *config_container,
                           CGeometry *geometry_container,
                           CElement  *element_container,
                           CSpatial  *spatial_container,
                           CSolver   *solver_container,
                           CInitial  *initial_container,
                           as3double  localTime);

  private:
    // Temporal signature.
    as3vector1d<as3double> TimeProcessed;
    // Processed data.
    as3vector2d<as3double> DataProcessed;
    // Sensor data for probing a location, if needed.
    as3data1d<CSensorData> sensor_container;
};


class CSensorData {

  public:
    // Constructor.
    CSensorData(CConfig               *config_container,
                CGeometry             *geometry_container,
                CElement              *element_container,
                as3vector1d<as3double> probe,
                unsigned short         iZone);
    // Destructor.
    ~CSensorData(void);

    // Getter: returns IndexZone.
    unsigned short GetIndexProbeZone(void)               const {return IndexProbeZone;}
    // Getter: returns IndexProbeElem.
    unsigned long GetIndexProbeElem(void)                const {return IndexProbeElem;}
    // Getter: returns ProbeLocation.
    const as3vector1d<as3double> &GetProbeLocation(void) const {return ProbeLocation;}
    // Getter: returns Interpolation.
    const as3vector1d<as3double> &GetInterpolation(void) const {return Interpolation;}

  private:
    // Zone ID of element holding current probe.
    unsigned short IndexProbeZone;
    // Global element index of the probe ownership.
    unsigned long IndexProbeElem;
    // Probe coordinates.
    as3vector1d<as3double> ProbeLocation;
    // Lagrange interpolation over this element.
    as3vector1d<as3double> Interpolation;

    // Function that computes the Lagrange interpolation entry in 1D for a given basis.
    as3double EvaluateEntryLagrangePolynomial(const as3double               xPoint,
                                              const as3vector1d<as3double> &rBasis,
                                              const unsigned short          iDegree);

};


