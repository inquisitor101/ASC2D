#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "riemann_structure.hpp"
#include "data_structure.hpp"
#include "initial_structure.hpp"


class CSpatial {

	public:
		// Constructor.
		CSpatial(CConfig   	   *config_container,
						 CGeometry 	   *geometry_container,
						 CElement  	  **element_container,
             CInitial      *initial_container,
						 unsigned short iZone);

		// Destructor.
		virtual ~CSpatial(void);

    // Function that applies a filter to the solution.
    void ComputeFilteredSolution(CConfig   *config_container,
                                 CGeometry *geometry_container,
                                 CElement  *element_container,
                                 CData     *data_container);

    // Pure virtual function that initializes the solution needed.
    // Note, this must be overriden by a derived class.
    virtual void InitializeSolution(CConfig                *config_container,
                                    CInitial               *initial_container,
                                    CElement               *element_container,
                                    const CGeometryElement *grid_element,
                                    CData                  *data_element,
                                    as3double               time) = 0;

    // Pure virtual function that computes the volume term contribution to
    // the residual. Note, this must be overriden by a derived class.
    virtual void ComputeVolumeResidual(CConfig                *config_container,
                                       CGeometry              *geometry_container,
                                       CElement               *element_container,
                                       CInitial               *initial_container,
                                       CData                  *data_container,
                                       const CGeometryElement *geometry_element,
                                       as3data1d<as3double>   &work_array,
                                       as3double               localTime,
                                       as3vector1d<as3double> &MonitoringData) = 0;

    // Pure virtual function that computes the surface term contribution to
    // the residual. Note, this must be overriden by a derived class.
    virtual void ComputeSurfaceResidual(CConfig              *config_container,
                                        CGeometry            *geometry_container,
                                        CElement             *element_container,
                                        as3element           &data_container,
                                        as3data1d<as3double> &work_array,
                                        as3double             localTime,
                                        unsigned long         iElem) = 0;

	protected:
    // Zone ID.
    unsigned short zoneID;

		// Number of solution DOFs per element (1D) in this zone.
		unsigned short nDOFsSol1D;
    // Number of solution DOFs per element (2D) in this zone.
    unsigned short nDOFsSol2D;
		// Number of integration DOFs per element (1D) in this zone.
		unsigned short nDOFsInt1D;
    // Number of integration DOFs per element (2D) in this zone.
    unsigned short nDOFsInt2D;

    // Filter matrix, if needed.
    as3double *FilterMatrix = nullptr;

    // Riemann container.
    CRiemann *riemann_container = nullptr;

    // Pure virtual function that preprocesses the riemann container.
    // Not, this must be overriden by a derived class.
    virtual void Riemann_Preprocessing(CConfig       *config_container,
                                       unsigned short iZone) = 0;

    // Function that computes the filter matrix: F = V*C*Vinv.
    void InitializeFilterMatrix(CConfig   *config_container,
                                CGeometry *geometry_container,
                                CElement  *element_container,
                                as3double *matF);

	private:

};


class CEESpatial : public CSpatial {

	public:
		// Constructor.
		CEESpatial(CConfig   	   *config_container,
						 	 CGeometry 	   *geometry_container,
						 	 CElement  	  **element_container,
               CInitial      *initial_container,
						 	 unsigned short iZone);

		// Destructor.
		~CEESpatial(void);

    // Function that initializes the solution needed.
    void InitializeSolution(CConfig                *config_container,
                            CInitial               *initial_container,
                            CElement               *element_container,
                            const CGeometryElement *grid_element,
                            CData                  *data_element,
                            as3double               time) override;

    // Function that computes the volume term contribution to the residual.
    void ComputeVolumeResidual(CConfig                *config_container,
                               CGeometry              *geometry_container,
                               CElement               *element_container,
                               CInitial               *initial_container,
                               CData                  *data_container,
                               const CGeometryElement *geometry_element,
                               as3data1d<as3double>   &work_array,
                               as3double               localTime,
                               as3vector1d<as3double> &MonitoringData) override;

    // Function that computes the surface term contribution to
    // the residual. Note, this must be overriden by a derived class.
    void ComputeSurfaceResidual(CConfig              *config_container,
                                CGeometry            *geometry_container,
                                CElement             *element_container,
                                as3element           &data_container,
                                as3data1d<as3double> &work_array,
                                as3double             localTime,
                                unsigned long         iElem) override;

	protected:
    // Function that preprocesses the riemann container.
    void Riemann_Preprocessing(CConfig       *config_container,
                               unsigned short iZone) final;

	private:

};



class CEESpongeSpatial : public CEESpatial {

	public:
		// Constructor.
		CEESpongeSpatial(CConfig   	    *config_container,
      						 	 CGeometry 	    *geometry_container,
      						 	 CElement   	 **element_container,
                     CInitial       *initial_container,
      						 	 unsigned short  iZone);

		// Destructor.
		~CEESpongeSpatial(void);

    // Function that initializes the solution needed.
    void InitializeSolution(CConfig                *config_container,
                            CInitial               *initial_container,
                            CElement               *element_container,
                            const CGeometryElement *grid_element,
                            CData                  *data_element,
                            as3double               time) final;

    // Function that computes the volume term contribution to the residual.
    void ComputeVolumeResidual(CConfig                *config_container,
                               CGeometry              *geometry_container,
                               CElement               *element_container,
                               CInitial               *initial_container,
                               CData                  *data_container,
                               const CGeometryElement *geometry_element,
                               as3data1d<as3double>   &work_array,
                               as3double               localTime,
                               as3vector1d<as3double> &MonitoringData) final;

    // Function that computes the surface term contribution to
    // the residual. Note, this must be overriden by a derived class.
    void ComputeSurfaceResidual(CConfig              *config_container,
                                CGeometry            *geometry_container,
                                CElement             *element_container,
                                as3element           &data_container,
                                as3data1d<as3double> &work_array,
                                as3double             localTime,
                                unsigned long         iElem) final;

	protected:

	private:
    // Whether or not an artificial-convection is needed on each face.
    as3vector1d<bool> ArtificialConvectionFace;
};


class CEEPMLSpatial : public CEESpatial {

	public:
		// Constructor.
		CEEPMLSpatial(CConfig   	   *config_container,
    						 	CGeometry 	   *geometry_container,
    						 	CElement  	  **element_container,
                  CInitial       *initial_container,
    						 	unsigned short  iZone);

		// Destructor.
		~CEEPMLSpatial(void);

    // Function that initializes the solution needed.
    void InitializeSolution(CConfig                *config_container,
                            CInitial               *initial_container,
                            CElement               *element_container,
                            const CGeometryElement *grid_element,
                            CData                  *data_element,
                            as3double               time) final;

    // Function that computes the volume term contribution to the residual.
    void ComputeVolumeResidual(CConfig                *config_container,
                               CGeometry              *geometry_container,
                               CElement               *element_container,
                               CInitial               *initial_container,
                               CData                  *data_container,
                               const CGeometryElement *geometry_element,
                               as3data1d<as3double>   &work_array,
                               as3double               localTime,
                               as3vector1d<as3double> &MonitoringData) final;

    // Function that computes the surface term contribution to
    // the residual. Note, this must be overriden by a derived class.
    void ComputeSurfaceResidual(CConfig              *config_container,
                                CGeometry            *geometry_container,
                                CElement             *element_container,
                                as3element           &data_container,
                                as3data1d<as3double> &work_array,
                                as3double             localTime,
                                unsigned long         iElem) final;

	protected:

	private:
    // Dispersion-relation correction coefficient.
    as3double DispersionCorrection;
    // Transverse background velocity (i.e. in y-direction).
    as3double VelocityTransverse;

};
