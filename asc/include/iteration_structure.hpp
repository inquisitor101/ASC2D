#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "spatial_structure.hpp"
#include "initial_structure.hpp"
#include "blas.hpp"



class CIteration {

	public:
		// Constructor.
		CIteration(CConfig  		 *config_container,
							 CGeometry		 *geometry_container,
							 CSolver  		**solver_container,
							 CElement 		**element_container,
							 unsigned short iZone);

		// Destructor.
		virtual ~CIteration(void);

		// Pure virtual function that preprocesses every iteration.
		// Must be overriden by a derived class.
		virtual void Preprocess(CConfig              *config_container,
														CGeometry            *geometry_container,
														CSolver             **solver_container,
														CElement            **element_container,
														CSpatial            **spatial_container,
                            as3data1d<as3double> &work_array,
														as3double             localTime,
                            unsigned long         iElem) = 0;

    // Pure virtual function that sweeps through a single element
    // in this zone and updates the residual accordingly.
    // Must be overriden by a derived class.
    virtual void ComputeResidual(CConfig                *config_container,
                                 CGeometry              *geometry_container,
                                 CSolver                *solver_container,
                                 CElement               *element_container,
                                 CSpatial               *spatial_container,
                                 CInitial               *initial_container,
                                 as3data1d<as3double>   &work_array,
                                 as3double               localTime,
                                 unsigned long           iElem,
                                 as3vector1d<as3double> &MonitoringData) = 0;

	protected:
		// Zone ID.
		unsigned short zoneID;
    // Number of elements in this zone.
    unsigned long  nElem;

	private:

};


class CEEIteration : public CIteration {

	public:
		// Constructor.
		CEEIteration(CConfig  		 *config_container,
							 	 CGeometry		 *geometry_container,
							 	 CSolver  		**solver_container,
							 	 CElement 		**element_container,
							 	 unsigned short iZone);

		// Destructor.
		~CEEIteration(void) override;

	protected:
		// Function that preprocesses every iteration.
		void Preprocess(CConfig              *config_container,
										CGeometry            *geometry_container,
										CSolver             **solver_container,
										CElement            **element_container,
										CSpatial            **spatial_container,
                    as3data1d<as3double> &work_array,
										as3double             localTime,
                    unsigned long         iElem) override;

    // Function that sweeps through a single element
    // in this zone and updates the residual accordingly.
    void ComputeResidual(CConfig                *config_container,
                         CGeometry              *geometry_container,
                         CSolver                *solver_container,
                         CElement               *element_container,
                         CSpatial               *spatial_container,
                         CInitial               *initial_container,
                         as3data1d<as3double>   &work_array,
                         as3double               localTime,
                         unsigned long           iElem,
                         as3vector1d<as3double> &MonitoringData) override;

	private:

};


