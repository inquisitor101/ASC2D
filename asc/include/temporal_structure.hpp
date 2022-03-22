#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "solver_structure.hpp"
#include "iteration_structure.hpp"
#include "spatial_structure.hpp"
#include "initial_structure.hpp"



class CTemporal {

	public:
		// Constructor.
		CTemporal(CConfig                    *config_container,
							CGeometry                  *geometry_container,
							CIteration                **iteration_container,
							CSolver                   **solver_container,
							CElement                  **element_container,
							CSpatial                  **spatial_container,
              as3vector2d<unsigned long> &input_MapGlobalToLocal);

		// Destructor.
		virtual ~CTemporal(void);

		// Pure virtual Function that performs an update over the entire simulation in time.
		// Must be overriden by one of the derived classes.
		virtual void TimeMarch(CConfig                *config_container,
									 				 CGeometry              *geometry_container,
									 				 CIteration            **iteration_container,
									 				 CSolver               **solver_container,
									 				 CElement              **element_container,
									 				 CSpatial              **spatial_container,
                           CInitial              **initial_container,
									 				 as3double               physicalTime,
													 as3double               dtTime,
                           as3vector1d<as3double> &MonitoringData) = 0;

	protected:
		// Number of zones.
		unsigned short nZone;
		// Number of elements, per zone.
		as3vector1d<unsigned long>  nElemZone;
		// Number of solution nodes in each element per zone.
		as3vector1d<unsigned short> nNodeZone;

    // Total number of elements in all zones combined.
    unsigned long nElemTotal;

    // Number of degrees-of-freedom in the working array.
    unsigned short nWorkingArrayDOFs;
    // Number of variables needed in the working array.
    unsigned short nWorkingArrayVar;
    // Number of data entries needed in the working array.
    unsigned short nWorkingArrayEntry;

    // Map data from [iZone][iElemZone] to [iElem]. This is use for
    // parallelization efficiency. Dimension: [iElem][iData], where:
    // [iData] is: [0]: iZone, [1]: iElemZone.
    as3vector2d<unsigned long> MapGlobalToLocal;

    // Function that initializes and defined the dimension parameters needed
    // in the working array.
    void InitializeWorkArrayDimension(CConfig   *config_container,
                                      CElement **element_container);

    // Function that initializes the working array by reserving its memory.
    void InitializeWorkArray(as3data1d<as3double> &work_array);

	private:

};


class CLSRK4Temporal : public CTemporal {

	public:
		// Constructor.
		CLSRK4Temporal(CConfig                    *config_container,
									 CGeometry                  *geometry_container,
									 CIteration                **iteration_container,
									 CSolver                   **solver_container,
									 CElement                  **element_container,
									 CSpatial                  **spatial_container,
                   as3vector2d<unsigned long> &input_MapGlobalToLocal);

		// Destructor.
		~CLSRK4Temporal(void) final;

		// Function that performs an update over the entire simulation in time.
		void TimeMarch(CConfig                *config_container,
					 				 CGeometry              *geometry_container,
					 				 CIteration            **iteration_container,
					 				 CSolver               **solver_container,
					 				 CElement              **element_container,
					 				 CSpatial              **spatial_container,
                   CInitial              **initial_container,
					 				 as3double               physicalTime,
									 as3double               dtTime,
                   as3vector1d<as3double> &MonitoringData) final;

	protected:

	private:
		// Number of RK stages: this is a 5-stage scheme.
		unsigned short nStageRK = LSRK4_N_STAGES;
		// Number of storage for tentative data.
		unsigned short nStorage = LSRK4_N_STORAGE;

		// LSRK4 coefficients.
		as3vector1d<as3double> rk4a;
		as3vector1d<as3double> rk4b;
		as3vector1d<as3double> rk4c;

		// Tentative/intermediate solution stored.
		// Dimension: [iZone][iElem][iVar][iNode].
		as3data3d<as3double> DataDOFsSolTentative;

    // Function that performs a single stage sweep of a LSRK4.
		void UpdateTime(CConfig                *config_container,
									  CGeometry              *geometry_container,
									  CIteration            **iteration_container,
									  CSolver               **solver_container,
									  CElement              **element_container,
									  CSpatial              **spatial_container,
                    CInitial              **initial_container,
									  as3double               localTime,
                    as3double               dtTime,
                    as3double               alpha,
                    as3double               beta,
                    as3vector1d<as3double> &MonitoringData);
};


class CSSPRK3Temporal : public CTemporal {

	public:
		// Constructor.
		CSSPRK3Temporal(CConfig                    *config_container,
									  CGeometry                  *geometry_container,
									  CIteration                **iteration_container,
									  CSolver                   **solver_container,
									  CElement                  **element_container,
									  CSpatial                  **spatial_container,
                    as3vector2d<unsigned long> &input_MapGlobalToLocal);

		// Destructor.
		~CSSPRK3Temporal(void) final;

		// Function that performs an update over the entire simulation in time.
		void TimeMarch(CConfig                *config_container,
					 				 CGeometry              *geometry_container,
					 				 CIteration            **iteration_container,
					 				 CSolver               **solver_container,
					 				 CElement              **element_container,
					 				 CSpatial              **spatial_container,
                   CInitial              **initial_container,
					 				 as3double               physicalTime,
									 as3double               dtTime,
                   as3vector1d<as3double> &MonitoringData) final;

	protected:

	private:
		// Number of RK stages: this is a 3-stage scheme.
		unsigned short nStageRK = SSPRK3_N_STAGES;
		// Number of storage for tentative data.
		unsigned short nStorage = SSPRK3_N_STORAGE;

		// LSRK4 coefficients.
		as3vector1d<as3double> rk4a;
		as3vector1d<as3double> rk4b;
		as3vector1d<as3double> rk4c;

		// Tentative/intermediate solution stored.
		// Dimension: [iZone][iElem][iVar][iNode].
		as3data3d<as3double> DataDOFsSolTentative;

    // Function that performs a single stage sweep of a SSPRK3.
		void UpdateTime(CConfig                *config_container,
									  CGeometry              *geometry_container,
									  CIteration            **iteration_container,
									  CSolver               **solver_container,
									  CElement              **element_container,
									  CSpatial              **spatial_container,
                    CInitial              **initial_container,
									  as3double               localTime,
                    as3double               dtTime,
                    as3double               alpha,
                    as3double               beta,
                    as3vector1d<as3double> &MonitoringData);
};

