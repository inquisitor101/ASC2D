#include "iteration_structure.hpp"



CIteration::CIteration
(
 CConfig  		 *config_container,
 CGeometry		 *geometry_container,
 CSolver  		**solver_container,
 CElement 		**element_container,
 unsigned short iZone
)
 /*
	* Constructor, used to initialize CIteration in zone: iZone.
	*/
{
	// Assign zone ID.
	zoneID = iZone;
  // Extract number of elements in this zone.
  nElem  = geometry_container->GetGeometryZone(iZone)->GetnElem();
}


CIteration::~CIteration
(
 void
)
 /*
	* Destructor for CIteration class, frees allocated memory.
	*/
{

}


CEEIteration::CEEIteration
(
 CConfig  		 *config_container,
 CGeometry		 *geometry_container,
 CSolver  		**solver_container,
 CElement 		**element_container,
 unsigned short iZone
)
	:
		CIteration
		(
		 config_container,
		 geometry_container,
		 solver_container,
		 element_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CEEIteration in zone: iZone.
	*/
{

}


CEEIteration::~CEEIteration
(
 void
)
 /*
	* Destructor for CEEIteration class, frees allocated memory.
	*/
{

}


void CEEIteration::Preprocess
(
 CConfig             *config_container,
 CGeometry           *geometry_container,
 CSolver            **solver_container,
 CElement           **element_container,
 CSpatial           **spatial_container,
 as3data1d<as3double> &work_array,
 as3double             localTime,
 unsigned long         iElem
)
 /*
	* Function that performs a preprocessing step for every iteration of the
  * EE solver.
	*/
{

}


void CEEIteration::ComputeResidual
(
 CConfig               *config_container,
 CGeometry             *geometry_container,
 CSolver               *solver_container,
 CElement              *element_container,
 CSpatial              *spatial_container,
 CInitial              *initial_container,
 as3data1d<as3double>  &work_array,
 as3double              localTime,
 unsigned long          iElem,
 as3vector1d<as3double> &MonitoringData
)
 /*
	* Function that performs a single element sweep across space to update
  * the residual in an EE solver.
	*/
{
  // Extract number of solution DOFs in 2D.
  unsigned short nDOFsSol2D = element_container->GetnDOFsSol2D();
  // Load inverse mass matrix.
  auto* InvMassMatrix  = solver_container->GetInvMassMatrix();
  // Load the grid of the element.
  auto* geometry_element = geometry_container->GetGeometryZone(zoneID)->GetGeometryElem(iElem);
  // Load overall data container for all elements.
  auto& data_container = solver_container->GetDataContainer();

  // Initialize residual to zero.
  data_container[iElem]->ResetResidual();

  // Step 1: Compute all volume terms contribution to the residual.
  spatial_container->ComputeVolumeResidual(config_container,
                                           geometry_container,
                                           element_container,
                                           initial_container,
                                           data_container[iElem],
                                           geometry_element,
                                           work_array,
                                           localTime,
                                           MonitoringData);

  // Step 2: Compute surface terms contribution to the residual.
  spatial_container->ComputeSurfaceResidual(config_container,
                                            geometry_container,
                                            element_container,
                                            data_container,
                                            work_array,
                                            localTime, iElem);

  // Extract current total residual.
  auto& res = data_container[iElem]->GetDataDOFsRes();

  // Step 3: Multiply residual with the inverse of the mass matrix.
  for(unsigned short iVar=0; iVar<res.size(); iVar++)
    gemv(nDOFsSol2D, nDOFsSol2D, InvMassMatrix, res[iVar], res[iVar]);
}



