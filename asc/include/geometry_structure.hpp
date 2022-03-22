#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "element_structure.hpp"

// Forward declaration to avoid compiler problems.
class CGeometryZone;
class CGeometryElement;


class CGeometry {

  public:
    // Constructor.
    CGeometry(CConfig   *config_container,
				      CElement **element_container);

    // Destructor.
    ~CGeometry(void);

		// Getter: returns nZone.
		unsigned short GetnZone(void)                                const {return nZone;}
    // Getter: returns nElemZone.
    const as3vector1d<unsigned long> &GetnElemZone(void)         const {return nElemZone;}
    // Getter: returns nPointSubElemP1.
		unsigned long GetnPointSubElemP1(void)                       const {return nPointSubElemP1;}
    // Getter: returns geometry_zone[iZone].
		const CGeometryZone *GetGeometryZone(unsigned short iZone)   const {return geometry_zone[iZone];}
    // GetteR: returns MatchingFace.
    const as3vector1d<unsigned short> &GetMatchingFace(void)     const {return MatchingFace;}
    // Getter: returns ElemBoundaryIndex.
    const as3vector3d<unsigned long> &GetElemBoundaryIndex(void) const {return ElemBoundaryIndex;}
    // Getter: returns ElemBoundaryIndex, per input zone.
    const as3vector2d<unsigned long> &GetElemBoundaryIndex(unsigned short iZone) const {
      return ElemBoundaryIndex[iZone];
    }
    // Getter: returns ElemBoundaryIndex, per zone, per boundary.
    const as3vector1d<unsigned long> &GetElemBoundaryIndex(unsigned short iZone,
                                                           unsigned short iBoundary) const {
      return ElemBoundaryIndex[iZone][iBoundary];
    }
    // Getter: returns UnitNormal, per input face.
    const as3vector1d<as3double> &GetUnitNormal(unsigned short iFace) const {return UnitNormal[iFace];}

	protected:

  private:
		// Number of zones.
		unsigned short nZone;
		// Total number of points in all zones, for nPoly=1 sub-elements.
		unsigned long nPointSubElemP1;

    // Polynomial order in every zone.
    as3vector1d<unsigned short> nPolyZone;
    // Number of elements in x, in every zone.
    as3vector1d<unsigned long> nxElemZone;
    // Number of elements in y, in every zone.
    as3vector1d<unsigned long> nyElemZone;
    // Total number of elements, in every zone.
    as3vector1d<unsigned long> nElemZone;

    // Element indices that share a boundary face, per zone, per boundary.
    // Dimension: [iZone][iBoundary][iElem].
    as3vector3d<unsigned long> ElemBoundaryIndex;

    // Face indices that link/connect to their opposite counterparts.
    as3vector1d<unsigned short> MatchingFace;
    // Unit-normal per face.
    // Dimension: [iFace][iDim].
    as3vector2d<as3double> UnitNormal;

		// Zone geometry.
		as3data1d<CGeometryZone> geometry_zone;

    // Function that computes the total number of grid points in all zones
    // based on the assumptioin of nPoly=1 sub-elements. This is used to
    // compute nPointSubElemP1 for the VTK output.
    void ComputePointsSubElementsP1(void);

    // Function that determines the global indices of each of the elements
    // sharing an external boundary.
    void IdentifyBoundaryElementIndices(void);

    // Function that maps every face to its matching counterpart.
    void IdentifyMatchingFace(void);

    // Function that computes the unit-normal used in every element face.
    void ComputeUnitNormal(void);
};


class CGeometryZone {

	public:
		// Constructor.
		CGeometryZone(CConfig 			*config_container,
									CElement      *element_container,
									unsigned short iZone);

		// Destructor.
		~CGeometryZone(void);

    // Function used to search if a set of coordinates belongs in it.
    bool SearchElementProbe(as3vector1d<as3double> probe,
                            unsigned long         &index,
                            bool                   unique) const;

		// Getter: returns zoneID.
	  unsigned short GetZoneID(void)                     const {return zoneID;}
		// Getter: returns nPolySol.
		unsigned short GetnPolySol(void)                   const {return nPolySol;}
		// Getter: returns nElem.
		unsigned long GetnElem(void)                       const {return nElem;}
    // Getter: returns nsElem.
    const as3vector1d<unsigned long> &GetnsElem(void)  const {return nsElem;}
    // Getter: returns nsElem, per input boundary.
    unsigned long GetnsElem(unsigned short iBoundary)  const {return nsElem[iBoundary];}
    // Getter: returns ZoneSize.
    const as3vector1d<as3double> &GetZoneSize()        const {return ZoneSize;}
    // Getter: returns InternalElemFace.
    const as3vector2d<bool> &GetInternalElemFace(void) const {return InternalElemFace;}
    // Getter: returns IndexNeighborInternalElement.
    const as3vector2d<unsigned long> &GetIndexNeighborInternalElement(void) const {return IndexNeighborInternalElement;}
    // Getter: returns geometry_element[elemID].
		const CGeometryElement *GetGeometryElem(unsigned long elemID)  const {return geometry_element[elemID];}

	protected:
		// Current grid zone ID.
		unsigned short zoneID;
    // Number of elements in x-direction.
    unsigned long  nxElem;
    // Number of elements in y-direction.
    unsigned long  nyElem;
		// Number of elements.
		unsigned long  nElem;
    // Number of elements on each boundary surface (line).
    // Dimension: [iBoundary].
    as3vector1d<unsigned long> nsElem;

    // Type of zone.
    unsigned short TypeZone;
    // Width of the zone.
    // Dimension: [0]: dx, [1]: dy.
    as3vector1d<as3double> ZoneSize;
    // Type of element face (internal or external).
    // Dimension: [iElem][iFace], true if internal.
    as3vector2d<bool> InternalElemFace;
    // Indices of internal elements and their neighbors that share a face.
    // Dimension: [iElem][iFace].
    as3vector2d<unsigned long> IndexNeighborInternalElement;

		// Solution polynomial order of elements.
		unsigned short nPolySol;
		// Number of solution nodes per line element (1D).
		unsigned short nDOFsSol1D;
		// Number of solution nodes per surface element (2D).
		unsigned short nDOFsSol2D;

		// Number of integration nodes per line element (1D).
		unsigned short nDOFsInt1D;
		// Number of integration nodes per surface element (2D).
		unsigned short nDOFsInt2D;

    // Function that generates the entire grid in this zone.
    void GenerateGridZone(CConfig  *config_container,
                          CElement *element_container);

    // Function that flags whether or not an element's face is
    // internal or boundary.
    void IdentifyElementFaceType(void);

    // Function that determines the indices of the internal elements that
    // share the same face.
    void IdentifyNeighborInternalElement(void);

    // Function that assembles the grid in this zone by
    // creating the needed physical elements.
    void AssembleGridZoneElements(CConfig                      *config_container,
                                  CElement                     *element_container,
                                  const as3vector1d<as3double> &BoundingBox,
                                  const as3vector2d<as3double> &ElementSizeZone);

	private:
		// Element geometry.
		as3data1d<CGeometryElement> geometry_element;
};


class CGeometryElement {

	public:
		// Constructor.
		CGeometryElement(CConfig                      *config_container,
                     CElement                     *element_container,
                     const as3vector1d<as3double> &rBasis,
                     const as3vector1d<as3double> &LocalBox);

		// Destructor.
		~CGeometryElement(void);

    // Element ID.
    unsigned long elemID;
    
    // Getter: returns ElemSize.
    const as3vector1d<as3double> &GetElemSize(void)   const {return ElemSize;}
		// Getter: returns CoordSolDOFs.
		const as3data1d<as3double> &GetCoordSolDOFs(void) const {return CoordSolDOFs;}
    // Getter: returns CoordIntDOFs.
    const as3data1d<as3double> &GetCoordIntDOFs(void) const {return CoordIntDOFs;}
    // Getter: returns CoordBoundary, per given face.
    as3double GetCoordBoundary(unsigned short iFace)  const {return CoordBoundary[iFace];}

	protected:

	private:
    // Element size.
    as3vector1d<as3double> ElemSize;

    // Nodal points at the physical element on its boundaries.
    // Dimension: [iFace].
    as3vector1d<as3double> CoordBoundary;

		// Nodal points at solution DOFs.
		// Dimension: [iDim][iNode].
		as3data1d<as3double> CoordSolDOFs;

    // Nodal-points at integration DOFs.
    // Dimension: [iDim][iNode].
    as3data1d<as3double> CoordIntDOFs;
};


