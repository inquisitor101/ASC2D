#pragma once

/*!
 * @file geometry_structure.hpp
 * @brief The file containing all the grid information.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "element_structure.hpp"

// Forward declaration to avoid compiler problems.
class CGeometryZone;
class CGeometryElement;


/*!
 * @brief A class used for generating a (multi-zone) grid.
 */
class CGeometry {

  public:
		/*!
		 * @brief Default constructor of CGeometry, which initializes the grid.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] element_container pointer to input standard element container.
		 */
    CGeometry(CConfig   *config_container,
				      CElement **element_container);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
    ~CGeometry(void);

		/*!
		 * @brief Getter function which returns the value of nZone.
		 *
		 * @return nZone
		 */
		unsigned short GetnZone(void)                                const {return nZone;}
    /*!
		 * @brief Getter function which returns the value of nElemZone.
		 *
		 * @return nElemZone
		 */
    const as3vector1d<unsigned long> &GetnElemZone(void)         const {return nElemZone;}
    /*!
		 * @brief Getter function which returns the value of nPointSubElemP1.
		 *
		 * @return nPointSubElemP1
		 */
		unsigned long GetnPointSubElemP1(void)                       const {return nPointSubElemP1;}
    /*!
		 * @brief Getter function which returns the value of geometry_zone[iZone].
		 *
		 * @param[in] iZone input zone ID.
		 *
		 * @return geometry_zone[iZone]
		 */
		const CGeometryZone *GetGeometryZone(unsigned short iZone)   const {return geometry_zone[iZone];}
    /*!
		 * @brief Getter function which returns the value of MatchingFace.
		 *
		 * @return MatchingFace
		 */
    const as3vector1d<unsigned short> &GetMatchingFace(void)     const {return MatchingFace;}
    /*!
		 * @brief Getter function which returns the value of ElemBoundaryIndex.
		 *
		 * @return ElemBoundaryIndex
		 */
    const as3vector3d<unsigned long> &GetElemBoundaryIndex(void) const {return ElemBoundaryIndex;}
    /*!
		 * @brief Getter function which returns the value of ElemBoundaryIndex[iZone].
		 *
		 * @param[in] iZone input zone ID.
		 *
		 * @return ElemBoundaryIndex[iZone]
		 */
    const as3vector2d<unsigned long> &GetElemBoundaryIndex(unsigned short iZone) const {
      return ElemBoundaryIndex[iZone];
    }
    /*!
		 * @brief Getter function which returns the value of ElemBoundaryIndex[iZone][iBoundary].
		 *
		 * @param[in] iZone input zone ID.
		 * @param[in] iBoundary input boundary ID.
		 *
		 * @return ElemBoundaryIndex[iZone][iBoundary]
		 */
    const as3vector1d<unsigned long> &GetElemBoundaryIndex(unsigned short iZone,
                                                           unsigned short iBoundary) const {
      return ElemBoundaryIndex[iZone][iBoundary];
    }
    /*!
		 * @brief Getter function which returns the value of UnitNormal[iFace].
		 *
		 * @param[in] iFace input face ID.
		 *
		 * @return UnitNormal[iFace]
		 */
    const as3vector1d<as3double> &GetUnitNormal(unsigned short iFace) const {return UnitNormal[iFace];}

	protected:

  private:
		unsigned short              nZone;             ///< Total number of zones.
		unsigned long               nPointSubElemP1;   ///< Total number of points in all zones, based on nPoly=1 sub-elements.
    as3vector1d<unsigned short> nPolyZone;         ///< Solution polynomial order in every zone.
    as3vector1d<unsigned long>  nxElemZone;        ///< Number of elements in x-direction in every zone.
    as3vector1d<unsigned long>  nyElemZone;        ///< Number of elements in y-direction in every zone.
    as3vector1d<unsigned long>  nElemZone;         ///< Total number of elemetns in every zone.
    as3vector3d<unsigned long>  ElemBoundaryIndex; ///< Element indices that share a boundary face, per zone, per boundary.
																									 ///< Dimension: [iZone][iBoundary][iElem].
    as3vector1d<unsigned short> MatchingFace;      ///< Face indices that are connected with their opposite counterparts.
    as3vector2d<as3double>      UnitNormal;        ///< Unit-normal per each face. Dimension: [iFace][iDim].
		as3data1d<CGeometryZone>    geometry_zone;     ///< Container that defines the grid in every zone.

    /*!
		 * @brief Function that computes the total number of grid points in all zones.
		 * Note, this is based on the assumptioin of nPoly=1 sub-elements. This is also 
		 * used to compute nPointSubElemP1 for the VTK output.
		 */
    void ComputePointsSubElementsP1(void);

    /*!
		 * @brief Function that determines the global indices of each of the elements sharing an external boundary.
		 */
    void IdentifyBoundaryElementIndices(void);

    /*!
		 * @brief Function that maps every face to its matching counterpart.
		 */
    void IdentifyMatchingFace(void);

    /*!
		 * @brief Function that computes the unit-normal used in every element face.
		 */
    void ComputeUnitNormal(void);
};


/*!
 * @brief A class used for generating a grid in a single zone.
 */
class CGeometryZone {

	public:
		/*!
		 * @brief Default constructor of CGeometryZone, which initializes the grid in a single zone.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] element_container pointer to input standard element container on this zone.
		 * @param[in] iZone current zone ID.
		 */
		CGeometryZone(CConfig 			*config_container,
									CElement      *element_container,
									unsigned short iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CGeometryZone(void);

    /*!
		 * @brief Function used to search if a set of coordinates belongs in it.
		 *
		 * @param[in] probe coordinates of the input probe location.
		 * @param[out] index reference to the index of the element containing the probe.
		 * @param[in] unique option whether or not the probe location is unique to a single element.
		 *
		 * @return whether or not the element owner of the probe has been found.
		 */
    bool SearchElementProbe(as3vector1d<as3double> probe,
                            unsigned long         &index,
                            bool                   unique) const;

		/*!
		 * @brief Getter function which returns the value of LayerOrientation.
		 *
		 * @return LayerOrientation
		 */
		as3vector1d<as3double> GetLayerOrientation(void)                        const {return LayerOrientation;}
		/*!
		 * @brief Getter function which returns the value of zoneID.
		 *
		 * @return zoneID
		 */
	  unsigned short GetZoneID(void)                                          const {return zoneID;}
		/*!
		 * @brief Getter function which returns the value of nPolySol.
		 *
		 * @return nPolySol
		 */
		unsigned short GetnPolySol(void)                                        const {return nPolySol;}
		/*!
		 * @brief Getter function which returns the value of nElem.
		 *
		 * @return nElem
		 */
		unsigned long GetnElem(void)                                            const {return nElem;}
		/*!
		 * @brief Getter function which returns the value of nDOFsSol2D.
		 *
		 * @return nDOFsSol2D
		 */
		unsigned short GetnDOFsSol2D(void)                                      const {return nDOFsSol2D;}
    /*!
		 * @brief Getter function which returns the value of nsElem.
		 *
		 * @return nsElem
		 */
    const as3vector1d<unsigned long> &GetnsElem(void)                       const {return nsElem;}
    /*!
		 * @brief Getter function which returns the value of nsElem[iBoundary].
		 *
		 * @param[in] iBoundary input boundary ID.
		 *
		 * @return nsElem[iBoundary]
		 */
    unsigned long GetnsElem(unsigned short iBoundary)                       const {return nsElem[iBoundary];}
    /*!
		 * @brief Getter function which returns the value of ZoneSize.
		 *
		 * @return ZoneSize
		 */
    const as3vector1d<as3double> &GetZoneSize(void)                         const {return ZoneSize;}
    /*!
		 * @brief Getter function which returns the value of InternalElemFace.
		 *
		 * @return InternalElemFace
		 */
    const as3vector2d<bool> &GetInternalElemFace(void)                      const {return InternalElemFace;}
    /*!
		 * @brief Getter function which returns the value of IndexNeighborInternalElement.
		 *
		 * @return IndexNeighborInternalElement
		 */
    const as3vector2d<unsigned long> &GetIndexNeighborInternalElement(void) const {return IndexNeighborInternalElement;}
    /*!
		 * @brief Getter function which returns the value of geometry_element[elemID].
		 *
		 * @param[in] elemID input element ID.
		 *
		 * @return geometry_element[elemID]
		 */
		const CGeometryElement *GetGeometryElem(unsigned long elemID)           const {return geometry_element[elemID];}

	protected:
		unsigned short             zoneID;                       ///< Current grid zone ID.
    unsigned long              nxElem;                       ///< Number of elements in x-direction in this zone.
    unsigned long              nyElem;                       ///< Number of elements in y-direction in this zone.
		unsigned long              nElem;                        ///< Total number of elements in this zone.
    as3vector1d<unsigned long> nsElem;                       ///< Number of elements on each boundary surface in this zone.
																														 ///< Dimension: [iBoundary].
		as3vector1d<as3double>     LayerOrientation;             ///< Current layer orientation w.r.t. main zone. 
																														 ///< Note, this is defined in non-main zone layers only.
																														 ///< Dimension: [iDim].

    unsigned short             TypeZone;                     ///< Type of current zone.
    as3vector1d<as3double>     ZoneSize;                     ///< Extrusion width of the current zone.
																														 ///< Dimension: [0]: dx, [1]: dy.
    as3vector2d<bool>          InternalElemFace;             ///< Type of element face, options: internal or external.
																														 ///< Dimension: [iElem][iFace], true if internal.
    as3vector2d<unsigned long> IndexNeighborInternalElement; ///< Indices of internal elements and their neighbors that share a face.
																														 ///< Dimension: [iElem][iFace].

		unsigned short             nPolySol;                     ///< Solution polynomial order of elements in current zone.
		unsigned short             nDOFsSol1D;                   ///< Number of solution DOFs in 1D in each element in current zone.
		unsigned short             nDOFsSol2D;                   ///< Number of solution DOFs in 2D in each element in current zone.
		unsigned short             nDOFsInt1D;                   ///< Number of integration points in 1D in each element in current zone.
		unsigned short             nDOFsInt2D;                   ///< Number of integration points in 2D in each element in current zone.

    /*!
		 * @brief Function that generates the entire grid in this zone.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] element_container pointer to the standard element container of this zone.
		 */
    void GenerateGridZone(CConfig  *config_container,
                          CElement *element_container);

    /*!
		 * @brief Function that flags whether or not an element's face is internal or boundary.
		 */
    void IdentifyElementFaceType(void);

    /*!
		 * @brief Function that determines the indices of the internal elements that share the same face.
		 */
    void IdentifyNeighborInternalElement(void);

    /*!
		 * @brief Function that assembles the grid in this zone by creating the needed physical elements.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] element_container pointer to the standard element container of this zone.
		 * @param[in] BoundingBox reference to dimensions of the extruded zone.
		 * @param[in] ElementSizeZone reference to element sizes in the extruded direction. 
		 */
    void AssembleGridZoneElements(CConfig                      *config_container,
                                  CElement                     *element_container,
                                  const as3vector1d<as3double> &BoundingBox,
                                  const as3vector2d<as3double> &ElementSizeZone);

	private:
		as3data1d<CGeometryElement> geometry_element;  ///< Grid geometry of each element in this zone.
};


/*!
 * @brief A class used for generating a grid in a element belonging to a specific zone.
 */
class CGeometryElement {

	public:
		/*!
		 * @brief Default constructor of CGeometryElement, which initializes the grid of a single element.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] element_container pointer to input standard element container on this zone.
		 * @param[in] rBasis reference to the solution DOFs in 1D in this zone/element.
		 * @param[in] LocalBox reference to dimension of the current element.
		 */
		CGeometryElement(CConfig                      *config_container,
                     CElement                     *element_container,
                     const as3vector1d<as3double> &rBasis,
                     const as3vector1d<as3double> &LocalBox);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CGeometryElement(void);

    unsigned long elemID;   ///< Current element ID.
    
    /*!
		 * @brief Getter function whih returns the value of ElemSize.
		 *
		 * @return ElemSize
		 */
    const as3vector1d<as3double> &GetElemSize(void)   const {return ElemSize;}
		/*!
		 * @brief Getter function which returns the value of CoordSolDOFs.
		 *
		 * @return CoordSolDOFs
		 */
		const as3data1d<as3double> &GetCoordSolDOFs(void) const {return CoordSolDOFs;}
    /*!
		 * @brief Getter function which returns the value of CoordIntDOFs.
		 *
		 * @return CoordIntDOFs
		 */
    const as3data1d<as3double> &GetCoordIntDOFs(void) const {return CoordIntDOFs;}
    /*!
		 * @brief Getter function which returns the value of CoordBoundary[iFace].
		 *
		 * @param[in] iFace input face ID.
		 *
		 * @return CoordBoundary[iFace]
		 */
    as3double GetCoordBoundary(unsigned short iFace)  const {return CoordBoundary[iFace];}

	protected:

	private:
    as3vector1d<as3double> ElemSize;      ///< Current element size in x- and y-dimension.
    as3vector1d<as3double> CoordBoundary; ///< Coordinates of solution DOFs at the element boundaries/faces in 1D.
																					///< Dimension: [iFace].
		as3data1d<as3double>   CoordSolDOFs;  ///< Coordinates of solution DOFs on the element in 2D.
																					///< Dimension: [iDim][iNode].
    as3data1d<as3double>   CoordIntDOFs;  ///< Coordinates of integration points on the element in 2D.
																					///< Dimension: [iDim][iNode].
};


