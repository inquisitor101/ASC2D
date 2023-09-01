#pragma once

/*!
 * @file boundary_structure.hpp
 * @brief The file containing all the boundary condition information.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "spatial_structure.hpp"
#include "solver_structure.hpp"
#include "initial_structure.hpp"

// Forward declaration to avoid compilor problems.
class CSolver;
class CInitial;

/*!
 * @brief A class used as an interface for different boundary conditions.
 */
class CBoundary {

	public:
		/*!
		 * @brief Default constructor of CBoundary, which initializes a generic boundary interface.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */
		CBoundary(CConfig       *config_container,
							CGeometry     *geometry_container,
              CInitial      *initial_container,
							CElement     **element_container,
							unsigned short iZone,
							unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~CBoundary(void);

		/*!
		 * @brief Getter function which returns the value of ElemIndexI.
		 *
		 * @return ElemIndexI
		 */
		const as3vector1d<unsigned long> &GetElemIndexI(void) const {return ElemIndexI;}

    /*!
		 * @brief Pure virtual function that applies the boundary condition. 
		 * Must be overriden by a derived class.
		 *
		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
    virtual void ImposeBoundaryCondition(CConfig    *config_container,
                                         CGeometry  *geometry_container,
                                         CSolver   **solver_container,
                                         CElement  **element_container,
                                         CSpatial  **spatial_container,
                                         as3double   localTime) = 0;

	protected:
		unsigned short             zoneID;      ///< Zone ID of this boundary.
		unsigned short             boundaryID;  ///< Current boundary ID.
    unsigned short             typeZone;    ///< Type of zone containing this boundary.
    unsigned short             nDOFsSol1D;  ///< Number of solution DOFs in 1D on this boundary.
		unsigned short             nDOFsInt1D;  ///< Number of integration points in 1D on this boundary.
    as3vector1d<as3double>     UnitNormal;  ///< Unit-normal of this boundary.

    as3double                  KronDelta11; ///< Kronecker-delta used in the formulation of the NSCBC normal/transverse terms.
																						///< This coefficient (11) implies: ||nx|| = 1,   ny   = 0.
    as3double                  KronDelta22; ///< Kronecker-delta used in the formulation of the NSCBC normal/transverse terms.
																						///< This coefficient (22) implies:   nx   = 0, ||ny|| = 1.

    as3vector1d<unsigned long> ElemIndexI;  ///< Vector containing all the element indices that share this boundary.

	private:

};


/*!
 * @brief A generic class used for implementing different Euler-equation (EE) boundary conditions.
 */
class CEEBoundary : public CBoundary {

	public:
		/*!
		 * @brief Default constructor of CEEBoundary, which initializes a generic Euler-equation (EE) boundary.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */	
		CEEBoundary(CConfig       *config_container,
								CGeometry     *geometry_container,
                CInitial      *initial_container,
								CElement     **element_container,
								unsigned short iZone,
								unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~CEEBoundary(void) override;

  protected:
  as3data2d<as3double> DataDOFsIntBoundary; ///< Prescribed data on this boundary, at the integration points (1D).
																						///< Dimension: [iElem][iVar][iNode].
  as3data2d<as3double> OrigDOFsIntBoundary; ///< Prescribed original data on this boundary, at integration points (1D).
																						///< Note, this is only used in case a modified boundary condition is specified.
																						///< Dimension: [iElem][iVar][iNode].
  as3data2d<as3double> GridDOFsIntBoundary; ///< Coordinates of integration points (1D) on this boundary.
																						///< Dimension: [iElem][iDim][iNode].

  /*!
	 * @brief Function that initializes and computes the target-state at the boundary.
	 *
	 * @param[in] config_container pointer to input configuration/dictionary file.
	 * @param[in] geometry_container pointer to the geometry container.
	 * @param[in] element_container pointer to current zone standard element container.
	 * @param[in] initial_container pointer to the current zone initial condition container.
	 * @param[in] iZone current zone ID.
	 * @param[in] iBoundary current boundary ID.
	 */
  void InitializePrescribedState(CConfig       *config_container,
                                 CGeometry     *geometry_container,
                                 CElement      *element_container,
                                 CInitial      *initial_container,
                                 unsigned short iZone,
                                 unsigned short iBoundary);

  /*!
	 * @brief Function that initializes the modified prescribed boundary condition data.
	 *
	 * @param[in] config_container pointer to input configuration/dictionary file.
	 * @param[in] geometry_container pointer to the geometry container.
	 * @param[in] element_container pointer to current zone standard element container.
	 * @param[in] initial_container pointer to the current zone initial condition container.
	 * @param[in] iZone current zone ID.
	 * @param[in] iBoundary current boundary ID.
	 */ 
	void InitializeModifiedBC(CConfig       *config_container,
                            CGeometry     *geometry_container,
                            CElement      *element_container,
                            CInitial      *initial_container,
                            unsigned short iZone,
                            unsigned short iBoundary);

  /*!
	 * @brief Function that modifies the boundary condition.
	 *
	 * @param[in] config_container pointer to input configuration/dictionary file.
	 * @param[in] geometry_container pointer to the geometry container.
	 * @param[in] solver_container pointer to the current zone solver container.
	 * @param[in] element_container pointer to the current zone standard element container.
	 * @param[in] spatial_container pointer to the currrent zone spatial container.
	 * @param[in] localTime current physical time.
	 */ 
	void ModifyBoundaryCondition(CConfig    *config_container,
                               CGeometry  *geometry_container,
                               CSolver    *solver_container,
                               CElement   *element_container,
                               CSpatial   *spatial_container,
                               as3double   localTime);

  private:

};


/*!
 * @brief A class used for implementing interface/periodic boundary conditions.
 */
class CEEInterfaceBoundary : public CEEBoundary {

	public:
		/*!
		 * @brief Default constructor of CEEInterfaceBoundary, which initializes an interface BC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */	
		CEEInterfaceBoundary(CConfig       *config_container,
												 CGeometry     *geometry_container,
                         CInitial      *initial_container,
												 CElement     **element_container,
												 unsigned short iZone,
												 unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CEEInterfaceBoundary(void) final;

    /*!
		 * @brief Function that applies an interface boundary condition.
		 *
		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:
		unsigned short             zoneMatchID;                         ///< Matching zone ID.
		unsigned short             boundaryMatchID;                     ///< Matching boundary ID.
    unsigned short             typeZoneMatch;                       ///< Type of zone that is coupled to this one through this interface boundary.
    unsigned short             nDOFsSol1DMatch;                     ///< Number of solution DOFs in 1D in the matching zone.
    as3vector1d<unsigned long> ElemIndexJ;                          ///< Element indices that share the matching boundary in the matching zone.
    as3double                 *lagrangeIntExt1DTranspose = nullptr; ///< Lagrange polynomial in 1D (transposed), whih interpolates from solution to integration points in 1D.
																																		///< Note, the solution DOFs belong to this zone and the integration points are that of the matching zone.

	private:
    /*!
		 * @brief Function that reports information on this boundary.
		 */
    void ReportOutput(void);
};


/*!
 * @brief A class used for implementing symmetry boundary conditions.
 */
class CEESymmetryBoundary : public CEEBoundary {

	public:
		/*!
		 * @brief Default constructor of CEESymmetryBoundary, which initializes a symmetry BC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */
		CEESymmetryBoundary(CConfig       *config_container,
												CGeometry     *geometry_container,
                        CInitial      *initial_container,
												CElement     **element_container,
												unsigned short iZone,
												unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CEESymmetryBoundary(void) final;

    /*!
		 * @brief Function that applies a symmetry boundary condition.
		 *
		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */   
		void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:

};


/*!
 * @brief A class used for implementing Riemann-extrapolation static outlet boundary conditions.
 */
class CEEStaticOutletBoundary : public CEEBoundary {

	public:
		/*!
		 * @brief Default constructor of CEEStaticOutletBoundary, which initializes a Riemann-extrapolation static outlet BC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */
		CEEStaticOutletBoundary(CConfig       *config_container,
    											  CGeometry     *geometry_container,
                            CInitial      *initial_container,
    											  CElement     **element_container,
    											  unsigned short iZone,
    											  unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CEEStaticOutletBoundary(void) final;

    /*!
		 * @brief Function that applies a subsonic Riemann-extrapolation static outlet boundary condition.
		 *
		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */    
		void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:
};


/*!
 * @brief A class used for implementing supersonic outlet boundary conditions.
 */
class CEESupersonicOutletBoundary : public CEEBoundary {

	public:
		/*!
		 * @brief Default constructor of CEESupersonicOutletBoundary, which initializes a supersonic outlet BC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */
		CEESupersonicOutletBoundary(CConfig       *config_container,
        											  CGeometry     *geometry_container,
                                CInitial      *initial_container,
        											  CElement     **element_container,
        											  unsigned short iZone,
        											  unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */	
		~CEESupersonicOutletBoundary(void) final;

    /*!
		 * @brief Function that applies a supersonic outlet boundary condition.
		 *
		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:

};


/*!
 * @brief A class used for implementing Riemann-extrapolation static inlet boundary conditions.
 */
class CEEStaticInletBoundary : public CEEBoundary {

	public:
		/*!
		 * @brief Default constructor of CEEStaticInletBoundary, which initializes a Riemann-extrapolation static inlet BC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */
		CEEStaticInletBoundary(CConfig       *config_container,
    											 CGeometry     *geometry_container,
                           CInitial      *initial_container,
    											 CElement     **element_container,
    											 unsigned short iZone,
    											 unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */		
		~CEEStaticInletBoundary(void) final;

    /*!
		 * @brief Function that applies a subsonic Riemann-extrapolation static inlet boundary condition.
		 *
		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:

};


/*!
 * @brief A class used for implementing Riemann-extrapolation total/stagnation inlet boundary conditions.
 */
class CEETotalInletBoundary : public CEEBoundary {

	public:
		/*!
		 * @brief Default constructor of CEETotalInletBoundary, which initializes a Riemann-extrapolation total/stagnation inlet BC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */
		CEETotalInletBoundary(CConfig       *config_container,
  											  CGeometry     *geometry_container,
                          CInitial      *initial_container,
  											  CElement     **element_container,
  											  unsigned short iZone,
  											  unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */		
		~CEETotalInletBoundary(void) final;

    /*!
		 * @brief Function that applies a subsonic Riemann-extrapolation total/stagnation inlet boundary condition.
		 *
		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:
    as3double udir;  ///< Flow direction, x-component.
    as3double vdir;  ///< Flow direction, y-component.
    as3double alpha; ///< Dot product between the normal and velocity direction.
};


/*!
 * @brief A class used for implementing supersonic inlet boundary conditions.
 */
class CEESupersonicInletBoundary : public CEEBoundary {

	public:
		/*!
		 * @brief Default constructor of CEESupersonicInletBoundary, which initializes a supersonic inlet BC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */	
		CEESupersonicInletBoundary(CConfig       *config_container,
        											 CGeometry     *geometry_container,
                               CInitial      *initial_container,
        											 CElement     **element_container,
        											 unsigned short iZone,
        											 unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */		
		~CEESupersonicInletBoundary(void) final;

    /*!
		 * @brief Function that applies a supersonic inlet boundary condition.
		 *
 		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
		void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:

};


/*!
 * @brief A generic class used for implementing characteristic boundary conditions.
 */
class CEECharacteristicBoundary : public CEEBoundary {

  public:
 		/*!
		 * @brief Default constructor of CEECharacteristicBoundary, which initializes a generic NSCBC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */	
		CEECharacteristicBoundary(CConfig       *config_container,
                              CGeometry     *geometry_container,
                              CInitial      *initial_container,
                              CElement     **element_container,
                              unsigned short iZone,
                              unsigned short iBoundary);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */		
    ~CEECharacteristicBoundary(void) override;

  protected:
    unsigned short       nData = 3;         ///< Number of data values in the working array. Data: sol, dSolDx, dSolDy.
    as3data2d<as3double> WorkingDataInt1D;  ///< Working data variables at integration points (1D).
																					  ///< Dimension: [iData][iVar][iInt1D].

    as3double           *MatrixC = nullptr; ///< Least-squares matrix used in the PC-formulation, referred to as C.
    as3double            Coef_dell;         ///< Coefficient of the normal Lagrange gradient in 1D.

    unsigned short       PsiIndex;          ///< Index of the incoming acoustic wave.
    unsigned short       IndexNormal;       ///< Indices used to identify the normal gradient on this boundary.
    unsigned short       IndexTransverse;   ///< Indices used to identify the transverse gradient on this boundary.

    as3double            Coef_eta;          ///< Safety-factor relaxation coefficient for the transverse terms.

    /*!
		 * @brief Function that computes the average of the Mach number on each element boundary.
		 *
		 * @param[in] weights reference to the integration weights in 1D.
		 * @param[in] Var pointer to the working variables on this boundary (1D).
		 *
		 * @return average of the Mach number on this element boundary surface.
		 */
    as3double ComputeAverageMachLocal(const as3vector1d<as3double> &weights,
                                      as3double                   **Var);

    /*!
		 * @brief Function that computes the average of the Mach number on entire boundary.
		 *
		 * @param[in] geometry_zone pointer to the grid geometry in the current zone.
		 * @param[in] data_container reference to the data container in this zone.
		 * @param[in] weights reference to the integration weights in 1D.
		 * @param[in] FaceIndexI reference to the indices of the solution DOFs on this boundary in 1D.
		 * @param[in] ellT pointer to the transpose of the Lagrange interpolating polynomial in 1D.
		 * @param[in] Var pointer to the working variables.
		 *
		 * @return average of the Mach number over this entire boundary.
		 */
    as3double ComputeAverageMachGlobal(const CGeometryZone               *geometry_zone,
                                       const as3element                  &data_container,
                                       const as3vector1d<as3double>      &weights,
                                       const as3vector1d<unsigned short> &FaceIndexI,
                                       const as3double                   *ellT,
                                       as3double                        **Var);


     /*!
			* @brief Function that computes the least-squares matrix needed in the NSCBC, if required.
			*
			* @param[in] lagrangeInt1D pointer to the Lagrange interpolating polynomial in 1D.
			* @param[in] lagrangeInt1DTranspose pointer to the transpose of the Lagrange interpolating polynomial in 1D.
			* @param[in] dell coefficient of the derivative of the Lagrange polynomial acting on a boundary.
			* @param[out] MatrixLeastSquares pointer to the least-squares matrix (C) that is computed. 
			*/
     void ComputeLeastSquaresMatrix(const as3double *lagrangeInt1D,
                                    const as3double *lagrangeInt1DTranspose,
                                    const as3double  dell,
                                    as3double       *MatrixLeastSquares);

  private:

};


/*!
 * @brief A class used for implementing characteristic static outlet boundary conditions.
 */
class CEEOutletCBC : public CEECharacteristicBoundary {

  public:
 		/*!
		 * @brief Default constructor of CEEOutletCBC, which initializes a static outlet NSCBC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */   
		CEEOutletCBC(CConfig       *config_container,
                 CGeometry     *geometry_container,
                 CInitial      *initial_container,
                 CElement     **element_container,
                 unsigned short iZone,
                 unsigned short iBoundary);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */		   
		~CEEOutletCBC(void) final;

    /*!
		 * @brief Function that applies an outlet CBC boundary condition.
		 *
 		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
  protected:

  private:
    bool      AdaptiveBeta_l; ///< Option for adaptive coupled transverse term relaxation.
    bool      AdaptiveBeta_t; ///< Option for adaptive uncoupled transverse term relaxation.
    as3double beta_l;         ///< Relaxation coefficient for coupled transverse terms.
    as3double beta_t;         ///< Relaxation coefficient for uncoupled transverse terms.
    as3double coefK;          ///< Relaxation coefficient for normal terms, such that: coefK = sigma/(2*len).
};


/*!
 * @brief An interface class used for implementing characteristic inlet boundary conditions.
 */
class CEEInletCBC : public CEECharacteristicBoundary {

  public:
 		/*!
		 * @brief Default constructor of CEEInletCBC, which initializes a generic inlet NSCBC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */   
		CEEInletCBC(CConfig       *config_container,
                 CGeometry     *geometry_container,
                 CInitial      *initial_container,
                 CElement     **element_container,
                 unsigned short iZone,
                 unsigned short iBoundary);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
    virtual ~CEEInletCBC(void) override;
};


/*!
 * @brief A class used for implementing characteristic static inlet boundary conditions.
 */
class CEEStaticInletCBC : public CEEInletCBC {

  public:
 		/*!
		 * @brief Default constructor of CEEStaticInletCBC, which initializes a static inlet NSCBC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */   
    CEEStaticInletCBC(CConfig       *config_container,
                      CGeometry     *geometry_container,
                      CInitial      *initial_container,
                      CElement     **element_container,
                      unsigned short iZone,
                      unsigned short iBoundary);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
    ~CEEStaticInletCBC(void) final;

    /*!
		 * @brief Function that applies a static inlet characteristic boundary condition.
		 *
 		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;

  protected:

  private:
    as3double coefS;  ///< Relaxation coefficient for normal acoustic terms.
    as3double coefE;  ///< Relaxation coefficient for normal vorticity/entropy terms.
};


/*!
 * @brief A class used for implementing characteristic total/stagnation inlet boundary conditions.
 */
class CEETotalInletCBC : public CEEInletCBC {

  public:
 		/*!
		 * @brief Default constructor of CEETotalInletCBC, which initializes a total/stagnation inlet NSCBC.
		 *
		 * @param[in] config_container pointer to input configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the current zone initial condition container.
		 * @param[in] element_container pointer to input standard element container.
		 * @param[in] iZone current zone ID.
		 * @param[in] iBoundary current boundary ID.
		 */ 
    CEETotalInletCBC(CConfig       *config_container,
                     CGeometry     *geometry_container,
                     CInitial      *initial_container,
                     CElement     **element_container,
                     unsigned short iZone,
                     unsigned short iBoundary);

 		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
    ~CEETotalInletCBC(void) final;

    /*!
		 * @brief Function that applies a total/stagnation inlet characteristic boundary condition.
		 *
 		 * @param[in] config_container pointer to the configuration/dictionary file.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container pointer to the solver container.
		 * @param[in] element_container pointer to the standard element container.
		 * @param[in] spatial_container pointer to the spatial container.
		 * @param[in] localTime current physical time.
		 */
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;

  protected:

  private:
    as3double      coefK;     ///< Relaxation coefficient for normal acoustic terms.
		unsigned short PhiIndex;  ///< Outgoing acoustic index.
};



