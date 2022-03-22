#pragma once

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


class CBoundary {

	public:
		// Constructor.
		CBoundary(CConfig       *config_container,
							CGeometry     *geometry_container,
              CInitial      *initial_container,
							CElement     **element_container,
							unsigned short iZone,
							unsigned short iBoundary);

		// Destructor.
		virtual ~CBoundary(void);

    // Pure virtual function that applies the boundary condition.
    // Must be overriden by a derived class.
    virtual void ImposeBoundaryCondition(CConfig    *config_container,
                                         CGeometry  *geometry_container,
                                         CSolver   **solver_container,
                                         CElement  **element_container,
                                         CSpatial  **spatial_container,
                                         as3double   localTime) = 0;

	protected:
		// Zone ID.
		unsigned short zoneID;
		// Boundary ID.
		unsigned short boundaryID;
    // Type of zone.
    unsigned short typeZone;
    // Number of solution nodes in 1D in this zone.
    unsigned short nDOFsSol1D;
    // Number of integration nodes in 1D in this zone.
		unsigned short nDOFsInt1D;
    // Unit-normal on this boundary.
    as3vector1d<as3double> UnitNormal;

    // Kronecker-delta used in the formulation of the normal/transverse terms.
    as3double KronDelta11; // i.e. ||nx|| = 1,   ny   = 0.
    as3double KronDelta22; // i.e.   nx   = 0, ||ny|| = 1.

    // Element indices that share this boundary.
    as3vector1d<unsigned long> ElemIndexI;

	private:

};


class CEEBoundary : public CBoundary {

	public:
		// Constructor.
		CEEBoundary(CConfig       *config_container,
								CGeometry     *geometry_container,
                CInitial      *initial_container,
								CElement     **element_container,
								unsigned short iZone,
								unsigned short iBoundary);

		// Destructor.
		virtual ~CEEBoundary(void) override;

  protected:
  // Prescribed data on this boundary at integration points.
  // Dimension: [iElem][iVar][iNode].
  as3data2d<as3double> DataDOFsIntBoundary;

  // Prescribed original data on this boundary at integration points.
  // Dimension: [iElem][iVar][iNode].
  as3data2d<as3double> OrigDOFsIntBoundary;

  // Coordinates on this boundary at integration points.
  // Dimension: [iElem][iDim][iNode].
  as3data2d<as3double> GridDOFsIntBoundary;

  // Function that initializes and computes the target-state at the boundary.
  void InitializePrescribedState(CConfig       *config_container,
                                 CGeometry     *geometry_container,
                                 CElement      *element_container,
                                 CInitial      *initial_container,
                                 unsigned short iZone,
                                 unsigned short iBoundary);

  // Function that initializes the modified prescribed boundary condition data.
  void InitializeModifiedBC(CConfig       *config_container,
                            CGeometry     *geometry_container,
                            CElement      *element_container,
                            CInitial      *initial_container,
                            unsigned short iZone,
                            unsigned short iBoundary);

  // Function that modifies the boundary condition.
  void ModifyBoundaryCondition(CConfig    *config_container,
                               CGeometry  *geometry_container,
                               CSolver    *solver_container,
                               CElement   *element_container,
                               CSpatial   *spatial_container,
                               as3double   localTime);

  private:

};


class CEEInterfaceBoundary : public CEEBoundary {

	public:
		// Constructor.
		CEEInterfaceBoundary(CConfig       *config_container,
												 CGeometry     *geometry_container,
                         CInitial      *initial_container,
												 CElement     **element_container,
												 unsigned short iZone,
												 unsigned short iBoundary);

		// Destructor.
		~CEEInterfaceBoundary(void) override;

    // Function that applies an interface boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) override;
	protected:
		// Matching zone ID.
		unsigned short zoneMatchID;
		// Boundary match ID.
		unsigned short boundaryMatchID;
    // Type of zone for matching zone ID.
    unsigned short typeZoneMatch;
    // Number of solution nodes in 1D in the matching zone.
    unsigned short nDOFsSol1DMatch;

    // Element indices that share the matching boundary.
    as3vector1d<unsigned long> ElemIndexJ;

    // Lagrange polynomials in 1D(transposed): interpolates from
    // SolDOFs of the matching zone to IntDOFs of this current zone.
    as3double *lagrangeIntExt1DTranspose = nullptr;

	private:
    // Function that reports information on this boundary.
    void ReportOutput(void);
};


class CEESymmetryBoundary : public CEEBoundary {

	public:
		// Constructor.
		CEESymmetryBoundary(CConfig       *config_container,
												CGeometry     *geometry_container,
                        CInitial      *initial_container,
												CElement     **element_container,
												unsigned short iZone,
												unsigned short iBoundary);

		// Destructor.
		~CEESymmetryBoundary(void) final;

    // Function that applies a symmetry boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:

};


class CEEStaticOutletBoundary : public CEEBoundary {

	public:
		// Constructor.
		CEEStaticOutletBoundary(CConfig       *config_container,
    											  CGeometry     *geometry_container,
                            CInitial      *initial_container,
    											  CElement     **element_container,
    											  unsigned short iZone,
    											  unsigned short iBoundary);

		// Destructor.
		~CEEStaticOutletBoundary(void) final;

    // Function that applies a subsonic static-based inlet boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:
};


class CEESupersonicOutletBoundary : public CEEBoundary {

	public:
		// Constructor.
		CEESupersonicOutletBoundary(CConfig       *config_container,
        											  CGeometry     *geometry_container,
                                CInitial      *initial_container,
        											  CElement     **element_container,
        											  unsigned short iZone,
        											  unsigned short iBoundary);

		// Destructor.
		~CEESupersonicOutletBoundary(void) final;

    // Function that applies a supersonic outlet boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:

};


class CEEStaticInletBoundary : public CEEBoundary {

	public:
		// Constructor.
		CEEStaticInletBoundary(CConfig       *config_container,
    											 CGeometry     *geometry_container,
                           CInitial      *initial_container,
    											 CElement     **element_container,
    											 unsigned short iZone,
    											 unsigned short iBoundary);

		// Destructor.
		~CEEStaticInletBoundary(void) final;

    // Function that applies a subsonic static-based inlet boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:

};


class CEETotalInletBoundary : public CEEBoundary {

	public:
		// Constructor.
		CEETotalInletBoundary(CConfig       *config_container,
  											  CGeometry     *geometry_container,
                          CInitial      *initial_container,
  											  CElement     **element_container,
  											  unsigned short iZone,
  											  unsigned short iBoundary);

		// Destructor.
		~CEETotalInletBoundary(void) final;

    // Function that applies a subsonic total-condition inlet boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:
    // Flow direction.
    as3double udir;
    as3double vdir;

    // Dot product between the normal and velocity direction.
    as3double alpha;
};


class CEESupersonicInletBoundary : public CEEBoundary {

	public:
		// Constructor.
		CEESupersonicInletBoundary(CConfig       *config_container,
        											 CGeometry     *geometry_container,
                               CInitial      *initial_container,
        											 CElement     **element_container,
        											 unsigned short iZone,
        											 unsigned short iBoundary);

		// Destructor.
		~CEESupersonicInletBoundary(void) final;

    // Function that applies a supersonic inlet boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:

};


class CEECharacteristicBoundary : public CEEBoundary {

  public:
    // Constructor.
    CEECharacteristicBoundary(CConfig       *config_container,
                              CGeometry     *geometry_container,
                              CInitial      *initial_container,
                              CElement     **element_container,
                              unsigned short iZone,
                              unsigned short iBoundary);

    // Destructor.
    ~CEECharacteristicBoundary(void) override;

  protected:
    // Number of data values in the working array.
    // Data: sol, dSolDx, dSolDy.
    unsigned short nData = 3;
    // Working data.
    // Dimension: [iData][iVar][iInt1D].
    as3data2d<as3double> WorkingDataInt1D;

    // Least-squares matrix used in the PC-formulation, referred to as C.
    as3double *MatrixC = nullptr;
    // Normal Lagrange gradient 1D coefficient.
    as3double Coef_dell;

    // Index of the incoming acoustic wave.
    unsigned short PsiIndex;
    // Indices used to identify the normal gradient on this boundary.
    unsigned short IndexNormal;
    // Indices used to identify the transverse gradient on this boundary.
    unsigned short IndexTransverse;

    // Safety-factor relaxation coefficient for the transverse terms.
    as3double Coef_eta;

    // Function that computes the average of the Mach number on each element boundary.
    as3double ComputeAverageMachLocal(const as3vector1d<as3double> &weights,
                                      as3double                   **Var);

    // Function that computes the average of the Mach number on entire boundary.
    as3double ComputeAverageMachGlobal(const CGeometryZone               *geometry_zone,
                                       const as3element                  &data_container,
                                       const as3vector1d<as3double>      &weights,
                                       const as3vector1d<unsigned short> &FaceIndexI,
                                       const as3double                   *ellT,
                                       as3double                        **Var);


     // Function that computes the least-squares matrix needed in the NSCBC, if required.
     void ComputeLeastSquaresMatrix(const as3double *lagrangeInt1D,
                                    const as3double *lagrangeInt1DTranspose,
                                    const as3double  dell,
                                    as3double       *MatrixLeastSquares);

  private:

};


class CEEOutletCBC : public CEECharacteristicBoundary {

  public:
    // Constructor.
    CEEOutletCBC(CConfig       *config_container,
                 CGeometry     *geometry_container,
                 CInitial      *initial_container,
                 CElement     **element_container,
                 unsigned short iZone,
                 unsigned short iBoundary);

    // Destructor.
    ~CEEOutletCBC(void) final;

    // Function that applies an outlet CBC boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
  protected:

  private:
    // Condition for adaptive coupled-transverse term relaxation.
    bool AdaptiveBeta_l;
    // Condition for adaptive uncoupled-transverse term relaxation.
    bool AdaptiveBeta_t;
    // Relaxation coefficient for coupled transverse terms.
    as3double beta_l; // eta*beta_l
    // Relaxation coefficient for uncoupled transverse terms.
    as3double beta_t; // eta*beta_t
    // Relaxation coefficient for normal terms.
    as3double coefK; // sigma/(2*len)

};


class CEEInletCBC : public CEECharacteristicBoundary {

  public:
    // Constructor.
    CEEInletCBC(CConfig       *config_container,
                 CGeometry     *geometry_container,
                 CInitial      *initial_container,
                 CElement     **element_container,
                 unsigned short iZone,
                 unsigned short iBoundary);

    // Destructor.
    ~CEEInletCBC(void) final;

    // Function that applies an inlet CBC boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;

  protected:

  private:
    // Condition for adaptive uncoupled-transverse term relaxation.
    bool AdaptiveBeta_t;
    // Relaxation coefficient for uncoupled transverse terms.
    as3double beta_t; // eta*beta_t
    // Relaxation coefficient for normal acoustic terms.
    as3double coefS; // sign(+,-)*sigma/(2*len)
    // Relaxation coefficient for normal vorticity/entropy terms.
    as3double coefE; // sigma/len
};



class CEEPMLInterfaceBoundary : public CEEInterfaceBoundary {

	public:
		// Constructor.
		CEEPMLInterfaceBoundary(CConfig       *config_container,
  												  CGeometry     *geometry_container,
                            CInitial      *initial_container,
  												  CElement     **element_container,
  												  unsigned short iZone,
  												  unsigned short iBoundary);

		// Destructor.
		~CEEPMLInterfaceBoundary(void) final;

    // Function that applies a PML interface boundary condition.
    void ImposeBoundaryCondition(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver   **solver_container,
                                 CElement  **element_container,
                                 CSpatial  **spatial_container,
                                 as3double   localTime) final;
	protected:

	private:
};
