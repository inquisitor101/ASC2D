#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"


class CRiemann {

  public:
    // Constructor.
    CRiemann(CConfig *config_container);

    // Destructor.
    virtual ~CRiemann(void);

		// Function that compute that determines the unique state of two variables,
		// purely based on upwinding.
		void ComputeVariableStateUpwinding(const as3vector1d<as3double>  &UnitNormal,
		                                   const as3vector1d<as3double>  &weights,
		                                   const as3vector1d<as3double>  &hElem,
		                                   as3double                     Velocity,
		                                   as3double                   **VarI,
		                                   as3double                   **VarJ,
		                                   as3double                   **VarF);

    // Pure virtual function that determines the unique state of the flux at a face.
    // Note, must be overriden by a derived class.
    virtual void ComputeFluxState(const as3vector1d<as3double>  &UnitNormal,
                                  const as3vector1d<as3double>  &weights,
                                  const as3vector1d<as3double>  &hElem,
                                  as3double                    **VarI,
                                  as3double                    **VarJ,
                                  as3double                    **Flux) = 0;

  protected:
    // Abbreivation: gamma minus one.
    as3double gm1;

  private:

};


class CRoeRiemann : public CRiemann {

  public:
    // Constructor.
    CRoeRiemann(CConfig *config_container);

    // Destructor.
    ~CRoeRiemann(void) final;

  protected:
    // Function that determines the unique state of the flux at a face.
    void ComputeFluxState(const as3vector1d<as3double>  &UnitNormal,
                          const as3vector1d<as3double>  &weights,
                          const as3vector1d<as3double>  &hElem,
                          as3double                    **VarI,
                          as3double                    **VarJ,
                          as3double                    **Flux) final;
  private:
    // Entropy fix.
    as3double Delta;
};


class CRusanovRiemann : public CRiemann {

  public:
    // Constructor.
    CRusanovRiemann(CConfig *config_container);

    // Destructor.
    ~CRusanovRiemann(void) final;

  protected:
    // Function that determines the unique state of the flux at a face.
    void ComputeFluxState(const as3vector1d<as3double>  &UnitNormal,
                          const as3vector1d<as3double>  &weights,
                          const as3vector1d<as3double>  &hElem,
                          as3double                    **VarI,
                          as3double                    **VarJ,
                          as3double                    **Flux) final;
  private:

};


class CRoeIsmailRiemann : public CRiemann {

  public:
    // Constructor.
    CRoeIsmailRiemann(CConfig *config_container);

    // Destructor.
    ~CRoeIsmailRiemann(void) final;

  protected:
    // Function that determines the unique state of the flux at a face.
    void ComputeFluxState(const as3vector1d<as3double>  &UnitNormal,
                          const as3vector1d<as3double>  &weights,
                          const as3vector1d<as3double>  &hElem,
                          as3double                    **VarI,
                          as3double                    **VarJ,
                          as3double                    **Flux) final;
  private:
    // Abbreviations involving gamma.
    as3double ovgm1;
    as3double gp1Ovg;
    as3double gm1Ovg;

    // Values to scale the acoustic eigenvalues to obtain an adequate amount
    // of dissipation to be entropy satisfying in the Ismail_Roe flux.
    as3double beta;
    as3double alphaMax;
};




