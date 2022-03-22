#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "initial_structure.hpp"



class CData {

	public:
		// Constructor.
		CData(CConfig   	  *config_container,
				  CGeometry 	  *geometry_container,
          CInitial      *initial_container,
					CElement  	  *element_container,
					unsigned short iZone,
					unsigned long  iElem);

		// Destructor.
		virtual ~CData(void);

    // Function that resets the residual.
    void ResetResidual(void){
      for(unsigned short i=0; i<DataDOFsRes.size(); i++)
        memset(DataDOFsRes[i], 0.0, nDOFsSol2D*sizeof(as3double));
    }

		// Getter: returns DataDOFsSol.
    as3data1d<as3double> &GetDataDOFsSol(void)                             {return DataDOFsSol;}
    // Getter: returns DataDOFsResSol.
    as3data1d<as3double> &GetDataDOFsRes(void)                             {return DataDOFsRes;}
    // Getter: returns DampDOFsInt.
    as3data1d<as3double> &GetDampDOFsInt(void)                             {return DampDOFsInt;}
    // Getter: returns DataDOFsIntMean.
    as3data1d<as3double> &GetDataDOFsIntMean(void)                         {return DataDOFsIntMean;}
    // Getter: returns derLagrangeGridStretching1D.
    as3data1d<as3double> &GetDerLagrangeGridStretching1D(void)             {return derLagrangeGridStretching1D;}
    // Getter: returns derLagrangeGridStretchingFace, per input face.
    as3double GetDerLagrangeGridStretchingFace(unsigned short iFace)       {return derLagrangeGridStretchingFace[iFace];}
    // Getter: returns derLagrangeArtificialConvection1D.
    as3data1d<as3double> &GetDerLagrangeArtificialConvection1D(void)       {return derLagrangeArtificialConvection1D;}
    // Getter: returns derLagrangeArtificialConvectionFace.
    as3vector1d<as3double> &GetDerLagrangeArtificialConvectionFace(void)   {return derLagrangeArtificialConvectionFace;}
    // Getter: returns derLagrangeArtificialConvectionFace, per input face.
    as3double GetDerLagrangeArtificialConvectionFace(unsigned short iFace) {return derLagrangeArtificialConvectionFace[iFace];}
    // Getter: returns DataDOFsSol.size().
    unsigned short GetnVarTotal(void)                                      {return DataDOFsSol.size();}
    // Getter: returns DataDOFsIntFace, per input face.
    as3data1d<as3double> &GetDataDOFsIntFace(unsigned short iFace)         {return DataDOFsIntFace[iFace];}
    // Getter: returns DataDOFsIntMeanFace, per input face.
    as3data1d<as3double> &GetDataDOFsIntMeanFace(unsigned short iFace)     {return DataDOFsIntMeanFace[iFace];}
    // Getter: returns DataDOFsIntAuxFace, per input face.
    as3data1d<as3double> &GetDataDOFsIntAuxFace(unsigned short iFace)      {return DataDOFsIntAuxFace[iFace];}

	protected:
    // Number of solution DOFs in 2D.
    unsigned short nDOFsSol2D;

		// Data DOFs at solution points on volume (surface) in 2D.
		// Dimension: [iVar][iNode2D].
		as3data1d<as3double> DataDOFsSol;

    // Data DOFs at integration points on face (line) in 1D.
    // Dimension: [iFace][iVar][iNode1D].
    as3data2d<as3double> DataDOFsIntFace;

    // Data residual DOFs at solution points on volume (surface) in 2D.
    // Dimension: [iVar][iNode2D].
    as3data1d<as3double> DataDOFsRes;

    // Damping function coefficients on the integration points in 2D.
    // Dimension: [iDim][iNode2D].
    as3data1d<as3double> DampDOFsInt;

    // Data of pseudo-mean flow or target state to damp against on the
    // integration points in 2D.
    // Dimension: [iVar][iNode2D].
    as3data1d<as3double> DataDOFsIntMean;

    // Data of pseudo-mean flow or target state to damp against on the
    // integration points in 1D, per each face.
    // Dimension: [iFace][iVar][iNode1D].
    as3data2d<as3double> DataDOFsIntMeanFace;

    // Data auxiliary DOFs at integration points in 1D, per each PML face.
    // Dimension: [iFace][iVar][iNode1D].
    as3data2d<as3double> DataDOFsIntAuxFace;

    // Grid-stretching lagrange differential operator 1D (solution-to-integration).
    // Dimension: [iDim][iData].
    as3data1d<as3double> derLagrangeGridStretching1D;

    // Grid-stretching coefficients on each face (0D).
    // Dimension: [iFace].
    as3vector1d<as3double> derLagrangeGridStretchingFace;

    // Artificial-convection lagrange differential operator 1D (solution-to-integration).
    // Dimension: [iDim][iData].
    as3data1d<as3double> derLagrangeArtificialConvection1D;

    // Artificial-convection coefficients on each face (0D).
    // Dimension: [iFace].
    as3vector1d<as3double> derLagrangeArtificialConvectionFace;

	private:

};


class CEEData : public CData {

	public:
		// Constructor.
		CEEData(CConfig   	  *config_container,
						CGeometry 	  *geometry_container,
            CInitial      *initial_container,
						CElement  	  *element_container,
						unsigned short iZone,
						unsigned long  iElem);

		// Destructor.
		~CEEData(void) override;

	protected:

	private:

};


class CEESpongeData : public CEEData {

	public:
		// Constructor.
		CEESpongeData(CConfig   	  *config_container,
    						  CGeometry 	  *geometry_container,
                  CInitial      *initial_container,
    						  CElement  	  *element_container,
    						  unsigned short iZone,
    						  unsigned long  iElem);

		// Destructor.
		~CEESpongeData(void) override;

	protected:

	private:
    // Function that initializes and defines the grid-stretching coefficients.
    void InitializeGridStretching(CConfig       *config_container,
                                  CGeometry     *geometry_container,
                                  CInitial      *initial_container,
                                  CElement      *element_container,
                                  unsigned short iZone,
                                  unsigned long  iElem);

    // Function that initializes and defines the artificial velocity coefficients.
    void InitializeArtificialConvection(CConfig       *config_container,
                                        CGeometry     *geometry_container,
                                        CInitial      *initial_container,
                                        CElement      *element_container,
                                        unsigned short iZone,
                                        unsigned long  iElem);

    // Function that initializes and defines the sponge damping coefficients.
    void InitializeSpongeDamping(CConfig       *config_container,
                                 CGeometry     *geometry_container,
                                 CElement      *element_container,
                                 unsigned short iZone,
                                 unsigned long  iElem);

};


class CEEPMLData : public CEESpongeData {

	public:
		// Constructor.
		CEEPMLData(CConfig   	   *config_container,
  						 CGeometry 	   *geometry_container,
               CInitial      *initial_container,
  						 CElement  	   *element_container,
  						 unsigned short iZone,
  						 unsigned long  iElem);

		// Destructor.
		~CEEPMLData(void) final;

	protected:

	private:
};

