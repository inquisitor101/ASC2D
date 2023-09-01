#pragma once

/*!
 * @file data_structure.hpp
 * @brief The file containing most of the element data.
 */

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "element_structure.hpp"
#include "initial_structure.hpp"


/*!
 * @brief A generic class used for storing most of the data on each element.
 */
class CData {

	public:
		/*!
		 * @brief Default constructor of CData, which initializes the class.
		 *
		 * @param[in] config_container pointer to the dictionary/configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the initial conditions container in this zone.
		 * @param[in] iZone input zone ID.
		 * @param[in] iElem input element ID.
		 */
		CData(CConfig   	  *config_container,
				  CGeometry 	  *geometry_container,
          CInitial      *initial_container,
					CElement  	  *element_container,
					unsigned short iZone,
					unsigned long  iElem);

		/*!
		 * @brief Destructor, which frees allocated memory.
		 */
		virtual ~CData(void);

    /*!
		 * @brief Function that resets the solution residual to zero.
		 */
    void ResetResidual(void){
      for(unsigned short i=0; i<DataDOFsRes.size(); i++)
        memset(DataDOFsRes[i], 0.0, nDOFsSol2D*sizeof(as3double));
    }

		/*!
		 * @brief Getter function which returns the value of MatchDOFsInt.
		 *
		 * @return MatchDOFsInt
		 */
		as3data1d<as3double> &GetMatchDOFsInt(void)                            {return MatchDOFsInt;}
		/*!
		 * @brief Getter function which returns the value of DataDOFsSol.
		 *
		 * @return DataDOFsSol
		 */
    as3data1d<as3double> &GetDataDOFsSol(void)                             {return DataDOFsSol;}
    /*!
		 * @brief Getter function which returns the value of DataDOFsResSol.
		 *
		 * @return DataDOFsResSol
		 */
    as3data1d<as3double> &GetDataDOFsRes(void)                             {return DataDOFsRes;}
    /*!
		 * @brief Getter function which returns the value of DampDOFsInt.
		 *
		 * @return DampDOFsInt
		 */
    as3data1d<as3double> &GetDampDOFsInt(void)                             {return DampDOFsInt;}
    /*!
		 * @brief Getter function which returns the value of DataDOFsIntMean.
		 *
		 * @return DataDOFsIntMean
		 */
    as3data1d<as3double> &GetDataDOFsIntMean(void)                         {return DataDOFsIntMean;}
    /*!
		 * @brief Getter function which returns the value of derLagrangeGridStretching1D.
		 *
		 * @return derLagrangeGridStretching1D
		 */
    as3data1d<as3double> &GetDerLagrangeGridStretching1D(void)             {return derLagrangeGridStretching1D;}
		/*!
		 * @brief Getter function which returns the value of derLagrangeGridStretchingFace[iFace].
		 *
		 * @param[in] iFace input face ID.
		 *
		 * @return derLagrangeGridStretchingFace[iFace]
		 */
    as3double GetDerLagrangeGridStretchingFace(unsigned short iFace)       {return derLagrangeGridStretchingFace[iFace];}
    /*!
		 * @brief Getter function which returns the value of derLagrangeArtificialConvection1D.
		 *
		 * @return derLagrangeArtificialConvection1D
		 */
    as3data1d<as3double> &GetDerLagrangeArtificialConvection1D(void)       {return derLagrangeArtificialConvection1D;}
    /*!
		 * @brief Getter function which returns the value of derLagrangeArtificialConvectionFace.
		 *
		 * @return derLagrangeArtificialConvectionFace.
		 */
    as3vector1d<as3double> &GetDerLagrangeArtificialConvectionFace(void)   {return derLagrangeArtificialConvectionFace;}
    /*!
		 * @brief Getter function which returns the value of derLagrangeArtificialConvectionFace[iFace].
		 *
		 * @param[in] iFace input face ID.
		 *
		 * @return derLagrangeArtificialConvectionFace[iFace]
		 */
    as3double GetDerLagrangeArtificialConvectionFace(unsigned short iFace) {return derLagrangeArtificialConvectionFace[iFace];}
    /*!
		 * @brief Getter function which returns the value of DataDOFsSol.size().
		 *
		 * @return DataDOFsSol.size()
		 */
    unsigned short GetnVarTotal(void)                                      {return DataDOFsSol.size();}
    /*!
		 * @brief Getter function which returns the value of DataDOFsIntFace[iFace].
		 *
		 * @param[in] iFace input face ID.
		 *
		 * @return DataDOFsIntFace[iFace]
		 */
    as3data1d<as3double> &GetDataDOFsIntFace(unsigned short iFace)         {return DataDOFsIntFace[iFace];}
    /*!
		 * @brief Getter function which returns the value of DataDOFsIntMeanFace, per input face.
		 *
		 * @return DataDOFsIntMeanFace[iFace]
		 */
    as3data1d<as3double> &GetDataDOFsIntMeanFace(unsigned short iFace)     {return DataDOFsIntMeanFace[iFace];}
    /*!
		 * @brief Getter function which returns the value of DataDOFsIntAuxFace, per input face.
		 *
		 * @return DataDOFsIntAuxFace[iFace]
		 */
    as3data1d<as3double> &GetDataDOFsIntAuxFace(unsigned short iFace)      {return DataDOFsIntAuxFace[iFace];}

	protected:
    unsigned short         nDOFsSol2D;                          ///< Number of solution DOFs in 2D.
		as3data1d<as3double>   DataDOFsSol;                         ///< Solution DOFs on the volume in 2D.
																																///< Dimension: [iVar][iNode2D].
    as3data2d<as3double>   DataDOFsIntFace;                     ///< Integration points on surface in 1D.
																																///< Dimension: [iFace][iVar][iNode1D].
    as3data1d<as3double>   DataDOFsRes;                         ///< Residual of solution DOFs on volume in 2D.
																																///< Dimension: [iVar][iNode2D].
    as3data1d<as3double>   DampDOFsInt;                         ///< Damping function on volume integration points in 2D.
																																///< Dimension: [iDim][iNode2D].
    as3data1d<as3double>   DataDOFsIntMean;                     ///< (pseudo-)Mean flow to damp against on volume integration points in 2D.
																																///< Dimension: [iVar][iNode2D].
    as3data2d<as3double>   DataDOFsIntMeanFace;                 ///< (pseudo-)Mean flow to damp against on surface integration points in 1D, per each face.
																																///< Dimension: [iFace][iVar][iNode1D].
    as3data2d<as3double>   DataDOFsIntAuxFace;                  ///< Auxiliary PML data on surface integration points in 1D, per each face.
																																///< Dimension: [iFace][iVar][iNode1D].
    as3data1d<as3double>   derLagrangeGridStretching1D;         ///< Modified Lagrange differential operator with grid-stretching in 1D (solution-to-integration).
																																///< Dimension: [iDim][iData].
    as3vector1d<as3double> derLagrangeGridStretchingFace;       ///< Grid-stretching coefficients on each face in 0D (i.e. points).
																																///< Dimension: [iFace].
    as3data1d<as3double>   derLagrangeArtificialConvection1D;   ///< Modified Lagrange differential opterator with artificial-convection in 1D (solution-to-integration).
																																///< Dimension: [iDim][iData].
    as3vector1d<as3double> derLagrangeArtificialConvectionFace; ///< Artificial-convection coefficients on each face in 0D (i.e. points).
																																///< Dimension: [iFace].
		as3data1d<as3double>   MatchDOFsInt;                        ///< Characteristic matching layer profile on volume integration points in 2D.
																																///< Dimension: [iDim][iNode2D].

	private:

};


/*!
 * @brief An interface class used for storing most of the data on each Euler-equations (EE) element.
 */
class CEEData : public CData {

	public:
		/*!
		 * @brief Default constructor of CEEData, which initializes the class.
		 *
		 * @param[in] config_container pointer to the dictionary/configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the initial conditions container in this zone.
		 * @param[in] iZone input zone ID.
		 * @param[in] iElem input element ID.
		 */
		CEEData(CConfig   	  *config_container,
						CGeometry 	  *geometry_container,
            CInitial      *initial_container,
						CElement  	  *element_container,
						unsigned short iZone,
						unsigned long  iElem);

		/*!
		 * @brief Destructor, which frees allocated memory.
		 */
		~CEEData(void) override;

	protected:

	private:

};

/*!
 * @brief A class used for storing most of the data on each sponge-layer element.
 */
class CEESpongeData : public CEEData {

	public:
		/*!
		 * @brief Default constructor of CEESpongeData, which initializes the class.
		 *
		 * @param[in] config_container pointer to the dictionary/configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the initial conditions container in this zone.
		 * @param[in] iZone input zone ID.
		 * @param[in] iElem input element ID.
		 */
		CEESpongeData(CConfig   	  *config_container,
    						  CGeometry 	  *geometry_container,
                  CInitial      *initial_container,
    						  CElement  	  *element_container,
    						  unsigned short iZone,
    						  unsigned long  iElem);

		/*!
		 * @brief Destructor, which frees allocated memory.
		 */
		~CEESpongeData(void) override;

	protected:

	private:
    /*!
		 * @brief Function that initializes and defines the grid-stretching coefficients.
		 */
    void InitializeGridStretching(CConfig       *config_container,
                                  CGeometry     *geometry_container,
                                  CInitial      *initial_container,
                                  CElement      *element_container,
                                  unsigned short iZone,
                                  unsigned long  iElem);

    /*!
		 * @brief Function that initializes and defines the artificial velocity coefficients.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the initial conditions container in this zone.
		 * @param[in] element_container pointer to the standard element container in this zone.
		 * @param[in] iZone input zone ID.
		 * @param[in] iElem input element ID.
		 */
    void InitializeArtificialConvection(CConfig       *config_container,
                                        CGeometry     *geometry_container,
                                        CInitial      *initial_container,
                                        CElement      *element_container,
                                        unsigned short iZone,
                                        unsigned long  iElem);

    /*!
		 * @brief Function that initializes and defines the sponge damping coefficients.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] element_container pointer to the standard element container in this zone.
		 * @param[in] iZone input zone ID.
		 * @param[in] iElem input element ID.
		 */
    void InitializeSpongeDamping(CConfig       *config_container,
                                 CGeometry     *geometry_container,
                                 CElement      *element_container,
                                 unsigned short iZone,
                                 unsigned long  iElem);


    /*!
		 * @brief Function that initializes and defines the characteristic matching coefficients.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] element_container pointer to the standard element container in this zone.
		 * @param[in] iZone input zone ID.
		 * @param[in] iElem input element ID.
		 */
    void InitializeCharacteristicMatching(CConfig       *config_container,
                                          CGeometry     *geometry_container,
                                          CElement      *element_container,
                                          unsigned short iZone,
                                          unsigned long  iElem);
};

/*!
 * @brief A class used for storing most of the data on each PML element.
 */
class CEEPMLData : public CEESpongeData {

	public:
		/*!
		 * @brief Default constructor of CEEPMLData, which initializes the class.
		 *
		 * @param[in] config_container pointer to the dictionary/configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] initial_container pointer to the initial conditions container in this zone.
		 * @param[in] iZone input zone ID.
		 * @param[in] iElem input element ID.
		 */
		CEEPMLData(CConfig   	   *config_container,
  						 CGeometry 	   *geometry_container,
               CInitial      *initial_container,
  						 CElement  	   *element_container,
  						 unsigned short iZone,
  						 unsigned long  iElem);

		/*!
		 * @brief Destructor, which frees allocated memory.
		 */
		~CEEPMLData(void) final;

	protected:

	private:
};

