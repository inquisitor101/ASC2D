#pragma once

#include "import_structure.hpp"
#include <cmath>



class CProcess
{
	public:
		// Default constructor.
		CProcess(const unsigned long nTime);
		
		// Destructor.
		~CProcess(void);

		// Getter: returns Time.
		const as3vector1d<as3double> &GetTime(void) const {return Time;}
		// Getter: returns Data.
		const as3vector3d<as3double> &GetData(void) const {return Data;}

		// Function which computes the necessary metrics between temporal
		// sampled I0 to I1.
		void ComputeMetrics(const CImport      *import_A_container,
				                const CImport      *import_B_container,
												const unsigned long I0,
												const unsigned long I1);

	protected:
		// Total number of elements.
		unsigned long  nElem;
		// Total number of nodes.
		unsigned short nNode;

		// Number of metrics considered.
		unsigned short nMetric;
		// Number of variables processed.
		unsigned short nVarProc;

		// Processed data.
		// Dimension: [iTime][iMetric][iVar].
		as3vector3d<as3double> Data;
		// Time stamps.
		as3vector1d<as3double> Time;

	private:
		// Function which determines the scaling factor used.
		inline as3vector1d<as3double> DetermineScalingFactor(const as3vector2d<as3double> &data)
		{
			const as3double rho   = data[0][0];
			const as3double ovrho = 1.0/rho;
			const as3double u     = ovrho*data[1][0];
			const as3double v     = ovrho*data[2][0];
			const as3double p     = GAMMA_MINUS_ONE*( data[3][0]
					                  - 0.5*( u*data[1][0] + v*data[2][0] ) );

			// Assemble scale metric.
			as3vector1d<as3double> scale = {1.0/rho, 1.0/u, 1.0/v, 1.0/p};

			// If any of the entries is greater than one, then just use a factor of unity.
			for(unsigned short i=0; i<scale.size(); i++) 
				if( fabs( scale[i] ) >= 1.0 ) scale[i] = 1.0; 

			return scale;
		}
};

