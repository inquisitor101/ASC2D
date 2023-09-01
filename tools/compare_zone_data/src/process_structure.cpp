#include "process_structure.hpp"



CProcess::CProcess
(
 const unsigned long nTime
)
 /*
	* Constructor, that initializes a process container. 
	*/
{
	// Specify the number of elements and nodes.
	nElem = nElemExpected;
	nNode = nNodeExpected;

	// Number of metrics considered.
	// These are: Linf, L1 and L2 norm.
	nMetric  = 3;

	// Number of variables considered in the processing of each metric.
	// Note, these are: rho, u, v, p.
	nVarProc = 4;

	// Reserve memory for the temporal stamps.
	Time.resize(nTime, 0.0);

	// Reserve memory for the processed data.
	Data.resize(nTime);
	
	// Reserve memory for each time sample.
	for(unsigned long iTime=0; iTime<Data.size(); iTime++)
	{
		// Reserve memory for each metric.
		// Note, here 3 are used to indicate: Linf, L1 and L2 norm.
		Data[iTime].resize(nMetric);
		
		// Reserve the memory for each variable per metric.
		for(unsigned short iMetric=0; iMetric<Data[iTime].size(); iMetric++)
		{
			// Reserve memory for each variable considered.
			// Note, here 4 are used to indicate: rho, u, v, p.
			Data[iTime][iMetric].resize(nVarProc, 0.0);
		}
	}
}


CProcess::~CProcess
(
 void
)
 /*
	* Destructor for CProcess, frees dynamic memory.
	*/
{

}


void CProcess::ComputeMetrics
(
 const CImport      *import_A_container,
 const CImport      *import_B_container,
 const unsigned long I0,
 const unsigned long I1
)
 /*
	* Function which computes the specified metrices between
	* indices: I0 and I1.
	*/
{
	// Report progress.
	std::cout << "   Processing data....";

	// Abbreviation involving gamma.
	const as3double gm1 = GAMMA_MINUS_ONE;
	
	// Obtain the temporal data, it does not matter from which container.
	auto& time   = import_A_container->GetSimTime();
	// Obtain the data, based on the container: A.
	auto& data_A = import_A_container->GetImportedData();
	// Obtain the data, based on the container: B.
	auto& data_B = import_B_container->GetImportedData();

	// Determine scaling factor, which is used to alleviate potential overflows.
	// Note, this is only needed in the L1 and L2 computations for the momentum
	// and pressure. For simplicity, take any (first?) entry in the reference 
	// simulation (B).
	as3vector1d<as3double> scale = DetermineScalingFactor( data_B[0][0] );

	// Consistency check.
	if( I1-I0 != data_A.size() || data_A.size() != data_B.size() ) ERROR("data sizes do not match.");

	// Loop over the selected files and process them.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
	for(unsigned long iTime=I0; iTime<I1; iTime++)
	{
		// Deduce actual local index, based on the chunk division.
		const unsigned long IDX = iTime-I0;

		// Initialize   L1-norm of variables.
		as3vector1d<as3double>   L1(nVarProc, 0.0);
		// Initialize   L2-norm of variables.
		as3vector1d<as3double>   L2(nVarProc, 0.0);
		// Initialize Linf-norm of variables.
		as3vector1d<as3double> Linf(nVarProc, 0.0);
		
		// Loop over each of the elements.
		for(unsigned long i=0; i<nElem; i++)
		{
			// Loop over the nodes owned by the element.
			for(unsigned short l=0; l<nNode; l++)
			{
				// Compute the primitive variables of container: A.
				const as3double A_r   = data_A[IDX][i][0][l];
				const as3double A_ovr = 1.0/A_r;
				const as3double A_u   = A_ovr*    data_A[IDX][i][1][l];
				const as3double A_v   = A_ovr*    data_A[IDX][i][2][l];
	    	const as3double A_p   = gm1*(     data_A[IDX][i][3][l]
    	                        - 0.5*( A_u*data_A[IDX][i][1][l] 
															      + A_v*data_A[IDX][i][2][l] ) );

				// Compute the primitive variables of container: B.
				const as3double B_r   = data_B[IDX][i][0][l];
				const as3double B_ovr = 1.0/B_r;
				const as3double B_u   = B_ovr*    data_B[IDX][i][1][l];
				const as3double B_v   = B_ovr*    data_B[IDX][i][2][l];
	    	const as3double B_p   = gm1*(     data_B[IDX][i][3][l]
    	                        - 0.5*( B_u*data_B[IDX][i][1][l] 
															      + B_v*data_B[IDX][i][2][l] ) );

				// Compute the difference between the primitive variables.
				// Note, it is assumed that container B is the reference data.
				const as3double dr = scale[0]*fabs( A_r - B_r );
				const as3double du = scale[1]*fabs( A_u - B_u );
				const as3double dv = scale[2]*fabs( A_v - B_v );
				const as3double dp = scale[3]*fabs( A_p - B_p );

				// Compute:   L1-norm.
				L1[0] += dr;
				L1[1] += du;
				L1[2] += dv;
				L1[3] += dp;

				// Compute:   L2-norm.
				L2[0] += dr*dr;
				L2[1] += du*du;
				L2[2] += dv*dv;
				L2[3] += dp*dp;

				// Compute: Linf-norm.
				Linf[0] = std::max( Linf[0], dr );
				Linf[1] = std::max( Linf[1], du );
				Linf[2] = std::max( Linf[2], dv );
				Linf[3] = std::max( Linf[3], dp );
			}
		}

		// Book-keep values and normalize L1 and L2 norms.
#pragma omp simd
		for(unsigned short k=0; k<nVarProc; k++)
		{
			// Divide by the scale, in order to obtain the correct
			// ratio. Note, the scale is used to alleviate potential
			// overflow problems.
			const as3double ovs = 1.0/scale[k];

			// L1-norm.
			Data[iTime][0][k] = ovs*L1[k];
			// L2-norm.
			Data[iTime][1][k] = ovs*sqrt( L2[k] );
			// Linf-norm.
			Data[iTime][2][k] = ovs*Linf[k];
		}

		// Copy temporal values.
		Time[iTime] = time[IDX]; 
	}

	// Report progress.
	std::cout << " Done." << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;
}






