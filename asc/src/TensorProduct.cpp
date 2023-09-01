#include "ASC2D.hpp"



void TensorProductSolAndGradVolume
(
 const unsigned short    M,
 const unsigned short    N,
 const unsigned short    K,
 const as3double        *A,
 const as3double        *ADer,
 const as3double *const *B,
 as3double             **C,
 as3double             **CDerR,
 as3double             **CDerS
)
 /*
	* Function that performs a tensor product multiplication.
	* Input,
	*   M    : Dimension of DOFs of output points in 1D.
	*   N    : Column dimension of both input and output vector.
	*   K    : Dimension of DOFs of input points in 1D.
	*   A    : Polynomial basis of interpolation in 1D (M-by-K).
	*   ADer : Derivative of polynomial basis of interpolation in 1D (M-by-K).
	*   B    : Input 2D data ( (K*K)-by-N ).
	*   C    : Output interpolated data ( (M*M)-by-N ).
	*   CDerR: Output x-derivative of interpolated data ( (M*M)-by-N ).
	*   CDerS: Output y-derivative of interpolated data ( (M*M)-by-N ).
	*   NOTE,
	*   A and ADer are transposed (for efficiency), since C++ uses a row-major format.
	*   Hence, in essense they are: ( K-by-M ).
	*/
{
	// Cast the arrays from 1D to 2D.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Loop over each N entry.
	for(unsigned short l=0; l<N; l++){

		// Create temporary storage.
		as3double tmpJ[K][M];
		as3double tmpI[M][M];

		// Cast array B from 2D to 3D. Note, 3D when including [l] index.
		const as3double (*b)[K] = (const as3double (*)[K]) B[l];

		for(unsigned short i=0; i<K; i++){
#pragma omp simd
			for(unsigned short s=0; s<M; s++) tmpJ[i][s] = 0.0;
			for(unsigned short jj=0; jj<K; jj++){
#pragma omp simd
        for(unsigned short j=0; j<M; j++)
					tmpJ[i][j] += a[jj][j] * b[jj][i];
			}
		}


		if( C ){

			// Cast array C from 2D to 3D. Note, 3D when including [l] index.
			as3double (*c)[M] = (as3double (*)[M]) C[l];

			for(unsigned short j=0; j<M; j++){
#pragma omp simd
        for(unsigned short s=0; s<M; s++) tmpI[j][s] = 0.0;
				for(unsigned short ii=0; ii<K; ii++){
#pragma omp simd
          for(unsigned short i=0; i<M; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to C.
			for(unsigned short j=0; j<M; j++)
#pragma omp simd
				for(unsigned short i=0; i<M; i++)
					c[j][i] = tmpI[j][i];
		}


		if( CDerR ){

			// Cast array CDerR from 2D to 3D. Note, 3D when including [l] index.
			as3double (*cDerR)[M] = (as3double (*)[M]) CDerR[l];

			for(unsigned short j=0; j<M; j++){
#pragma omp simd
				for(unsigned short s=0; s<M; s++) tmpI[j][s] = 0.0;
				for(unsigned short ii=0; ii<K; ii++){
#pragma omp simd
          for(unsigned short i=0; i<M; i++)
						tmpI[j][i] += aDer[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerR.
			for(unsigned short j=0; j<M; j++)
#pragma omp simd
				for(unsigned short i=0; i<M; i++)
					cDerR[j][i] = tmpI[j][i];
		}


		if( CDerS ){

			// Cast array CDerS from 2D to 3D. Note, 3D when including [l] index.
			as3double (*cDerS)[M] = (as3double (*)[M]) CDerS[l];

			for(unsigned short i=0; i<K; i++){
#pragma omp simd
        for(unsigned short s=0; s<M; s++) tmpJ[i][s] = 0.0;
				for(unsigned short jj=0; jj<K; jj++){
#pragma omp simd
          for(unsigned short j=0; j<M; j++)
						tmpJ[i][j] += aDer[jj][j] * b[jj][i];
				}
			}

			for(unsigned short j=0; j<M; j++){
#pragma omp simd
        for(unsigned short s=0; s<M; s++) tmpI[j][s] = 0.0;
				for(unsigned short ii=0; ii<K; ii++){
#pragma omp simd
          for(unsigned short i=0; i<M; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerS.
			for(unsigned short j=0; j<M; j++)
#pragma omp simd
				for(unsigned short i=0; i<M; i++)
					cDerS[j][i] = tmpI[j][i];
		}
	}
}


void TensorProductSolAndGradFace
(
 const unsigned short    iFace,
 const unsigned short    M,
 const unsigned short    N,
 const unsigned short    K,
 const unsigned short   *IDX,
 const as3double        *A,
 const as3double        *ADer,
 const as3double        *ADerFace,
 const as3double *const *B,
 as3double             **C,
 as3double             **CDerR,
 as3double             **CDerS
)
 /*
	* Function that serves as an interface to the relevant direction of the
	* face tensor-product functions.
	*/
{
	// Call J-direction tensor product, if the face is in the s-direction.
	// Otherwise, use the I-direction variant, since it must be in the r-direction.
	if( iFace == IDX_SOUTH || iFace == IDX_NORTH )
		TensorProductSolAndGradFaceJ(M, N, K, IDX,
																 A, ADer, ADerFace,
																 B, C, CDerR, CDerS);
	else
		TensorProductSolAndGradFaceI(M, N, K, IDX,
																 A, ADer, ADerFace,
																 B, C, CDerR, CDerS);
}


void TensorProductSolAndGradFaceI
(
 const unsigned short    M,
 const unsigned short    N,
 const unsigned short    K,
 const unsigned short   *IDX,
 const as3double        *A,
 const as3double        *ADer,
 const as3double        *ADerFace,
 const as3double *const *B,
 as3double             **C,
 as3double             **CDerR,
 as3double             **CDerS
)
 /*
  * Function that computes the tensor product on a face in the I-direction.
	* Input,
	*   M    	  : Dimension of DOFs of output points in 1D.
	*   N    	  : Column dimension of both input and output vector.
	*   K    	  : Dimension of DOFs of input points in 1D.
	*   IDX     : Indices on B that return the surface data in 1D of current face (K-by-1).
	*   A    	  : Polynomial basis of interpolation in 1D (M-by-K).
	*   ADer 	  : Derivative of polynomial basis of interpolation in 1D (M-by-K).
	*   ADerFace: Derivative of polynomial basis of interpolation at face in 1D (K-by-1).
	*   B    		: Input 2D data ( (K*K)-by-N ).
	*   C    		: Output 1D interpolated data ( M-by-N ).
	*   CDerR		: Output 1D x-derivative of interpolated data ( M-by-N ).
	*   CDerS		: Output 1D y-derivative of interpolated data ( M-by-N ).
	*   NOTE,
	*   A and ADer are transposed (for efficiency), since C++ uses a row-major format.
	*   Hence, in essense they are: ( K-by-M ).
  */
{
	// Cast the arrays from 1D to 2D.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Rename ADerFace in order to be consistent with other arrays.
	const as3double *aDerFace = ADerFace;

	// Loop over each N entry.
	for(unsigned short l=0; l<N; l++){

		// Create temporary storage.
		as3double tmpI[M], tmpJ[K];

		// Cast array B from 2D to 1D.
		const as3double *bb = B[l];
		// Cast array B from 2D to 3D. Note, 3D when including [l] index.
		const as3double (*b)[K] = (const as3double (*)[K]) B[l];

		// Store boundary surface info in tmpJ.
// #pragma omp simd safelen(K)
    for(unsigned short s=0; s<K; s++) tmpJ[s] = bb[IDX[s]];

		if( C ){

			// Cast array C from 2D to 1D.
			as3double *c = C[l];

#pragma omp simd
			for(unsigned short s=0; s<M; s++) tmpI[s] = 0.0;
			for(unsigned short k=0; k<M; k++){
// #pragma omp simd safelen(K)
				for(unsigned short jj=0; jj<K; jj++)
					tmpI[k] += a[jj][k] * tmpJ[jj];
			}

			// Copy values to C.
#pragma omp simd
      for(unsigned short j=0; j<M; j++)
				c[j] = tmpI[j];
		}

		if( CDerS ){

			// Cast array CDerS from 2D to 1D.
			as3double *cDerS = CDerS[l];

#pragma omp simd
			for(unsigned short s=0; s<M; s++) tmpI[s] = 0.0;
			for(unsigned short k=0; k<M; k++){
// #pragma omp simd safelen(K)
        for(unsigned short jj=0; jj<K; jj++)
					tmpI[k] += aDer[jj][k] * tmpJ[jj];
			}

			// Copy values to CDerS.
#pragma omp simd
      for(unsigned short j=0; j<M; j++)
				cDerS[j] = tmpI[j];
		}

		if( CDerR ){

			// Cast array CDerR from 2D to 1D.
			as3double *cDerR = CDerR[l];

// #pragma omp simd safelen(K)
			for(unsigned short s=0; s<K; s++) tmpJ[s] = 0.0;
			for(unsigned short j=0; j<K; j++){
// #pragma omp simd safelen(K)
        for(unsigned short ii=0; ii<K; ii++)
					tmpJ[j] += aDerFace[ii] * b[j][ii];
			}

#pragma omp simd
			for(unsigned short s=0; s<M; s++) tmpI[s] = 0.0;
			for(unsigned short k=0; k<M; k++){
// #pragma omp simd safelen(K)
        for(unsigned short jj=0; jj<K; jj++)
					tmpI[k] += a[jj][k] * tmpJ[jj];
			}

			// Copy values to CDerR.
#pragma omp simd
			for(unsigned short j=0; j<M; j++)
				cDerR[j] = tmpI[j];
		}
	}
}


void TensorProductSolAndGradFaceJ
(
 const unsigned short    M,
 const unsigned short    N,
 const unsigned short    K,
 const unsigned short   *IDX,
 const as3double        *A,
 const as3double        *ADer,
 const as3double        *ADerFace,
 const as3double *const *B,
 as3double             **C,
 as3double             **CDerR,
 as3double             **CDerS
)
 /*
  * Function that computes the tensor product on a face in the J-direction.
	* Input,
	*   M    	  : Dimension of DOFs of output points in 1D.
	*   N    	  : Column dimension of both input and output vector.
	*   K    	  : Dimension of DOFs of input points in 1D.
	*   IDX     : Indices on B that return the surface data in 1D of current face (K-by-1).
	*   A    	  : Polynomial basis of interpolation in 1D (M-by-K).
	*   ADer 	  : Derivative of polynomial basis of interpolation in 1D (M-by-K).
	*   ADerFace: Derivative of polynomial basis of interpolation at face in 1D (K-by-1).
	*   B    		: Input 2D data ( (K*K)-by-N ).
	*   C    		: Output 1D interpolated data ( M-by-N ).
	*   CDerR		: Output 1D x-derivative of interpolated data ( M-by-N ).
	*   CDerS		: Output 1D y-derivative of interpolated data ( M-by-N ).
	*   NOTE,
	*   A and ADer are transposed (for efficiency), since C++ uses a row-major format.
	*   Hence, in essense they are: ( K-by-M ).
  */
{
	// Cast the arrays from 1D to 2D.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Rename ADerFace in order to be consistent with other arrays.
	const as3double *aDerFace = ADerFace;

	// Loop over each N entry.
	for(unsigned short l=0; l<N; l++){

		// Create temporary storage.
		as3double tmpI[M], tmpJ[K];

		// Cast array B from 2D to 1D.
		const as3double *bb = B[l];
		// Cast array B from 2D to 3D. Note, 3D when including [l] index.
		const as3double (*b)[K] = (const as3double (*)[K]) B[l];

		// Store boundary surface info in tmpJ.
// #pragma omp simd safelen(K)
    for(unsigned short s=0; s<K; s++) tmpJ[s] = bb[IDX[s]];

		if( C ){

			// Cast array C from 2D to 1D.
			as3double *c = C[l];

#pragma omp simd
			for(unsigned short s=0; s<M; s++) tmpI[s] = 0.0;
			for(unsigned short k=0; k<M; k++){
// #pragma omp simd safelen(K)
        for(unsigned short ii=0; ii<K; ii++)
					tmpI[k] += a[ii][k] * tmpJ[ii];
			}

			// Copy values to C.
#pragma omp simd
			for(unsigned short i=0; i<M; i++)
				c[i] = tmpI[i];
		}

		if( CDerR ){

			// Cast array CDerR from 2D to 1D.
			as3double *cDerR = CDerR[l];

#pragma omp simd
			for(unsigned short s=0; s<M; s++) tmpI[s] = 0.0;
			for(unsigned short k=0; k<M; k++){
// #pragma omp simd safelen(K)
        for(unsigned short ii=0; ii<K; ii++)
					tmpI[k] += aDer[ii][k] * tmpJ[ii];
			}

			// Copy values to CDerR.
#pragma omp simd
			for(unsigned short i=0; i<M; i++)
				cDerR[i] = tmpI[i];
		}

		if( CDerS ){

			// Cast array CDerS from 2D to 1D.
			as3double *cDerS = CDerS[l];

// #pragma omp simd safelen(K)
			for(unsigned short s=0; s<K; s++) tmpJ[s] = 0.0;
			for(unsigned short i=0; i<K; i++){
// #pragma omp simd safelen(K)
        for(unsigned short jj=0; jj<K; jj++)
					tmpJ[i] += aDerFace[jj] * b[jj][i];
			}

#pragma omp simd
			for(unsigned short s=0; s<M; s++) tmpI[s] = 0.0;
			for(unsigned short k=0; k<M; k++){
// #pragma omp simd safelen(K)
        for(unsigned short ii=0; ii<K; ii++)
					tmpI[k] += a[ii][k] * tmpJ[ii];
			}

			// Copy values to CDerS.
#pragma omp simd
			for(unsigned short i=0; i<M; i++)
				cDerS[i] = tmpI[i];
		}
	}
}


void TensorProductVolumeResidual
(
 const unsigned short M,
 const unsigned short N,
 const unsigned short K,
 const as3double     *Ar,
 const as3double     *As,
 as3double          **B,
 as3double          **C
)
 /*
	* Function that performs a tensor product multiplication and adds the
  * volume terms contribution to the residual.
	* Input,
	*   M  : Dimension of DOFs of input points in 1D.
	*   N  : Column dimension of both input and output vector.
	*   K  : Dimension of DOFs of output points in 1D.
	*   Ar : Polynomial basis of interpolation in 1D in r-direction (K-by-M).
	*   As : Polynomial basis of interpolation in 1D in s-direction (K-by-M).
	*   B  : Input 2D data ( (M*M)-by-N ).
	*   C  : Output residual data ( (K*K)-by-N ).
	*   NOTE,
	*   Ar and As are not transposed, unlike in the TensorProductVolumeResidual.
	*/
{
	// Cast the arrays from 1D to 2D.
	const as3double (*ar)[K] = (const as3double (*)[K]) Ar;
	const as3double (*as)[K] = (const as3double (*)[K]) As;

	// Loop over each N entry.
	for(unsigned short l=0; l<N; l++){

		// Create temporary storage.
		as3double tmpJ[M][K];
		as3double tmpI[K][K];

		// Cast array B from 2D to 3D. Note, 3D when including [l] index.
		const as3double (*b)[M] = (const as3double (*)[M]) B[l];

    // Cast array C from 2D to 3D. Note, 3D when including [l] index.
    as3double (*c)[K] = (as3double (*)[K]) C[l];


    // Tensor product in s-direction.
		for(unsigned short i=0; i<M; i++){
// #pragma omp simd safelen(K)
  		for(unsigned short s=0; s<K; s++) tmpJ[i][s] = 0.0;
			for(unsigned short jj=0; jj<M; jj++){
// #pragma omp simd safelen(K)
        for(unsigned short j=0; j<K; j++)
					tmpJ[i][j] += as[jj][j] * b[jj][i];
			}
		}

    // Tensor product in r-direction.
    for(unsigned short j=0; j<K; j++){
// #pragma omp simd safelen(K)
      for(unsigned short s=0; s<K; s++) tmpI[j][s] = 0.0;
			for(unsigned short ii=0; ii<M; ii++){
// #pragma omp simd safelen(K)
        for(unsigned short i=0; i<K; i++)
					tmpI[j][i] += ar[ii][i] * tmpJ[ii][j];
			}
		}

		// Copy values to C.
		for(unsigned short j=0; j<K; j++)
#pragma omp simd
			for(unsigned short i=0; i<K; i++)
				c[j][i] += tmpI[j][i];
  }
}


void TensorProductSurfaceResidual
(
 const unsigned short  M,
 const unsigned short  N,
 const unsigned short  K,
 const unsigned short *IDX,
 const as3double      *A,
 as3double           **B,
 as3double           **C
)
 /*
	* Function that performs a tensor product multiplication and adds the
  * surface term contribution to the residual.
	* Input,
	*   M  : Dimension of DOFs of input points in 1D.
	*   N  : Column dimension of both input and output vector.
	*   K  : Dimension of DOFs of output points in 1D.
  *   IDX: Indices on B that return the surface data in 1D of current face (K-by-1).
	*   A  : Polynomial basis of interpolation in 1D (K-by-M).
	*   B  : Input 1D data ( (M-by-N ).
	*   C  : Output residual data ( (K*K)-by-N ).
	*   NOTE,
	*   A is not transposed, unlike in the TensorProductVolumeResidual.
	*/
{
  // Cast the arrays from 1D to 2D.
	const as3double (*a)[K] = (const as3double (*)[K]) A;

  // Loop over each N entry.
  for(unsigned short l=0; l<N; l++){

    // Create temporary storage.
		as3double tmpI[K];

		// Cast array B from 2D to 1D.
		const as3double *b = B[l];
    // Cast array C from 2D to 1D.
    as3double *c = C[l];

    // Tensor product on a face.
// #pragma omp simd safelen(K)
    for(unsigned short s=0; s<K; s++) tmpI[s] = 0.0;
    for(unsigned short k=0; k<K; k++){
#pragma omp simd
      for(unsigned short ii=0; ii<M; ii++)
        tmpI[k] += a[ii][k] * b[ii];
    }

    // Copy values to C.
// #pragma omp simd safelen(K)
    for(unsigned short s=0; s<K; s++) c[IDX[s]] += tmpI[s];
  }
}




