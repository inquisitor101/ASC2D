#pragma once





// NOTE,
// the rational behind these implementations of the tensor product
// is to breakdown the outer loop, in order to reduce computations.
// That is,
// instead of doing the following:
//
// for( iRowJ < M )
//   for( iRowI < M )
//     for( iColJ < K )
//       for( iColI < K )
//         ComputeTensorProduct()
//
// ... which is: M2*K2 operations.
//
// We do the below:
//
// for( i < K )
//   for( jj < K )
//     for( j < M )
//       ComputeTensorProduct()
//
// ... then follow up with:
//
// for( j < M )
//   for( ii < K )
//     for( i < M )
//       ComputeTensorProduct()
//
// ... these yield: K2*M + M2*K operations.
//
// Thus, this is cheaper than a standard M2*K2 approach.
//

// Function that computes the tensor product on a volume.
void TensorProductSolAndGradVolume(const unsigned short    M,
																	 const unsigned short    N,
																	 const unsigned short    K,
																	 const as3double        *A,
																	 const as3double        *ADer,
																	 const as3double *const *B,
																	 as3double             **C,
																	 as3double             **CDerR,
																	 as3double             **CDerS);


// Function that serves as an interface to select which direction the
// tensor products on a face are used: I or J, depending on the input iFace.
void TensorProductSolAndGradFace(const unsigned short    iFace,
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
																 as3double             **CDerS);


// Function that computes the tensor product on a face in the I-direction.
void TensorProductSolAndGradFaceI(const unsigned short    M,
																	const unsigned short    N,
																	const unsigned short    K,
																	const unsigned short   *IDX,
																	const as3double        *A,
																	const as3double        *ADer,
																	const as3double        *ADerFace,
																	const as3double *const *B,
																	as3double             **C,
																	as3double             **CDerR,
																	as3double             **CDerS);


// Function that computes the tensor product on a face in the J-direction.
void TensorProductSolAndGradFaceJ(const unsigned short    M,
																	const unsigned short    N,
																	const unsigned short    K,
																	const unsigned short   *IDX,
																	const as3double        *A,
																	const as3double        *ADer,
																	const as3double        *ADerFace,
																	const as3double *const *B,
																	as3double             **C,
																	as3double             **CDerR,
																	as3double             **CDerS);


// Function that computes the tensor product on a volume and adds the
// contribution to the residual accordingly.
void TensorProductVolumeResidual(const unsigned short M,
                                 const unsigned short N,
                                 const unsigned short K,
                                 const as3double     *Ar,
                                 const as3double     *As,
                                 as3double          **B,
                                 as3double          **C);

// Function that computes the tensor product on a surface and adds the
// contribution to the residual accordingly.
void TensorProductSurfaceResidual(const unsigned short  M,
                                  const unsigned short  N,
                                  const unsigned short  K,
                                  const unsigned short *IDX,
                                  const as3double      *A,
                                  as3double           **B,
                                  as3double           **C);
