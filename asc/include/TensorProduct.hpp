#pragma once

/*!
 * @file TensorProduct.hpp
 * @brief The file containing all the tensor product implementations.
 *
 * Note,
 * The rational behind these implementations of the tensor product
 * is to breakdown the outer loop, in order to reduce computations.
 * That is, instead of doing the following:
 *
 * for( iRowJ < M )
 *   for( iRowI < M )
 *     for( iColJ < K )
 *       for( iColI < K )
 *         ComputeTensorProduct()
 *
 * ... which is: M2*K2 operations.
 *
 * We do the below:
 *
 * for( i < K )
 *   for( jj < K )
 *     for( j < M )
 *       ComputeTensorProduct()
 *
 * ... then follow up with:
 *
 * for( j < M )
 *   for( ii < K )
 *     for( i < M )
 *       ComputeTensorProduct()
 *
 * ... these yield: K2*M + M2*K operations.
 *
 * Thus, this is cheaper than a standard M2*K2 approach.
 */


/*!
 * @brief Function that computes the tensor product on a volume.
 *
 * @param[in] M dimension of DOFs of output points in 1D.
 * @param[in] N column dimension of both input and output vector.
 * @param[in] K dimension of DOFs of input points in 1D.
 * @param[in] A pointer to the interpolation polynomial in 1D (M-by-K).
 * @param[in] ADer pointer to the derivative of the interpolation polynomial in 1D (M-by-K).
 * @param[in] B pointer to input 2D data ( (K*K)-by-N ).
 * @param[out] C pointer to output interpolated data ( (M*M)-by-N ).
 * @param[out] CDerR pointer to output x-derivative of interpolated data ( (M*M)-by-N ).
 * @param[out] CDerS pointer to output y-derivative of interpolated data ( (M*M)-by-N ).
 */
void TensorProductSolAndGradVolume(const unsigned short    M,
																	 const unsigned short    N,
																	 const unsigned short    K,
																	 const as3double        *A,
																	 const as3double        *ADer,
																	 const as3double *const *B,
																	 as3double             **C,
																	 as3double             **CDerR,
																	 as3double             **CDerS);


/*!
 * @brief Function that serves as an interface to tensor product on a face.
 * This selects which direction the tensor products on a face are used.
 * For instance, could be: I or J, depending on the input iFace.
 *
 * @param[in] iFace input face ID.
 * @param[in] M dimension of DOFs of output points in 1D.
 * @param[in] N column dimension of both input and output vector.
 * @param[in] K dimension of DOFs of input points in 1D.
 * @param[in] IDX pointer to the indices of a 2D element on the current face in 1D.
 * @param[in] A pointer to the interpolation polynomial in 1D (M-by-K).
 * @param[in] ADer pointer to the derivative of the interpolation polynomial in 1D (M-by-K).
 * @param[in] ADerFace pointer to the derivative of the interpolation polynomial at a face in 1D (K-by-1).
 * @param[in] B pointer to input 2D data ( (K*K)-by-N ).
 * @param[out] C pointer to output interpolated data in 1D ( M-by-N ).
 * @param[out] CDerR pointer to output x-derivative of interpolated data in 1D ( M-by-N ).
 * @param[out] CDerS pointer to output y-derivative of interpolated data in 1D ( M-by-N ).
 */
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


/*!
 * @brief Function that computes the tensor product on a face in the I(x-)direction.
 *
 * @param[in] M dimension of DOFs of output points in 1D.
 * @param[in] N column dimension of both input and output vector.
 * @param[in] K dimension of DOFs of input points in 1D.
 * @param[in] IDX pointer to the indices of a 2D element on the current face in 1D.
 * @param[in] A pointer to the interpolation polynomial in 1D (M-by-K).
 * @param[in] ADer pointer to the derivative of the interpolation polynomial in 1D (M-by-K).
 * @param[in] ADerFace pointer to the derivative of the interpolation polynomial at a face in 1D (K-by-1).
 * @param[in] B pointer to input 2D data ( (K*K)-by-N ).
 * @param[out] C pointer to output interpolated data in 1D ( M-by-N ).
 * @param[out] CDerR pointer to output x-derivative of interpolated data in 1D ( M-by-N ).
 * @param[out] CDerS pointer to output y-derivative of interpolated data in 1D ( M-by-N ).
 */
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


/*!
 * @brief Function that computes the tensor product on a face in the J(y-)direction.
 *
 * @param[in] M dimension of DOFs of output points in 1D.
 * @param[in] N column dimension of both input and output vector.
 * @param[in] K dimension of DOFs of input points in 1D.
 * @param[in] IDX pointer to the indices of a 2D element on the current face in 1D.
 * @param[in] A pointer to the interpolation polynomial in 1D (M-by-K).
 * @param[in] ADer pointer to the derivative of the interpolation polynomial in 1D (M-by-K).
 * @param[in] ADerFace pointer to the derivative of the interpolation polynomial at a face in 1D (K-by-1).
 * @param[in] B pointer to input 2D data ( (K*K)-by-N ).
 * @param[out] C pointer to output interpolated data in 1D ( M-by-N ).
 * @param[out] CDerR pointer to output x-derivative of interpolated data in 1D ( M-by-N ).
 * @param[out] CDerS pointer to output y-derivative of interpolated data in 1D ( M-by-N ).
 */
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


/*!
 * @brief Function that computes the tensor product on a volume.
 * Also, this adds the contribution to the residual accordingly.
 *
 * @param[in] M dimension of DOFs of input points in 1D.
 * @param[in] N column dimension of both input and output vector.
 * @param[in] K dimension of DOFs of output points in 1D.
 * @param[in] Ar pointer to the interpolation polynomial in 1D in r-direction (K-by-M).
 * @param[in] As pointer to the interpolation polynomial in 1D in s-direction (K-by-M).
 * @param[in] B pointer to input 2D data ( (M*M)-by-N ).
 * @param[out] C pointer to output residual data ( (K*K)-by-N ).
 */
void TensorProductVolumeResidual(const unsigned short M,
                                 const unsigned short N,
                                 const unsigned short K,
                                 const as3double     *Ar,
                                 const as3double     *As,
                                 as3double          **B,
                                 as3double          **C);

/*!
 * @brief Function that computes the tensor product on a surface.
 * Also, this adds the contribution to the residual accordingly.
 *
 * @param[in] M dimension of DOFs of input points in 1D.
 * @param[in] N column dimension of both input and output vector.
 * @param[in] K dimension of DOFs of output points in 1D.
 * @param[in] IDX pointer to the indices of a 2D element on the current face in 1D (K-by-1).
 * @param[in] A pointer to the interpolation polynomial in 1D (K-by-M).
 * @param[in] B pointer to input 1D data ( M-by-N ).
 * @param[out] C pointer to output residual data ( (K*K)-by-N ).
 */
void TensorProductSurfaceResidual(const unsigned short  M,
                                  const unsigned short  N,
                                  const unsigned short  K,
                                  const unsigned short *IDX,
                                  const as3double      *A,
                                  as3double           **B,
                                  as3double           **C);
