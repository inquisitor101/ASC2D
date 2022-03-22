#include "blas.hpp"



void ComputeInverseMatrix
(
  const unsigned short  M,
  as3double            *A
)
 /*
  * Function that computes the inverse of a given matrix A.
  */
{
	// Cast A into an eigen matrix.
	Eigen::Map<Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matA(A, M, M);

  // Compute the inverse.
  matA = matA.inverse();
}


void ComputeInverseMatrix
(
  const unsigned short  M,
  const as3double      *A,
  as3double            *B
)
 /*
  * Function that computes the inverse of a given matrix A.
  */
{
  // Cast A into an eigen matrix.
  const Eigen::Map<const Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matA(A, M, M);
	// Cast B into an eigen matrix.
	Eigen::Map<Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matB(B, M, M);

  // Compute the inverse.
  matB = matA.inverse();
}


void gemv
(
  const unsigned short  M,
  const unsigned short  K,
  const as3double      *A,
  const as3double      *x,
  as3double            *y
)
 /*
  * Function that computes: y = A*x.
  */
{
  // Dimension of matrix A: M-by-K.
	const Eigen::Map<const Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matA(A, K, M);
  // Dimension of vector x: K-by-1.
	const Eigen::Map<const Eigen::Vector<as3double, Eigen::Dynamic> > vecX(x, K);
  // Dimension of vector y: M-by-1.
	Eigen::Map<Eigen::Vector<as3double, Eigen::Dynamic> > vecY(y, M);

  // Matrix-vector multiplication: y = A*x.
  vecY = vecX.transpose()*matA;
}


void lgemv
(
  const unsigned short  M,
  const unsigned short  K,
  const as3double      *A,
  const as3double      *x,
  as3double            *y
)
 /*
  * Function that computes: y = A*x using a lazy evaluation.
  */
{
  // Dimension of matrix A: M-by-K.
	const Eigen::Map<const Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matA(A, K, M);
  // Dimension of vector x: K-by-1.
	const Eigen::Map<const Eigen::Vector<as3double, Eigen::Dynamic> > vecX(x, K);
  // Dimension of vector y: M-by-1.
	Eigen::Map<Eigen::Vector<as3double, Eigen::Dynamic> > vecY(y, M);

  // Matrix-vector multiplication: y = A*x.
  vecY.noalias() = vecX.transpose()*matA;
}


void lgemv
(
  const unsigned short  M,
  const unsigned short  K,
  const as3double       scalar,
  const as3double      *A,
  const as3double      *x,
  as3double            *y
)
 /*
  * Function that computes: y = scalar*(A*x) using a lazy evaluation.
  */
{
  // Dimension of matrix A: M-by-K.
	const Eigen::Map<const Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matA(A, K, M);
  // Dimension of vector x: K-by-1.
	const Eigen::Map<const Eigen::Vector<as3double, Eigen::Dynamic> > vecX(x, K);
  // Dimension of vector y: M-by-1.
	Eigen::Map<Eigen::Vector<as3double, Eigen::Dynamic> > vecY(y, M);

  // Matrix-vector multiplication: y = A*x.
  vecY.noalias() = scalar*vecX.transpose()*matA;
}


void lgemm
(
  const unsigned short  M,
  const unsigned short  N,
  const unsigned short  K,
  const as3double      *A,
  const as3double      *B,
  as3double            *C
)
 /*
  * Function that computes: C = A*B using a lazy evaluation.
  */
{
  // Dimension of matrix A: M-by-K.
  const Eigen::Map<const Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matA(A, K, M);
  // Dimension of matrix B: K-by-N.
  const Eigen::Map<const Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matB(B, N, K);
  // Dimension of matrix C: M-by-N.
  Eigen::Map<Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matC(C, N, M);

  // Matrix-matrix multiplication: C = A*B.
  matC.noalias() = matB*matA;
}


void lgemm
(
  const unsigned short  M,
  const unsigned short  N,
  const unsigned short  K,
  const as3double       scalar,
  const as3double      *A,
  const as3double      *B,
  as3double            *C
)
 /*
  * Function that computes: C = scalar*(A*B) using a lazy evaluation.
  */
{
  // Dimension of matrix A: M-by-K.
  const Eigen::Map<const Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matA(A, K, M);
  // Dimension of matrix B: K-by-N.
  const Eigen::Map<const Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matB(B, N, K);
  // Dimension of matrix C: M-by-N.
  Eigen::Map<Eigen::Matrix<as3double, Eigen::Dynamic, Eigen::Dynamic> > matC(C, N, M);

  // Matrix-matrix multiplication: C = A*B.
  matC.noalias() = scalar*matB*matA;
}






