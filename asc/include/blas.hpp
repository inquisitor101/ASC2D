#pragma once

/*!
 * @file blas.hpp
 * @brief An interface to the Eigen library.
 */

#include "option_structure.hpp"
#include <../external/eigen/Eigen/Dense>


/*!
 * @brief Function that computes the inverse matrix of A.
 * 
 * @param[in] M dimension of the square matrix.
 * @param[in,out] A pointer to the square matrix.
 */
void ComputeInverseMatrix(const unsigned short  M,
                          as3double            *A);

/*!
 * @brief (overloaded) Function that computes the inverse matrix of A.
 * 
 * @param[in] M dimension of the square matrix.
 * @param[in] A pointer to the square matrix.
 * @param[out] B pointer to the inverse of the input matrix.
 */
void ComputeInverseMatrix(const unsigned short  M,
                          const as3double      *A,
                          as3double            *B);

/*!
 * @brief Function that computes a generalized matrix-vector multiplication.
 * 
 * @param[in] M number of rows of the matrix.
 * @param[in] K number of columns of the matrix or dimension of the vector.
 * @param[in] A pointer to the input matrix.
 * @param[in] x pointer to the input vector.
 * @param[out] y pointer to the output of the matrix-vector multiplication: A*x.
 */
void gemv(const unsigned short  M,
          const unsigned short  K,
          const as3double      *A,
          const as3double      *x,
          as3double            *y);


/*!
 * @brief Function that computes a generalized matrix-vector multiplication using lazy evaluation. Note, be careful of aliasing errors.
 * 
 * @param[in] M number of rows of the matrix.
 * @param[in] K number of columns of the matrix or dimension of the vector.
 * @param[in] A pointer to the input matrix.
 * @param[in] x pointer to the input vector.
 * @param[out] y pointer to the output of the matrix-vector multiplication: A*x.
 */
void lgemv(const unsigned short  M,
           const unsigned short  K,
           const as3double      *A,
           const as3double      *x,
           as3double            *y);

/*!
 * @brief (overloaded) Function that computes a scaled generalized matrix-vector multiplication using lazy evaluation. Note, be careful of aliasing errors.
 * 
 * @param[in] M number of rows of the matrix.
 * @param[in] K number of columns of the matrix or dimension of the vector.
 * @param[in] scalar scalar that is used to multiply the output.
 * @param[in] A pointer to the input matrix.
 * @param[in] x pointer to the input vector.
 * @param[out] y pointer to the output of the scaled matrix-vector multiplication: scalar*(A*x).
 */
void lgemv(const unsigned short  M,
           const unsigned short  K,
           const as3double       scalar,
           const as3double      *A,
           const as3double      *x,
           as3double            *y);

/*!
 * @brief Function that computes a generalized matrix-matrix multiplication using lazy evaluation. Note, be careful of aliasing errors.
 * 
 * @param[in] M number of rows of the matrix A.
 * @param[in] N number of columns of the matrix A or number of rows of matrix B.
 * @param[in] K number of columns of matrix B. 
 * @param[in] A pointer to the first input matrix.
 * @param[in] B pointer to the second input matrix.
 * @param[out] C pointer to the output of the matrix-matrix multiplication: A*B.
 */
void lgemm(const unsigned short  M,
           const unsigned short  N,
           const unsigned short  K,
           const as3double      *A,
           const as3double      *B,
           as3double            *C);

/*!
 * @brief (overloaded) Function that computes a scaled generalized matrix-matrix multiplication using lazy evaluation. Note, be careful of aliasing errors.
 * 
 * @param[in] M number of rows of the matrix A.
 * @param[in] N number of columns of the matrix A or number of rows of matrix B.
 * @param[in] K number of columns of matrix B.
 * @param[in] scalar scalar that is used to multiply the output.
 * @param[in] A pointer to the first input matrix.
 * @param[in] B pointer to the second input matrix.
 * @param[out] C pointer to the output of the scaled matrix-matrix multiplication: scalar*(A*B).
 */
void lgemm(const unsigned short  M,
           const unsigned short  N,
           const unsigned short  K,
           const as3double       scalar,
           const as3double      *A,
           const as3double      *B,
           as3double            *C);






