#pragma once

#include "option_structure.hpp"
#include <../external/eigen/Eigen/Dense>



// Function that computes the inverse matrix of A.
void ComputeInverseMatrix(const unsigned short  M,
                          as3double            *A);

// Function (overloaded) that computes the inverse matrix of A.
void ComputeInverseMatrix(const unsigned short  M,
                          const as3double      *A,
                          as3double            *B);

// Function that computes a generalized matrix-vector multiplication.
void gemv(const unsigned short  M,
          const unsigned short  K,
          const as3double      *A,
          const as3double      *x,
          as3double            *y);

// Function that computes a generalized matrix-vector multiplication.
// Note, this uses a lazy evaluation. Careful of aliasing errors.
void lgemv(const unsigned short  M,
           const unsigned short  K,
           const as3double      *A,
           const as3double      *x,
           as3double            *y);

// Function that computes a generalized matrix-vector multiplication, multiplied by a scalar.
// Note, this uses a lazy evaluation. Careful of aliasing errors.
void lgemv(const unsigned short  M,
           const unsigned short  K,
           const as3double       scalar,
           const as3double      *A,
           const as3double      *x,
           as3double            *y);

// Function that computes a generalized matrix-matrix multiplication.
// Note, this uses a lazy evaluation.
void lgemm(const unsigned short  M,
           const unsigned short  N,
           const unsigned short  K,
           const as3double      *A,
           const as3double      *B,
           as3double            *C);

// Function that computes a generalized matrix-matrix multiplication, multiplied by scalar.
// Note, this uses a lazy evaluation.
void lgemm(const unsigned short  M,
           const unsigned short  N,
           const unsigned short  K,
           const as3double       scalar,
           const as3double      *A,
           const as3double      *B,
           as3double            *C);






