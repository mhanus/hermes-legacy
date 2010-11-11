// This file is part of Hermes
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#include "umfpack_solver.h"

#ifdef WITH_UMFPACK
  extern "C" {
    #include <umfpack.h>
  }
#endif

#include "../trace.h"
#include "../error.h"
#include "../utils.h"
#include "../callstack.h"

static int find_position(int *Ai, int Alen, int idx) {
  _F_
  assert (idx >= 0);
  
  register int lo = 0, hi = Alen - 1, mid;
  
  while (1) 
  {
    mid = (lo + hi) >> 1;
    
    if (idx < Ai[mid]) hi = mid - 1;
    else if (idx > Ai[mid]) lo = mid + 1;
    else break;
    
    // Sparse matrix entry not found (raise an error when trying to add 
    // value to this position, return 0 when obtaining value there).
    if (lo > hi) 
    {
      mid = -1;
      break;
    }
  }
  return mid;
}

UMFPackMatrix::UMFPackMatrix() {
  _F_
  size = 0; nnz = 0;
  Ap = NULL;
  Ai = NULL;
  Ax = NULL;
}

UMFPackMatrix::UMFPackMatrix(int size) {
  _F_
  this->size = size;
  this->alloc();
}

UMFPackMatrix::~UMFPackMatrix() {
  _F_
  free();
}

void UMFPackMatrix::alloc() {
  _F_
  assert(pages != NULL);
  assert(size > 0);

  // initialize the arrays Ap and Ai
  Ap = new int [size + 1];
  MEM_CHECK(Ap);
  int aisize = get_num_indices();
  Ai = new int [aisize];
  MEM_CHECK(Ai);

  // sort the indices and remove duplicities, insert into Ai
  int i, pos = 0;
  for (i = 0; i < size; i++) {
    Ap[i] = pos;
    pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
  }
  Ap[i] = pos;

  delete [] pages;
  pages = NULL;

  nnz = Ap[size];
  
  Ax = new scalar [nnz];
  MEM_CHECK(Ax);
  memset(Ax, 0, sizeof(scalar) * nnz);
}

void UMFPackMatrix::free() {
  _F_
  nnz = 0;
  delete [] Ap; Ap = NULL;
  delete [] Ai; Ai = NULL;
  delete [] Ax; Ax = NULL;
}

scalar UMFPackMatrix::get(int m, int n)
{
  _F_
  // Find m-th row in the n-th column.
  int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);

  if (mid < 0) // if the entry has not been found
    return 0.0;   
  else 
    return Ax[Ap[n]+mid];
}

void UMFPackMatrix::zero() {
  _F_
  memset(Ax, 0, sizeof(scalar) * nnz);
}

void UMFPackMatrix::add(int m, int n, scalar v) {
  _F_
  if (v != 0.0 && m >= 0 && n >= 0)   // ignore dirichlet DOFs
  {
    // Find m-th row in the n-th column.
    int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
    // Make sure we are adding to an existing non-zero entry.
    if (pos < 0) 
      error("Sparse matrix entry not found");
    
    Ax[Ap[n]+pos] += v;
  }
}

void UMFPackMatrix::add(int m, int n, scalar **mat, int *rows, int *cols) {
  _F_
  for (int i = 0; i < m; i++)       // rows
    for (int j = 0; j < n; j++)     // cols
      add(rows[i], cols[j], mat[i][j]);
}

/// dumping matrix and right-hand side
///
bool UMFPackMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
  _F_
  switch (fmt) 
  {
    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", size, size, nnz, nnz);
      for (int j = 0; j < size; j++)
        for (int i = Ap[j]; i < Ap[j + 1]; i++)
          fprintf(file, "%d %d " SCALAR_FMT "\n", Ai[i] + 1, j + 1, SCALAR(Ax[i]));
      fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

      return true;

    case DF_HERMES_BIN: 
    {
      hermes_fwrite("H3DX\001\000\000\000", 1, 8, file);
      int ssize = sizeof(scalar);
      hermes_fwrite(&ssize, sizeof(int), 1, file);
      hermes_fwrite(&size, sizeof(int), 1, file);
      hermes_fwrite(&nnz, sizeof(int), 1, file);
      hermes_fwrite(Ap, sizeof(int), size + 1, file);
      hermes_fwrite(Ai, sizeof(int), nnz, file);
      hermes_fwrite(Ax, sizeof(scalar), nnz, file);
      return true;
    }

    case DF_PLAIN_ASCII:
      EXIT(HERMES_ERR_NOT_IMPLEMENTED);
      return false;

    default:
      return false;
  }
}

int UMFPackMatrix::get_matrix_size() const {
  _F_
  assert(Ap != NULL);
  /*          Ai             Ax                      Ap                     nnz     */    
  return (sizeof(int) + sizeof(scalar)) * nnz + sizeof(int)*(size+1) + sizeof(int);
}

double UMFPackMatrix::get_fill_in() const {
  _F_
  return nnz / (double) (size * size);
}


// UMFPackVector ///////

UMFPackVector::UMFPackVector() {
  _F_
  v = NULL;
  size = 0;
}

UMFPackVector::UMFPackVector(int size) {
  _F_
  this->size = size;
        this->alloc(size);
}

UMFPackVector::~UMFPackVector() {
  _F_
  free();
}

void UMFPackVector::alloc(int n) {
  _F_
  free();
  this->size = n;
  v = new scalar [n];
  MEM_CHECK(v);
  this->zero();
}

void UMFPackVector::zero() {
  _F_
  memset(v, 0, size * sizeof(scalar));
}

void UMFPackVector::free() {
  _F_
  delete [] v;
  v = NULL;
  size = 0;
}

void UMFPackVector::set(int idx, scalar y) {
  _F_
  if (idx >= 0) v[idx] = y;
}

void UMFPackVector::add(int idx, scalar y) {
  _F_
  if (idx >= 0) v[idx] += y;
}

void UMFPackVector::add(int n, int *idx, scalar *y) {
  _F_
  for (int i = 0; i < n; i++)
    if (idx[i] >= 0) v[idx[i]] += y[i];
}

bool UMFPackVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
  _F_
  switch (fmt) 
  {
    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx1\n%s = [\n", size, var_name);
      for (int i = 0; i < size; i++)
        fprintf(file, SCALAR_FMT "\n", SCALAR(v[i]));
      fprintf(file, " ];\n");
      return true;

    case DF_HERMES_BIN: 
    {
      hermes_fwrite("H3DR\001\000\000\000", 1, 8, file);
      int ssize = sizeof(scalar);
      hermes_fwrite(&ssize, sizeof(int), 1, file);
      hermes_fwrite(&size, sizeof(int), 1, file);
      hermes_fwrite(v, sizeof(scalar), size, file);
      return true;
    }

    case DF_PLAIN_ASCII:
      EXIT(HERMES_ERR_NOT_IMPLEMENTED);
      return false;

    default:
      return false;
  }
}

// UMFPack solver //////

#if !defined (H2D_COMPLEX) && !defined (H3D_COMPLEX)
  // real case
  #define umfpack_symbolic(m, n, Ap, Ai, Ax, S, C, I)   umfpack_di_symbolic(m, n, Ap, Ai, Ax, S, C, I)
  #define umfpack_numeric(Ap, Ai, Ax, S, N, C, I)       umfpack_di_numeric(Ap, Ai, Ax, S, N, C, I)
  #define umfpack_solve(sys, Ap, Ai, Ax, X, B, N, C, I) umfpack_di_solve(sys, Ap, Ai, Ax, X, B, N, C, I)
  #define umfpack_free_symbolic                         umfpack_di_free_symbolic
  #define umfpack_free_numeric                          umfpack_di_free_numeric
  #define umfpack_defaults                              umfpack_di_defaults
#else
  // macros for calling complex UMFPACK in packed-complex mode
  #define umfpack_symbolic(m, n, Ap, Ai, Ax, S, C, I)   umfpack_zi_symbolic(m, n, Ap, Ai, (double *) (Ax), NULL, S, C, I)
  #define umfpack_numeric(Ap, Ai, Ax, S, N, C, I)       umfpack_zi_numeric(Ap, Ai, (double *) (Ax), NULL, S, N, C, I)
  #define umfpack_solve(sys, Ap, Ai, Ax, X, B, N, C, I) umfpack_zi_solve(sys, Ap, Ai, (double *) (Ax), NULL, (double *) (X), NULL, (double *) (B), NULL, N, C, I)
  #define umfpack_free_symbolic                         umfpack_di_free_symbolic
  #define umfpack_free_numeric                          umfpack_zi_free_numeric
  #define umfpack_defaults                              umfpack_zi_defaults
#endif


UMFPackLinearSolver::UMFPackLinearSolver(UMFPackMatrix *m, UMFPackVector *rhs)
  : LinearSolver(HERMES_FACTORIZE_FROM_SCRATCH), m(m), rhs(rhs), symbolic(NULL), numeric(NULL)
{
  _F_
#ifdef WITH_UMFPACK
#else
  error(UMFPACK_NOT_COMPILED);
#endif
}


UMFPackLinearSolver::~UMFPackLinearSolver() {
  _F_
  free_factorization_data();
}

#ifdef WITH_UMFPACK

static void check_status(const char *fn_name, int status) {
  _F_
  switch (status) {
    case UMFPACK_OK: break;
    case UMFPACK_WARNING_singular_matrix:       warning("%s: singular matrix!", fn_name); break;
    case UMFPACK_ERROR_out_of_memory:           warning("%s: out of memory!", fn_name); break;
    case UMFPACK_ERROR_argument_missing:        warning("%s: argument missing", fn_name); break;
    case UMFPACK_ERROR_invalid_Symbolic_object: warning("%s: invalid Symbolic object", fn_name); break;
    case UMFPACK_ERROR_invalid_Numeric_object:  warning("%s: invalid Numeric object", fn_name); break;
    case UMFPACK_ERROR_different_pattern:       warning("%s: different pattern", fn_name); break;
    case UMFPACK_ERROR_invalid_system:          warning("%s: invalid system", fn_name); break;
    case UMFPACK_ERROR_n_nonpositive:           warning("%s: n nonpositive", fn_name); break;
    case UMFPACK_ERROR_invalid_matrix:          warning("%s: invalid matrix", fn_name); break;
    case UMFPACK_ERROR_internal_error:          warning("%s: internal error", fn_name); break;
    default:                                    warning("%s: unknown error (%d)", fn_name, status); break;
  }
}

#endif

bool UMFPackLinearSolver::solve() {
  _F_
#ifdef WITH_UMFPACK
  assert(m != NULL);
  assert(rhs != NULL);

  assert(m->size == rhs->size);

  TimePeriod tmr;

  int status;

  if ( !setup_factorization() )
  {
    warning("LU factorization could not be completed.");
    return false;
  }

  if(sln)
    delete [] sln;
  sln = new scalar[m->size];
  MEM_CHECK(sln);
  memset(sln, 0, m->size * sizeof(scalar));

  status = umfpack_solve(UMFPACK_A, m->Ap, m->Ai, m->Ax, sln, rhs->v, numeric, NULL, NULL);
  if (status != UMFPACK_OK) {
    check_status("umfpack_di_solve", status);
    return false;
  }

  tmr.tick();
  time = tmr.accumulated();
  
  return true;
#else
  return false;
#endif
}

bool UMFPackLinearSolver::setup_factorization()
{
  _F_
#ifdef WITH_UMFPACK
  int status;
  switch(factorization_scheme)
  {
    case HERMES_FACTORIZE_FROM_SCRATCH:
      if (symbolic != NULL) umfpack_free_symbolic(&symbolic);
      
      //debug_log("Factorizing symbolically.");
      status = umfpack_symbolic(m->size, m->size, m->Ap, m->Ai, m->Ax, &symbolic, NULL, NULL);
      if (status != UMFPACK_OK) {
        check_status("umfpack_di_symbolic", status);
        return false;
      }
      if (symbolic == NULL) EXIT("umfpack_di_symbolic error: symbolic == NULL");
      
    case HERMES_REUSE_MATRIX_REORDERING:
    case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING:
      if (numeric != NULL) umfpack_free_numeric(&numeric);
      
      //debug_log("Factorizing numerically.");
      status = umfpack_numeric(m->Ap, m->Ai, m->Ax, symbolic, &numeric, NULL, NULL);
      if (status != UMFPACK_OK) {
        check_status("umfpack_di_numeric", status);
        return false;
      }
      if (numeric == NULL) EXIT("umfpack_di_numeric error: numeric == NULL");
  }
  
  return true;
#else
  return false;
#endif
}

void UMFPackLinearSolver::free_factorization_data()
{ 
  _F_
#ifdef WITH_UMFPACK
  if (symbolic != NULL) umfpack_free_symbolic(&symbolic);
  symbolic = NULL;
  if (numeric != NULL) umfpack_free_numeric(&numeric);
  numeric = NULL;
#endif
}


