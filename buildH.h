// contains useful functions to build the Hamiltonian

#include <slepceps.h>

#ifndef MY_HEADER_FILE_H
#define MY_HEADER_FILE_H

int isOffdiag(int a, int b);

PetscScalar uniformRnd();

PetscScalar gaussRnd(PetscScalar mu, PetscScalar sig);

PetscInt getNZoffDiagColVal(PetscInt *col, PetscScalar *val, PetscInt i, PetscInt n, PetscScalar J);

#endif /* MY_HEADER_FILE_H */
