// contains useful functions to build the Hamiltonian

#include <math.h>
#include <stdlib.h>
// #include <stdio.h>
#include <time.h>
#include <slepceps.h>

// #define PI = 4*atan(1)

int isOffdiag(int a, int b){

	int xor_result = a ^ b;  // XOR the two numbers
	return (xor_result != 0) && ((xor_result & (xor_result - 1)) == 0);

}

PetscScalar uniformRnd(){
    return (double)rand() / RAND_MAX;
}

PetscScalar gaussRnd(PetscScalar mu, PetscScalar sig){

    double u1, u2;
    const long double pi = acosl(-1.0L);

    // srand(time(0));
    u1 = uniformRnd();
    // srand(time(0) + 1);
    u2 = uniformRnd();

    // printf("u1 = %f, u2 = %f \n", u1, u2);

    return (sig*cos(2*pi*u2)*sqrt(-2.*log(u1)) + mu);
}

PetscInt getNZoffDiagColVal(PetscInt *col, PetscScalar *val, PetscInt i, PetscInt n, PetscScalar G){
    int c = 0;
    for(int j=i+1; j < n; j++){
	    if(isOffdiag(i, j)){
                col[c] = j;
                // val[c] = -G;
		val[c] = gaussRnd(0., sqrt(1 / (1 + pow( (j-i)/G, 2 ))));
                c++;
            }
    }
    return c;
}

