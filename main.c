#include <time.h>
#include "buildH.h"
#include <math.h>
#include <slepceps.h>

int main(int argc, char **argv){

    Mat A;
    // Vec D;
    EPS eps;
    EPS eps_si;
    Vec xr;
    ST st;
    KSP ksp;
    PC pc;
    // int fd;
    // PetscReal tol=1000*PETSC_MACHINE_EPSILON;;
    PetscErrorCode ierr;
    PetscInt L=10, n, *col, c, start, end, nev=10, nconv, itr=1, i;
    PetscInt off=0;
    PetscScalar *val, b = 1.0, d, a = 0., g=1, kr;
    PetscScalar emin, emax, sig, espec=0.;

    char fvecname[100], fvalname[100];  
    // char vecnum[20];
    // int fd;

    ierr = SlepcInitialize(&argc, &argv, NULL, NULL); if (ierr) return ierr;

    // Get the process rank and size
    PetscMPIInt petscRank, petscSize;
    MPI_Comm_rank(PETSC_COMM_WORLD, &petscRank);
    MPI_Comm_size(PETSC_COMM_WORLD, &petscSize);

    PetscOptionsGetInt(NULL,NULL,"-l",&L,NULL);
    PetscOptionsGetInt(NULL,NULL,"-nev",&nev,NULL);
    PetscOptionsGetInt(NULL,NULL,"-itr",&itr,NULL);
    PetscOptionsGetReal(NULL,NULL,"-g",&g,NULL);
    PetscOptionsGetReal(NULL,NULL,"-e",&espec,NULL);
    PetscOptionsGetReal(NULL,NULL,"-b",&b,NULL);
    PetscOptionsGetReal(NULL,NULL,"-a",&a,NULL);

    n = round(pow(2, L));
     
    col = (int*)malloc(n * sizeof(int));
    val = (double*)malloc(n * sizeof(double));

    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
    MatSetFromOptions(A);
    MatSetUp(A);
    MatGetOwnershipRange(A, &start, &end);

    srand(time(0));

    for(i = start; i <= end-1; i++){

        c = getNZoffDiagColVal(col, val, i, n, b);
	off += c;
        // dig[i] = getDiagNR(Delta, ind[i], L);
        MatSetValues(A, 1, &i, c, col, val, INSERT_VALUES);
        MatSetValues(A, c, col, 1, &i, val, INSERT_VALUES);
        
    }
    
    free(col);
    free(val);
    // PetscPrintf(PETSC_COMM_WORLD, "tot off diag = %d\n", off);
    // PetscPrintf(PETSC_COMM_WORLD, "a = %f\n", a);
    srand(time(0));

    for(i = start; i <= end-1; i++){

        d = gaussRnd(0., sqrt(L/2));
	// d = gaussRnd(0., sqrt((L/2) / (1 + pow(abs(i - (n/2)), a))));
	// d = gaussRnd(g*pow(i+1, a), sqrt(L/2));
        MatSetValue(A, i, i, d, INSERT_VALUES);
    
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    MatCreateVecs(A,NULL,&xr);

    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);		

    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, A, NULL);
    EPSSetProblemType(eps, EPS_HEP);
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    EPSSetFromOptions(eps);
    EPSSolve(eps);
    EPSGetConverged(eps, &nconv);

    EPSGetEigenvalue(eps, 0, &kr, NULL);
    emin = (double)(kr);
    PetscPrintf(PETSC_COMM_WORLD, "min = %f\n", emin);

    EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
    EPSSetFromOptions(eps);
    EPSSolve(eps);
    EPSGetConverged(eps, &nconv);

    EPSGetEigenvalue(eps, 0, &kr, NULL);
    emax = (double)(kr);
    PetscPrintf(PETSC_COMM_WORLD, "max = %f\n", emax);

    // sig = (emin + emax) / 2;
    // sig = espec*(emax - emin) + emin;
    sig = espec;

    // PetscPrintf(PETSC_COMM_WORLD, "sigma = %f\n", sig);
    EPSDestroy(&eps);

    EPSCreate(PETSC_COMM_WORLD, &eps_si);
    EPSSetOperators(eps_si, A, NULL);
    EPSSetProblemType(eps_si, EPS_HEP);
    EPSSetFromOptions(eps_si);
    EPSSetWhichEigenpairs(eps_si, EPS_TARGET_REAL);
    EPSSetTarget(eps_si, sig);

    EPSSetDimensions(eps_si, nev, PETSC_DEFAULT, PETSC_DEFAULT);
    EPSGetST(eps_si, &st);
    STSetType(st, STSINVERT);
    STGetKSP(st,&ksp);
    KSPSetType(ksp, KSPPREONLY);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);
    // PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
    PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

    EPSSolve(eps_si);
    EPSGetConverged(eps_si, &nconv);

    PetscPrintf(PETSC_COMM_WORLD, "#converged = %d\n", nconv);

    // snprintf(fvalname, 100, "./dataBinQREM_L=%d/eigval_L=%d_g=%0.3f_itr_%d_e=%0.1f.dat", L, L, G, itr, espec);
    // snprintf(fvecname, 100, "./dataBinQREM_L=%d/eigvec_L=%d_g=%0.3f_itr_%d_e=%0.1f.dat", L, L, G, itr, espec);

    snprintf(fvalname, 100, "./dataBinQREM_L=%d/eigval_L=%d_b=%0.1f_itr_%d_e=%0.1f.dat", L, L, b, itr, espec);
    snprintf(fvecname, 100, "./dataBinQREM_L=%d/eigvec_L=%d_b=%0.1f_itr_%d_e=%0.1f.dat", L, L, b, itr, espec);

    // snprintf(fvalname, 100, "./dataBinQREM_L=%d/eigval_L=%d_b=%0.1f_a=%0.1f_itr_%d_e=%0.1f.dat", L, L, b, a, itr, espec);
    // snprintf(fvecname, 100, "./dataBinQREM_L=%d/eigvec_L=%d_b=%0.1f_a=%0.1f_itr_%d_e=%0.1f.dat", L, L, b, a, itr, espec);

    // snprintf(fvalname, 100, "./dataBinQREM_L=%d/eigval_L=%d_b=%0.1f_a=%0.1f_g=%0.1f_itr_%d_e=%0.1f.dat", L, L, b, a, g, itr, espec);
    // snprintf(fvecname, 100, "./dataBinQREM_L=%d/eigvec_L=%d_b=%0.1f_a=%0.1f_g=%0.1f_itr_%d_e=%0.1f.dat", L, L, b, a, g, itr, espec);

    FILE* file1 = fopen(fvecname, "wb");
    FILE* file2 = fopen(fvalname, "wb");

    for (i=0; i<nconv; i++){
        
        EPSGetEigenpair(eps_si,i,&kr,PETSC_NULLPTR,xr,PETSC_NULLPTR);
        // PetscPrintf(PETSC_COMM_WORLD, "e = %f\n", kr);
        
        Vec Vec_local;
        if (petscSize==1) { VecCreateSeq(PETSC_COMM_SELF,n,&Vec_local); VecCopy(xr,Vec_local); }
        else
            {
                VecScatter ctx; VecScatterCreateToZero(xr,&ctx,&Vec_local);
                VecScatterBegin(ctx,xr,Vec_local,INSERT_VALUES,SCATTER_FORWARD); 
                VecScatterEnd(ctx,xr,Vec_local,INSERT_VALUES,SCATTER_FORWARD);
                VecScatterDestroy(&ctx);
            }

        if (petscRank==0)
            {
                fwrite(&kr, sizeof(double), 1, file2);
                PetscScalar * state; VecGetArray(Vec_local, &state);
                
                // fwrite(&state, sizeof(double), n, file);
                for (int i = 0; i < n; i++) {
                fwrite(&state[i], sizeof(double), 1, file1);
                }	
                
            }
        VecDestroy(&Vec_local);
    }

    fclose(file1);
    fclose(file2);
    
    EPSDestroy(&eps_si);
    VecDestroy(&xr);

    PetscPrintf(MPI_COMM_WORLD, "itr = %d\n", itr);
    
    MatDestroy(&A);
    SlepcFinalize();

    return 0;
}

