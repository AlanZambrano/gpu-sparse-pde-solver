#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <cusparse_v2.h>

#include "magma_lapack.h"

#include "magma_v2.h"
#include "magmasparse.h"


// Esto resuelve Ax = B para matrices densas en doble precision

int main (int argc, char **argv){



// otro código de test (tomado de "Sparse Overview" de la doc de Magma):

// initialize MAGMA-sparse and create some LA objects
    magma_dopts dopts; //"dopts"?
    magma_queue_t queue;
    magma_queue_create(0, &queue);
    magma_d_matrix A, x, b, A_d, x_d, b_d;

    // pass linear system to MAGMA-sparse
    magma_dcsrset( m, n, row, col, val, &A, queue ); // {triada} (host) -> &A (device)
    magma_dvset( m, n, valb, &b, queue ); // valb (host) -> &b (device)
    
    // copy the linear system to the device
    magma_d_vtransfer( b, &b_d, Magma_CPU, Magma_DEV, queue ); // b -> b_d
    magma_d_mtransfer( A, &A_d, Magma_CPU, Magma_DEV, queue ); // A -> A_d   
    // allocate solution vector - on device
    magma_d_vinit( &x_d, Magma_DEV, A.num_cols, one, queue ); //"init"? "one"?
    
    // configure solver
    dopts.solver_par.solver = Magma_PGMRES; dopts.solver_par.restart = 30; 
    dopts.solver_par.rtol = 1e-10;
    magma_dsolverinfo_init( &dopts.solver_par, &dopts.precond_par, queue );
    
    // configure the preconditioner
    dopts.precond_par.solver = Magma_ILU; dopts.precond_par.levels = 0; 
    dopts.precond_par.trisolver = Magma_CUSOLVE;
    magma_d_precondsetup( A, b, &dopts.solver_par, &dopts.precond_par, queue );

    // solve the linear system
    magma_d_solver( B_d, b_d, &x_d, &dopts, queue );

    // copy the solution vector back to the host and pass it back to the application
    magma_d_mtransfer( x_d, &x, Magma_DEV, Magma_CPU, queue );
    magma_dvget( x, &m, &n, &valx, queue );
    
    // clean up the memory
    magma_dsolverinfo_free( &dopts.solver_par, &dopts.precond_par, queue );
    magma_d_mfree(&A_d, queue );
    magma_d_mfree(&A, queue );
    magma_d_mfree(&x_d, queue );
    magma_d_mfree(&b_d, queue );
    magma_d_mfree(&b, queue );

    // finalize MAGMA
    magma_queue_destroy( queue );
    magma_finalize();




return 0;
}


// All iterative solvers and eigensolvers included in the MAGMA-sparse package work on the device
magma_int_t magma_dfgmres (		magma_d_matrix	A_d,
		magma_d_matrix  		 b_d,
		magma_d_matrix 			*x_d,
		magma_d_solver_par  	*solver_par,
		magma_d_preconditioner  *precond_par,
		magma_queue_t  			queue 
	) 	


magma_zcsrset_gpu(
    magma_int_t m,
    magma_int_t n,
    magmaIndex_ptr row_d,
    magmaIndex_ptr col_d,
    magmaDoubleComplex_ptr val_d,
    magma_z_matrix *A, /* este es el A en MAGMA (donde vive?) */
    magma_queue_t queue );

//Similarly, matrices handled in MAGMA can be returned via the functions
// (pero qué retorno? la A de magma o la tríada? cuál es el input,
// cuál es el output y donde están alocados?)
magma_zcsrget(
    magma_z_matrix A,
    magma_int_t *m,
    magma_int_t *n,
    magma_index_t **row,
    magma_index_t **col,
    magmaDoubleComplex **val,
    magma_queue_t queue );
				// cuál es la diferencia entre ambos?
magma_zcsrget_gpu(
    magma_z_matrix A,
    magma_int_t *m,
    magma_int_t *n,
    magmaIndex_ptr *row,
    magmaIndex_ptr *col,
    magmaDoubleComplex_ptr *val,
    magma_queue_t queue );


//Vectors are handled as dense matrices (Magma_DENSE) and can be initialized inside MAGMA via
magma_z_vinit(
    magma_z_matrix *x,
    magma_location_t memory_location, /*Magma_CPU or Magma_DEV*/
    magma_int_t num_rows,
    magma_int_t num_cols, /*no debería ser siempre 1?*/
    magmaDoubleComplex values, /*por eso "_z_" ?*/
    magma_queue_t queue );

//or passed from/to main memory:
magma_zvset(
    magma_int_t m,
    magma_int_t n,
    magmaDoubleComplex *val,
    magma_z_matrix *b,
    magma_queue_t queue );

magma_zvget(
    magma_z_matrix b
    magma_int_t *m,
    magma_int_t *n,
    magmaDoubleComplex **val,
    magma_queue_t queue );



// "Copy" example (tomado de "magma-tutorial2.pdf)
int device;
magma_queue_t queue;
magma_getdevice(&device);
magma_queue_create(device, &queue);

// copu A, X -> dA, dX
magma_dsetmatrix( n, n,
				A, lda,
				dA, ldda, queue	);

magma_dsetmatrix( n, nrhs,
				A, lda,
				dA, ldda, queue	);

// ... solve AX = B

// copu result dX to X
magma_dgetmatrix( n, nrhs,
				  dX, lddx,
				  X, ldx, queue );

// print result

magma_queue_destroy ( queue );

// ... cleanup




//////////////////////
// ejemplo tomado de "/magma.2.7.1/example/example_sparse.cpp"

 // The software does e.g. discretization of a PDE,
    // ends up with a sparse linear system in CSR format and a RHS.
    // Let's assume this system is a diagonal system of size m.
    
    int i, m=200, n=1;
    int *col, *row;
    double *val, *rhs, *sol;
    
    row = (int*) calloc(m+1, sizeof(int));
    col = (int*) calloc(m,   sizeof(int));
    
    val = (double*) calloc(m, sizeof(double));
    rhs = (double*) calloc(m, sizeof(double));
    sol = (double*) calloc(m, sizeof(double));
    
    for (i = 0; i < m; ++i) {
        col[i] = i;
        row[i] = i;
        val[i] = 55.0;
        rhs[i] = 3.0;
        sol[i] = 0.0;
    }
    row[m] = m;
    
    // Initialize MAGMA and create some LA structures.
    magma_init();
    magma_dopts opts; // creo que dopts es como un descriptor double pero del solver
    magma_queue_t queue;
    magma_queue_create( 0, &queue );
    
    memset(&opts, 0, sizeof(magma_dopts));
    
    magma_d_matrix A={Magma_CSR}, dA={Magma_CSR};
    magma_d_matrix b={Magma_CSR}, db={Magma_CSR};
    magma_d_matrix x={Magma_CSR}, dx={Magma_CSR};
    
    // Pass the system to MAGMA.
    magma_dcsrset( m, m, row, col, val, &A, queue );
    magma_dvset( m, 1, rhs, &b, queue );
    
    // Choose a solver, preconditioner, etc. - see documentation for options.
    opts.solver_par.solver     = Magma_PIDRMERGE;
    opts.solver_par.restart    = 8;
    opts.solver_par.maxiter    = 1000;
    opts.solver_par.rtol       = 1e-10;
    opts.solver_par.maxiter    = 1000;
    opts.precond_par.solver    = Magma_ILU;
    opts.precond_par.levels    = 0;
    opts.precond_par.trisolver = Magma_CUSOLVE;
    
    // Initialize the solver.
    magma_dsolverinfo_init( &opts.solver_par, &opts.precond_par, queue );
    
    // Copy the system to the device
    magma_dmtransfer( A, &dA, Magma_CPU, Magma_DEV, queue );
    magma_dmtransfer( b, &db, Magma_CPU, Magma_DEV, queue );
    // initialize an initial guess for the iteration vector
    magma_dvinit( &dx, Magma_DEV, b.num_rows, b.num_cols, 0.0, queue );

    // Generate the preconditioner.
    magma_d_precondsetup( dA, db, &opts.solver_par, &opts.precond_par, queue );

    // In case we only wanted to generate a preconditioner, we are done.
    // The preconditioner in the opts.precond_par structure - in this case an ILU.
    // The lower ILU(0) factor is in opts.precond_par.L and
    // the upper ILU(0) factor is in opts.precond_par.U (in this case on the device).
    
    // If we want to solve the problem, run:
    magma_d_solver( dA, db, &dx, &opts, queue );
    
    // Then copy the solution back to the host...
    magma_dmtransfer( dx, &x, Magma_DEV, Magma_CPU, queue );
    // and back to the application code
    magma_dvcopy( x, &m, &n, sol, queue );
    
    // Free the allocated memory...
    magma_dmfree( &dx, queue );
    magma_dmfree( &db, queue ); 
    magma_dmfree( &dA, queue );
    magma_dmfree( &b, queue );  // won't do anything as MAGMA does not own the data. 
    magma_dmfree( &A, queue );  // won't do anything as MAGMA does not own the data. 
    // and finalize MAGMA.
    magma_queue_destroy( queue );
    magma_finalize();
    
    // From here on, the application code may continue with the solution in sol...
    for (i = 0; i < 20; ++i) {
        printf("%.4f\n", sol[i]);
    }
    free(val);
    free(col);
    free(row);
    free(sol);
    free(rhs);



