#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <cusparse_v2.h>
#include "magma_v2.h"
#include "magma_lapack.h"


// Esto resuelve Ax = B para matrices densas en doble precision

int main (int argc, char **argv){

//	std::cout << "hello world" << std::endl;


	magma_init ();
	magma_queue_t queue = NULL;
	
	magma_int_t dev =0;
	magma_queue_create (dev, &queue );
	real_Double_t gpu_time;

	magma_int_t * piv, info ; // piv - array of indices of inter -
	magma_int_t m = 8192; // changed rows ; a , d_a - mxm matrices
	magma_int_t n = 100; // b ,c , d_c - mxn matrices
	magma_int_t mm = m * m ; // size of a , d_a
	magma_int_t mn = m * n ; // size of b ,c , d_c

	double *a; // a - mxm matrix on the host
	double *b; // b - mxn matrix on the host
	double *c; // c - mxn matrix on the host
	double *d_a; // d_a - mxm matrix a on the device
	double *d_c; // d_c - mxn matrix c on the device
	magma_int_t ione = 1;
	magma_int_t ISEED [4] = { 0 ,0 ,0 ,1 }; // seed
	magma_int_t err;
	const double alpha = 1.0; // alpha =1
	const double beta = 0.0; // beta =0

	// allocate matrices
	err = magma_dmalloc_cpu (&a, mm ); // host memory for a
	err = magma_dmalloc_cpu (&b, mn ); // host memory for b
	err = magma_dmalloc_cpu (&c, mn ); // host memory for c
	err = magma_dmalloc (&d_a, mm ); // device memory for a
	err = magma_dmalloc (&d_c, mn ); // device memory for c
	piv =( magma_int_t *) malloc ( m * sizeof ( magma_int_t )); // host mem .

	// generate matrices // for piv
	// (LAPACK ESTÁ DANDO PROBLEMAS)
	//lapackf77_dlarnv(&ione, ISEED, &mm, a ); // randomize a
	//lapackf77_dlarnv(&ione, ISEED, &mn, b ); // randomize b


	// modificación de Alan
	for(magma_int_t i=0; i < mm; i++ ){
		a[i] = (double) rand()/RAND_MAX;
	}
	for(magma_int_t i=0; i < mn; i++ ){
		b[i] = (double) rand()/RAND_MAX;
	}
	// fin de modificacion de Alan


	printf( " expected solution :\ n " );
	magma_dprint( 4 , 4 , b , m ); // part of the expected solution

	// (ERROR: no reconoce "dgemm") right hand side c = a * b
	//blasf77_dgemm ( " N " ," N " ,&m ,& n ,& m ,& alpha ,a ,& m ,b ,& m ,& beta ,c ,& m );

	// (aquí copié el b -> d_c en lugar de c -> d_c)
	magma_dsetmatrix(m, m, a, m, d_a, m, queue); // copy a -> d_a
	//magma_dsetmatrix(m, n, c, m, d_c, m, queue); // copy c -> d_c
	magma_dsetmatrix(m, n, b, m, d_c, m, queue); // copy c -> d_c


	// MAGMA
	// solve the linear system d_a * x = d_c , d_a - mxm matrix ,
	// d_c - mxn matrix , d_c is overwritten by the solution ;
	// LU decomposition with partial pivoting and row
	// interchanges is used , row i is interchanged with row piv ( i )
	gpu_time = magma_sync_wtime ( NULL );

	magma_dgesv_gpu(m, n, d_a, m, piv, d_c, m, &info);
	gpu_time = magma_sync_wtime ( NULL ) - gpu_time;

	printf ("magma_dgesv_gpu time : %7.5 f sec .\ n ", gpu_time );
	magma_dgetmatrix (m, n, d_c, m, c, m, queue);
	printf ( " solution :\ n " );
	magma_dprint ( 4 , 4 , b , m ); // part of Magma solution
	free ( a ); // free host memory
	free ( b ); // free host memory
	free ( c ); // free host memory
	free ( piv ); // free host memory
	magma_free ( d_a ); // free device memory
	magma_free ( d_c ); // free device memory
	magma_queue_destroy ( queue ); // destroy queue


	magma_finalize (); // finalize Magma


return 0;
}

/*
magma_init (); // initialize Magma
	real_Double_t gpu_time;

	magma_int_t *piv, info;

	magma_int_t m = 8192; // a - mxm matrix
	magma_int_t n = 100; // c - mxn matrix
	magma_int_t mm = m * m; // size of a
	magma_int_t mn = m * n; // size of c

	double *a; // a - mxm matrix on the host
	double *b; // b - mxn matrix on the host
	double *c; // c - mxn matrix on the host

	magma_int_t ione = 1; // random uniform distr . in (0 ,1)
	magma_int_t ISEED [4] = { 0 ,0 ,0 ,1 }; // seed

	magma_int_t err;

	const double alpha = 1.0; // alpha =1
	const double beta  = 0.0; // beta =0

	// allocate matrices on the host
	err = magma_dmalloc_pinned (&a, mm); // host memory for a
	err = magma_dmalloc_pinned (&b, mn); // host memory for b
	err = magma_dmalloc_pinned (&c, mn); // host memory for c
	piv = (magma_int_t *) malloc (m * sizeof(magma_int_t ) );

	// generate random matrices a , b ;
	lapackf77_dlarnv (&ione, ISEED, &mm, a ); // randomize a
	lapackf77_dlarnv (&ione, ISEED, &mn, b ); // randomize b

	printf(" expected solution :\ n " );
	magma_dprint ( 4 , 4 , b , m );

	// right hand side c = a * b
	blasf77_dgemm ( " N " ," N " ,&m , &n, &n, &alpha, a, &m, b, &m, &beta, c, &m );

	// solve the linear system a * x = c
	// c - mxn matrix , a - mxm matrix ;
	// c is overwritten by the solution
	gpu_time = magma_sync_wtime ( NULL );

	magma_dgesv(m, n, a, m, piv, c, m, &info);

	gpu_time = magma_sync_wtime ( NULL ) - gpu_time ;
	printf ( " magma_dgesv time : %7.5 f sec .\ n " , gpu_time ); // time
	printf ( " solution :\ n " );
	magma_dprint ( 4 , 4 , c , m ); // part of the solution

	magma_free_pinned ( a ); // free host memory
	magma_free_pinned ( b ); // free host memory
	magma_free_pinned ( c ); // free host memory
	free ( piv ); // free host memory

	magma_finalize (); // finalize Magma
*/






// otro código de test:
/*
// initialize MAGMA-sparse and create some LA objects
    magma_dopts dopts;
    magma_queue_t queue;
    magma_queue_create( 0, &queue );
    magma_d_matrix A, A_d, x, x_d, b, b_d;

    // pass linear system to MAGMA-sparse
    magma_dcsrset( m, n, row, col, val, &A, queue );
    magma_dvset( m, n, valb, &b, queue );
    
    // copy the linear system to the device
    magma_d_vtransfer( b, &b_d, Magma_CPU, Magma_DEV, queue );    
    magma_d_mtransfer( A, &A_d, Magma_CPU, Magma_DEV, queue );    
    // allocate solution vector - on device
    magma_d_vinit( &x_d, Magma_DEV, A.num_cols, one, queue );
    
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
*/



