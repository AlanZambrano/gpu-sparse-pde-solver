#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Sparse>
#include <chrono>

#include "../include/funciones.h"
#include "../include/estado.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
	
#include "magma_v2.h"
#include "magmasparse.h"

#include <cuda_runtime.h>


// parámetros espacio-temporales
double L = 50.0;		// Espacio 
double dx = 0.1;

#ifndef TAMANIO
int N = (int)L/dx;
#else
int N = TAMANIO;
#endif


double Tmax = 5.0;	// Tiempo
double dt = dx*dx/4.0;
int T = 4000;// (int)(Tmax/dt + 1);

// parámetros difusivos
double rho = 0.5; 			   // herding parameter
double D0 = 1.0/(4.0*rho*rho); // D0 tal que ncrít = 0.5 
							   //

// parámetros de persecución
double e0   = 0.0;//10.0;	// espantativity
double chi0 = 0.0;//1.5;	// chemotaxis 
					// duda: ¢cómo hago para que chi0
					// sea más importante en la dinámica?

// parámetros depredativos 
double p  = 0.0;//1.0;	// predation rate
double mu = 0.0;//0.7;	// puma mortality

// parametros de cond. iniciales
double ancho  = N/4;
double altura = 0.1;

double radio  = 2;//N/2;
int i_cent = N/2, j_cent = N/2;

//////////////////////////////////
// definiciones variables Magma //
//////////////////////////////////
magma_dopts opts;		// objeto de opciones del solver
magma_queue_t queue;	// cola

// declaro global los objetos del Magma (matriz y vectores)

magma_d_matrix A={Magma_CSR}, A_d={Magma_CSR};
magma_d_matrix b={Magma_CSR}, b_d={Magma_CSR};
magma_d_matrix x={Magma_CSR}, x_d={Magma_CSR};


// para compilar:
// g++ -I ./path/to/eigen/ main.cpp -o main -std=c++11
int main(int arg, char **argv){
		// clock
    auto startTime = std::chrono::high_resolution_clock::now();

	// space-time grid
	//std::cout << "N = " << N << std::endl;
	//std::cout << "T = " << T << std::endl;
	int sets = 100; // cada cuantos pasos escribo en file

	// poblaciones (Oveja, Puma, difusividad)
	Eigen::VectorXd n_2d((N+1)*(N+1)), m_2d((N+1)*(N+1)), D_2d((N+1)*(N+1)); // 2D
	
	//cond. iniciales // 1D
		double c = 0.0;

		//constante2d(n2d, c);
		//recta2d(n2d);
		//parabola2d(n2d);
		sombrero2d(n_2d, i_cent, j_cent, radio, altura);

		//constante2d(n2d, altura);
		constante2d(m_2d, c);
		//recta2d(m2d);
		//parabola2d(m2d);		
		//sombrero2d(m2d, i_cent, j_cent, radio, altura);


	// Init del estado del sistema
	Estado estado(n_2d, m_2d);

	// file
	std::fstream file;
	//file.open("oveja_puma_difusividad.dat", std::ios::out);
	file.open("oveja_puma_difusividad_N_" + std::to_string(N) + ".dat", std::ios::out);
	
	// Matriz del sistema
	Eigen::SparseMatrix<double,Eigen::RowMajor> M_2d(2*(N+1)*(N+1),2*(N+1)*(N+1));
		
	// vector rhs del sistema y vector de incógnitas
	Eigen::VectorXd H_2d(2*n_2d.size()), du_2d(2*n_2d.size());


	#ifndef MAGMA
		std::cout << "CPU:" << std::endl;
		// instanciar el solver
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	#else
		//std::cout << "MAGMA: ON" << std::endl;
		// Initialize MAGMA and create some LA structures.
		magma_init(); 			// inicialización del Magma		
		magma_queue_create( 0, &queue );		
		memset(&opts, 0, sizeof(magma_dopts));
			
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

	#endif	
	
	int printeos_counter = -1;
	// main loop
	for(int t=0; t<=T; t++){	

		// printeo
		if( t%sets==0 ){
			//printf("%d/%d\n",t, 4000);
			printeos_counter++;
			//printf("\nn_2d(%d): \n", t);
			n_2d = estado.get_n_2d();
			m_2d = estado.get_m_2d();
			D_2d = estado.get_D_2d();
			//print_pob(n_2d);

			writeToFile(file, n_2d, m_2d, D_2d);

		}
		
		// ensamblaje
		estado.setM(M_2d);
		estado.setrhs(H_2d);
			
		// solución	
		du_2d = estado.eigenSolve(M_2d, H_2d);

		// actualización
		estado.actualizar(du_2d);
	
		// memory usage
		if(t==T){
			
			int device = 0;
			cudaSetDevice(device);

		    size_t free_bytes;
		    size_t total_bytes;	
			
			cudaMemGetInfo( &free_bytes, &total_bytes );

			//cudaError_t cuda_status ;
			//​cudaError_t cuda_status = cudaMemGetInfo( &free_bytes, &total_bytes );
			
			
			double free_db  = (double)free_bytes ;
		    double total_db = (double)total_bytes ;
		    double used_db = total_db - free_db ;

			printf("Size, UsedMem(MiB), FreeMem(MiB), TotMem(MiB), TotTime(mseg):\n");
			printf("N = %d, %.2f, %.2f, %.2f, ", 
				N, used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

		}

	}
	file.close();

#ifdef MAGMA
	// Free the allocated memory...
		magma_dmfree( &x_d, queue );
		magma_dmfree( &b_d, queue ); 
		magma_dmfree( &A_d, queue );
		magma_dmfree( &b,   queue );  // won't do anything as MAGMA does not own the data. 
		magma_dmfree( &A,   queue );  // won't do anything as MAGMA does not own the data. 
		
	// and finalize MAGMA.
		magma_queue_destroy( queue );
		magma_finalize();
#endif

    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
	

	//std::cout << "N = " << N << ", " << duration.count() << " mseg"  << std::endl;

		std::cout << duration.count() << std::endl;
	std::cout << "Nº pasos de tiempo escritos = " << printeos_counter << std::endl;
return 0;
}











