#include <iostream>
#include <Eigen/Sparse>
#include "../include/parametros.h"
#include "../include/funciones.h"
#include "../include/estado.h"

//definiciones
#include "estado_operadores.cpp"
#include "estado_ensamblaje.cpp"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
	
#include "magma_v2.h"
#include "magmasparse.h"

// declaraciones externas de variables Magma
extern magma_dopts opts;		// objeto de opciones del solver
extern magma_queue_t queue;		// cola

// Instancio los objetos del Magma (matriz y vectores)
extern magma_d_matrix A, A_d;
extern magma_d_matrix b, b_d;
extern magma_d_matrix x, x_d;

	// constructor
	Estado::Estado(
				   const Eigen::VectorXd &n_2d_0, const Eigen::VectorXd &m_2d_0
				  ){

		n_2d_ = n_2d_0;
		m_2d_ = m_2d_0;

		derivadas();
		difusividades(); 		
	}
	
	// destructor
	Estado::~Estado(){}

	// derivadas
	void Estado::derivadas(){
				
		// 2D
		dndx_2d_ = partialx_2d(n_2d_); d2ndx2_2d_ = partial2x_2d(n_2d_); 
		dmdx_2d_ = partialx_2d(m_2d_); d2mdx2_2d_ = partial2x_2d(m_2d_);

		dndy_2d_ = partialy_2d(n_2d_); d2ndy2_2d_ = partial2y_2d(n_2d_); 
		dmdy_2d_ = partialy_2d(m_2d_); d2mdy2_2d_ = partial2y_2d(m_2d_);

	}

	// difusividades
	void Estado::difusividades(){
		
		// 2D
		D_2d_  = D (n_2d_);
		D1_2d_ = D1(n_2d_);			
	}

	// actualización del sistema
	void Estado::actualizar(Eigen::VectorXd &du_2d){

		// 2D
		int I = 0;
		for(int i=0; i<=N; i++){
	
			for(int j=0; j<=N; j++){
				
				n_2d_(i*(N+1) + j) = n_2d_(i*(N+1) + j) + du_2d(I); I++;
				m_2d_(i*(N+1) + j) = m_2d_(i*(N+1) + j) + du_2d(I); I++;	
			}		
		}

		derivadas();
		difusividades();
	}

	// getters 2D
		// pobs
	Eigen::VectorXd Estado::get_n_2d(){return n_2d_;}
	Eigen::VectorXd Estado::get_m_2d(){return m_2d_;}
		// derivada x
	Eigen::VectorXd Estado::get_dndx_2d(){return dndx_2d_;}
	Eigen::VectorXd Estado::get_dmdx_2d(){return dmdx_2d_;}
		// derivada y
	Eigen::VectorXd Estado::get_dndy_2d(){return dndy_2d_;}
	Eigen::VectorXd Estado::get_dmdy_2d(){return dmdy_2d_;}
		// 2da derivada x
	Eigen::VectorXd Estado::get_d2ndx2_2d(){return d2ndx2_2d_;}
	Eigen::VectorXd Estado::get_d2mdx2_2d(){return d2mdx2_2d_;}
		// 2da derivada y
	Eigen::VectorXd Estado::get_d2ndy2_2d(){return d2ndy2_2d_;}
	Eigen::VectorXd Estado::get_d2mdy2_2d(){return d2mdy2_2d_;}
		// difusividades
	Eigen::VectorXd Estado::get_D_2d(){return D_2d_;}
	Eigen::VectorXd Estado::get_D1_2d(){return D1_2d_;}


	// du = eigenSolve(M, H)
	Eigen::VectorXd 
	Estado::eigenSolve(
		Eigen::SparseMatrix<double,Eigen::RowMajor> &M,
		Eigen::VectorXd &H
		){

	#ifdef MAGMA
		int m = 2*(N+1)*(N+1);// M_(m*m)		
		int n = m;		
		M.makeCompressed();

		// extraigo la triada sparse de Eigen en host
		int    *row = M.outerIndexPtr(); // acumulador
		int    *col = M.innerIndexPtr(); // índice de columna
		double *val = M.valuePtr();	 	 // valores/datos

		// el rhs y aloco donde estará la solución
		double *rhs = H.data();

		Eigen::VectorXd _du(m);
		double *sol = _du.data();

		// Pass the system to MAGMA.
	    magma_dcsrset( m, m, row, col, val, &A, queue ); // {row, col, val} -> A
	    magma_dvset( m, 1, rhs, &b, queue );			 // rhs -> b

		// Copy the system to the device
		magma_dmtransfer( A, &A_d, Magma_CPU, Magma_DEV, queue );
		magma_dmtransfer( b, &b_d, Magma_CPU, Magma_DEV, queue );

		// initialize an initial guess for the iteration vector
		magma_dvinit( &x_d, Magma_DEV, b.num_rows, b.num_cols, 0.0, queue );

		// Generate the preconditioner.
		magma_d_precondsetup( A_d, b_d, &opts.solver_par, &opts.precond_par, queue );

		// *** SOLVER ***
		magma_d_solver( A_d, b_d, &x_d, &opts, queue );

		// Then copy the solution back to the host...
		magma_dmtransfer( x_d, &x, Magma_DEV, Magma_CPU, queue );
		// and back to the application code
		magma_dvcopy( x, &m, &n, sol, queue );
		
		_du = Eigen::Map<Eigen::VectorXd>( sol, m);

		return _du;
	#else

		// declaro vector de incógnitas del sistema
		Eigen::VectorXd _du(H.size());

		// Calculo la inversa de M 
		solver.compute(M);
		
			// si la descomposición falla -> print error
			if(solver.info() != Eigen::Success) {
				std::cout << "decomposition failed" << std::endl;		
			}
		
		// Alimento con el rhs y retorno la solución
		_du = solver.solve(H);
			// si la solución falla -> print error
			if(solver.info() != Eigen::Success) {
				std::cout << "Solving  failed" << std::endl;		
			}
		return _du;	
	#endif

	}
