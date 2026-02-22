#include <iostream>
#include <fstream>
#include <Eigen/Sparse>
#include "../include/parametros.h"

//=================================================================
//=========================  2D   =================================

//## funciones de test ##


//(constante)
void constante2d(Eigen::VectorXd &pob, double &c){
	int I;
	for(int i=0; i<=N; i++){

		for(int j=0; j<=N; j++){
			I = i*(N+1)+j;
			pob(I) = c;	
		}
	}
return;
}



//(recta)
void recta2d(Eigen::VectorXd &pob){
	for(int i=0; i<=N; i++){
		
		for(int j=0; j<=N; j++){

			pob(i*(N+1)+j) = i*dx;// + j*dx;	
		}	
	}
return;
}


//(parabola)
void parabola2d(Eigen::VectorXd &pob){
	for(int i=0; i<=N; i++){

		for(int j=0; j<=N; j++){

			pob(i*(N+1)+j) = /*(i*dx)*(i*dx) +*/ (j*dx)*(j*dx);
		
		}

	}
return;
}

// poblacion 2D
void print_pob(Eigen::VectorXd &pob){

	for(int i=0; i<=N; i++){

		for(int j=0; j<=N; j++){

			std::cout << pob(i*(N+1)+j) << "  ";
		}

		std::cout<< "\n" << std::endl;
	}

return;
}

// condicion inicial sombrero 2D
void sombrero2d
(Eigen::VectorXd &pob, int i_cent, int j_cent, int radio, double altura){
	
	int imin, jmin, imax, jmax;
	imin = (i_cent - radio + N + 1)%(N+1);
	jmin = (j_cent - radio + N + 1)%(N+1);

	imax = (i_cent + radio )%(N+1);
	jmax = (j_cent + radio )%(N+1);

	bool i_in, j_in;
	for(int i=0; i <=N; i++){   

		i_in = (imin <= i) && (i <= imax);         

		for(int j=0; j<=N; j++){

			j_in = (jmin <= j) && (j <= jmax);
	
			pob(i*(N+1)+j) = (i_in && j_in) ? altura : pob(i*(N+1)+j); 
		}       
    }
return;
}


// ## derivadas ## //
// primera derivada respecto x (diferencia central) 2D
Eigen::VectorXd partialx_2d(const Eigen::VectorXd &pob){
	
	Eigen::VectorXd dpobdx((N+1)*(N+1));

	int inext, iprev;
    for(int i = 0; i <= N; i++){
		inext = (i+1)%(N+1);
		iprev = (i+N)%(N+1);

		for(int j = 0; j<= N; j++){

			dpobdx(i*(N+1)+j) = 
				(pob(inext*(N+1)+j)-pob(iprev*(N+1)+j))/(2.0*dx);
		}
    }				

return dpobdx;
}

// primera derivada respecto y (diferencia central) 2D
Eigen::VectorXd partialy_2d(const Eigen::VectorXd &pob){
	
	Eigen::VectorXd dpobdy((N+1)*(N+1));

	int jnext, jprev;
    for(int i = 0; i <= N; i++){

		for(int j = 0; j<= N; j++){
			jnext = (j+1)%(N+1);
			jprev = (j+N)%(N+1);		
		
			dpobdy(i*(N+1)+j) = 
				(pob(i*(N+1)+jnext) - pob(i*(N+1)+jprev))/(2.0*dx);			
		}
    }

return dpobdy;
}

// segunda derivada respecto x (diferencia central) 2D
Eigen::VectorXd partial2x_2d(const Eigen::VectorXd &pob){
	
	Eigen::VectorXd d2pobdx2((N+1)*(N+1));

	int inext, iprev;
    for(int i = 0; i <= N; i++){

		inext = (i+1)%(N+1);
		iprev = (i+N)%(N+1);

		for(int j = 0; j<= N; j++){

			d2pobdx2(i*(N+1)+j) = 
				(     pob(inext*(N+1)+j)

				 -2.0*pob(i*(N+1)+j) 

				 	+ pob(iprev*(N+1)+j)

				)/(dx*dx);
		}
    }				

return d2pobdx2;
}

// segunda derivada respecto y (diferencia central) 2D
Eigen::VectorXd partial2y_2d(const Eigen::VectorXd &pob){
	
	Eigen::VectorXd d2pobdy2((N+1)*(N+1));

	int jnext, jprev;
    for(int i = 0; i <= N; i++){

		for(int j = 0; j<= N; j++){

			jnext = (j+1)%(N+1);
			jprev = (j+N)%(N+1);		
		
			d2pobdy2(i*(N+1)+j) = 
				(	  pob(i*(N+1)+jnext) 

				 -2.0*pob(i*(N+1)+j) 

				    + pob(i*(N+1)+jprev)

				)/(dx*dx);			
		}
    }

return d2pobdy2;
}

// ## difusividades ##
// D(n)
Eigen::VectorXd D(const Eigen::VectorXd &n) {
	int size = n.size();
	Eigen::VectorXd D_(size);

	for(int i=0; i<size; i++){

		D_(i) = D0 - n(i) + n(i)*n(i);
	}	
	
return D_;	
}

Eigen::VectorXd D1(const Eigen::VectorXd &n) {
	int size = n.size();
	Eigen::VectorXd D1_(size);

	for(int i=0; i<size; i++){
		D1_(i) = - 1.0 + 2.0*n(i);
	}	
	
return D1_;	
}

// writeToFile(n,m,D)
void writeToFile(std::fstream &file, 
Eigen::VectorXd &n_2d, 
Eigen::VectorXd &m_2d,
Eigen::VectorXd &D_2d){
	
	file << "#i" << ",\t" 
		 << "j" << ",\t"
		 << "n(i,j)" << ",\t" 
		 << "m(i,j)" << ",\t" 
		 << "D(i,j)" << std::endl;	

	file << "\n"; // una segunda línea en blanco para cada llamada

	for(int i=0; i<=N; i++){

		for(int j=0; j<=N; j++){
			
			file << i << ",\t" 
				 << j << ",\t"
				 << n_2d(i*(N+1) + j) << ",\t" 
				 << m_2d(i*(N+1) + j) << ",\t" 
				 << D_2d(i*(N+1) + j) << std::endl;
			
		}
		file << "\n"; // una línea en blanco para cada cambio de i
	}

return;
}



















