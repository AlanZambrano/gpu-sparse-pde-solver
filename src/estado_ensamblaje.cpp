//estado_ensamblaje.cpp. Definición de los métodos dentro del Class "Estado"
#include <vector>
#include <Eigen/Sparse>
#include "../include/estado.h"
#include "../include/parametros.h"

#include <stdio.h>
//------Ensamblaje------------//
//	dado el sistema       //
//		dndt = F      //
//		dmdt = G      //
//	construye el sistema  //
//		Mx = H        //
//----------------------------//

void Estado::setM(
				  Eigen::SparseMatrix<double,Eigen::RowMajor> &M_2d
				 ){
	
	typedef Eigen::Triplet<double> Tr;
	std::vector<Tr> coefs2;

	coefs2.reserve( 20*(N+1)*(N+1) ); // caso 2D
	

	int Mrows = 2*(N+1)*(N+1);
	int iprev, inext, jprev, jnext;
	int I,J;
	for(int i = 0; i<= N; i++){

		iprev = (N+i)%(N+1);
		inext = (i+1)%(N+1);

		for(int j = 0; j <= N; j++){
			I = i*(N+1) + j;
				
			jprev = (N+j)%(N+1);
			jnext = (j+1)%(N+1);	

						// even rows (sheeps) //
			// up
			J = (2*(I-2) + Mrows)%Mrows;
			coefs2.push_back(Tr(2*I, J  , -(dt/2.0)*dFdn(i, j, iprev, j) ));
			coefs2.push_back(Tr(2*I, J+1, -(dt/2.0)*dFdm(i, j, iprev, j) ));

			// left
			J = (2*(I-1) + Mrows)%Mrows;
			coefs2.push_back(Tr(2*I, J  , -(dt/2.0)*dFdn(i, j, i, jprev) ));
			coefs2.push_back(Tr(2*I, J+1, -(dt/2.0)*dFdm(i, j, i, jprev) ));

			// center
			J = 2*I;;
			coefs2.push_back(Tr(2*I, J  ,1-(dt/2.0)*dFdn(i, j, i, j) ));
			coefs2.push_back(Tr(2*I, J+1, -(dt/2.0)*dFdm(i, j, i, j) ));

			// right
			J = (2*(I+1))%Mrows;
			coefs2.push_back(Tr(2*I, J  , -(dt/2.0)*dFdn(i, j, i, jnext) ));
			coefs2.push_back(Tr(2*I, J+1, -(dt/2.0)*dFdm(i, j, i, jnext) ));

			// down
			J = (2*(I+2))%Mrows;
			coefs2.push_back(Tr(2*I, J  , -(dt/2.0)*dFdn(i, j, inext, j) ));
			coefs2.push_back(Tr(2*I, J+1, -(dt/2.0)*dFdm(i, j, inext, j) ));


			// odd rows (pumas) //
			// up
			J = (2*(I-2) + Mrows)%Mrows;
			coefs2.push_back(Tr(2*I+1, J  , -(dt/2.0)*dGdn(i, j, iprev, j) ));
			coefs2.push_back(Tr(2*I+1, J+1, -(dt/2.0)*dGdm(i, j, iprev, j) ));

			// left
			J = (2*(I-1) + Mrows)%Mrows;
			coefs2.push_back(Tr(2*I+1, J  , -(dt/2.0)*dGdn(i, j, i, jprev) ));
			coefs2.push_back(Tr(2*I+1, J+1, -(dt/2.0)*dGdm(i, j, i, jprev) ));

			// center 
			J = 2*I;
			coefs2.push_back(Tr(2*I+1, J  ,  -(dt/2.0)*dGdn(i, j, i, j) ));
			coefs2.push_back(Tr(2*I+1, J+1, 1-(dt/2.0)*dGdm(i, j, i, j) ));
	
			// right
			J = (2*(I+1))%Mrows;
			coefs2.push_back(Tr(2*I+1, J  , -(dt/2.0)*dGdn(i, j, i, jnext) ));
			coefs2.push_back(Tr(2*I+1, J+1, -(dt/2.0)*dGdm(i, j, i, jnext) ));

			// down
			J = (2*(I+2))%Mrows;
			coefs2.push_back(Tr(2*I+1, J  , -(dt/2.0)*dGdn(i, j, inext, j) ));
			coefs2.push_back(Tr(2*I+1, J+1, -(dt/2.0)*dGdm(i, j, inext, j) ));
		}
	}

	M_2d.setFromTriplets(coefs2.begin(), coefs2.end()); 
return;
}

// construye H a partir de F, G
void Estado::setrhs(Eigen::VectorXd &H_2d){

	int I = 0;
	for(int i=0; i<=N; i++){
	
		for(int j=0; j<=N; j++){
			H_2d(I) = dt * F(i,j);			
			I++;

			H_2d(I) = dt * G(i,j);
			I++;
		}		
	}

return;
}
