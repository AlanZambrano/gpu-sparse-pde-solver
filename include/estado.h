#pragma once

class Estado{
public:

	// Atributos 2D
	Eigen::VectorXd
	n_2d_, m_2d_, dndx_2d_, dmdx_2d_, d2ndx2_2d_, d2mdx2_2d_, D_2d_, D1_2d_, 					  dndy_2d_, dmdy_2d_, d2ndy2_2d_, d2mdy2_2d_;
	
	// Constructor	
	Estado(
		   const Eigen::VectorXd &n_2d_0, const Eigen::VectorXd &m_2d_0		
		  );
	
	// Destructor
	~Estado();

	// Métodos
	void derivadas();
	
	void difusividades();

	void actualizar(Eigen::VectorXd &du_2d);

	// getters 2D
		// pobs
	Eigen::VectorXd get_n_2d();
	Eigen::VectorXd get_m_2d();
		// derivada x
	Eigen::VectorXd get_dndx_2d();	 
	Eigen::VectorXd get_dmdx_2d();
		// derivada y
	Eigen::VectorXd get_dndy_2d();	 
	Eigen::VectorXd get_dmdy_2d();
		// 2da derivada x
	Eigen::VectorXd get_d2ndx2_2d(); 
	Eigen::VectorXd get_d2mdx2_2d();
		// 2da derivada y
	Eigen::VectorXd get_d2ndy2_2d(); 
	Eigen::VectorXd get_d2mdy2_2d();
		// difusividades
	Eigen::VectorXd get_D_2d();
	Eigen::VectorXd get_D1_2d();

	// operadores como métodos para que conozcan los atributos

	// 2D
	double F(int i, int j);
	double G(int i, int j);

	// elementos del Jacobiano 2D
	double dFdn(int i, int j, int l, int m);
	double dFdm(int i, int j, int l, int m);
	double dGdn(int i, int j, int l, int m);
	double dGdm(int i, int j, int l, int m);
	
	// Matriz del sistema a partir del vector de tripletes
	void setM(
	     	  Eigen::SparseMatrix<double,Eigen::RowMajor> &M_2d
			 );

	// rhs del sistema matricial (H)
	void setrhs(Eigen::VectorXd &H_2d);

	// solver de tipo SparseLU
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	
	// resolver sistema
	Eigen::VectorXd eigenSolve(
		Eigen::SparseMatrix<double,Eigen::RowMajor> &M_2d,
		Eigen::VectorXd &H_2d
		);

};
