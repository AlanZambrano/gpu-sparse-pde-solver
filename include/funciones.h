#pragma once
#include <Eigen/Sparse>
// funciones test 1D
void constante (Eigen::VectorXd &pob, double &c);
void recta     (Eigen::VectorXd &pob);
void parabola  (Eigen::VectorXd &pob);

// write to file, function overloading
void writeToFile(std::fstream &file, 
	Eigen::VectorXd &n);

void writeToFile(std::fstream &file, 
	Eigen::VectorXd &n, 
	Eigen::VectorXd &m);

void writeToFile(std::fstream &file, 
	Eigen::VectorXd &n_2d, 
	Eigen::VectorXd &m_2d,
	Eigen::VectorXd &D_2d);

// condiciones iniciales
void sombreroDoble(Eigen::VectorXd &pob, double alt);

// derivadas
Eigen::VectorXd partialx (const Eigen::VectorXd &pob);
Eigen::VectorXd partial2x(const Eigen::VectorXd &pob);

// difusividades
Eigen::VectorXd D (const Eigen::VectorXd &n);
Eigen::VectorXd D1(const Eigen::VectorXd &n);



////////////////////  2D  //////////////////////////

// funciones test
void constante2d(Eigen::VectorXd &pob, double &c);
void recta2d    (Eigen::VectorXd &pob);
void parabola2d (Eigen::VectorXd &pob);
/*
// write to file, function overloading
void writeToFile(std::fstream &file, 
	Eigen::VectorXd &n);

void writeToFile(std::fstream &file, 
	Eigen::VectorXd &n, 
	Eigen::VectorXd &m);

void writeToFile(std::fstream &file, 
	Eigen::VectorXd &n, 
	Eigen::VectorXd &m,
	Eigen::VectorXd &d);
*/

void print_pob(Eigen::VectorXd &pob);

// condiciones iniciales
void sombrero2d
(Eigen::VectorXd &pob, int i_cent, int j_cent, int radio, double altura);

// derivadas
Eigen::VectorXd partialx_2d (const Eigen::VectorXd &pob);
Eigen::VectorXd partialy_2d (const Eigen::VectorXd &pob);
Eigen::VectorXd partial2x_2d(const Eigen::VectorXd &pob);
Eigen::VectorXd partial2y_2d(const Eigen::VectorXd &pob);

/*
// difusividades
Eigen::VectorXd D (const Eigen::VectorXd &n);
Eigen::VectorXd D1(const Eigen::VectorXd &n);
*/

