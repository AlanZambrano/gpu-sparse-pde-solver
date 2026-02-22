#include "common.h" // para usar los "CHECK"
#include <iostream>
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cuda.h>

//#define EIGEN_NO_CUDA // para quitar los warnings de Eigen
#include <Eigen/Sparse>



// g++ -I ./path/to/eigen/ main.cpp -o main -std=c++11
int main(int arg, char **argv){
	
	int N = 4;

	Eigen::SparseMatrix<double> Mat(N,N);
	
	Mat.insert(0, 0) = 1.0;
    Mat.insert(1, 1) = 1.0;
	Mat.insert(2, 2) = 1.0;
	Mat.insert(3, 3) = 1.0;
	Mat.insert(4, 4) = 1.0;
	Mat.makeCompressed();

return 0;
}
