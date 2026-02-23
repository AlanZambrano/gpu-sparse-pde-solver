GPU-Accelerated Nonlinear PDE Solver

C/C++ · Eigen · MAGMA · CUDA · HPC (SGE)

------------------------------------------------------------------------

Overview

High-performance numerical solver for a coupled nonlinear
reaction–diffusion PDE system with:

-   Density-dependent diffusion
-   Predator–prey pursuit–avoidance dynamics
-   2D spatial discretization
-   Hybrid CPU/GPU acceleration

The nonlinear PDE system is discretized into large-scale sparse linear
systems and solved iteratively at each timestep.

The implementation targets Linux-based HPC environments with GPU
resources.

------------------------------------------------------------------------

Numerical Model

Model characteristics:

-   Coupled nonlinear PDE system (reaction–diffusion type)
-   Spatially explicit 2D discretization
-   Density-dependent diffusion
-   Iterative time-stepping scheme
-   Sparse matrix assembly at each timestep

For grid resolution N:

System dimension: dim = 2 × (N + 1)^2

With default parameters (N = 500), the solver handles approximately:

~5×10^5 degrees of freedom per timestep

Each iteration performs:

1.  Sparse matrix assembly
2.  RHS vector construction
3.  Iterative sparse linear solve
4.  State update

------------------------------------------------------------------------

Technical Stack

Language: - C/C++

Sparse Linear Algebra (CPU): - Eigen (SparseLU)

GPU Acceleration: - MAGMA 2.7.1 - Solver: PIDRMERGE - Preconditioner:
ILU (level 0) - Backend: CUDA 11.x

Build System: - Makefile-based compilation

Execution Environment: - Linux HPC Cluster

Scheduler: - SGE (qsub)

------------------------------------------------------------------------

Computational Features

-   Modular C++ architecture
-   Operator-based sparse matrix assembly
-   Iterative sparse linear solvers
-   ILU preconditioning
-   Hybrid CPU/GPU execution paths
-   GPU memory monitoring
-   Runtime profiling
-   CPU vs GPU benchmarking
-   Reproducible simulation workflows

------------------------------------------------------------------------

Performance Benchmark

Benchmarked across increasing grid sizes.

Results show:

-   Near-linear computational complexity (~O(N^0.99))
-   Up to ~6× GPU speedup vs CPU
-   Improved scalability for larger problem sizes
-   Stable hybrid execution across cluster queues

Benchmark figures and datasets are included in the repository.

------------------------------------------------------------------------

Build & Execution

1.  Load required modules (example):

module load magma-2.7.1 module load cuda-11.x

2.  Compile:

cd src make

3.  Submit to cluster:

make submit

The job script requests GPU resources via SGE and configures MAGMA and
CUDA runtime libraries accordingly.

------------------------------------------------------------------------

Repository Structure

include/ Header files src/ Core solver implementation scripts/ Cluster
job submission scripts figures/ Output visualizations

------------------------------------------------------------------------

Relevance

This project demonstrates:

-   Large-scale sparse system engineering
-   GPU-accelerated scientific computing
-   Hybrid CPU/GPU workflows
-   HPC cluster deployment and scheduling
-   Performance benchmarking and scalability analysis
-   Engineering-oriented numerical modeling

------------------------------------------------------------------------
