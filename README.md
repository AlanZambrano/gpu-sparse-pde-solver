
GPU-ACCELERATED NONLINEAR PDE SOLVER
C/C++ | Eigen | MAGMA | CUDA | HPC Cluster (SGE)

============================================================
1. OVERVIEW
============================================================

This project implements a high-performance numerical solver for a coupled 
nonlinear reaction–diffusion PDE system with density-dependent diffusion 
and pursuit–avoidance dynamics.

The nonlinear PDE system is discretized into large-scale sparse linear 
systems and solved iteratively at each timestep using hybrid CPU/GPU 
acceleration.

The implementation was designed for execution in a Linux-based HPC 
cluster environment using GPU resources.

============================================================
2. NUMERICAL MODEL
============================================================

• Coupled nonlinear PDE system (reaction–diffusion type)
• Spatially explicit 2D discretization
• Density-dependent diffusion
• Time-stepping iterative scheme
• Sparse matrix assembly at each timestep

For a grid resolution N:

    System dimension:
        dim = 2 * (N+1)^2

With default parameters (N = 500), the solver handles approximately:

    ~5×10^5 degrees of freedom per timestep.

At each iteration:
1. Sparse matrix assembly
2. RHS vector construction
3. Iterative sparse linear solve
4. State update

============================================================
3. TECHNICAL STACK
============================================================

Language:
    C/C++

Sparse Linear Algebra (CPU):
    Eigen (SparseLU)

GPU Acceleration:
    MAGMA 2.7.1
    Solver: PIDRMERGE
    Preconditioner: ILU (level 0)
    Backend: CUDA 11.x

Build System:
    Makefile-based compilation

Execution Environment:
    Linux HPC Cluster

Scheduler:
    SGE (qsub)

============================================================
4. COMPUTATIONAL FEATURES
============================================================

• Modular C++ architecture
• Operator-based sparse matrix assembly
• Iterative sparse linear solvers
• ILU preconditioning
• Hybrid CPU/GPU execution paths
• GPU memory monitoring
• Runtime profiling
• Benchmarking CPU vs GPU execution
• Reproducible simulation workflows

============================================================
5. PERFORMANCE BENCHMARK
============================================================

The solver was benchmarked across increasing grid sizes.

Results demonstrate:

• Near-linear computational complexity (~O(N^0.99))
• Up to ~6× GPU speedup compared to CPU
• Improved scalability for larger problem sizes
• Stable hybrid execution behavior across cluster queues

Benchmark figures and data are included in the repository.

============================================================
6. BUILD & EXECUTION
============================================================

1. Load required modules (example):

    module load magma-2.7.1
    module load cuda-11.x

2. Compile:

    cd src
    make

3. Submit to cluster:

    make submit

The job script requests GPU resources via SGE and configures 
MAGMA and CUDA runtime libraries accordingly.

============================================================
7. RELEVANCE
============================================================

This project demonstrates:

• Large-scale sparse system engineering
• GPU-accelerated scientific computing
• Hybrid CPU/GPU workflows
• HPC cluster deployment and scheduling
• Performance benchmarking and scalability analysis
• Engineering-oriented numerical modeling

============================================================
END OF README
============================================================
