# YHAMG
## 一 软件概述
### 1.1 软件功能
YHAMG是由国防科技大学袁璠等研发的并行代数多重网格算法库，主要用于求解形如Ax = b的稀疏线性系统，其中A是一个分布式大型稀疏矩阵，x, b是实向量。YHAMG的目的是在大规模并行计算机上为用户提供高性能可扩展的并行线性求解器。求解方案是预条件子和子空间迭代法，包含了以下算法：
- 稀疏线性代数基本运算；
- 子空间迭代法，包括CG、GMRES和BiCGStab；
- 经典代数多重网格算法（AMG）；
- Jacobi、SOR、SSOR基本迭代法；
- 不完全LU分解。

### 1.2 软件结构
YHAMG主要包括以下模块：
- 线性代数：提供向量，矩阵/算子和线性代数运算；
- 线性求解器：提供并行子空间迭代法；
- 预条件子：提供代数预条件子。

### 1.3 软件特点
YHAMG是用C++编写的，具有以下特点：
- 面向对象设计，支持MPI和OpenMP；
- 完整独立的代码，不依赖其它库； 
- 良好的算法可扩展性和并行可扩展性，适用于大型稀疏线性系统求解。

## 二 软件安装
### 2.1 系统配置要求
硬件要求：
- 内存4G以上，CPU主频2G以上。

软件环境要求：
- 操作系统：Ubuntu，版本16.04；
- MPI编译器：MPICH，版本3.2.1；
- C++编译器：GCC，版本7.1。

### 2.2 安装YHAMG
1. 解压YHAMG包到安装目录。
1. 执行make。

### 2.3 链接YHAMG
用户需要添加下面的编译命令来链接YHAMG库。
- 添加库路径： -L*YHAMG_DIR*/lib；
- 添加库链接： -lyhamg；
- 添加头文件路径： -I*YHAMG_DIR*/include；
其中 *YHAMG_DIR* 表示YHAMG主目录路径。并且在代码的标头中包含：
```
#include <yhamg.h>
```

## 三 应用程序接口
参考用户手册。

## 四 示例
YHAMG提供了一个示例程序，于test目录下。运行前，请设置系统环境变量OMP_NUM_THREADS控制OpenMP线程数。
```
export OMP_NUM_THREADS=1
```

### 4.1 运行
```
> mpirun -n 8 ./test -laplace -n 100 100 100 -P 2 2 2 -pre 4 -solver 0

# Test Problem
--------------------------------------------------
3D Laplace Problem on a Cube
Global Problem Size: (200, 200, 200)
Local Domain Size: (100, 100, 100)
MPI Process Size: (2, 2, 2)

# AMG Setup Phase
--------------------------------------------------
Multigrid Hierarchy:
Level	Rows	Nz
--------------------------------------------------
0	8000000	55760000
1	482611	18098598
2	72128	4481012
3	9449	670789
4	936	45482
5	130	4268
6	19	196
7	4	15
--------------------------------------------------
Number of Levels: 8
Grid Complexity: 1.07066
Operator Complexity: 1.41787
Setup Time: 2.76183

# CG Solve Phase
--------------------------------------------------
Iter	Residual
--------------------------------------------------
0	2828.43
1	16795.2
2	2702.6
3	491.845
4	89.0701
5	14.9227
6	2.49919
7	0.417837
8	0.0692848
9	0.0115655
10	0.00192553
11	0.000324108
12	5.44729e-05
13	9.07726e-06
--------------------------------------------------
Iterations: 13
Final Relative Residual: 3.2093e-09
Solve Time: 1.91141
Time/Iteration: 0.147032

```

### 4.2 参数列表
```
-laplace         3D Laplace problem on a cube(default)
-27pt            3D Laplace problem with 27-point stencil
-n <nx ny nz>    local domain size
-b <bx by bz>    nonzero block size
-P <Px Py Pz>    MPI process size

-pre <ID>        preconditioner ID
                 0(default): Jacobi
                 1: SOR
                 2: ILU(0)
				 3: Chebyshev
                 4: AMG

-solver <ID>     solver ID
                 0(default): CG
                 1: GMRES(10)
                 2: BiCGStab
```
