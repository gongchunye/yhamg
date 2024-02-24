# YHAMG
## һ �������
### 1.1 �������
YHAMG���ɹ����Ƽ���ѧ�з��Ĳ��д������������㷨�⣬��Ҫ�����������Ax = b��ϡ������ϵͳ������A��һ���ֲ�ʽ����ϡ�����x, b��ʵ������YHAMG��Ŀ�����ڴ��ģ���м������Ϊ�û��ṩ�����ܿ���չ�Ĳ����������������ⷽ����Ԥ�����Ӻ��ӿռ�������������������㷨��
- ϡ�����Դ����������㣻
- �ӿռ������������CG��GMRES��BiCGStab��
- ����������������㷨��AMG����
- Jacobi��SOR��SSOR������������
- ����ȫLU�ֽ⡣

### 1.2 ����ṹ
YHAMG��Ҫ��������ģ�飺
- ���Դ������ṩ����������/���Ӻ����Դ������㣻
- ������������ṩ�����ӿռ��������
- Ԥ�����ӣ��ṩ����Ԥ�����ӡ�

### 1.3 ����ص�
YHAMG����C++��д�ģ����������ص㣺
- ���������ƣ�֧��MPI��OpenMP��
- ���������Ĵ��룬�����������⣻ 
- ���õ��㷨����չ�ԺͲ��п���չ�ԣ������ڴ���ϡ������ϵͳ��⡣

## �� �����װ
### 2.1 ϵͳ����Ҫ��
Ӳ��Ҫ��
- �ڴ�4G���ϣ�CPU��Ƶ2G���ϡ�

�������Ҫ��
- ����ϵͳ��Ubuntu���汾16.04��
- MPI��������MPICH���汾3.2.1��
- C++��������GCC���汾7.1��

### 2.2 ��װYHAMG
1. ��ѹYHAMG������װĿ¼��
1. ִ��make��

### 2.3 ����YHAMG
�û���Ҫ�������ı�������������YHAMG�⡣
- ��ӿ�·���� -L*YHAMG_DIR*/lib��
- ��ӿ����ӣ� -lyhamg��
- ���ͷ�ļ�·���� -I*YHAMG_DIR*/include��
���� *YHAMG_DIR* ��ʾYHAMG��Ŀ¼·���������ڴ���ı�ͷ�а�����
```
#include <yhamg.h>
```

## �� Ӧ�ó���ӿ�
�ο��û��ֲᡣ

## �� ʾ��
YHAMG�ṩ��һ��ʾ��������testĿ¼�¡�����ǰ��������ϵͳ��������OMP_NUM_THREADS����OpenMP�߳�����
```
export OMP_NUM_THREADS=1
```

### 4.1 ����
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

### 4.2 �����б�
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