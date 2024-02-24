clear

n = 100;
I = speye(n,n);
E = sparse(2:n,1:n-1,1,n,n);
D = 2*I - E - E';
A = kron(D,I) + kron(I,D);

b = ones(size(A,1),1);

opts.lmax = 20;           % maximum number of levels
opts.cmin = 100;          % minimum size of coarse level
opts.theta = 0.25;        % strength threshold
opts.interp = 'AGG';      % type of interpolation
                          % 'D1', 'D2' or 'AGG'
opts.smooth = 'ILU';      % type of smoother
                          % 'Jacobi', 'SOR', 'SSOR' or 'ILU'
opts.omega = 1.0;         % relaxation factor
opts.pre = 1;             % number of pre-smooth iterations
opts.post = 1;            % number of post-smooth iterations
cycle = 1;                % number of iterations at each level
                          % 1 for V-cycle, 2 for W-cycle
maxit = 100;
tol = 1.0e-10;

tic;
M = MultiGridSetup(A, opts);
tset = toc;

tic;
% [x, flag, relres, iter, resvec] = MultiGridSolve(b, M, cycle, tol, maxit);
[x, flag, relres, iter, resvec] = pcg(A, b, tol, maxit, @MultiGridSolve, [], [], M, cycle, eps, 1);
tsol = toc;
