%
% Solve the multigrid system
%
% x = MultiGridSolve(b, M)
% x = MultiGridSolve(b, M, cycle)
% x = MultiGridSolve(b, M, cycle, tol)
% x = MultiGridSolve(b, M, cycle, tol, maxit)
% x = MultiGridSolve(b, M, cycle, tol, maxit, x0)
% [x, flag] = MultiGridSolve(b, M, ...)
% [x, flag, relres] = MultiGridSolve(b, M, ...)
% [x, flag, relres, iter] = MultiGridSolve(b, M, ...)
% [x, flag, relres, iter, resvec] = MultiGridSolve(b, M, ...)
%
% INPUT        b:  rhs vector
%              M:  multigrid hierarchy
%          cycle:  number of iterations at each level,
%                  1 for V-cycle, 2 for W-cycle
%            tol:  relative tolerance
%          maxit:  maximum number of iterations
%             x0:  initial guess
%
% OUTPUT       x:  solution vector
%           flag:  convergence flag
%                   0: converged to within the prescribed tolerance.
%                   1: reached the maximum number of iterations but did not converge.
%                   2: ill-conditioned smoother.
%                   3: stagnated after two consecutive iterations were the same.
%         relres:  relative residual error
%           iter:  number of iterations
%         resvec:  vector of residual error at each iteration
%
% VERSION     1.0
% DATE        2020.8.27
% EMAIL       fyuan.nudt@foxmail.com
%

function [x, flag, relres, iter, resvec] = MultiGridSolve(b, M, cycle, tol, maxit, x0)
if nargin < 6 || isempty(x0); x0 = zeros(size(b)); end
if nargin < 5 || isempty(maxit); maxit = 500; end
if nargin < 4 || isempty(tol); tol = 1.0e-08; end
if nargin < 3 || isempty(cycle); cycle = 1; end
resvec(1,1) = norm(b - M.A{1}*x0);
relres = resvec(1,1)/norm(b);
if relres < tol
	x = x0;
	flag = 0;
	iter = 0;
	return
end
for iter = 1:maxit
	x = MultiGridCycle(b, 1, M.A, M.P, M.smooth, M.S, M.pre, M.post, cycle, x0);
	resvec(iter + 1,1) = norm(b - M.A{1}*x);
    relres = resvec(iter + 1,1)/norm(b);
    if relres < tol 
		flag = 0;
		return
    end
    if ~all(isfinite(x))
		flag = 2;
		return
    end
	if norm(x-x0) < eps*norm(x)
		flag = 3;
		return
	end
	x0 = x;
end
flag = 1;
end

function x = MultiGridCycle(b, l, A, P, smooth, S, pre, post, cycle, x)
if l == length(A)
    x = A{l}\b;
else
    for k = 1:pre
        x = smooth(b, A{l}, S{l}, x);
    end
    r = b - A{l}*x;
    g = P{l}'*r;
    e = zeros(size(g));
    for k = 1:cycle
        e = MultiGridCycle(g, l + 1, A, P, smooth, S, pre, post, cycle, e);
        if (l + 1 == length(A))
            break
        end
    end
    x = x + P{l}*e;
    for k = 1:post
        x = smooth(b, A{l}, S{l}, x);
    end
end
end
