%
% Construct a multigrid hierarchy
%
% M = MultiGridSetup(A)
% M = MultiGridSetup(A, opts)
%
% INPUT        A:  input matrix
%           opts:  structure with up to eight fields below
%                   lmax   -  maximum number of levels
%                   cmin   -  minimum size of coarse level
%                   theta  -  strength threshold
%                   interp -  type of interpolation
%                             'D1', 'D2' or 'AGG'
%                   smooth -  type of smoother,
%                             'Jacobi', 'SOR', 'SSOR' or 'ILU'
%                   omega  -  relaxation factor
%                   pre    -  number of pre-smooth iterations
%                   post   -  number of post-smooth iterations
%
% OUTPUT       M:  multigrid hierarchy
%
% VERSION     1.0
% DATE        2020.8.27
% EMAIL       fyuan.nudt@foxmail.com
%

function M = MultiGridSetup(A, opts)
if nargin == 1 || isempty(opts); opts = struct(); end
if ~isfield(opts,'lmax'); opts.lmax = 20; end
if ~isfield(opts,'cmin'); opts.cmin = 100; end
if ~isfield(opts,'theta'); opts.theta = 0.25; end
if ~isfield(opts,'interp'); opts.interp = 'AGG'; end
if ~isfield(opts,'smooth'); opts.smooth = 'SSOR'; end
if ~isfield(opts,'omega'); opts.omega = 1.0; end
if ~isfield(opts,'pre'); opts.pre = 1; end
if ~isfield(opts,'post'); opts.post = 1; end
M.A = {A};
M.P = {};
M.S = {};
M.pre = opts.pre;
M.post = opts.post;
switch opts.smooth
    case 'Jacobi'
        M.smooth = @(b, A, S, x0) x0 + S.d.*(b - A*x0);
    case 'SOR'
        M.smooth = @(b, A, S, x0) S.M\(b + S.N*x0);
	case {'SSOR', 'ILU'}
        M.smooth = @(b, A, S, x0) x0 + S.U\(S.L\(b - A*x0));
end
for l = 1:opts.lmax-1
    n = size(M.A{l}, 1);
	switch opts.smooth
		case 'Jacobi'
			M.S{l}.d = opts.omega./spdiags(M.A{l},0);
		case 'SOR'
			D = spdiags(spdiags(M.A{l},0),0,n,n);
			M.S{l}.M = D/opts.omega + tril(M.A{l},-1);
			M.S{l}.N = (1/opts.omega - 1)*D - triu(M.A{l},1);
		case 'SSOR'
			D = spdiags(spdiags(M.A{l},0),0,n,n);
			M.S{l}.L = (D + opts.omega*tril(M.A{l},-1))/D;
			M.S{l}.U = (D/opts.omega + triu(M.A{l},1))/(2-opts.omega);
		case 'ILU'
			[M.S{l}.L,M.S{l}.U] = ilu(M.A{l});
	end
    s = MultiGridStrength(M.A{l}, opts.theta);
    [c, f] = MultiGridCoarsening(s);
    M.P{l} = MultiGridInterpolation(M.A{l}, s, c, f, opts.interp);
    M.A{l+1} = M.P{l}'*M.A{l}*M.P{l};
    if length(c) <= opts.cmin; break; end
end
end

function s = MultiGridStrength(A, theta)
n = size(A,1);
B = tril(A,-1) + triu(A,1);
[i,j,x] = find(B);
b = theta*full(max(abs(B),[],2));
x = abs(x);
k = x > b(i);
s = sparse(i(k),j(k),x(k),n,n);
end

function [c, f] = MultiGridCoarsening(s)
s = spones(s);
n = size(s,1);
w = full(sum(s,1))' + rand(n,1);
[w,v] = sort(w,'descend');
c = [];
k = w < 1;
f = v(k);
v(k) = [];
s = s(v,v);
while ~isempty(v)
    p = ~full(any(tril(s),2));
    q = ~full(any(triu(s),1))';
    d = p & q;
    e = full(any(s(:,d),2));
    fnew = v(e); 
    cnew = v(d);
    c = [c;cnew];
    f = [f;fnew];
    k = d | e;
    v(k) = [];
    s(k,:) = [];
    s(:,k) = [];
end
c = sort(c);
f = sort(f);
end

function P = MultiGridInterpolation(A, s, c, f, interp)
switch interp
    case 'D1'
        P = MultiGridInterpolationDistance1(A, s, c, f);
    case 'D2'
        P = MultiGridInterpolationDistance2(A, s, c, f);
    case 'AGG'
        P = MultiGridSmoothedAggregation(A, s, c, f);
end
end

function P = MultiGridInterpolationDistance1(A, s, c, f)
n = size(A,1);
m = length(c);
b = 1 - full(sum(A(f,:),2))./diag(A(f,f));
z = full(sum(s(f,c),2));
t = z ~= 0;
h = f(t);
z = b(t)./z(t);
[i,j,x] = find(s(h,c));
P = sparse(c,1:m,1,n,m) + sparse(h(i),j,x.*z(i),n,m);
end

function P = MultiGridInterpolationDistance2(A, s, c, f)
n = size(A,1);
m = length(c);
b = 1 - full(sum(A(f,:),2))./diag(A(f,f));
z = full(sum(s(f,c),2));
t = z ~= 0;
h = f(t);
z = b(t)./z(t);
[i,j,x] = find(s(h,c));
P = sparse(c,1:m,1,n,m) + sparse(h(i),j,x.*z(i),n,m);
s(:,f(~t)) = 0;
z = full(sum(s(f,:),2));
t = z ~= 0; 
h = f(t);
z = b(t)./z(t);
[i,j,x] = find(s(h,:));
Q = sparse(c,c,1,n,n) + sparse(h(i),j,x.*z(i),n,n);
P = Q * P;
end

function P = MultiGridSmoothedAggregation(A, s, c, f)
n = size(A,1);
m = length(c);
D = full(diag(A));
y = rand(n,1);
lambda = norm(y);
for i = 1:10
    y = (A*(y/lambda))./D;
    lambda = norm(y);
end
[z, j] = max(s(f,c),[],2);
t = z ~= 0;
P = sparse(c,1:m,1,n,m) + sparse(f(t),j(t),1,n,m);
[i,j,x] = find(A*P);
z = (4/(3*lambda))./D;
P = P - sparse(i,j,x.*z(i),n,m);
end
