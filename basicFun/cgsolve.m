function [ x, res, iter ] = cgsolve( A, b, tol, maxiter, verbose )
%CGSOLVE Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 5), verbose = 1; end

implicit = isa(A, 'function_handle');

x = zeros(length(b),1);
r = b;
d = r;
delta = r'*r;
delta0 = b'*b;
numiter = 0;
bestx = x;
bestres = sqrt(delta/delta0);
while ((numiter < maxiter) && (delta > tol^2*delta0))
    if (implicit), q = A(d); else q = A*d; end
    
    alpha = delta/(d'*q);
    x = x + alpha*d;
    
    if (mod(numiter+1,50) == 0)
        if (implicit), r = b - A(x); else r = b - A*x; end
    else
        r = r - alpha * q;
    end
    
    deltaold = delta;
    delta = r'*r;
    beta = delta/deltaold;
    d = r + beta*d;
    numiter = numiter + 1;
    if (sqrt(delta/delta0) < bestres)
        bestx = x;
        bestres = sqrt(delta/delta0);
    end
    
    if ((verbose) && (mod(numiter, verbose)==0))
        disp(sprintf('cg: Iter = %d, Best residual = %8.3e, Current residual = %8.3e',...
            numiter, bestres, sqrt(delta/delta0)));
    end
    
end

if (verbose)
    disp(sprintf('cg:Iterations = %d, best residual = %14.8e', numiter, bestres));
end
x = bestx;
res = bestres;
iter = numiter;

end

