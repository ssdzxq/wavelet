function xp = l1eq_pd( x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter )
% Solve
% min_x ||x||_1  s.t.  Ax = b

% Recast as linear program
% min_{x,u} sum(u)  s.t.  -u <= x <= u,  Ax=b
% and use primal-dual interior point method
%
% Usage: xp = l1eq_pd(x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - Either a handle to a function that takes a N vector and returns a K
%     vector , or a KxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes a K vector and returns an N vector.
%      If A is a KxN matrix, At is ignored.
%
% b - Kx1 vector of observations.
%
% pdtol - Tolerance for primal-dual algorithm (algorithm terminates if
%     the duality gap is less than pdtol). 
%     Default = 1e-3.
%
% pdmaxiter - Maximum number of primal-dual iterations. 
%     Default = 50.
%
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.
%     Default = 1e-8.
%
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored
%     if A is a matrix.
%     Default = 200.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%
 
largescale = isa(A, 'function_handle');

if (nargin < 5), pdtol = 1e-3;end
if (nargin < 6), pdmaxiter = 50;end
if (nargin < 7), cgtol = 1e-8;end
if (nargin < 8), cgmaxiter = 200;end

N = length(x0);

alpha = 0.01; % the user-specified parameter, for delta-value selection
beta = 0.5; % update s: s'=beta*s, s is the step length.
mu = 10; % meaning??

gradf0 = [zeros(N,1); ones(N,1)]; % c0，是随便选取的吗？？
x = x0; 
u = (0.95)*abs(x0) + (0.10)*max(abs(x0)); % 后面会更新，更新的式子是什么？？

% from the condition -u <= x <= u,there are two inequality constraits
fu1 = x - u;
fu2 = -x - u;

%initial lamu1 and lamu2 when k = 0 in the eq(2) in page 4
lamu1 = -1./fu1;
lamu2 = -1./fu2;

if (largescale)
    v = -A(lamu1 - lamu2);
    Atv = At(v);
    rpri = A(x) - b;
else
    v = -A*(lamu1-lamu2); % v是怎么计算的？？
    Atv = A'*v;
    rpri = A*x - b;
end

sdg = -(fu1'*lamu1 + fu2'*lamu2); % fu1, fu2都是负的，lamu1, lamu2都是正的，所以sdg>0
tau = mu*2*N/sdg; % 分子是如何确定的？？

rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)]; % 后面两项补充的一半是怎么得来的？
resnorm = norm([rdual; rcent; rpri]);

pditer = 0;
done = (sdg < pdtol) | (pditer >= pdmaxiter);
while (~done)
    pditer = pditer + 1;
    
    % 设置这6个变量的是为了求解微分量方程组
    % [w1;w2]是一个2N*1的向量，为了保持维数一致。自变量z,v,lamu的维数都是N*1
    w1 = -1/tau*(-1./fu1 + 1./fu2) - Atv; % 对应page4 (5)式右边第一项的前N行，c0置零？
    w2 = -1 - 1/tau*(1./fu1 + 1./fu2); % 对应page4 (5)式右边第一项的后N行，c0置1？
    w3 = -rpri; % 对应page4 (5)式右边第二项
    
    sig1 = -lamu1./fu1 - lamu2./fu2; % 对应page4 (5)式左边第一个矩阵第一项的前N行
    sig2 = lamu1./fu1 - lamu2./fu2; % 对应page4 (5)式左边第一个矩阵第一项的后N行
    sigx = sig1 - sig2.^2./sig1; % 求解方程组(5)过程中的一步吗？？
    
    if (largescale)
        w1p = w3 - A(w1./sigx - w2.*sig2./(sigx.*sig1));
        h11pfun = @(z) -A(1./sigx.*At(z));
        [dv, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);
        if (cgres > 1/2)
            disp('Primal-dual: Cannot solve system. Returning previous iterate.');
            xp = x;
            return;
        end
        dx = (w1 - w2.*sig2./sig1 - At(dv))./sigx;
        Adx = A(dx);
        Atdv = At(dv);
    else
        H11p = -A*diag(1./sigx)*A';
        w1p = w3 - A*(w1./sigx - w2.*sig2./(sigx.*sig1));
        % 为什么不能直接把前面得到的系数矩阵和值矩阵作为输入,这样不是可以直接得到[Δz Δv]吗？？
        % 为什么要构造H11p和w1p，而且求解的只有Δv？？
        [dv, hcond] = linsolve(H11p,w1p); 
        if (hcond < 1e-14)
            disp('Primal-dual: Matrix ill-conditioned. Returning previous iterate.');
            xp = x;
            return
        end
        dx = (w1 - w2.*sig2./sig1 - A'*dv)./sigx;
        Adx = A*dx;
        Atdv = A'*dv;
    end
    
    du = (w2 - sig2.*dx)./sig1;
    
    dlamu1 = (lamu1./fu1).*(-dx+du) - lamu1 - (1/tau)*1./fu1;
    dlamu2 = (lamu2./fu2).*(dx+du) - lamu2 - 1/tau*1./fu2;
    
    % make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
    indp = find(dlamu1 < 0); indn = find(dlamu2 < 0);
    s = min([1; -lamu1(indp)./dlamu1(indp); -lamu2(indn)./dlamu2(indn)]);
    indp = find((dx-du) > 0); indn = find((-dx-du) > 0);
    s = (0.99)*min([s; -fu1(indp)./(dx(indp)-du(indp)); -fu2(indn)./(-dx(indn)-du(indn))]);
    
    % backtracking line search
    backiter = 0;
    xp = x + s*dx; up = u + s*du;
    vp = v + s*dv; Atvp = Atv + s*Atdv;
    lamu1p = lamu1 + s*dlamu1; lamu2p = lamu2 + s*dlamu2;
    fu1p = xp - up; fu2p = -xp - up;
    rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
    rpp = rpri + s*Adx;
    
    % 迭代
    while (norm([rdp; rcp; rpp]) > (1-alpha*s)*resnorm)
        s = beta*s;
        xp = x + s*dx; up = u + s*du;
        vp = v + s*dv; Atvp = Atv + s*Atdv;
        lamu1p = lamu1 + s*dlamu1; lamu2p = lamu2 + s*dlamu2;
        fu1p = xp - up; fu2p = -xp - up;
        rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
        rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
        rpp = rpri + s*Adx;
        backiter = backiter + 1;
        if (backiter > 32)
            disp('Stuck backtracking, returning last iterate.');
            xp = x;
            return
        end
    end
    
    % update values
    x = xp; u = up;
    v = vp; Atv = Atvp;
    lamu1 = lamu1p; lamu2 = lamu2p;
    fu1 = fu1p; fu2 = fu2p;
    
    sdg = -(fu1'*lamu1 + fu2'*lamu2);
    tau = mu*2*N/sdg;
    rpri = rpp;
    rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
    rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
    resnorm = norm([rdual; rcent; rpri]);
    
    done = (sdg < pdtol) | (pditer >= pdmaxiter);
    
    disp(sprintf('Iteration = %d, tau = %8.3e, Primal = %8.3e, PDGap = %8.3e, Dual res = %8.3e, Primal res = %8.3e',...
        pditer, tau, sum(u), sdg, norm(rdual), norm(rpri)));
    if (largescale)
        disp(sprintf('                  CG Res = %8.3e, CG Iter = %d', cgres, cgiter));
    else
        disp(sprintf('                  H11p condition number = %8.3e', hcond));
    end

end

end

