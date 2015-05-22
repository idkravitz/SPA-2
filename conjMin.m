function [fmin, xmin] = conjMin(f, fdot, x0, version)
%CONJMIN minimizer based on approximation of conjugate function at 0
% f - function, that is subject to minimization
% fdot - function gradient, approximation is fine too
% x0 - some value of f argument, used to determine shape of X space
% maxit = 1e6;
maxit = 500;

ERR = 1e-7;
N = size(x0, 1);
wmin = -inf;
wmax = inf;
xk = x0;

% Coefficients for LP target functions
c = [zeros(N, 1); 1];
lb = ones(N + 1, 1) * -inf;
ub = ones(N + 1, 1) * inf;

eps = 1;
P = [
fdot(x0), zeros(N, 1);
x0' * fdot(x0) - f(x0), 500
];

% rand('seed', 13569);
% P = [
% fdot(x0), eye(N) * eps, eye(N) * -eps;
% x0' * fdot(x0) - f(x0), 500 + rand(1, 2*N)*3
% ];
converged = false;
for k = 1:maxit
    wk = -f(xk);
    wmin = max(wmin, wk);
    xv = glpk(c, [P(1:end-1,:)', ones(size(P, 2), 1) * -1], P(end, :)', lb, ub, repmat('U', 1, size(P, 2)), repmat('C', 1, N + 1));
    xtilde = xv(1:end-1);
    fmin = -xv(end);
    wmax = min(wmax, fmin);
    % wmin
    % wmax
    if fmin > wmax
        printf('alert');
    end

    err = abs(wmax - wmin);
    fprintf('%3d. err %g\n', k, err);
    if abs(wmax - wmin) < ERR
        converged = true;
        break;
    end

    w = [zeros(N, 1); wmin];

    % Projection
    [z, reps, iter, lmb, kvec, R1, info] = SimPro(P - repmat(w, 1, size(P, 2)), 1e8, 1e-16, [-1, -1], [], []);
    
    dirZ = z / -z(end);
    
    % Search for next polytope point
    xbar = dirZ(1:N);

    lmDot = @(lm) [
        f((xbar + lm(2)*xtilde)/(1+lm(1)+lm(2))) + wmax - fdot((xbar + lm(2)*xtilde)/(1+lm(1)+lm(2)))'*(xbar + lm(2)*xtilde)/(1+lm(1)+lm(2)),
        f((xbar + lm(2)*xtilde)/(1+lm(1)+lm(2))) + wmax - fdot((xbar + lm(2)*xtilde)/(1+lm(1)+lm(2)))'*(xbar - (1 + lm(1))*xtilde)/(1+lm(1)+lm(2))
    ];
    
    % P
    % xbar
    % gbar_old = fdot(xbar);
    % fbar_old = xbar' * gbar_old - f(xbar);
    % xtilde
    % fmin

    % if ( index(lastwarn(), 'division by zero') > 0 )
     % error('division by zero before cmMin')
    % endif

    % Version selection, possible variants: 0 for SPA, 1 for SPA-1, 2 for SPA-2
    if version == 2
        lm = cmMin(lmDot, [0, 100, 0; 0, 0, 100], 0);
    % if ( index(lastwarn(), 'division by zero') > 0 )
     % error('division by zero after cmMin')
    % endif

        xnew = (xbar + lm(2)*xtilde) / (1+lm(1)+lm(2));
    elseif version == 1
        lmFun = @(l) (1+l)*f(xbar/(1+l)) + l*wmax;
        l = fminbnd(lmFun, 0, 10000, struct('TolX', 1e-12, 'TolFun', 1e-12));
        xnew = xbar/(1+l);
    elseif version == 0
        xnew = xbar;
    end
    gnew = fdot(xnew);
    fnew = xnew' * gnew - f(xnew);
    % fnew

    pk = [gnew; fnew];
    P = [P, pk];

    xk = xnew;
end

fmin = -(wmax + wmin) * 0.5;
xmin = xk;

end