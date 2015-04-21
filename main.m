% f = @(x) (abs(x) + 5)'*ones(size(x));
% fdot = @(x) sign(x);

% f = @(x) (abs(x + 5));
% fdot = @(x) sign(x + 5);

%f = @(x) ((x + 5) ** 2 + 5);
%fdot = @(x) 2 * (x + 5);

%f = @(x) (abs(x) + 7);
%fdot = @(x) (sign(x));

% f = @(x) (norm(x - 1) ** 2 + 5);
% fdot = @(x) (2 * (x - 1));

% f = @(x) ((x-1)^2 + 5);
% fdot = @(x) (2 * (x - 1));
global B b m;
n = 10; m = 5;

B = zeros(n, n, m);
b = zeros(m, n);
en = [ 1:n ];
a1 = en'*ones(1,n);
a2 = exp(min(a1, a1') ./ max(a1, a1'));
a3 = cos(en'*en);
for k = 1:m
        A = a2 .* a3 * sin(k);
        A = A - diag(diag(A));
        B(:,:,k) = A + diag(sum(abs(A)) + abs(sin(k))*en/n);
        b(k, :) = exp(en/k) .* sin(en*k);
endfor


function [ f g ] = pwq(x)
global B b m;
ff = zeros(1,m);
for k = 1:m
        ff(k) =  x'*B(:,:,k)*x + b(k, :)*x;
endfor
[ f indx ] = max(ff);
g = 2*x'*B(:,:,indx) + b(indx,:);
g = g';
endfunction

function [ f g ] = pwqfinal(x)
global B b m;
ff = zeros(1,m);
for k = 1:m
        ff(k) =  x'*B(:,:,k)*x + b(k, :)*x;
endfor
[ f indx ] = max(ff);
g = vec(2*x'*B(:,:,indx) + b(indx,:));
endfunction

function y = maxQuadF(x)
	[y g] = pwqfinal(x);
endfunction

function g = maxQuadFDot(x)
	[y g] = pwqfinal(x);
endfunction	

[fmin, xmin] = conjMin(@maxQuadF, @maxQuadFDot, ones(10, 1) * 5)