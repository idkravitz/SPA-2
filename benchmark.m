function benchmark(fun, fundot, dim=10)
	x0 = ones(dim, 1) * 5;
	printf('benchmark of SPA\n');
	[fmin, xmin] = conjMin(fun, fundot, x0, 0);
	printf('benchmark of SPA-1\n');
	[fmin, xmin] = conjMin(fun, fundot, x0, 1);
	printf('benchmark of SPA-2\n');
	[fmin, xmin] = conjMin(fun, fundot, x0, 2);
	printf('benchmark of R-ALG\n');
endfunction

global a = rand();
global b = rand();

function y = linf(x)
	global a;
	global b;
	szx = size(x);
	if !isequal(size(a), szx)
		a = (rand(szx) * 1) * 3;
	end
	if !isequal(size(b), szx)
		b = (rand(szx) * 1 - .5) * 10;
	end
	y = a' * abs(x - b);
endfunction

function dy = linfdot(x)
	global a;
	global b;
	szx = size(x);
	if !isequal(size(a), szx)
		a = (rand(szx) * 1) * 3;
	end
	if !isequal(size(b), szx)
		b = (rand(szx) * 1 - 0.5) * 10;
	end
	dy = a .* sign(x - b);
endfunction

dim = 10;
printf('Simplest cases, dim = %d\n', dim);
benchmark(@linf, @linfdot, dim);
dim = 100;
printf('Simplest cases, dim = %d\n', dim);
benchmark(@linf, @linfdot, dim);