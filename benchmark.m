function benchmark(fun, fundot)
	x0 = ones(10, 1) * 5;
	printf('benchmark of SPA\n');
	[fmin, xmin] = conjMin(fun, fundot, x0, 0);
	printf('benchmark of SPA-1\n');
	[fmin, xmin] = conjMin(fun, fundot, x0, 1);
	printf('benchmark of SPA-2\n');
	[fmin, xmin] = conjMin(fun, fundot, x0, 2);
	printf('benchmark of R-ALG\n');
endfunction

printf('Simplest cases\n');

global a = rand();
global b = rand();

function y = linf(x)
	global a;
	global b;
	szx = size(x);
	if !isequal(size(a), szx)
		a = rand(szx) * 1;
	end
	if !isequal(size(b), szx)
		b = rand(szx) * 1 - .5;
	end
	y = a' * abs(x - b);
endfunction

function dy = linfdot(x)
	global a;
	global b;
	szx = size(x);
	if !isequal(size(a), szx)
		a = rand(szx) * 1;
	end
	if !isequal(size(b), szx)
		b = rand(szx) * 1 - 0.5;
	end
	dy = a .* sign(x - b);
endfunction

benchmark(@linf, @linfdot)