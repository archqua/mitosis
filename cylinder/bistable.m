global A = 1.128e-04
%A
function res = B(a0)
	res = 1.00204 - 5.88e-05 .* a0;
endfunction
function res = C(a0)
	res = -2.0028 .* a0;
endfunction

function res = D(a0)
	global A;
	%A
	%B(a0)
	%C(a0)
	res = B(a0).*B(a0) - 4 .* A .* C(a0);
endfunction

function res = b1(a0)
	global A;
	res = (-B(a0) - sqrt(D(a0))) ./ (2 * A);
endfunction
function res = b2(a0)
	global A;
	res = (-B(a0) + sqrt(D(a0))) ./ (2 * A);
endfunction

figure(1);
xx = linspace(1e+03, 1e+05, 101);
y1 = b1(xx);
y2 = b2(xx);
clf;
plot( xx, y1
    , xx, y2
    , xx, xx
    );
grid on;
legend([ "b_1"
       ; "b_2"
       ; "a_0"]
       , "location", "northwest"
       , "fontsize", 14);
xlabel("total aur (a_0)", "fontsize", 14);
ylabel("active aur", "fontsize", 14);

global gama = 1.4e-03;
global phi = 0.042;
function res = a1(a0)
	global A; global phi; global gama;
	res = (A/phi * b1(a0) - gama * a0) / (2 - gama);
endfunction
function res = a2(a0)
	global A; global phi; global gama;
	res = (A/phi * b2(a0) - gama * a0) / (2 - gama);
endfunction

figure(2);
yy1 = a1(xx);
yy2 = a2(xx);
clf;
plot( xx, yy1
    , xx, yy2
    , xx, xx
    );
grid on;
legend([ "a_1"
       ; "a_2"
       ; "a_0"]
       , "location", "northwest"
       , "fontsize", 14);
xlabel("total aur (a_0)", "fontsize", 14);
ylabel("inactive aur", "fontsize", 14);
