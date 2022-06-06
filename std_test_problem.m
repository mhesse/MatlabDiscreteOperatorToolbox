function f = std_test_problem(x) % MDOT repo
% The test profile all over LeVeque's book
f = exp(-40*(x+.5).^2) + (abs(x-.5)<.25);% + eps*rand(size(x));
