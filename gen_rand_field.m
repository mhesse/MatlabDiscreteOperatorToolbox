function [K,X,Y] = gen_rand_field(Grid,X,Y,corr_length,amp,Kmean,type,rng_state) % MDOT repo
% [X,Y] = meshgrid(Grid.xc,Grid.yc);
x = X(:); y = Y(:);
if strcmp(type,'exp')
    sig = -log(.1)/corr_length;
else
    error('Unknown covariance model')
end
Cov = zeros(Grid.N,Grid.N);
for i = 1:Grid.N
    dist = sqrt((x(i) - x).^2 + (y(i) - y).^2);
    if strcmp(type,'exp')
        Cov(i,:) = exp(-sig * dist);
    else
        error('Unknown covariance model')
    end
end
rng(rng_state);
% Cholesky factorization is equivalent to square root of a matrix 
% Cov = L*L' <=> L ~ sqrt(Cov)
L = chol(Cov,'lower');
Kpert = reshape(L*randn(Grid.N,1),Grid.Ny,Grid.Nx);
K = Kmean + amp*Kpert;

