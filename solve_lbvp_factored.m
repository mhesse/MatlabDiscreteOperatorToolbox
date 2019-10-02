function [u] = solve_lbvp_factored(Lred_L,Lred_U,L,f,B,g,N)
% author: Marc Hesse
% date: 26 Sept 2014, 27 Apr 2016
% Description
% Computes the solution $u$ to the linear differential problem given by
%
% $$\mathcal{L}(u)=f \quad x\in \Omega $$
%
% with boundary conditions
%
% $$\mathcal{B}(u)=g \quad x\in\partial\Omega$$.
% To improve the performance with multiple r.h.s.'s the LU-factorization
% of the linear operator is used
%
% Input:
% Lred_lo = lower diagonal matrix from the LU factorization of Lred
% Lred_up = upper diagonal matrix from the LU factorization of Lred
%           where Lred = N'*L*N
% L = matrix representing the discretized linear operator of size N by N, 
%     where N is the number of degrees of fredom
% f = column vector representing the discretized r.h.s. and contributions
%     due non-homogeneous Neumann BC's of size N by 1
% B = matrix representing the constraints arising from Dirichlet BC's of
%     size Nc by N
% g = column vector representing the non-homogeneous Dirichlet BC's of size
%     Nc by 1.
% N = matrix representing a orthonormal basis for the null-space of B and
%     of size N by (N-Nc).
% Output:
% u = column vector of the solution of size N by 1
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> L = -D*G; fs = ones(Grid.Nx,1);
% >> dof_dir = 1;
% >> B = I(dof_dir,:); g = 1; 
% >> N = I; N(:,dof_dir) = [];
% >> h = solve_lbvp(L,fs,B,g,N);

if ~istril(Lred_L)
    error('Matrix Lred_L is not lower triangular - likely permutation.')
end
if  ~istriu(Lred_U)
    error('Matrix Lred_U is not upper triangular - likely permutation.')
end
if isempty(B) % no constraints
    u = Lred_U\(Lred_L\f);
else
    up = spalloc(length(f),1,length(g));
    up = B'*(B*B'\g);
    u0 = N*( Lred_U\(Lred_L\(N'*(f-L*up))) );
    u = u0 + up;
end