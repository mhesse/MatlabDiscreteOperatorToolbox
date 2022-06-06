function [A] = flux_LaxWend(q,Grid,G,dt,Param) % MDOT repo
% author: Marc Hesse, Jake Jordan
% date: 26 Jul 2016, 5 Aug 2016
% Description: Operator implementation of the 2nd-order Lax-Wendroff flux
%              See LeVeque (red book) Chp. 6.1, p. 100
%              My notes from 19 Apr 2016
%              The Lax-Wendroff flux is a combination of a upwind flux
%              and an anti-diffusion term
% Input:
% q    =  Nf by 1 flux vector
% Grid = sturcture containing grid information
% G    = discrete gradient operator
% dt   = timestep
% dof_out = column vector with the dof's of the outflow cells
%
% Output:
% A    = Nf by Nf matrix containing the Lax-Wendroff flux 
%
% Example call:
% BTD

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;
Nfx = Grid.Nfx; % # of x faces
Nfy = Grid.Nfy; % # of y faces
Nf = length(q);
dx = Grid.dx;

%% First-order upwind (Godunov) fluxes
A = flux_upwind(q,Grid);

%% Second-order Lax-Wendroff modification
% Build vector with cell sizes
if     Nx>1  && Ny==1
    dxy = ones(Grid.Nfx,1)*Grid.dx;
elseif Nx==1 && Ny>1 % 1D
    dxy = ones(Grid.Nfy,1)*Grid.dy;
    
else
    dxy = [ones(Grid.Nfx,1)*Grid.dx;ones(Grid.Nfy,1)*Grid.dy];
end

% Save upwind fluxes at outflow 
Aup = A(Param.dof_out,:);

%% 2nd-order Lax-Wendroff flux
A = A + spdiags(dxy.*abs(q)/2,0,Nf,Nf)*(speye(Nf)-spdiags(dt*abs(q)./dxy,0,Nf,Nf))*G;

% Reduce flux to upwind at outflow boundaries
A(Param.dof_out,:) = Aup;