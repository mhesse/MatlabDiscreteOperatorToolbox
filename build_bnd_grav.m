function [B,N,fn] = build_bnd_grav(BC,Grid,I,Lam,Rho,grav_vec) % messing with this
% author: Marc Hesse
% date: 06/09/2015, 13 Jan 2023
% Description:
% This function computes the operators and r.h.s vectors for both Dirichlet
% and Neumann boundary conditions. 
% Note:
% N is not created from I the same way B is created from I, because 
% the vector dof_dir contains the columns that must be eliminated rather
% then the colums that are retained in N. If you wanted to do it this way
% you would have to create a new vector
% dof_non_dir = setdiff(dof,dof_dir)
% I suspect that the set-operators are expensive on large vectors, hence
% we simply eliminate the rows.
%
% Input:
% Grid  = structure containing all pertinent information about the grid.
% BC    = structure containing all information about the physical problem
%         in particular this function needs the fields
%         BC.dof_dir = Nc by 1 column vector containing 
%                      the dof's of the Dirichlet boundary.
%         BC.dof_neu = Nn by 1 column vector containing 
%                      the dof's of the Neumann boundary.
%         BC.qb      = Nn by 1 column vector of prescribed fluxes on Neuman bnd.
% I    = identity matrix in the full space
% Lam  = Nf by Nf diagonal matrix of averaged mobilities (lam = k/mu).
% Rho  = Nf by Nf diagonal matrix of averaged fluid densities
% grav_vec = Nf by 1 vector of gravitational accelerations -
%            This is the vector that would be obtained from solving 
%            for the gravitational field, i.e., typically negative! 
%
% Output:
% B = Nc by N matrix of the Dirichlet constraints
% N = (N-Nc) by (N-Nc) matrix of the nullspace of B
% fn = N by 1 r.h.s. vector of Neuman contributions
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);
% >> Param.dof_dir   = Grid.dof_xmin;    % identify cells on Dirichlet bnd
% >> Param.dof_f_dir = Grid.dof_f_xmin;  % identify faces on Dirichlet bnd
% >> Param.dof_neu   = Grid.dof_xmax;    % identify cells on Neumann bnd
% >> Param.dof_f_neu = Grid.dof_f_xmax;  % identify cells on Neumann bnd
% >> Param.qb = 1;                   % set bnd flux
% >> [B,N,fn] = build_bnd(Param,Grid,I);

%% Check input format
if isrow(BC.dof_dir)   && length(BC.dof_dir)>1;               error('BC.dof_dir is not a column vector');   end
if isrow(BC.dof_neu)   && length(BC.dof_neu)>1;               error('BC.dof_neu is not a column vector');   end
if isrow(BC.dof_f_dir) && length(BC.dof_f_dir)>1;             error('BC.dof_f_dir is a not column vector'); end
if isrow(BC.dof_f_neu) && length(BC.dof_f_neu)>1;             error('BC.dof_f_neu is a not column vector'); end
if isfield(BC,'qb')    && isrow(BC.qb) && length(BC.qb)>1; error('BC.qb is not a column vector');        end

%% Check density input
% Worth checking because comp_flux_grav and build_bnd_grav need different
% rho's
[row,col] = size(Rho);
% if row == col; error('build_bnd_grav.m needs column vector containing cell densities.'); end

%% Dirichlet boundary conditions
% It seems we don't need this if-statement (Marc, 6 June 2017)
% if isempty(BC.dof_dir)
%     B = [];
%     N = [];
% else
%     B = I(BC.dof_dir,:);
%     N = I; N(:,BC.dof_dir) = [];
% end

B = I(BC.dof_dir,:);
N = I; N(:,BC.dof_dir) = [];

%% Neumann boundary conditions
% Note: Here it IS important that bnd entries of Kd ~= 0! (Feb 8, 2023 - I think this comment was because grav_vec wasn't zero on boundary)
if isempty(BC.dof_neu)
    fn = spalloc(Grid.N,1,0);                     % allocate sparse zero vector
else
    fn = spalloc(Grid.N,1,length(BC.dof_neu)); % allocate sparse vector
%     fn(BC.dof_neu) = ( BC.qb - diag(Kd(BC.dof_f_neu,BC.dof_f_neu)).*rho(BC.dof_neu).*grav_vec(BC.dof_f_neu) ).*Grid.A(BC.dof_f_neu)./Grid.V(BC.dof_neu);
    % q_grav = diag(Lam(BC.dof_f_neu)).*diag(Rho(BC.dof_f_neu))*grav_vec(BC.dof_f_neu);
    q_grav = Lam*Rho*grav_vec; 
    fn(BC.dof_neu) = (BC.qb - q_grav(BC.dof_f_neu)).*Grid.A(BC.dof_f_neu)./Grid.V(BC.dof_neu);
end

