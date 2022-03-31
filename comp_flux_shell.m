function [q] = comp_flux_shell(G,Kd,u,Grid)
% author: Marc Hesse
% date: Mar 30, 2022
% Description: Computes the flux of a potential field on a spherical shell.
% An extra function is needed, because two entire boundaries are mapped to 
% the poles. The fluxes at these faces on the poles need to to be computed
% separately.
%
% Input:
% G = Nf by N discrete gradient matrix with periodice BC's in the azimuthal
% direction and natural BC's in the longitudinal direction
% Kd = Nf by Nf diagonal matrix cntaining the appropriate average of the
% conductivity on the faces
% u = N by 1 column vector of the potential
% Grid = structure containing useful information about the grid

%% Standard gradient flux
% This requires a discrete gradient with periodic rather than natural
% boundary conditions on the azimuthal boundaries.
q = -Kd*G*u;

%% Compute flux at south pole
% At the pole
% [dof_f_bnd,dof_f] = find_faces(Grid.dof_xmin,Dref,Grid);
% dof_f_pole = dof_f_bnd(~ismember(dof_f_bnd,[Grid.dof_f_xmin;Grid.dof_f_ymin;Grid.dof_f_ymax]));
[u_max_pole,imax] = max(u(Grid.dof_xmin));
[u_min_pole,imin] = min(u(Grid.dof_xmin));
[r,c] = size(Kd);

if r == 1 && c == 1
    K_pole = Kd;
else 
    K = diag(Kd); % 
    K_pole = mean(K(Grid.dof_xmin)); % have not thought about this in detail - arithmetic mean should be fine
end
q_max_pole = -K_pole*(u_max_pole-u_min_pole)/Grid.dx; % no K because dimensionless
dir_max = Grid.yc(imax);
q_pole = q_max_pole*cos(dir_max-Grid.yc);
q(Grid.dof_f_xmin) = q_pole;

%% Flux at North pole is not computed - because we don't need it for Mars
