function [q] = comp_flux_ade(flux,res,u,uold,v,theta,dt,Grid,BC)
% author: Marc Hesse
% date: 29 Apr 2020
% Description:
% Computes the fuxes on the interior from discrete ADE operator and 
% reconstructs the fluxes on the boundary faces from the residuals in 
% the adjacent boundary cells, res(u). 
%
% Note to self: with more experience it seems impossible to write a truely
% generic comp_flux function, so it may be better to be specific

% Input:
% u = unknown at time level n+1
% uold = unknown at time level n
% v = velocity field
% flux = anonymous function computing the flux (correct in the interior)
% res = anonymous function computing the residual 
% Grid = structure containing pertinent information about the grid
% BC = structure containing pertinent information about BC's
% 
% Output:
% q = correct flux everywhere
%
% Example call:


%% Compute interior fluxes
q = flux(dt,theta,v,u,uold);

%% Residual 
r = res(dt,theta,v,u,uold);

%% Compute boundary fluxes
% 1) Identify the faces and cells on the boundary
dof_cell = [BC.dof_dir;BC.dof_neu];
dof_face = [BC.dof_f_dir;BC.dof_f_neu];
% 2) Determine sign of flux: Convention is that flux is positive in
%    coordinate direction. So the boundary flux, qb is not equal to q*n,
%    were n is the outward normal! 
sign = ismember(dof_face,[Grid.dof_f_xmin;Grid.dof_f_ymin])...
      -ismember(dof_face,[Grid.dof_f_xmax;Grid.dof_f_ymax]);
% 3) Compute residuals and convert them to bnd fluxes

q(dof_face) =  sign.*r(dof_cell).*Grid.V(dof_cell)./Grid.A(dof_face);