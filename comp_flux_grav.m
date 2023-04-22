function [q] = comp_flux_grav(D,Lam,G,p,fs,Grid,BC,Rhod,grav_vec) % messing with this
% author: Marc Hesse
% date: 25 Nov 2014, 10 Jul 2015, 19 Jun 2017, 13 Jan 2023
% Description:
% Computes the mass conservative fluxes across all boundaries from the 
% residual of the compatability condition over the boundary cells.
% Note: Current implmentation works for all cases where one face 
%       is assigend to each bnd cell. So conrner cells must have
%       natural BC's on all but one face.
%
% Input:
% D = N by Nf discrete divergence matrix.
% Lam  = Nf by Nf diagonal matrix of averaged mobilities (lam = k/mu).
% Rho  = Nf by Nf diagonal matrix of averaged fluid densities
% G = Nf by N discrete gradient matrix.
% p = N by 1 vector of pressures in cell centers.
% fs = N by 1 right hand side vector containing only source terms.
% Grid = structure containing grid information.
% BC = structure contaning problem paramters and information about BC's
% Rhod = Nf by Nf diagonal density matrix 
% grav_vec = Nf by 1 vector of gravitational accelerations -
%            This is the vector that would be obtained from solving 
%            for the gravitational field, i.e., typically negative! 

%
% Output:
% q = Nf by 1 vector of fluxes across all cell faces
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I] = build_ops(Grid);
% >> L = -D*G; fs = ones(Grid.Nx,1);
% >> BC.dof_dir = Grid.dof_xmin;
% >> BC.dof_f_dir = Grid.dof_f_xmin;
% >> g = 0;
% >> BC.dof_neu = []; BC.dof_f_neu =[];
% >> [B,N,fn] = build_bnd(BC,Grid);
% >> h = solve_lbvp(L,fs+fn,B,g,N);
% >> q = comp_flux(D,1,G,h,fs,Grid,BC);

%% Check density input
% Worth checking because comp_flux_grav and build_bnd_grav need different
% rho's
[row,col] = size(Rhod);
if row ~= col; error('comp_flux_grav needs diagonal matrix containing densities averaged to faces.'); end

%% Identify boundary cells
dof_cell = [BC.dof_dir;BC.dof_neu];
dof_face = [BC.dof_f_dir;BC.dof_f_neu];

%% Compute interior fluxes
% Note: Recovery of bnd fluxes from the residual requires q = 0 on bnd.
%       This is automatic without gravity but has to be
%       set explicitly with gravity,because Rhod*grav_vec~=0 on bnd.
%       if the gravity vector is constructed to be zero on bnd we don't
%       need to set q(dof_face) = 0.
q = -Lam*(G*p-Rhod*grav_vec);
% q(dof_face) = 0;
res = D*q-fs;

%% Compute boundary fluxes
sign = ismember(dof_face,[Grid.dof_f_xmin;Grid.dof_f_ymin])...
      -ismember(dof_face,[Grid.dof_f_xmax;Grid.dof_f_ymax]);
q(dof_face) =  sign.*res(dof_cell).*Grid.V(dof_cell)./Grid.A(dof_face);  
% q(dof_face) =  sign.*( D(dof_cell,:) * q - fs(dof_cell)).*Grid.V(dof_cell)./Grid.A(dof_face);