function [Q,q] = comp_flux_dupuit(D,Kd,G,h,fs,Grid,BC) % MDOT repo
% author: Marc Hesse
% date: 25 Nov 2014, 10 Jul 2015, 23 Mar 2021
% Description:
% Computes the mass conservative flow rate across all boundaries from the 
% residual of the compatability condition over the boundary cells.
% Note: Current implmentation works for all cases where one face 
%       is assigend to each bnd cell. So conrner cells must have
%       natural BC's on all but one face.
%
% Input:
% D = N by Nf discrete divergence matrix.
% Kd = Nf by Nf conductivity matrix.
% G = Nf by N discrete gradient matrix.
% h = N by 1 vector of flow potential in cell centers.
% fs = N by 1 right hand side vector containing only source terms.
% Grid = structure containing grid information.
% BC = structure contaning problem paramters and information about BC's
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

%% Compute interior fluxes
M = Grid.dx/2*abs(G); % matrix to the arithmetic mean
h_mean = M*h; h(h==0) = nan;
Hd = spdiags(h_mean,0,Grid.Nfx,Grid.Nfx);
Q = -Kd*Hd*G*h;

%% Compute boundary fluxes
% note: check if o.k. for homoeneous Neumann problem
dof_cell = [BC.dof_dir;BC.dof_neu];
dof_face = [BC.dof_f_dir;BC.dof_f_neu];
sign = ismember(dof_face,[Grid.dof_f_xmin;Grid.dof_f_ymin])...
      -ismember(dof_face,[Grid.dof_f_xmax;Grid.dof_f_ymax]);

Q(dof_face) =  sign.*( D(dof_cell,:) * Q - fs(dof_cell)).*Grid.V(dof_cell)./Grid.A(dof_face);
q = Q./h_mean; 