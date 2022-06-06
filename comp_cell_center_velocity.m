function [Vx_c,Vy_c] = comp_cell_center_velocity(v,Xc,Yc,Grid) % MDOT repo
% author: Marc Hesse 
% date: 29 May 2020
% Description:
% This functions averages the face velocities to the cell centers for a
% standard tensor-product staggered mesh.

% Interploate x-velocities
Vx = reshape(v(1:Grid.Nfx),Grid.Ny,Grid.Nx+1);
[Xx,Yx] = meshgrid(Grid.xf,Grid.yc);
Vx_c = interp2(Xx,Yx,Vx,Xc,Yc);

% Interpolate y-velocities
Vy = reshape(v(Grid.Nfx+1:end),Grid.Ny+1,Grid.Nx);
[Xy,Yy] = meshgrid(Grid.xc,Grid.yf);
Vy_c = interp2(Xy,Yy,Vy,Xc,Yc);
end