function [Qxc,Qyc] = interp_f2c(q,Grid,Xc,Yc) % MDOT repo
% author: Marc Hesse
% date: 4 Mar 2020
% Description: Interpolate quantities from the faces to the cell centers
% Reshape x and y fluxes
Qx = reshape(q(1:Grid.Nfx)       ,Grid.Ny  ,Grid.Nx+1); 
Qy = reshape(q(Grid.Nfx+1,Grid.Nf,Grid.Ny+1,Grid.Nx  );

% Coordinates of the x and y face centers
[Xx,Yx] = meshgrid(Grid.xf,Grid.yc);
[Xy,Yy] = meshgrid(Grid.xc,Grid.yf);

Qxc = interp2(Xx,Yx,Qx,Xc,Yc);
Qyc = interp2(Xy,Yy,Qy,Xc,Yc);
end
