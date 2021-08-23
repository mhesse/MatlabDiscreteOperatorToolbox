function [D,Edot,Dp,Gp,Z,I]=build_stokes_ops(Grid) % repo
% authors: Marc Hesse, Evan Carnahan
% date: 4 Oct. 2019

%% Build operators for each variable
[Dp,Gp,Ip] = build_ops(Grid.p);
[Dx,Gx,Ix] = build_ops(Grid.x);
[Dy,Gy,Iy] = build_ops(Grid.y);

%% Extract x and y components of the velocity operators
Gxx = Gx(1:Grid.x.Nfx,:); Gxy = Gx(Grid.x.Nfx+1:Grid.x.Nf,:);
Gyx = Gy(1:Grid.y.Nfx,:); Gyy = Gy(Grid.y.Nfx+1:Grid.y.Nf,:);

Dxx = Dx(:,1:Grid.x.Nfx); Dxy = Dx(:,Grid.x.Nfx+1:Grid.x.Nf);
Dyx = Dy(:,1:Grid.y.Nfx); Dyy = Dy(:,Grid.y.Nfx+1:Grid.y.Nf);

% % % To test symmetry
% Gxx = -Dxx'; Gxy = -Dxy';
% Gyx = -Dyx'; Gyy = -Dyy';
% Gp = -Dp';

%% Assemble Stokes operators
% Symmetric derivative
Zxy = spalloc(Grid.x.Nfx,Grid.y.N,0);
Zyx = spalloc(Grid.y.Nfy,Grid.x.N,0);
Edot = [Gxx, Zxy;...
        Zyx, Gyy;...
     .5*Gxy, .5*Gyx];

 % Divergence of deviatoric stess tensor
D = [Dxx ,Zyx',Dxy;...
     Zxy',Dyy ,Dyx];
 
 % Zero block for the system matrix L
 Z = spalloc(Grid.p.N,Grid.p.N,0);

 % Identity for all dof's
 I = speye(Grid.p.Nf+Grid.p.N);