function [] = plot_flownet(Nh,Ns,h,PSI,head,stream,Grid) % MDOT repo
% author: Marc Hesse
% date: 24 Oct 2015
% Description: Plots a flownet with an equal number of equally spaced head 
%              contours and streamlines.        
% Input: N = numer of contours & streamlines
%        h = Grid.N by 1 column vector of heads
%        PSI = Ny+1 by Nx+1 matrix containing the stream function
%        head = string specifying the linestyle for the head contours
%        stream = string specifying the linestyle for the streamlines
%        Grid = structure containing all information about the grid.
% Output: none
hmin = min(h); hmax = max(h);
psi_min = min(PSI(:)); psi_max = max(PSI(:));

h_cont = linspace(hmin,hmax,Nh);
psi_cont = linspace(psi_min,psi_max,Ns);

[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
contour(Xc,Yc,reshape(h,Grid.Ny,Grid.Nx),h_cont,head), hold on
[Xp,Yp] = meshgrid(Grid.xf,Grid.yf);
contour(Xp,Yp,PSI,psi_cont,stream)
xlim([Grid.xmin Grid.xmax]), ylim([Grid.ymin Grid.ymax])
xlabel 'x', ylabel 'y'