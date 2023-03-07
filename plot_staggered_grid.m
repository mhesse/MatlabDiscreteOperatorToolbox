function [fig1] = plot_staggered_grid(Grid) % MDOT repo
% author: Marc Hesse
% date: 15 Apr 2020
% Description: Plots staggered grid with unknowns
x = Grid.xf; y = Grid.yf; Nx = Grid.Nx+1; Ny = Grid.Ny+1;
[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);
[Xx,Yx] = meshgrid(Grid.xf,Grid.yc);
[Xy,Yy] = meshgrid(Grid.xc,Grid.yf);

fig1 = figure;
plot([x';x'],[Grid.xmin*ones(1,Nx);Grid.xmax*ones(1,Nx)],'k-'), hold on
plot([Grid.ymin*ones(1,Ny);Grid.ymax*ones(1,Ny)],[y';y'],'k-')
plot(Xc,Yc,'ko')
plot(Xx,Yx,'ro','markerfacecolor','w')
plot(Xy,Yy,'bo','markerfacecolor','w')
end