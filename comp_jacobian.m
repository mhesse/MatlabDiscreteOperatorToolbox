function [J] = comp_jacobian(r,u,eps)
% This code was taken from the following lecture notes
% http://www.maths.lth.se/na/courses/FMN081/FMN081-06/lecture7.pdf
% by Claus Fuehrer at the University of Lund - thanks.
n=length(u);
u_perturb=u;
for i=1:n
   u_perturb(i)=u_perturb(i)+eps;
   J(:,i)=(r(u_perturb)-r(u))/eps;
   u_perturb(i)=u(i);
end
end