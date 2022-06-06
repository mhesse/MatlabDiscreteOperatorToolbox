function [mm] = minmod(a,b) % MDOT repo
% name: minmod.m
% author: Marc Hesse - actually got it from web 
%         https://github.com/wme7/aero-matlab/blob/master/CFD/minmod.m
% date: 27 Jul 2016
% Description: compares two vaues and chooses the smaller one if they have
% the same sign or returns a zero if they have opposite sign. 

mm = zeros(size(a));
mm =    ((abs(a)<=abs(b)).*(a.*b>0)).*a;
mm = mm+((abs(b)<abs(a)).*(a.*b>0)).*b;