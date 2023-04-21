% Copyright or © or Copr. 
% Authors: Christos Papoutsellis, Matthieu Mercier, Nicolas Grisouard
% Initial version 2021, Current Version 2023
% 
% cpapoutsellis@gmail.com
% matthieu.mercier@imft.fr
% nicolas.grisouard@utoronto.ca
% 
% This software is a computer program whose purpose is to solve numerically a
% a mathematical model describing the generation of internal tides.
% 
% This software is governed by the [CeCILL|CeCILL-B|CeCILL-C] license under
% French law (see details in ITCMS.m)



function f = gradientp(y,dx)
%%  numerical differentiation of a scalar function with 4th-order finite differences
N = length(y);
f = zeros(size(y));
% one-sided 
f(1) = ( (-25/12)*y(1) + 4*y(2) - 3*y(3) +(4/3)*y(4) -(1/4)*y(5) )/dx;
% assymetric
f(2) = ( -(1/4)*y(1) - (5/6)*y(2) + (3/2)*y(3) -(1/2)*y(4) +(1/12)*y(5) )/dx;
% central
for i=3:N-2
f(i) = ( y(i-2)-8*y(i-1)+8*y(i+1)-y(i+2) )/12/dx;
end
% assymetric
f(N-1) = ( (1/4)*y(N) +(5/6)*y(N-1) - (3/2)*y(N-2) + (1/2)*y(N-3) - (1/12)*y(N-4) )/dx;
% one-sided 
f(N) = ( (25/12)*y(N) - 4*y(N-1) + 3*y(N-2) - (4/3)*y(N-3) + (1/4)*y(N-4) )/dx;



