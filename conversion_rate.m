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



function [C,Cleft,Cright,Cint] = conversion_rate(ph,phx,x,h,hx,xLi,xRi,Q,N,fc,h0,omega,rho)
% Calculation of the conversion rates C_{\pm} and C_{int}

mu = sqrt(N^2-omega^2)/sqrt(omega^2-fc^2);
Nx = length(h);
Nm = size(ph,1);
Nmm = Nm;

% left boundary point
i = xLi; 
Cleft = -rho*(N^2-omega^2)/(4*omega*mu)*pi* [1:Nmm]*( ph(1:Nmm,i).*conj(ph(1:Nmm,i)) ) ;
%right boundary point
i = xRi; 
Cright = rho*(N^2-omega^2)/(4*omega*mu)*pi*[1:Nmm]*( ph(1:Nmm,i).*conj(ph(1:Nmm,i)) ) ;
% total conversion rate
C = Cright-Cleft;
% alternative calculation of C (total) through integration of the rhs of
nn = 1:Nmm;
I1 = zeros(1,Nx);
I2 = I1;
for i = 1:Nx
    I1(i) = -2*(-1).^nn./nn/pi * hx(i)^2/h(i)*imag(ph(1:Nmm,i))*Q;
    I2(i) =  - (-1).^nn./nn/pi * hx(i)*imag(phx(1:Nmm,i))*Q;
end
Cint = -N^2/2/omega*(1-omega^2/N^2)*rho*trapz(x(xLi:xRi),I1(xLi:xRi)+I2(xLi:xRi));




