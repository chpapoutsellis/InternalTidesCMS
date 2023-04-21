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


function [pphi,pphix,pphiz,pphih,pphihx,pphihz,pphir1, zz  ] = ...
                            phis(ph,phx,h,hx,hxx,Nv,epsilon,epsilonc,omega,fc,Q)
                        
% Modal reconstruction of the daggered spatial fields
% INPUTS
%   ph      : x-dependent modal coefficients in vectors (the output of CMS_solution)
% h,hx,hxx  : the topography and its derivatives
% Nv        : Number of vertical points
% omega     : tidal frequency
%   fc      : Coriolis frequency
%   Q       : flux U00*h0
% epsilon   : criticality
% epsilonc  : if epsilon>epsilonc then a filter is applied in order to
% attenuate the Gibbs phenomenon present when the baroclinic solution is singular

% OUTPUT:
%  pphi,pphix, pphiz      : phi^{\dagger} and its derivatives
%  pphih,pphihx,pphihz    :  Phi^(0) and its derivatives
%  pphir1                 : first-order non-hydrostatic response
%  zz                     : the discretised interval [-h(x),0] for all x

hwait = waitbar(0, 'Modal reconstruction...');



Nx = length(h);
Nmm = size(ph,1);

zz = zeros(Nx,Nv);
pphi = zeros(Nv,Nx); 
pphiz = pphi;
pphix = pphi;

sigman = sin(pi*[1:Nmm]/(Nmm))./(pi*[1:Nmm]/(Nmm))*(epsilon>=epsilonc) + ones(1,Nmm)*(epsilon<epsilonc) ;
sigman = sigman.^2;                        

for j=1:Nx
    zz(j,:) = linspace(-h(j),0,Nv);
    ZZ = funZ(zz(j,:)',h(j) , Nmm);
    ZZZ = funZz(zz(j,:)',h(j) , Nmm);
    ZZX = funZx(zz(j,:)',h(j) , hx(j), Nmm);
    for i=1:Nv
        pphi(i,j) = ZZ(i,:)*(ph(1:Nmm,j).*sigman');
        pphiz(i,j) = ZZZ(i,:)*(ph(1:Nmm,j).*sigman');
        pphix(i,j)= ZZX(i,:)*(ph(1:Nmm,j).*sigman') +ZZ(i,:)*(phx(1:Nmm,j).*sigman');
    end
end


pphih = zeros(Nv,Nx); 
pphihz = pphih;
pphihx = pphih;
pphir1 = pphih;

S = 2*hx.^2./h.^3 - hxx./h.^2;
for j=1:Nx
    for i=1:Nv
        pphih(i,j)   = zz(j,i)*Q/(-h(j));
        pphihz(i,j)  = Q/(-h(j));
        pphihx(i,j)  = -zz(j,i)*Q*hx(j)/(-h(j)).^2;
        pphir1(i,j)  = Q*(sqrt(1-fc^2/omega^2))^2* S(j)*1/6*(zz(j,i)^3-h(j)^2*zz(j,i)); 
    end
end


delete(hwait)


function y = funZ(z,h,Nm)
y = zeros(length(z),Nm);
for n=1:Nm
    %y(:,n) = sin(n*pi*z/(h0-h));
    y(:,n) = sin(n*pi*z/h);
end

function y = funZz(z,h,Nm)
y = zeros(length(z),Nm);
for n=1:Nm
    y(:,n) = cos(n*pi*z/h).*(n*pi./h);
end

function y = funZx(z,h,hx,Nm)
y = zeros(length(z),Nm);
for n=1:Nm
    y(:,n) = -cos(n*pi*z/h).*(n*pi*z./h.^2)*hx;
end






