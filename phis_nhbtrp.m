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

function [pphinhCMS,pphixnhCMS,pphiznhCMS,forcing, zz  ] = ...
                            phis_nhbtrp(phbar_nhCMS,phxbar_nhCMS,h,hx,hxx,Nv)

% Modal reconstruction of spatial fields corresponding to non-hydrostatic barotropic response 
% INPUTS
%  phbar_nhCMS      : x-dependent modal coefficients in vectors (the output of CMS_solution_brtrp)
% h,hx,hxx  : the topography and its derivatives
% Nv        : Number of vertical points

% OUTPUT:
%  pphinhCMS,pphixnhCMS, pphiznhCMS : Phi^{r} and its derivatives
%  forcing                          :  the body-forcing term -(Phi^(0))_xx 
%  zz                     : the discretised interval [-h(x),0] for all x

hwait = waitbar(0, 'Modal reconstruction (non-hydr. bartrp)...');

Nx = length(h);
Nbar = size(phbar_nhCMS,1);

pphinhCMS =  zeros(Nv,Nx);
pphiznhCMS = pphinhCMS;
pphixnhCMS = pphinhCMS;
forcing = pphinhCMS;
S = 2*hx.^2./h.^3 - hxx./h.^2;
for j=1:Nx
    zz(j,:) = linspace(-h(j),0,Nv);
    ZZ = funZ(zz(j,:)',h(j) , Nbar);
    ZZZ = funZz(zz(j,:)',h(j) , Nbar);
    ZZX = funZx(zz(j,:)',h(j) , hx(j), Nbar);
    for i=1:Nv
        pphinhCMS(i,j) = ZZ(i,:)*(phbar_nhCMS(1:Nbar,j));
        pphiznhCMS(i,j) = ZZZ(i,:)*(phbar_nhCMS(1:Nbar,j));
        pphixnhCMS(i,j)= ZZX(i,:)*(phbar_nhCMS(1:Nbar,j)) +ZZ(i,:)*(phxbar_nhCMS(1:Nbar,j));
    end
    forcing(:,j) = zz(j,:)*S(j);
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
