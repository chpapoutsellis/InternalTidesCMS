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



function [ph] = CMS_solution_brtrp( NT , x ,h,hx,hxx, hL , hR , mu0 ,h0,U0) 

% Solution of the CMS with a 4th-order finite-difference method (non-hydrostatic barotropic case)
% INPUTS
%   Nm      : Number of modes
%    x      : a vector of uniformly discretised grid points
% h,hx,hxx  : the topography and its derivatives
%  hL,hR    : The depths at the left and right
%   mu0     : The parameter mu0 = mu0 = sqrt(1-fc^2/omega^2);    
%   U0      : horizontal barotropic velocity at -\infty
%   h0      : depth at infinity

% OUTPUT:
%   ph      : modal coefficients

hwait = waitbar(0, 'Solving the CMS (non-hydr. barotropic case)...');


Nx = length(x); % the number of grid points, including boundaries 
Nm = NT;
% right-hand side
% right-hand side
Q = U0*(h0);
S = 2*hx.^2./h.^3 - hxx./h.^2;
f = zeros(Nm*Nx,1);
for n = 1:Nm
    gn = (-1)^(n+1)/n/pi;
    f((n-1)*Nx+1:(n)*Nx) = 2*Q*gn*h.*S;
end



% Set-up of the sparse fourth-order finite-difference matrix
dx = abs(x(2)-x(1));
NZ = 5*Nx+2;
J = zeros(NZ,1) ;
K = J ;
X = J ;
NZb = NT*(NT)*NZ;
JJ = zeros(NZb,1) ;
KK = JJ ;
XX = JJ;

J ( 1 ) = 1 ; K( 1 ) = 1 ;
J ( 2 ) = 2 ; K( 2 ) = 1 ;
J ( 3 ) = 3 ; K( 3 ) = 1 ;

J ( 4 ) = 1 ; K( 4 ) = 2 ;
J ( 5 ) = 2 ; K( 5 ) = 2 ;
J ( 6 ) = 3 ; K( 6 ) = 2 ;
J ( 7 ) = 4 ; K( 7 ) = 2 ;

J ( 8 )  = 1 ;  K( 8 ) = 3 ;
J ( 9 )  = 2 ;  K( 9 ) = 3 ;
J ( 10 ) = 3 ; K( 10 ) = 3 ;
J ( 11 ) = 4 ; K( 11 ) = 3 ;
J ( 12 ) = 5 ; K( 12 ) = 3 ;

J ( 13 ) = 1 ; K( 13 ) = 4 ;
J ( 14 ) = 2 ; K( 14 ) = 4 ;
J ( 15 ) = 3 ; K( 15 ) = 4 ;
J ( 16 ) = 4 ; K( 16 ) = 4 ;
J ( 17 ) = 5 ; K( 17 ) = 4 ;
J ( 18 ) = 6 ; K( 18 ) = 4 ;

J ( 19 ) = 1 ; K( 19 ) = 5 ;
J ( 20 ) = 2 ; K( 20 ) = 5 ;
J ( 21 ) = 3 ; K( 21 ) = 5 ;
J ( 22 ) = 4 ; K( 22 ) = 5 ;
J ( 23 ) = 5 ; K( 23 ) = 5 ;
J ( 24 ) = 6 ; K( 24 ) = 5 ;
J ( 25 ) = 7 ; K( 25 ) = 5 ;

J( 26 ) = 2 ;  K( 26 )  = 6 ;

k = 27; i = 6;

J( k )   = 4;  K( k )   = i ;
J( k+1 ) = 5;  K( k+1 ) = i ;
J( k+2 ) = 6 ; K( k+2 ) = i ;
J( k+3 ) = 7;  K( k+3 ) = i ;
J( k+4 ) = 8;  K( k+4 ) = i ;

k = 32 ;
for  i = 7 :Nx-6
    J( k )  = i-2;   K( k ) = i ;
    J( k+1 ) = i-1;  K( k+1 ) = i ;
    J( k+2 ) = i ;   K( k+2 ) = i ;
    J( k+3 ) = i+1;  K( k+3 ) = i ;
    J( k+4 ) = i+2;  K( k+4 ) = i ;
    k = k + 5 ;
end

k = 30; i = Nx-5;
J( NZ-k )  = i-2;        K(NZ - k ) = i ;
J( NZ-(k-1) ) = i-1;     K(NZ - (k-1) ) = i ;
J( NZ - (k-2) ) = i ;    K( NZ -(k-2) ) = i ;
J( NZ -(k-3) ) = i+1;    K( NZ -(k-3) ) = i ;
J( NZ -(k-4) ) = i+2;    K( NZ -(k-4) ) = i ;
J( NZ - 25 ) =   Nx-1 ;  K( NZ - 25 ) = i ;

J ( NZ - 24 ) = Nx-6;   K( NZ-24 ) = Nx-4;
J ( NZ - 23 ) = Nx-5;   K( NZ-23 ) = Nx-4;
J ( NZ - 22 ) = Nx-4 ;  K( NZ -22 ) = Nx-4 ;
J ( NZ - 21 ) = Nx-3 ;  K( NZ - 21 ) = Nx-4 ;
J ( NZ - 20 ) = Nx-2 ;  K( NZ -20  ) = Nx-4 ;
J ( NZ - 19 ) = Nx-1;   K( NZ - 19 ) = Nx-4;
J ( NZ - 18 ) = Nx;     K( NZ - 18 ) = Nx-4;

J ( NZ-17 ) = Nx-5;   K( NZ-17 ) = Nx-3;
J ( NZ-16 ) = Nx-4;   K( NZ-16 ) = Nx-3;
J ( NZ-15 ) = Nx-3;   K( NZ-15 ) = Nx-3;
J ( NZ-14 ) = Nx-2;   K( NZ-14 ) = Nx-3;
J ( NZ-13 ) = Nx-1;   K( NZ-13 ) = Nx-3;
J ( NZ-12 ) = Nx;     K( NZ-12 ) = Nx-3;

J ( NZ-11 ) = Nx-4;   K( NZ-11 ) = Nx-2;
J ( NZ-10 ) = Nx-3;   K( NZ-10 ) = Nx-2;
J ( NZ-9 ) = Nx-2;    K( NZ-9 )  = Nx-2;
J ( NZ-8 ) = Nx-1;    K( NZ-8 )  = Nx-2;
J ( NZ-7 ) = Nx;      K( NZ-7 )  = Nx-2;

J ( NZ-6 ) = Nx-3;   K( NZ-6 ) = Nx-1;
J ( NZ-5 ) = Nx-2;   K( NZ-5 ) = Nx-1;
J ( NZ-4 ) = Nx-1;   K( NZ-4 ) = Nx-1;
J ( NZ-3 ) = Nx;     K( NZ-3 ) = Nx-1;

J ( NZ-2 ) = Nx-2;   K( NZ-2 ) = Nx;
J ( NZ-1 ) = Nx-1;   K( NZ-1 ) = Nx;
J ( NZ ) = Nx;       K( NZ )   = Nx;

As = zeros(Nm,Nm,Nx);
for n=1:Nm
    As(n,n,:) = ones(1,Nx);
end
Bs = zeros(Nm,Nm,Nx);
for n=1:Nm
    for m = 1:Nm
        if m==n
            Bs(n,n,:) =  hx./h;
        else
            Bs(m,n,:) = 4*(-1)^(m+n)*m*n/(m^2-n^2) * hx./h;
        end
    end
end
Cs = zeros(Nm,Nm,Nx);
for n=1:Nm
    for m = 1:Nm
        if m==n
            c = -0.5 -1/3*n^2*pi^2;
            d = 0.5;
            Cs(n,n,:) = c*hx.^2./h.^2 + d*hxx./h - (1/mu0^2)*n^2*pi^2./h.^2;
        else
            c = -4*(-1)^(m+n)*m*n*(m^2+n^2)/(m^2-n^2)^2;
            d =  2*(-1)^(m+n)*m*n/(m^2-n^2);
            Cs(m,n,:) = c*hx.^2./h.^2 + d*hxx./h;
        end
    end
end

pos = 1;
tic
for m = 1:NT
    for n=1:NT
        Ind = (m==n);
        Ass = squeeze(As(m,n,:));
        Bss = squeeze(Bs(m,n,:));
        Css = squeeze(Cs(m,n,:));
        X( 1 ) = ((-25/12)/dx + 0*sqrt(-1)*(n*pi/hL/mu0))*Ind;
        X( 2 ) = (5/6)*Ass(2)/dx/dx - (1/4)*Bss(2)/dx ;
        X( 3 ) = (-(Ass(3))/dx/(dx*12)) + ((Bss(3))/12/dx);

        X( 4 ) = 4/dx*Ind;
        X( 5 ) = (-5/4)*Ass(2)/dx/dx - (5/6)*Bss(2)/dx + Css(2);
        X( 6 ) = ((16*Ass(3))/dx/(dx*12)) - (8*(Bss(3))/12/dx) ;
        X( 7 ) = (-(Ass(4))/dx/(dx*12)) + ((Bss(4))/12/dx);

        X( 8 ) = (-3)/dx*Ind;
        X( 9 ) = -(1/3)*Ass(2)/dx/dx + (3/2)*Bss(2)/dx ;
        X( 10 ) = (-30*(Ass(3))/dx/(dx*12)) + (Css(3));
        X( 11 ) = ((16*Ass(4))/dx/(dx*12)) - (8*(Bss(4))/12/dx) ;
        X( 12 ) = -(Ass(5))/dx/(dx*12) + ((Bss(5))/12/dx);

        X( 13 ) = (4/3)/dx*Ind;
        X( 14 ) = (7/6)*Ass(2)/dx/dx - (1/2)*Bss(2)/dx ;
        X( 15 ) =  16*Ass(3)/dx/(dx*12) + 8*Bss(3)/12/dx  ;
        X( 16 ) = -30*(Ass(4)/dx/(dx*12)) + Css(4);
        X( 17 ) = 16*Ass(5)/dx/(dx*12) - 8*Bss(5)/12/dx ;
        X( 18 ) = -Ass(6)/dx/(dx*12) + Bss(6)/12/dx;

        X( 19 ) = (-1/4)/dx*Ind;
        X( 20 ) = (-1/2)*Ass(2)/dx/dx + (1/12)*Bss(2)/dx ;
        X( 21 ) =  -Ass(3)/(12*dx^2) - Bss(3)/12/dx ;
        X( 22 ) =   16*Ass(4)/(12*dx^2) + 8*(Bss(4))/12/dx  ;
        X( 23 ) =  -30*(Ass(5))/(12*dx^2) + Css(5);
        X( 24 ) =     16*Ass(6)/(12*dx^2) - 8*Bss(6)/12/dx ;
        X( 25 ) =     -Ass(7)/(12*dx^2) + Bss(7)/12/dx;

        X( 26 ) = (1/12)*Ass(2)/dx/dx;

        k = 27; i = 6;

        X( k ) = (-(Ass(i-2))/dx/(dx*12)) - ((Bss(i-2))/12/dx)  ;
        X( k+1 ) = ((16*Ass(i-1))/dx/(dx*12)) + (8*(Bss(i-1))/12/dx)  ;
        X( k+2) = (-30*(Ass(i))/dx/(dx*12)) + (Css(i));
        X( k+3) = ((16*Ass(i+1))/dx/(dx*12)) - (8*(Bss(i+1))/12/dx)  ;
        X( k+4) = (-(Ass(i+2))/dx/(dx*12)) + ((Bss(i+2))/12/dx)  ;

        k = 32 ;
        for  i = 7 :Nx-6
            X( k ) = -Ass(i-2)/dx/(dx*12) - Bss(i-2)/12/dx  ;
            X( k+1 ) = 16*Ass(i-1)/dx/(dx*12) + 8*Bss(i-1)/12/dx  ;
            X( k+2 ) = -30*Ass(i)/dx/(dx*12) + Css(i);
            X( k+3 ) = 16*Ass(i+1)/dx/(dx*12) - (8*(Bss(i+1))/12/dx)  ;
            X( k+4 ) = -(Ass(i+2))/dx/(dx*12) + ((Bss(i+2))/12/dx)  ;
            k = k + 5 ;
        end
        
        k = 30; i = Nx-5;
        X(end- k ) = (-(Ass(i-2))/dx/(dx*12)) - ((Bss(i-2))/12/dx)  ;
        X( end-(k-1) ) = ((16*Ass(i-1))/dx/(dx*12)) + (8*(Bss(i-1))/12/dx)  ;
        X( end - (k-2)) = (-30*(Ass(i))/dx/(dx*12)) + (Css(i));
        X( end -(k-3)) = ((16*Ass(i+1))/dx/(dx*12)) - (8*(Bss(i+1))/12/dx)  ;
        X( end -(k-4)) = (-(Ass(i+2))/dx/(dx*12)) + ((Bss(i+2))/12/dx)  ;
        X( end-25 ) = (1/12)*Ass(Nx-1)/dx/dx;

        X( end-24) = -Ass(Nx-6)/(12*dx^2) - Bss(Nx-6)/12/dx;
        X( end-23 ) = (16*(Ass(Nx-5))/dx/(dx*12)) + (8*(Bss(Nx-5))/12/dx);
        X( end -22 ) =  -30*(Ass(Nx-4))/(12*dx^2) + Css(Nx-4);
        X( end - 21 ) =   16*Ass(Nx-3)/(12*dx^2) - 8*(Bss(Nx-3))/12/dx  ;
        X( end -20  ) =  -Ass(Nx-2)/(12*dx^2) + Bss(Nx-2)/12/dx ;
        X( end-19 ) =  -(1/2)*Ass(Nx-1)/dx/dx - (1/12)*Bss(Nx-1)/dx ;
        X( end-18 ) = (1/4)/dx*Ind;

        X( end-17 ) = (-(Ass(Nx-5))/dx/(dx*12)) - ((Bss(Nx-5))/12/dx);
        X( end-16 ) = ((16*Ass(Nx-4))/dx/(dx*12)) + (8*(Bss(Nx-4))/12/dx) ;
        X( end-15 ) = (-30*(Ass(Nx-3))/dx/(dx*12)) + (Css(Nx-3));
        X( end-14 ) = ((16*Ass(Nx-2))/dx/(dx*12)) - (8*(Bss(Nx-2))/12/dx)  ;
        X( end-13 ) = (7/6)*Ass(Nx-1)/dx/dx + (1/2)*Bss(Nx-1)/dx ;
        X( end-12 ) = (-(4/3)/dx)*Ind;

        X( end-11 ) = (-(Ass(Nx-4))/dx/(dx*12)) - ((Bss(Nx-4))/12/dx);
        X( end-10 ) = ((16*Ass(Nx-3))/dx/(dx*12)) + (8*(Bss(Nx-3))/12/dx) ;
        X( end-9 ) = (-30*(Ass(Nx-2))/dx/(dx*12)) + (Css(Nx-2));
        X( end-8 ) = -(1/3)*Ass(Nx-1)/dx/dx - (3/2)*Bss(Nx-1)/dx ;
        X( end-7 ) = ((3)/dx)*Ind;

        X( end-6 ) = -(Ass(Nx-3)/dx/(dx*12)) - (Bss(Nx-3)/12/dx);
        X( end-5 ) = 16*Ass(Nx-2)/dx/(dx*12) + (8*(Bss(Nx-2))/12/dx) ;
        X( end-4 ) = -(5/4)*Ass(Nx-1)/dx/dx + (5/6)*Bss(Nx-1)/dx + Css(Nx-1);
        X( end-3 ) = (-4)/dx*Ind;

        X( end-2 ) = -(Ass(Nx-2))/dx/(dx*12) - Bss(Nx-2)/12/dx  ;
        X( end-1 ) = (5/6)*Ass(Nx-1)/dx/dx + (1/4)*Bss(Nx-1)/dx ;
        X( end ) = ((25/12)/dx - 0*sqrt(-1)*(n*pi/hR/mu0))*Ind;
        
        JJ( (pos-1)*NZ + 1 : pos*NZ  ) = J +(m-1)*Nx   ;
        KK( (pos-1)*NZ + 1 : pos*NZ   ) = K +(n-1)*Nx ;
        XX( (pos-1)*NZ + 1 : pos*NZ   ) = X;
        pos = pos+1;
    end
end
Am = sparse( JJ ,KK,XX) ;
Bv = f(1:NT*Nx);
for n=1:NT
    Bv((n-1)*Nx+1) = 0;
    Bv((n-1)*Nx+Nx) = 0;
end

% Solve the linear system
Bvs = sparse(Bv);
Xm = Am\Bvs;

% Arrange modal coefficients in vectors
ph = zeros(NT,Nx);
for n=1:NT
    ph(n,:)=Xm(Nx*(n-1)+1:Nx*(n))';
end


delete(hwait)
