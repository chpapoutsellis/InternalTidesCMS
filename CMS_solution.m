function [ph] = CMS_solution(  Nm , x , h, hx, hxx, hL , hR , mu , U0, h0) 
% Solution of the CMS with a 4th-order finite-difference method
% INPUTS
%   Nm      : Number of modes
%    x      : a vector of uniformly discretised grid points
% h,hx,hxx  : the topography and its derivatives
%  hL,hR    : The depths at the left and right
%   mu      : The parameter \mu = sqrt(N^2-omega^2)/sqrt(omega^2-fc^2);
%   U0      : horizontal barotropic velocity at -\infty
%   h0      : depth at infinity

% OUTPUT:
%   ph      : modal coefficients

hwait = waitbar(0, 'Solving the CMS...');

Nx = length(x);  

% Variable Coefficients of the CMS 
% The coefficients are given in the form:
% A phi_xx + B phi_x + C phi = F
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
            Cs(n,n,:) = c*hx.^2./h.^2 + d*hxx./h + (1/mu^2)*n^2*pi^2./h.^2;
        else
            c = -4*(-1)^(m+n)*m*n*(m^2+n^2)/(m^2-n^2)^2;
            d =  2*(-1)^(m+n)*m*n/(m^2-n^2);
            Cs(m,n,:) = c*hx.^2./h.^2 + d*hxx./h;
        end
    end
end


% right-hand side
Q = U0*(h0);
S = 2*hx.^2./h.^3 - hxx./h.^2;
f = zeros(Nm*Nx,1);
for n = 1:Nm
    gn = (-1)^(n+1)/n/pi;
    f((n-1)*Nx+1:(n)*Nx) = 2*Q*gn*h.*S;
end

% lateral boundary conditions
for n=1:Nm 
    f((n-1)*Nx+1)  = 0; 
    f((n-1)*Nx+Nx) = 0; 
end

% Set-up of the sparse fourth-order finite-difference matrix
dx = abs(x(2)-x(1));
NZ = 5*Nx+2;
J = zeros(NZ,1) ;
K = J ;
X = J ;
NZb = Nm*(Nm)*NZ;
JJ = zeros(NZb,1) ;
KK = JJ;
XX = JJ;
pos = 1;

for m = 1:Nm
    for n=1:Nm
        Ind = (m==n);
        Ass = squeeze(As(m,n,:));
        Bss = squeeze(Bs(m,n,:));
        Css = squeeze(Cs(m,n,:));
        J( 1 ) = 1 ; K( 1 ) = 1 ; X( 1 ) = (-25/12/dx + sqrt(-1)*n*pi/hL/mu)*Ind;
        J( 2 ) = 2 ; K( 2 ) = 1 ; X( 2 ) = 5/6*Ass(2)/dx/dx - 1/4*Bss(2)/dx ;
        J( 3 ) = 3 ; K( 3 ) = 1 ; X( 3 ) = -Ass(3)/dx/(dx*12) + Bss(3)/12/dx;
         
        J( 4 ) = 1 ; K( 4 ) = 2 ; X( 4 ) = 4/dx*Ind;
        J( 5 ) = 2 ; K( 5 ) = 2 ; X( 5 ) = -5/4*Ass(2)/dx/dx - 5/6*Bss(2)/dx + Css(2);
        J( 6 ) = 3 ; K( 6 ) = 2 ; X( 6 ) = 16*Ass(3)/dx/(dx*12) - 8*Bss(3)/12/dx ;
        J( 7 ) = 4 ; K( 7 ) = 2 ; X( 7 ) = -Ass(4)/dx/(dx*12) + Bss(4)/12/dx;
        
        J( 8 )  = 1 ;  K( 8 ) = 3 ; X( 8 )  = -3/dx*Ind;
        J( 9 )  = 2 ;  K( 9 ) = 3 ; X( 9 )  = -1/3*Ass(2)/dx/dx + 3/2*Bss(2)/dx ;
        J( 10 ) = 3 ; K( 10 ) = 3 ; X( 10 ) = -30*Ass(3)/dx/(dx*12) + Css(3);
        J( 11 ) = 4 ; K( 11 ) = 3 ; X( 11 ) = 16*Ass(4)/dx/(dx*12) - 8*Bss(4)/12/dx ;
        J( 12 ) = 5 ; K( 12 ) = 3 ; X( 12 ) = -Ass(5)/dx/(dx*12) + Bss(5)/12/dx;
        
        J( 13 ) = 1 ; K( 13 ) = 4 ; X( 13 ) = 4/3/dx*Ind;
        J( 14 ) = 2 ; K( 14 ) = 4 ; X( 14 ) = 7/6*Ass(2)/dx/dx - 1/2*Bss(2)/dx ;
        J( 15 ) = 3 ; K( 15 ) = 4 ; X( 15 ) =  16*Ass(3)/dx/(dx*12) + 8*Bss(3)/12/dx  ;   
        J( 16 ) = 4 ; K( 16 ) = 4 ; X( 16 ) = -30*Ass(4)/dx/(dx*12) + Css(4);
        J( 17 ) = 5 ; K( 17 ) = 4 ; X( 17 ) =  16*Ass(5)/dx/(dx*12) - 8*Bss(5)/12/dx ;  
        J( 18 ) = 6 ; K( 18 ) = 4 ; X( 18 ) =    -Ass(6)/dx/(dx*12) + Bss(6)/12/dx;
        
        J( 19 ) = 1 ; K( 19 ) = 5 ; X( 19 ) = -1/4/dx*Ind;
        J( 20 ) = 2 ; K( 20 ) = 5 ; X( 20 ) = -1/2*Ass(2)/dx/dx + 1/12*Bss(2)/dx ;
        J( 21 ) = 3 ; K( 21 ) = 5 ; X( 21 ) =     -Ass(3)/(12*dx^2) - Bss(3)/12/dx ;
        J( 22 ) = 4 ; K( 22 ) = 5 ; X( 22 ) =   16*Ass(4)/(12*dx^2) + 8*Bss(4)/12/dx  ; 
        J( 23 ) = 5 ; K( 23 ) = 5 ; X( 23 ) =  -30*Ass(5)/(12*dx^2) + Css(5);
        J( 24 ) = 6 ; K( 24 ) = 5 ; X( 24 ) =   16*Ass(6)/(12*dx^2) - 8*Bss(6)/12/dx ;
        J( 25 ) = 7 ; K( 25 ) = 5 ; X( 25 ) =     -Ass(7)/(12*dx^2) + Bss(7)/12/dx; 
        
        J( 26 ) = 2 ;  K( 26 )  = 6 ;   X( 26 ) = 1/12*Ass(2)/dx/dx;
        
        k = 27; i = 6;
        
        J( k )   = 4;  K( k )   = i ;  X( k )   =   -Ass(i-2)/dx/(dx*12) - Bss(i-2)/12/dx  ;   
        J( k+1 ) = 5;  K( k+1 ) = i ;  X( k+1 ) = 16*Ass(i-1)/dx/(dx*12) + 8*Bss(i-1)/12/dx  ;  
        J( k+2 ) = 6 ; K( k+2 ) = i ;  X( k+2 ) = -30*Ass(i)/dx/(dx*12) + Css(i);
        J( k+3 ) = 7;  K( k+3 ) = i ;  X( k+3 ) =  16*Ass(i+1)/dx/(dx*12) - 8*Bss(i+1)/12/dx  ;  
        J( k+4 ) = 8;  K( k+4 ) = i ;  X( k+4 ) =    -Ass(i+2)/dx/(dx*12) + Bss(i+2)/12/dx  ;   
       
        k = 32 ;
        for  i = 7 :Nx-6
            J( k )   = i-2;  K( k )   = i ; X( k )   =      -Ass(i-2)/dx/(dx*12) - Bss(i-2)/12/dx  ;   
            J( k+1 ) = i-1;  K( k+1 ) = i ; X( k+1 ) =  16*Ass(i-1)/dx/(dx*12) + 8*Bss(i-1)/12/dx  ;  
            J( k+2 ) = i  ;  K( k+2 ) = i ; X( k+2 ) = -30*Ass(i)/dx/(dx*12)   + Css(i);
            J( k+3 ) = i+1;  K( k+3 ) = i ; X( k+3 ) =  16*Ass(i+1)/dx/(dx*12) - 8*Bss(i+1)/12/dx  ;  
            J( k+4 ) = i+2;  K( k+4 ) = i ; X( k+4 ) =    -Ass(i+2)/dx/(dx*12) + Bss(i+2)/12/dx  ;   
            k = k + 5 ;
        end

        k = 30; i = Nx-5;
        J( end-k )     = i-2;   K(end-k )      = i ;  X(end-k )      =    -Ass(i-2)/dx/(dx*12)  -   Bss(i-2)/12/dx  ;   
        J( end-(k-1) ) = i-1;   K(end-(k-1) )  = i ;  X( end-(k-1) ) =  16*Ass(i-1)/dx/(dx*12)  + 8*Bss(i-1)/12/dx  ;  
        J( end-(k-2) ) = i;     K( end-(k-2) ) = i ;  X( end-(k-2) ) = -30*Ass(i)/dx/(dx*12)    +    Cs(m,n,i);
        J( end-(k-3) ) = i+1;   K( end-(k-3) ) = i ;  X( end-(k-3) ) =  16*Ass(i+1)/dx/(dx*12)  -  8*Bss(i+1)/12/dx  ;  
        J( end-(k-4) ) = i+2;   K( end-(k-4) ) = i ;  X( end-(k-4) ) =     -Ass(i+2)/dx/(dx*12) +    Bss(i+2)/12/dx  ;  
        J( end-25 )    = Nx-1;  K( end- 25 )   = i ;  X( end-25 )    = 1/12*Ass(Nx-1)/dx/dx;
         
        J( end-24 ) = Nx-6;   K( end-24 ) = Nx-4;   X( end-24)  = -Ass(Nx-6)/(12*dx^2) - Bs(m,n,Nx-6)/12/dx;
        J( end-23 ) = Nx-5;   K( end-23 ) = Nx-4;   X( end-23 ) = 16*Ass(Nx-5)/dx/(dx*12) + 8*Bss(Nx-5)/12/dx;
        J( end-22 ) = Nx-4 ;  K( end-22 ) = Nx-4;   X( end -22 ) =  -30*Ass(Nx-4)/(12*dx^2) + Css(Nx-4);
        J( end-21 ) = Nx-3 ;  K( end-21 ) = Nx-4;   X( end - 21 ) =   16*Ass(Nx-3)/(12*dx^2) - 8*Bss(Nx-3)/12/dx  ; 
        J( end-20 ) = Nx-2 ;  K( end-20 ) = Nx-4;   X( end -20  ) =  -Ass(Nx-2)/(12*dx^2) + Bss(Nx-2)/12/dx ;
        J( end-19 ) = Nx-1;   K( end-19 ) = Nx-4;   X( end-19 ) =  -1/2*Ass(Nx-1)/dx/dx - 1/12*Bss(Nx-1)/dx ;
        J( end-18 ) = Nx;     K( end-18 ) = Nx-4;   X( end-18 ) = 1/4/dx*Ind;

        J( end-17 ) = Nx-5;   K( end-17 ) = Nx-3;  X( end-17 ) = -Ass(Nx-5)/dx/(dx*12) - Bss(Nx-5)/12/dx;
        J( end-16 ) = Nx-4;   K( end-16 ) = Nx-3;  X( end-16 ) = 16*Ass(Nx-4)/dx/(dx*12) + 8*(Bss(Nx-4))/12/dx ;  
        J( end-15 ) = Nx-3;   K( end-15 ) = Nx-3;  X( end-15 ) = -30*(Ass(Nx-3))/dx/(dx*12) + (Css(Nx-3));
        J( end-14 ) = Nx-2;   K( end-14 ) = Nx-3;  X( end-14 ) = 16*Ass(Nx-2)/dx/(dx*12) - 8*(Bss(Nx-2))/12/dx  ;   
        J( end-13 ) = Nx-1;   K( end-13 ) = Nx-3;  X( end-13 ) = 7/6*Ass(Nx-1)/dx/dx + 1/2*Bss(Nx-1)/dx ;
        J( end-12 ) = Nx;     K( end-12 ) = Nx-3;  X( end-12 ) = -(4/3)/dx*Ind;

        J( end-11 ) = Nx-4;   K( end-11 ) = Nx-2;  X( end-11 ) = -Ass(Nx-4)/dx/(dx*12) - Bss(Nx-4)/12/dx;
        J( end-10 ) = Nx-3;   K( end-10 ) = Nx-2;  X( end-10 ) = 16*Ass(Nx-3)/dx/(dx*12) + 8*(Bss(Nx-3))/12/dx ;
        J( end-9 ) = Nx-2;    K( end-9 )  = Nx-2;  X( end-9 ) = -30*Ass(Nx-2)/dx/(dx*12) + Css(Nx-2);
        J( end-8 ) = Nx-1;    K( end-8 )  = Nx-2;  X( end-8 ) = -1/3*Ass(Nx-1)/dx/dx - 3/2*Bss(Nx-1)/dx ;
        J( end-7 ) = Nx;      K( end-7 )  = Nx-2;  X( end-7 ) = 3/dx*Ind;
         
        J( end-6 ) = Nx-3;   K( end-6 ) = Nx-1;  X( end-6 ) = -Ass(Nx-3)/dx/(dx*12) - Bss(Nx-3)/12/dx;
        J( end-5 ) = Nx-2;   K( end-5 ) = Nx-1;  X( end-5 ) = 16*Ass(Nx-2)/dx/(dx*12) + 8*(Bss(Nx-2))/12/dx ;
        J( end-4 ) = Nx-1;   K( end-4 ) = Nx-1;  X( end-4 ) = -5/4*Ass(Nx-1)/dx/dx + 5/6*Bss(Nx-1)/dx + Cs(m,n,Nx-1);
        J( end-3 ) = Nx;     K( end-3 ) = Nx-1;  X( end-3 ) = -4/dx*Ind;
         
        J( end-2 ) = Nx-2;   K( end-2 ) = Nx;  X( end-2 ) = -Ass(Nx-2)/dx/(dx*12) - Bss(Nx-2)/12/dx  ;  
        J( end-1 ) = Nx-1;   K( end-1 ) = Nx;  X( end-1 ) = 5/6*Ass(Nx-1)/dx/dx + 1/4*Bss(Nx-1)/dx ;
        J( end ) = Nx;       K( end )   = Nx;  X( end ) = (25/12/dx - sqrt(-1)*(n*pi/hR/mu))*Ind;
        
        JJ( (pos-1)*NZ + 1 : pos*NZ ) = J +(m-1)*Nx ; 
        KK( (pos-1)*NZ + 1 : pos*NZ ) = K +(n-1)*Nx ;
        XX( (pos-1)*NZ + 1 : pos*NZ ) = X;
        pos = pos+1;
    end
end
Am = sparse( JJ ,KK,XX) ;
fs = sparse(f);

% Solve the linear system
Xm = Am\fs;

% Arrange modal coefficients in vectors
ph = zeros(Nm,Nx);
for n=1:Nm
    ph(n,:) = Xm(Nx*(n-1)+1:Nx*(n))';
end
delete(hwait)

end

