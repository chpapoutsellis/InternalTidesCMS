%*************************************************************************%
%  Simulation of 2D Internal-Tides using the Couple-Mode System obtained  %
%               from the body-forcing formulation                         %
%*************************************************************************%

% In this example, we calculate the IT generated by a Gaussian topography
%  
% Reference: Ch. Papoutsellis, Matthieu Mercier, Nicolas Grisouard,
% Internal tide generation from non-uniform barotropic body forcing
% 
% Author : Christos Papoutsellis 
% E-mail : cpapoutsellis@gmail.com                
% Web    : https://cpapoutsellis.wixsite.com/chpapoutsellis
% GitHub : https://github.com/ChPapoutsellis                


function [C] = ITCMS(epsilon,delta,M,s,varargin)
% Syntax:  IT_CMS(epsilon,delta,M,s  ,plots,rh0,N,fc,omega,U00,h0,par,Nv,trunc)
%
% Calculation of the energy conversion rate for a Gaussian ridge and plot
% of some fields
%
% Inputs:
%   epsilon - criticality
%   delta   - relative height
%   M       - number of modes
%   s       - horizontal discretisation parameter
%   (optional)
%   plots  -  plot fields if plots=1, make a video if plots=2 (default = 0)
%   rho0  - reference density
%   N     - constant stratification
%   fc    - Coriolis frequency
%   omega - tidal frequency
%   U00   - barotropic current at infinity
%   h0    - depth at infinity
%   par   - controls the depth at infinity after truncation of the Gaussian (default = 0.0001)
%   Nv    - number of vertical points for the visualisation (should be odd)
%   trunc - in the singular case, for the visulisation of the solution trunc*M modes are used (default =1) 
%
%
% Examples:
%    subcritical ridge
%    ITCMS(0.8, 0.5, 64, 6)
%    ITCMS(0.8, 0.5, 64, 6,'plots',1) (with plots)
%    ITCMS(0.8, 0.5, 64, 6,'plots',2) (with plots and video)
%    supercitical ridge
%    ITCMS(1.2, 0.5, 128, 12)
%    ITCMS(1.2, 0.5, 128, 12, 'plots',1) (with plots)
%    ITCMS(1.2, 0.5, 128, 12, 'plots',2) (with plots and video)
%
close all

% Set default values for optional input arguments
plots = 0;
rho0  = 1000;
N     = 0.0015; 
fc    = 1*10^(-4);
omega = 2*pi/12.4/3600;
U00   = 0.04;
h0    = 3000;
par   = 0.0001;
Nv    = 101;
trunc = 1;
% Parse optional input arguments
parser = inputParser;
addOptional(parser, 'plots', plots);
addOptional(parser, 'rho0', rho0);
addOptional(parser, 'N', N);
addOptional(parser, 'fc', fc);
addOptional(parser, 'omega', omega);
addOptional(parser, 'U00', U00);
addOptional(parser, 'h0', h0);
addOptional(parser, 'par', par);
addOptional(parser, 'Nv', Nv);
addOptional(parser, 'trunc', trunc);

parse(parser, varargin{:});
plots = parser.Results.plots;
rho0  = parser.Results.rho0;
N     = parser.Results.N;
fc    = parser.Results.fc;
omega = parser.Results.omega;
U00   = parser.Results.U00;
h0    = parser.Results.h0;
par   = parser.Results.par;
Nv    = parser.Results.Nv;

% In case the user does not provide any argument
if nargin==0
    delta = 0.5;
    epsilon = 0.8;
    M = 64;
    s = 6;
    plots = 2;
end

if N<omega
    error('N should be greater than omega');
end
if omega<fc
    error('omega should be greater than fc')
end
%--------------------------------------------------------------------------



% parameter mu
mu     = sqrt(N^2-omega^2)/sqrt(omega^2-fc^2);

% Create the topography and the numerical grid
Lambda = delta*h0; % max height of the topography
L = sqrt(1/exp(1))*Lambda*mu/epsilon;
LM = 2*pi/(M*pi/(h0-Lambda)/mu); % wave length of the last mode over the summit
dx = LM/s;              
prefL = 0;
xLL = -sqrt(2)*L*sqrt(log(Lambda/par)); xRR = -xLL;  
hLL = Lambda.*exp(-xLL.^2./2/L^2); hRR = hLL;
hL = h0-hLL;
hR = h0-hRR;
L1 = 2*pi/(1*pi/(h0)/mu);
x =  xLL-prefL*L1:dx:xRR+prefL*L1; 
Nx = length(x);
% derivatives of the Gaussian profile
hG = (Lambda.*exp(-x.^2/2/L^2));
hGx =   -hG.*x./L^2;
hGxx =  hG.*(-L^2+x.^2)/L^4;
h = h0 - hG;
hx = -hGx; % could be replaced by gradientp(h,dx) if an analytical expression of h is not availiable
hxx = -hGxx; % or gradientp(gradientp(h,dx),dx)

figure('color','w','position',[490  567  554  194])
box on
line(x,-h,'linewidth',1.5)
line(x,zeros(1,Nx),'color','k','linewidth',1.5)
ylim([-h0 0.1*h0])
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$z$ (m)','Interpreter','latex')
legend('seabed', 'rigid lid')

% Solve the CMS
[ph] = CMS_solution( M , x , h , hx , hxx , h0 ,h0 , mu , U00 , h0 ) ;

% take the horizontal derivative of the modal coefficients
phx = zeros(M,Nx);
for n=1:M
    phx(n,:) = gradientp(ph(n,:),dx);
end
% Calculate the conversion rate
Q = U00*h0;
[C,Cleft,Cright,Cint] = conversion_rate(ph(1:M,:),phx(1:M,:),x,h,hx,1,Nx,Q,N,fc,h0,omega,rho0);

% Print results
clc
fprintf('*** Gaussian ridge ***.\n');
fprintf('depth at infinity                h0 = %0.2f (m)\n', h0)
fprintf('Barotropic current at infinity   U0 = %0.2f (m/s)\n', U00)
fprintf('height Lambda = %0.2f m  (relative height %0.2f)\n', Lambda, delta)
fprintf('RMS width   L = %0.2f m  (criticality %0.2f) \n' , L, epsilon)
fprintf('Brunt–Väisälä frequency  N = %0.5f (1/s)\n', N)
fprintf('Tidal frequency      omega = %0.5f (1/s)\n', omega)
fprintf('Coriolis frequency      fc = %0.5f (1/s)\n', fc)
fprintf('------------------------------------------------------------------ \n');
fprintf('Number of modes      M = %d \n', M);
fprintf('Hor. discretisation  dx = %0.2f (m)\n', dx);
fprintf('Total conversion rate C = %f (W/m per unit ridge length) \n', C);
fprintf('Conservation rel. error (C+ - C- - Cint)/C = %f \n', (Cright-Cleft-Cint)/C);
save OUTPUT/C.mat C

%%
if plots>0
    %----------------------------------------------------------------------
    % modal reconstruction of the spatial part of the mixed (daggered) fields
    % from the solution of the CMS
    % and the hydrostatic barotropic fields
    if epsilon>1
        Mt = round(trunc*M);
    else
        Mt = M;
    end
    [phid,phidx,phidz, Phi0,Phi0x,Phi0z, Phir1, ~  ] = ...
        phis(ph(1:Mt,:),phx(1:Mt,:),h,hx,hxx,Nv,epsilon,1,omega,fc,Q);
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Solution of the Non-hydrostatic Barotropic CMS
    mu0 = sqrt(1-fc^2/omega^2);
    Nbar = 20;
    phbar_nhCMS = CMS_solution_brtrp( Nbar , x ,h,hx,hxx, h0 , h0 , mu0 , h0, U00) ;
    
    % take the horizontal derivative of the modal coefficients
    phxbar_nhCMS = 0*phbar_nhCMS;
    for n=1:Nbar
        phxbar_nhCMS(n,:) = gradientp(phbar_nhCMS(n,:),dx);
    end
    
    % modal reconstruction of the spatial part of the non-hydrostatic barotropic fields
    [pphinhCMS,pphixnhCMS,pphiznhCMS,forcing, zz  ] = ...
        phis_nhbtrp(phbar_nhCMS(1:Nbar,:),phxbar_nhCMS(1:Nbar,:),h,hx,hxx,Nv);
    %----------------------------------------------------------------------

    % Define various fields and plot some snapshots at time t
    T = 2*pi/omega;
    tind = 0;
    t = tind*T;
    im = sqrt(-1);

    % mixed (daggered) fields (psi^{\dagger}, u^^{\dagger}, w^^{\dagger} etc.)
    psid =  real(phid*exp(-im*omega*t));
    ud   = -real(phidz*exp(-im*omega*t));
    wd   =  real(phidx*exp(-im*omega*t));
    vd   =  fc/omega*real(im*phidz*exp(-im*omega*t));
    bd   =  (N^2/omega)*(imag(phidx*exp(-im*omega*t)));

    % hydrostatic barotropic fields  (Psi^(0), U^(0), W^(0) etc.)
    Psi0 =  real(Phi0*exp(-im*omega*t));
    U0   = -real(Phi0z*exp(-im*omega*t));
    W0   =  real(Phi0x*exp(-im*omega*t));
    V0   =  fc/omega*real(im*Phi0z*exp(-im*omega*t));
    B0   =  (N^2/omega)*(imag(Phi0x*exp(-im*omega*t)));

    % non-hydrostatic barotropic fields  (Psi^r, U^r, W^r etc.)
    Psir =  real(pphinhCMS*exp(-im*omega*t));
    Ur   = -real(pphiznhCMS*exp(-im*omega*t));
    Wr   =  real(pphixnhCMS*exp(-im*omega*t));
    Vr   =  fc/omega*real(im*pphiznhCMS*exp(-im*omega*t));
    Br   =  (N^2/omega)*(imag(pphixnhCMS*exp(-im*omega*t)));

    % purely baroclinic fields  (psi^#, u^#, w^# etc.)
    psi = psid - Psir;
    u   = ud   - Ur;
    v   = vd   - Vr;
    w   = wd   - Wr;
    b   = bd   - Br;



    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                            PLOTS                                       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fontsz = 12; % fontsize
    XX = repmat(x,Nv,1);
    ZZ = zz';


    %                     Baroclinic Stream function
    close all

    figure('color','w','position',[ 10  500  480  240],'Name','Baroclinic Stream function');
    box on
    pcolor(XX,ZZ,psi);
    hold on
    line(x,-h,'color','k')
    line(x,-0*h,'color','k')
    shading interp
    colormap(cmocean('curl', 'pivot',0))
    colorbar
    ylabel('$z$ (m)','interpreter','latex','fontsize',fontsz)
    xlabel('$x$ (m)','interpreter','latex','fontsize',fontsz)
    title( '$\psi^{\#}$ (m$^2$/s)','interpreter','latex','fontsize',fontsz);
    set(gcf,'Paperpositionmode','auto');
    print(gcf,'OUTPUT/psi.png','-dpng','-r600');


    %                      Baroclinic Horizontal velocity
    %--------------------------------------------------------------------------
    figure('color','w','position',[ 520  500  480  240],'Name','Horizontal Velocity');
    box on
    pcolor(XX,ZZ,u);
    hold on
    line(x,-h,'color','k')
    line(x,-0*h,'color','k')
    shading interp
    verm = ver('MATLAB');
    % Check if MATLAB version is R2022a or later
    if str2double(verm.Version) <= 9.12
        % Code for MATLAB R2021a or later
        caxis([-U00 U00])
    else
        clim([-U00 U00])
    end
    colormap(cmocean('curl', 'pivot',0))
    colorbar
    ylabel('$z$ (m)','interpreter','latex','fontsize',fontsz)
    xlabel('$x$ (m)','interpreter','latex','fontsize',fontsz)
    title( '$u^{\#}$ (m/s)','interpreter','latex','fontsize',fontsz);
    set(gcf,'Paperpositionmode','auto');
    print(gcf,'OUTPUT/u.png','-dpng','-r600');

    %                        Baroclinic Vertical velocity
    %--------------------------------------------------------------------------
    figure('color','w','position',[ 1030  500  480  240],'Name','Vertical Velocity');
    box on
    pcolor(XX,ZZ,w);
    hold on
    line(x,-h,'color','k')
    line(x,-0*h,'color','k')
    shading interp
    colormap(cmocean('curl', 'pivot',0))
    colorbar
    ylabel('$z$ (m)','interpreter','latex','fontsize',fontsz)
    xlabel('$x$ (m)','interpreter','latex','fontsize',fontsz)
    title( '$w^{\#}$ (m/s)','interpreter','latex','fontsize',fontsz);

    %                              Non-hydrostatic brtrp response
    %--------------------------------------------------------------------------
    figure('color','w','position', [ 10  170  480  240],'Name','non-hydrost. brtrp stream function');
    box on
    pcolor(XX,ZZ,Psir);
    hold on
    line(x,-h,'color','k')
    line(x,-0*h,'color','k')
    shading interp
    colormap(cmocean('curl', 'pivot',0))
    colorbar
    ylabel('$z$ (m)','interpreter','latex','fontsize',fontsz)
    xlabel('$x$ (m)','interpreter','latex','fontsize',fontsz)
    title( '$\Psi^{r}$ (m$^2$/s)','interpreter','latex','fontsize',fontsz);
    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------------
    %                  Nonhydrostatic Barotropic field                        %

    figure('color','w','position', [ 520  170  480  240],'Name','non-hydrost. brtrp horizontal velocity ');
    box on
    pcolor(XX,ZZ,Ur);
    hold on
    line(x,-h,'color','k')
    line(x,-0*h,'color','k')
    shading interp
    % clim([min(min(ur1CMS)) max(max(ur1CMS))])
    colormap(cmocean('curl', 'pivot',0))
    colorbar
    ylabel('$z$ (m)','interpreter','latex','fontsize',fontsz)
    xlabel('$x$ (m)','interpreter','latex','fontsize',fontsz)
    title( '$U^{r}$ (m/s)','interpreter','latex','fontsize',fontsz);


    figure('color','w','position', [ 1030  170  480  240],'Name','non-hydrost. brtrp vertical velocity ');
    box on
    pcolor(XX,ZZ,Wr);
    hold on
    line(x,-h,'color','k')
    line(x,-0*h,'color','k')
    shading interp
    colormap(cmocean('curl', 'pivot',0))
    colorbar
    ylabel('$z$ (m)','interpreter','latex','fontsize',fontsz)
    xlabel('$x$ (m)','interpreter','latex','fontsize',fontsz)
    title( '$W^{r}$ (m/s)','interpreter','latex','fontsize',fontsz);
end
if plots>1
    figure('color','w','position',[371.4000  336.2000  800  400],'Name','Horizontal Velocity');
    box on
    hv = pcolor(XX,ZZ,u);
    hold on
    line(x,-h,'color','k')
    line(x,-0*h,'color','k')
    shading interp
    colorbar
    ylabel('$z$ (m)','interpreter','latex','fontsize',fontsz)
    xlabel('$x$ (m)','interpreter','latex','fontsize',fontsz)
    title( '$u^{\#}$ (m/s)','interpreter','latex','fontsize',fontsz);
    verm = ver('MATLAB');
    % Check if MATLAB version is R2022a or later
    if str2double(verm.Version) <= 9.12
        % Code for MATLAB R2021a or later
        caxis([-U00 U00])
    else
        clim([-U00 U00])
    end
    colormap(cmocean('curl', 'pivot',0))
    avi_col = VideoWriter('OUTPUT/VIDEO.mp4','MPEG-4');
    avi_col.Quality = 100;
    time = 0:T/40:4*T;
    avi_col.FrameRate = round(length(time)/30);
    open(avi_col);
    for it = 1:1:length(time)
        t = time(it);
        ud   = -real(phidz*exp(-im*omega*t));
        %wd   =  real(phidx*exp(-im*omega*t));
        Ur   = -real(pphiznhCMS*exp(-im*omega*t));
        %Wr   =  real(pphixnhCMS*exp(-im*omega*t));
        u   = ud   - Ur;
        %w   = wd   - Wr;
        set(hv,'Cdata',u)
        drawnow
        frame = getframe(gcf);
        writeVideo(avi_col,frame);

    end
    close(avi_col);
end
% 
% 
% 
