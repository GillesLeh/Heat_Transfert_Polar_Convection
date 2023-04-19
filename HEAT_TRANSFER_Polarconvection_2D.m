clear variables,close all,clc,colormap(jet)
% load('T_init_res.mat')
Animation   =    0;  % 1 makes movie, 0 makes no movie
% physics==================================================================
rin         =    1;                                          % internal radius, earth core radius
rout        =    0.8387*rin;                                 % external radius, mantle radius
disc        =    50;%151;                                        % 155 is cell value of 660km discontinuity in mantle
rho         =    1;                                          % density [kg/m3]
rhog0_up    =    1.0;%0.95;                                          % upper mantle density*gravital
rhog0_dw    =    1.0;%1.05;                                          % lower mantle density*gravital, variation between 1 and 1.1
eta         =    1;                                          % viscosity [Pa*s]
eta_up      =    1;%0.1;                                          % upper mantle viscosity
eta_dw      =    1;                                          % lower mantle viscosity, variation between 1 and 1.0241
alph        =    1;                                          % thermal expansion coefficient [1/K]
beta        =    0;                                          % compressibility [1/Pa]
Ttop        =   -1;                                          % surface mantle temperature
Tbot        =    1;                                          % bottom mantle temperature
rhoCp       =    1;                                          % density times specific heat capacity
lam         =    1e-6;                                       % thermal diffusivity, k=lam/rhoCp--> lam=k*rhoCp
ttime       =    1;                                          % total time
% numerics=================================================================
% angular numerics---------------------------------------------------------
np          = 181;                                           % 101 number of angular cells/section
Lp          = 2*pi + 2*pi/(np-1);                            % angular length plus one cell(for polar bounduaries condition)
dphi        = Lp/(np-1);                                     % angular spacing
% radial numerics----------------------------------------------------------
nr          = 101;                                           % number of radial cells
Lr          = rout;                                          % radial length
drad        = Lr/(nr-1);                                     % radial spacing
% density numerics---------------------------------------------------------
rhog_v                     = ones(nr,np);                    % initial density matrix
rhog_v(round(1:disc),:)    = rhog0_dw;                       % lower half with different density
rhog_v(round(disc:end),:)  = rhog0_up;                       % upper half with 1 density
% vicosity numerics--------------------------------------------------------
% radial and angular stresses----------------------------------------------
eta_M                       = ones(nr,np);                   % initial viscosity matrix, noeuds
eta_M(round(1:disc),:)      = eta_dw;                        % lower half with different density
eta_M(round(disc:end),:)    = eta_up;                        % upper half with 1 density
% shear stress-------------------------------------------------------------
eta_Mrp                     = zeros(nr-1,np-1);              % initial shear matrix, center
% Staggered grid in 2 dimension rad and phi--------------------------------
[radn,phin]= ndgrid((       -Lr/2:drad:(Lr     )/2 ) + (Lr/2 + rin),       -Lp/2:dphi:Lp/2       ); % matrix for noeuds variables calculation  (?T,P,Stressrr,Stresspp,Err,Epp)
[radc,phic]= ndgrid((-(Lr-drad)/2:drad:(Lr-drad)/2 ) + (Lr/2 + rin),-(Lp-dphi)/2:dphi:(Lp-dphi)/2); % matrix for center variables calculation  (div V, Stressrp)
[radr,phir]= ndgrid((-(Lr+drad)/2:drad:(Lr+drad)/2 ) + (Lr/2 + rin),-(Lp     )/2:dphi:(Lp     )/2); % matrix for radial variables calculation  (Vr)
[radp,phip]= ndgrid((-(Lr     )/2:drad:(Lr     )/2 ) + (Lr/2 + rin),-(Lp+dphi)/2:dphi:(Lp+dphi)/2); % matrix for angular variables calculation (Vp)
% initial==================================================================

P      = -rhog_v.*radn;                                                                             % lithostatic pressure (with variable density), noeuds:same number of cells in just in angular direction
taurr  = zeros(nr  ,np  );                                                                          % stressrr, noeuds:same number of cells in the two directions
taupp  = zeros(nr  ,np  );                                                                          % stresspp, noeuds:same number of cells in the two directions
taurp  = zeros(nr-1,np-1);                                                                          % stressrp, center: 1 cell short in both two direction
Vr     = zeros(nr+1,np  );                                                                          % initial radial velocity, radial: 1 cell larger in radial direction and same number for angular
Vp     = zeros(nr  ,np+1);                                                                          % initial angular velocity, angular:  1 cell larger in angular direction and same number for radial
% T      = T_init;                                                                                    % Loading same random T
%T      = 0.1*exp(-(radn - 3/2*rin).^2*50-phin.^2*0.*radn.^2).*sin(phin)  + 1*(rand(nr,np)-0.5);      % initial random temperature: noeuds:same number of cells in the two directions
T      = 0.1*exp(-(radn - 3/2*rin).^2*50-phin.^2*0.*radn.^2).*sin(phin);
T_init =    T;                                          % initial random temperature(same in all simulation)
%save('C:\Users\cinga\Desktop\Università\Travaille de Bachelor\Documenti lavoro\Matlab\images mod\T_init_res','T');
% preprocessing============================================================
dtit   = min(drad,dphi)^2/(max(eta_M(:))/max(rhog_v(:)))/4/2e0;            % timesteps, mechanic part
dtdiff = min(drad,dphi)^2/(lam/rhoCp)/4/2e0;                               % timesteps, initial difference between cells, ? lam/rhoCp = k, thermal part
disp(dtit)
ntime  = 7000;%fix(ttime/dtit);                                                  % round towards zero, timesteps definition
% for mechanics conditions-------------------------------------------------
bet    = 1e-2;                                                             % numerical compressibility not physical
tol    = 0.05;                                                             % tolerance for mechanics solutions
% iteration & postprocessing===============================================
niter  = 1000;                                                             % number of iteration, not time
nout   = 3;                                                                % variable for postprocessing
k      = 1;                                                                % for getframe
finish = 4000;                                                             % for exit in shorter time
% for earth's structure visualisation---------------------------------------
p = nsidedpoly(1000,'Center', [0 0], 'Radius', rin+rout);                  % earth model without crust
q = nsidedpoly(1000,'Center', [0 0], 'Radius', rin);                       % nucleus model
% action===================================================================
if Animation == 1
    vidObj = VideoWriter('Convection_1e7.avi');
    set(vidObj,'FrameRate',8)
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end
for it = 1:ntime                                                           % physical time/real time
    % Buoyancy-------------------------------------------------------------
    rhofg                   = rhog_v.*(1 - alph*T + beta*P);               % state equation, with variable density
    % Mechanics, must converg to0----------------------------------------->
    error                   = 2*tol;
    iter                    = 0;
    while error>tol | iter<300 %iter = 1:niter                                                     % not a time but iteration
        iter                = iter+1;
        Srr                 = -P(:,2:end-1) + taurr(:,2:end-1);            % total stress in radial direction, sigma
        Spp                 = -P(2:end-1,:) + taupp(2:end-1,:);            % total stress in angular direction, sigma
        deltaS              = avr(taurr(:,2:end-1) - taupp(:,2:end-1));
        dVrdt               =  1/rho*( ...
            diff(Srr,1,1)./diff(radn(:,2:end-1),1,1) +...
            (diff(taurp,1,2)./diff(phic,1,2) + deltaS)./radr(2:end-1,2:end-1) ...
            -                 (rhofg(1:end-1,2:end-1) + rhofg(2:end  ,2:end-1))/2 );
        dVpdt               =  1/rho*( ...
            diff(taurp,1,1)./diff(radc,1,1) +...
            (diff(Spp,1,2)./diff(phin(2:end-1,:),1,2) + 2*avr(taurp))./radp(2:end-1,2:end-1) );      %pour coordonnée cylindrique partie avec 2*avr(taurp...
        Vr(2:end-1,2:end-1) = Vr(2:end-1,2:end-1) + dVrdt*dtit;
        Vp(2:end-1,2:end-1) = Vp(2:end-1,2:end-1) + dVpdt*dtit;
        Vr(:,end)           = Vr(:,2);
        Vr(:,1)             = Vr(:,end-1);
        Vp(:,end)           = Vp(:,3);
        Vp(:,1)             = Vp(:,end-2);
        % deviatoric strain rate-------------------------------------------
        Err                 =   diff(Vr           ,1,1)./diff(radr           ,1,1);
        Epp                 =  (diff(Vp           ,1,2)./diff(phip           ,1,2) + avr(Vr           ))./radn;
        Erp                 = ((diff(Vr(2:end-1,:),1,2)./diff(phir(2:end-1,:),1,2) - avr(Vp(:,2:end-1)))./radc ...
            +                   diff(Vp(:,2:end-1),1,1)./diff(radp(:,2:end-1),1,1))/2;
        %------------------------------------------------------------------
        divV                = Err + Epp;                                    % divergence of velocities 
        dPdt                = -1/bet*divV;                                  % (-)1/bet --> bulk modulus
        P                   = P + dPdt*dtit;                                % pression update
        % not incompressible, compressible---------------------------------
        taurr               = 2.*eta_M.* (Err - 1/3*divV);                  % deviatoric normal stresses + normal deformation/compression
        taupp               = 2.*eta_M.* (Epp - 1/3*divV);                  % deviatoric normal stresses + normal deformation/compression

        eta_Mrp             = av_N2C(eta_M);                                % noeuds to center, viscosity matrix for deviatoric shear stress
        taurp               = 2.*eta_Mrp.*Erp;                              % deviatoric shear stresses
        % exit criteria----------------------------------------------------
        %         error   = (max(abs(dVrdt(:))) + max(abs(dVpdt(:))) + max(abs(dPdt(:))));
        error   = max([ max(abs(dVrdt(:))) , max(abs(dVpdt(:))) , max(abs(dPdt(:))) ]);
        %         if mod(iter,500)==0
        %             error
        %         end
    end
    

    % Thermics------------------------------------------------------------->we use the results from mechanicss to calculate thermics
    dt         = min(drad/max(abs(Vr(:)))/2/1e0,dtdiff);                                                           % change of dt value if dt calculated is smaller than standard
    
    % Thermo---------------------------------------------------------------
    dTdt                = -max(0,Vr(2:end-2,2:end-1)).*diff(T(1:end-1,2:end-1),1,1)/drad ...                       % maximum T (not negative parameters)in radial direction (positive)
        -                  min(0,Vr(3:end-1,2:end-1)).*diff(T(2:end  ,2:end-1),1,1)/drad ...                       % minimum T (not negative parameters)in radial direction (negative)
        -                  max(0,Vp(2:end-1,2:end-2)).*diff(T(2:end-1,1:end-1),1,2)/dphi./radn(2:end-1,2:end-1) ...% maximum T (not negative parameters)in angular direction (positive)
        -                  min(0,Vp(2:end-1,3:end-1)).*diff(T(2:end-1,2:end  ),1,2)/dphi./radn(2:end-1,2:end-1) ...% minimum T (not negative parameters)in angular direction (negative)
        +                 (diff(lam*diff(T(:,2:end-1),1,1)/drad,1,1)/drad              ...                         % advection
        +                  (lam*diff(avr(T(:,2:end-1)),1,1)/drad./radn(2:end-1,2:end-1)) ...
        +                  diff(lam*diff(T(2:end-1,:),1,2)/dphi,1,2)/dphi./radn(2:end-1,2:end-1).^(2))/rhoCp;      % diffusion
    
    T(2:end-1,2:end-1)  = T(2:end-1,2:end-1) + dTdt*dt;                                                            % ? update of T values, we change the values of T, variation of T over dt
    % Polar boundaries conditions-------------------------------------------
    T(end,:)            = Ttop;                                                 % reset of the mantle surface T values
    T(1,  :)            = Tbot;                                                 % reset of the mantle bottom T values
    T(:,end)            = T(:,2);                                               % reset of the last angular cell with the second cell (they are superimposed)
    T(:,1)              = T(:,end-1);                                           % reset of the second-to-last angular cell with the first cell (they are superimposed)
    % postprocessing=======================================================
    % polar to cartesian converion-----------------------------------------
    x = radn.*cos(phin);                                                        % conversion to cartesian coordinates x
    y = radn.*sin(phin);                                                        % conversion to cartesian coordinates y
    Vx = Vr(2:end,:).*cos(phir(2:end,:)) - Vp(:,2:end).*sin(phip(:,2:end));     % conversion to cartesian coordinates Vx
    Vy = Vr(2:end,:).*sin(phir(2:end,:)) + Vp(:,2:end).*cos(phip(:,2:end));     % conversion to cartesian coordinates Vy
    % convection visualisation=============================================
    % control visualisation------------------------------------------------
    disc_x = x(disc,:);                                                         % for 660km discordance visualisation
    disc_y = y(disc,:);                                                         % for 660km discordance visualisation
    if it == 1e9                                                                  % to have the control values plotted just one time
        figure(2)
        set(gcf,'position',[140.0000 72.0000 1.3e+03 715.0000])
        % structure visualisation------------------------------------------
        subplot(1,3,1)
        hold on
        plot(p, 'FaceColor', 'b')
        plot(q, 'FaceColor', 'r')
        plot(disc_x,disc_y,'--', 'color', 'k')
        axis image, title('Structure')
        legend('Mantle section   2890 [km]', 'Nucleus section 3446 [km]','Discontinuity      660   [km]','Location','southoutside')
        % density matix visualisation--------------------------------------
        subplot(1,3,2)
        hold on
        pcolor(x,y,rhog_v)
        plot(disc_x,disc_y,'--', 'color', 'w','LineWidth',2)
        axis image,colorbar,shading interp,title('Density[kg/m^{3}]')
        % viscosity matrix visualisation-----------------------------------
        subplot(1,3,3)
        hold on
        pcolor(x,y,eta_M)
        plot(disc_x,disc_y,'--', 'color', 'w', 'LineWidth',2)
        axis image,colorbar,shading interp,title('Viscosity[Pa*s]')
    end
    if mod(it,nout)==0                                                          % modulus after division to avoid rounding error
        % principal visualisation==========================================
        %         figure(1)
        pcolor(x,y,T)
        hold on,st=max(1,fix(nr/20));st2=max(1,fix(nr/80));
        %         quiver( x(1:st:end,1:st2:end), y(1:st:end,1:st2:end) ...                  % plots velocity vectors as arrows with components (u,v) at the points (x,y)
        %             ,  Vx(1:st:end,1:st2:end),Vy(1:st:end,1:st2:end),'w'),
        %         plot(disc_x,disc_y,'--', 'color', 'w', 'linewidth',1),
        hold off
        axis equal,axis([-1 1 -1 1]*1.9),colorbar,shading interp,
        title(['Rayleigh number: ',num2str(1/lam,2),'; Time ',num2str(it)],['Iteration ',num2str(iter)])
        caxis([-1 1 ]/10)
        drawnow
        % for getframe-----------------------------------------------------
        %         M(k) = getframe;                                                        % to get actual frame for rapid visualisation
        %         k    = k+1;                                                             % upload of k for visualisation

        % ANIMATION
        if Animation == 1
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        end

    end
end

if Animation == 1
    close(vidObj);
end
% save MOV_1e7 M

% save('C:\Users\cinga\Desktop\Università\Travaille de Bachelor\Documenti lavoro\Matlab\Images mod\Using T_init Realistic\Res_301X301\Modellizazione_risoluzione','M') % to save visualisation
% shortcut function for average============================================
function U = avr(U)
U = 0.5*(U(1:end-1,:) + U(2:end,:));
end
function U = avp(U)
U = 0.5*(U(:,1:end-1) + U(:,2:end));
end
% shortcut function for average, noeuds to center==========================
function U = av_N2C(U)
U = avp(avr(U));
end