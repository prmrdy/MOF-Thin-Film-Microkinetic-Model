
clc; clear all; close all; tic

%% Reaction Conditions
T      = 35;                       % Temperature   
c_Zr   = 14.4e-3;                  % Concentration of Zirconium metal cluster (moles/L)
c_BDC  = 14.4e-3;                  % Concentration of H2BDC linker (moles/L)

%% Parameters
k1     = (5.66709e+08*exp(-6.36828e+03/(273.15+T)));        % Rate Constant for initiation  (previous work)
k2     = (1.06928e+25*exp(-1.42615e+04/(273.15+T)));        % Rate Constant for aggregation (previous work)                       
k3     = k2;                                                % Rate Constant for aggregation on film 
ksf    = (1.06928e+25*exp(-1.42615e+04/(273.15+30)));       % Forward rate constant for adsorption on substrate
Ksi    = 15/ksf;                                            % Rate constant for adsorption on substrate

kssf   = 25*(1.06928e+25*exp(-1.42615e+04/(273.15+30)));    % Forward rate constant for adsorption on film
Kssi   = 85/kssf;                                           % Rate constant for adsorption on substrate

alpha  = 9e3;                                               % Autocatalytic constant 
beta   = 10;                                                % Autocatalytic constant
gamma  = 5;                                                 % Autocatalytic constant

n      = 1400;                                              % Maximum possible cluster size
nk     = 700;                                               % Maximum aggregating cluster size                            
i_crit = 4;                                                 % Critical size of nucleus

tlim   = 39; b_new  = 1.2;                                  % Reduced linker order

v_par  = [0.3, 4, 0.07, 0.2, 0.008];                        % Evaporation related constants
stp_sz = 0.1;                                               % Maximum step size for the ODE run
tf     = 400;                                               % Simulation end time
tol_Val = 1e-16;                                            % Tolerance Value for ODE run

%% Initial Values
m      = 6;l = 6;a = 1.2;b = 1.8;
c_DMF  = 10.5446;                  % Concentration of DMF   (moles/L)
vi     = 60e-6;                    % Volume of mother solution (L)
vl     = (60-0.2673-0.0931)*1e-6;  % Volume of DMF+AA in the mothter solution (L)
vs     = vi - vl;                  % Volume of solutes in the mother solution (L)
ti     = 0;

%% Volume calculations
cm3toL = 1e-3;            
dx     = 1;dy=0.01;dz=0.02;dv=dx*dy*dz*cm3toL; % Volume of the differential element (L)
dvl    = dv*vl/vi;                             % Volume of DMSO in the differential element (L)
dvs    = dv - dvl;                             % Volume of solutes in the differential element (L)

M0     = c_Zr *dv;
L0     = c_BDC*dv;

bvf    = v_par(1); 
va     = v_par(2);
vb     = v_par(3);
vc     = v_par(4);
vd     = v_par(5);

%% Other
lc_moles   = min([c_Zr/m c_BDC/l])  *  dv;
V_unit_latt= 20.978^3;                            % Angstorm^3
Avog_No    = 6.02214e23;
no_of_p1   = 4;
A          = dx*dy;
Ap1        = (V_unit_latt/no_of_p1)^(2/3);        % Angstorm^2
Ns0        = dx*dy*1e16/(Avog_No*Ap1);            % nano-moles
all_i      = [1:n, 1:n, 1:n, 1:n, 1:n, 1, 1, 1]';
Np_inp     = [zeros(5*n,1); Ns0; M0; L0];

%%  Store parameters
params = struct('T',T,'k1',k1,'k2',k2,'ksf',ksf,'Ksi',Ksi,'kssf',kssf,'Kssi',Kssi,'alpha',alpha,'beta',beta,'gamma',gamma,'n',n,'nk',nk,...
                'i_crit',i_crit,'tf',tf,'stp_sz',stp_sz, 'm',m,'l',l,'a',a,'b',b,'b_new',b_new,'c_Zr',c_Zr,...
                'c_BDC',c_BDC,'c_DMF',c_DMF,'vi',vi,'vl',vl,'vs',vs,'dx',dx,'dy',dy,'dz',dz,...
                'dv',dv,'dvl',dvl,'dvs',dvs,'M0',M0,'L0',L0,'lc_moles',lc_moles,'Ns0',Ns0,...
                'all_i',all_i,'k3',k3,'Np_inp',Np_inp,'bvf',bvf,...
                'va',va,'vb',vb,'vc',vc,'vd',vd, 'tlim',tlim,'A',A);

%% Loading Experimental Values
cf_exp = load("../../CFvsT_v2.mat").cf35; 
t_exp  = load("../../CFvsT_v2.mat").t35;

opts = odeset('nonnegative',1:5*n+3,'MaxStep',stp_sz,'RelTol',tol_Val,'AbsTol',tol_Val,'OutputFcn',@odeprog);
                                                
%% Solving

[t_out, solN] = ode45(@(t,Np) MOF_rxnschm(t,Np,params,n,nk,k1,k2,k3,ksf,Ksi,kssf,Kssi,a,b,all_i,Ns0,alpha,beta, gamma), ...
                                          ti:1:tf,Np_inp,opts);
Npmoles      = solN(:,1:end-3);
subsN        = solN(:,end-2);
M            = solN(:,end-1);
L            = solN(:,end);


%% Post-Solving
clearvars -except V lc_moles t_out solN subsN Npmoles params sol start_time tic all_i n i_crit cf_exp t_exp dvt M L conF
crit_siz   = [n+i_crit:2*n,2*n+i_crit:3*n, 3*n+i_crit:4*n,4*n+i_crit:5*n];
Np1        = repmat(params.all_i(1:end-3)',[size(Npmoles,1), 1]).*Npmoles;
mof_ads    = sum(Np1(:,crit_siz),2);
cf_ads     = mof_ads/lc_moles;                       %Crystal fraction due to adsorbed growth units   
mof_adsblk = sum(Np1(:,[i_crit:n,crit_siz]),2);
cf_adsblk  = mof_adsblk/lc_moles;                      %Crystal fraction due to bulk and adsorbed growth units 

%% Grain Size 
gSiz = sum(Np1,2)./sum(Npmoles,2);
bta  = 4.3328e-4;
as   = ((3*gSiz)./(bta*sqrt(2))).^(1/3);

%% Crystalline Fraction Plot

fig1 =figure();
p = semilogx(t_out,cf_adsblk,'+',t_out,cf_ads,t_exp,cf_exp,'*','LineWidth',2);
xlim([0 3000])
ylim([0 1.1])
title("T = "+params.T+"^oC")
xlabel("Time (s)",'FontWeight','Bold')
ylabel("Crystalline Fraction",'FontWeight','Bold')
grid on
yyaxis right
semilogx(t_out, as/10)
ylabel("Grain Size (nm)",'FontWeight','Bold')
legend("Theoretical Ads+Blk","Theoretical Ads", "Experimental","Grain Size",'Location','best')


%% Finish
elapsedTime = toc/60;