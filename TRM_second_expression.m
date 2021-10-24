function [output] = TRM_second_expression()
syms SWin LWin Ta q albedo ra rs emis rhoa Ps G

% constants
globalconstant= getConstants();

%% deal with ground heat flux and its derivatives with respect to albedo, ra and rs
Rn_star = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*(Ta).^4; % potential to acutual temperature

%% important intermediate variables
delta   = desdT(Ta);
gamma = (globalconstant.cp.*Ps)./(globalconstant.epsilon.*globalconstant.Lv);
qsat    = qs(Ta, Ps);
lambda_o_prime = 1./(rhoa*globalconstant.cp);
ro =  rhoa*globalconstant.cp./(4*emis.*globalconstant.sb.*Ta.^3);
f_TRM = 1./ro+(1./ra).*(1 + (delta./gamma).*(ra./(ra + rs)));
AA = rhoa.*globalconstant.Lv.*(qsat-q);

%% computation of Ts, H, LE, and energy balance based on second-order Taylor expansion SEB (this is from Paw U 1987, J. Therm. Biol) 
a = lambda_o_prime.*6.*emis.*globalconstant.sb.*Ta.^2 + 1/2*(des2dT2(Ta)./gamma).*(1./(ra+rs)); 
b = f_TRM;
c = -lambda_o_prime.*(Rn_star - G -AA./(ra+rs));

Ts = Ta + (-b + sqrt(b.^2 - 4*a.*c))./(2*a);
H = rhoa.*globalconstant.cp.*(Ts - Ta)./ra;
LE = rhoa.*globalconstant.Lv.*(qs(Ts, Ps)-q)./(ra+rs); 
energy_balance = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*Ts.^4  - H - LE - G;
Rn = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*Ts.^4;

%% sensitivities 
output.dTs_dALBEDO = diff(Ts,albedo,1);
output.dTs_dRA = diff(Ts,ra,1);
output.dTs_dRS = diff(Ts,rs,1);
output.dTs_dEMIS = diff(Ts,emis,1); 
output.dTs_dSWin = diff(Ts,SWin,1);
output.dTs_dLWin = diff(Ts,LWin,1);
output.dTs_dQA = diff(Ts,q,1);
output.dTs_dTA = diff(Ts,Ta,1);
output.dTs_dG = diff(Ts,G,1);
output.dTs_dPs = diff(Ts,Ps,1);
output.Ts=Ts;
output.H=H;
output.LE=LE;
output.energy_balance=energy_balance;
output.Rn=Rn;
end

function [ qsat ] = qs( T, P )
globalconstant= getConstants();
esat = 611*exp(17.27*(T-273.15)./(T-273.15+237.3));  % Eq. 3.9a from Dingman in Pa
qsat = globalconstant.epsilon*esat./P; % not account for vapor pressure 
end

function delta= desdT(T1)
% T in K, delta is Pa/K
delta = 611.*exp(((1727*T1)/100 - 9434601/2000)./(T1 - 717/20)).*(1727./(100*(T1 - 717/20)) - ((1727*T1)/100 - 9434601/2000)./(T1 - 717/20).^2);
end

function delta2= des2dT2(T)
delta2 = 611*exp(((1727*T)./100 - 9434601/2000)./(T - 717/20)).*(1727./(100*(T - 717/20)) - ...
    ((1727*T)/100 - 9434601/2000)./(T - 717/20).^2).^2 - 611*exp(((1727*T)./100 - ...
    9434601/2000)./(T - 717/20)).*(1727./(50*(T - 717/20).^2) - (2*((1727*T)./100 - 9434601/2000))./(T - 717/20).^3);
end

function globalconstant= getConstants()
globalconstant.cp  = 1004.64;      % specific heat at constant pressure, J/kg/K
globalconstant.Lv  = 2.4665*10^6; % latent heat of vaparization
globalconstant.R   = 287.058;       % dry air gas constant J/kg/K
globalconstant.Rv  = 461.5;       % water vapor gas constant J/kg/K
globalconstant.sb         = 5.670367*10^(-8); % stephan-boltzman constant, W/(m^2 K^4)
globalconstant.epsilon    = globalconstant.R/globalconstant.Rv; 
end






