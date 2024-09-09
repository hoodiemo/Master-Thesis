function [nfx, nfy, nfxp, nfyp, nvarshock, nETASHOCK,tby] = financial_friction_ss(param_estim);
%[nfx, nfy, nfxp, nfyp, nvarshock, nETASHOCK,tby] = financial_friction_ss(param_estim) introduces calibrated parameters, computes the steady state, and evaluates the first derivatives of the equilibrium conditions of the DSGE model with financial frictions presented in ``Real Business Cycles in emerging Countries?,'' by J. Garcia-Cicco, R. Pancrazi, and M. Uribe (AER forthcoming).
%
%Calls: financial_friction_num_eval. This  file is produced by running financial_friction.m
%
%(c)Martin Uribe
%
%Date September 2009

SIGMAG = param_estim(1); %0;
RHOG = param_estim(2);
SIGMAA = param_estim(3); %0;
RHOA = param_estim(4);
SIGMANU = param_estim(5); %0;
RHONU = param_estim(6);
SIGMAS = 0;%param_estim(7); %0;
RHOS = param_estim(8);
SIGMAMU = param_estim(9); %0;
RHOMU = param_estim(10);
PHI = param_estim(11);
PSSI = param_estim(12); %0.001;
ETA = param_estim(13);
STDmey = param_estim(14);
STDmec = param_estim(15);
STDmeiv = param_estim(16);
STDmetby = param_estim(17);


%PSSI = 0.001, %ssg

GAMA = 2; %intertemporal elasticity of substitution

G = 1.01011; %mean growth rate of output in Argentina for the period 1901 to 2023. 

DELTA = 1.03^4-1; %Depreciation rate

ALFA = 0.32; %Capital elasticity of the production function

OMEGA = 1.6; %exponent of labor in utility function

RSTAR = 1.1; %parameter of the interest rate function

SHARE_S = 0.10; %Share of public spending in GDP

tby = 0.005; %trade balance to output ratio 

BETTA = 1/RSTAR * G^GAMA; %World interest rate

r = RSTAR;  %Country interest rate

k_over_gh = ((G^GAMA/BETTA - 1 + DELTA) / ALFA)^(1/(ALFA-1)); %K/(G*H)

h = ((1-ALFA)* G * k_over_gh^ALFA / (1+ETA*(r-1)/r))^(1/(OMEGA-1)); %hours

k = k_over_gh * G * h; %capital

ivv = (G-1+DELTA) * k; %investment

yy = k^ALFA * (h*G)^(1-ALFA); %output

YY = yy; %parameter of the premium function

s = SHARE_S *yy; %Level of detrended government spending

tb = tby * yy; %detrended trade balance

S = s; %parameter of law of motion of the share of government spending

d = -tb/(G/r-1); %net external debt detrended

DY = d/YY; %parameter of interest-rate function, equal to the steady state level of net foreign debt

c = (G/r-1) * d + yy - SHARE_S*yy - ivv; %Consumption

la = (c - 1/OMEGA*h^OMEGA)^(-GAMA); %marginal utility of wealth

k1 = k; %Auxiliary variable for capital

a = 1; %productivity shock 

g = G; %Growth rate of nonstationary productivity shock

%set StSt levels
gc = g;
giv = g;
yyback = yy;
cback = c;
ivvback = ivv;
gback = g;
gy = g;

nu = 1;
mu = 1;

financial_friction_num_eval; %The file financial_friction_num_eval.m is produced by running financial_friction.m
 

nvarshock = nETASHOCK*nETASHOCK';