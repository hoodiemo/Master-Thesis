function [GAMA, DELTA, ALFA, PSSI, OMEGA, G, SIGMAG, RHOG, SIGMAA, RHOA, PHI, RSTAR, DY, BETTA, c, cp, h, hp, k, kp, k1, k1p, d, dp, iv, ivp, tb, tbp, la, lap, a, ap, tby, tbyp, g, gp, gy, gyp, gc, gcp, giv, givp, yy, YY, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, gback, gbackp]=rbc_ss(b)
%produces the deep structural parameters and computes the steady state of a 
%small open economy described in the paper ``Real Business Cycles in Emerging Countries?,'' 
%by Javier Garcia-Cicco, Roberto Pancrazi, and Martin Uribe (forthcoming AER). 
%
%Input: vector b of estimated structural parameters
%
%(c)Martin Uribe
%
%September 2009

SIGMAG = b(1); %STD of innovation in permanent technology shock
RHOG = b(2); %Serial correlation of innovation in permanent technology shock
SIGMAA = b(3); %STD of innovation in transitory technology shock
RHOA = b(4); %Serial correlation of transitory technology shock
PHI = b(5); %Adjustment cost parameter

G = 1.01011; %mean growth rate of output in Argentina for the period 1901 to 2023. 

GAMA = 2; %intertemporal elasticity of substitution

DELTA = 1.03^4-1;%Depreciation rate of 12.55%

ALFA = 0.32; %Capital elasticity of the production function

PSSI = 0.001;%parameter governing the debt elasticity of the interest rate. 

OMEGA = 1.6; %exponent of labor in utility function

% in paper we derive the relationship of beta and r* like this:
    % BETTA = 0.98^4;%discount factor 0.9223
    % 
    % RSTAR = 1/BETTA * G^GAMA; %World interest rate 1.107
% here I follow the approach taken in the handbook

BETTA = 0.98^4; %FROM PAPER

RSTAR = 1/BETTA * G^GAMA; %TAKEN FROM PAPER

% RSTAR = 1.1; %parameter of the interest rate function
% 
% BETTA = 1/RSTAR * G^GAMA; %World interest rate (gives 0.9286)

r=RSTAR;  %Country interest rate

k_over_gh = ((G^GAMA/BETTA - 1 + DELTA) / ALFA)^(1/(ALFA-1)); %K/(G*H)

h = ((1-ALFA)* G * k_over_gh^ALFA)^(1/(OMEGA-1)); %hours

k = k_over_gh * G * h; %capital

iv = (G-1+DELTA) * k; %investment

yy = k^ALFA * (h*G)^(1-ALFA); %output

YY = yy; %parameter of the premium function

DY = 0.05; %parameter of interest-rate function, equal to the steady state level of net foreign debt.

d = DY * YY; %net external debt detrended

c = (G/r-1) * d + yy - iv; %Consumption

tb = yy - c - iv; %trade balance

tby = tb / yy; %trade balance-to-output ratio

la = (c - 1/OMEGA*h^OMEGA)^(-GAMA); %marginal utility of wealth

k1 = k; %Auxiliary variable

a = 1; %productivity shock 

g = G; %Growth rate of nonstationary productivity shock


%Log variables
c = log(c);
k = log(k);
k1 = log(k1);
iv = log(iv);
h = log(h);
d = log(d);
la = log(la);
a = log(a);
g = log(g);
yy = log(yy);

%in steady state all growth rates are governed by the non-stationary growth
%rate
gy = g;
gc = g;
giv = g;

%eliminate past and future variables
yyback = yy;
cback = c;
ivback = iv;
gback = g;
cp=c;
kp=k;
k1p=k;
ivp= iv;
hp=h;
dp=d; 
lap=la;
ap=a;
gp = g;
tbp = tb;
tbyp = tby;
gyp=gy;
gcp=gc;
givp=giv;
yyp = yy;
yybackp=yy;
cbackp=cp;
ivbackp=ivp;
gbackp = g;