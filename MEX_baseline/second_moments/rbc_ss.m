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

G =    1.0116; %mean growth rate of output in Mexico for the period 1901 to 2023.

% ETA = 0 as there's no working capital constraint 

GAMA = 2; %intertemporal elasticity of substitution

DELTA = 0.082;%Depreciation rate

ALFA = 0.32; %Capital elasticity of the production function

PSSI = 0.001;%parameter governing the debt elasticity of the interest rate. 

OMEGA = 1.6; %exponent of labor in utility function

BETTA = 0.98^4; %discount factor

RSTAR = 1/BETTA * G^GAMA; %world interest rate

r=RSTAR;  %Country interest rate

k_over_gh = ((G^GAMA/BETTA - 1 + DELTA) / ALFA)^(1/(ALFA-1)); %K/(G*H)

h = ((1-ALFA)* G * k_over_gh^ALFA)^(1/(OMEGA-1)); %hours

k = k_over_gh * G * h; %capital

iv = (G-1+DELTA) * k; %investment

yy = k^ALFA * (h*G)^(1-ALFA); %output

YY = yy; %parameter of the premium function

DY = 0.05; %PAPER %parameter of interest-rate function, equal to the steady state level of net foreign debt.

d = DY * YY; %PAPER  %net external debt detrended

c = (G/r-1) * d + yy - iv; %Consumption

la = (c - 1/OMEGA*h^OMEGA)^(-GAMA); %marginal utility of wealth

tb = yy - c - iv; %trade balance

tby = tb / yy; %trade balance to output ratio

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
gy = g;
gc = g;
giv = g;
yy = log(yy);
yyback = yy;
cback = c;
ivback = iv;
gback = g;

%Next-period variables
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