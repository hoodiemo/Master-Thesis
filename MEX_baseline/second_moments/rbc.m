function [fx,fxp,fy,fyp,f] = rbc
%[fx,fxp,fy,fyp,f] = rbc computes a log-linear approximation to the  function f for a small 
%open economy with a debt-elastic interest-rate premium. The model is described in the paper 
%``Real Business Cycles in Emerging Countries?,'' by Javier Garcia-Cicco, Roberto 
%Pancrazi, and Martin Uribe. The function f defines  the DSGE model (a p denotes next-
%period variables) : 
%  E_t f(yp,y,xp,x) =0. 
%
%Inputs: none
%
%Output: Numerical first derivative of f
%
%Calls: anal_deriv.m. This file is  taken from the suite of programs for 
%first- and second-order approximation of DSGE models prepared by S. Schmitt-Grohe and Martin 
%Uribe and available online at www.columbia.edu/~mu2166/2nd_order.htm
%
%(c) Martin Uribe
%
%Date September 2009

%Define parameters
syms BETTA GAMA DELTA ALFA PSSI RHOA RHOG OMEGA PHI RSTAR DY G SIGMAG SIGMAA	

%Define variables 
syms c cp k kp k1 k1p a ap h hp d dp  yy yyp iv ivp tb tbp la lap tby tbyp gy gc giv yyback cback ivback gyp gcp givp yybackp cbackp ivbackp g gp gback gbackp

%Give functional form for utility, production, and interest-rate premium functions

%argument in expressions related to capital adjustment costs
ac = kp/k*g-G;
acp = k1p/k1*gp-G;

%Trade balance -> needs to be multiplied by k
tb = yy - c - iv - PHI/2 * (kp/k*g -G)^2*k;

%Interest Rate
syms YY
r = RSTAR + PSSI * (exp(dp/YY-DY) - 1); 

%Martinal utility of consumption
la = (c - 1/OMEGA * h^OMEGA)^(-GAMA);
lap = (cp - 1/OMEGA * hp^OMEGA)^(-GAMA);

%Write equations (e1, e2,...en)
%Note: we take a linear, rather than log-linear, approximation with respect to tb, the trade balance)
e1 = d - tb - dp*g/r;

e2 = -yy + a * k^ALFA *  (g*h)^(1-ALFA);

e3 = -iv + kp*g - (1-DELTA) *k;

e4 = - la  + BETTA / g^GAMA * r * lap;

e5 = -h^(OMEGA-1) + (1-ALFA) * a * g^(1-ALFA) * (k/h)^ALFA;

e6 = -la * (1+ PHI * ac) + BETTA/ g^GAMA * lap * (1 - DELTA + ALFA * ap * (gp*hp/kp)^(1-ALFA) + PHI * k1p/k1*gp * acp - PHI/2 * acp^2 );

e7 = -k1 + kp;

e8 = -log(ap) + RHOA * log(a); 

e9 = -log(gp/G) + RHOG * log(g/G); 

e10 = -tby + tb / yy; 

e11 = -gy + yy/yyback*gback;

e12 = -gc + c/cback*gback;

e13 = -giv + iv/ivback*gback;

e14 = -yybackp + yy;

e15 = -cbackp + c;

e16 = -ivbackp + iv;

e17 = -gbackp + g;

%Create function f
f = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17];

% Define the vector of controls, y, and states, x
x = [yyback cback ivback gback k d g a];
xp = [yybackp cbackp ivbackp gbackp kp dp gp ap];

y = [gy gc giv tby h yy c iv k1];
yp = [gyp gcp givp tbyp hp yyp cp ivp k1p];

%variables to substitute from levels to logs
logvar = [gy gc giv h yy c iv k1 yyback cback ivback gback k d g a];
logvarp = [gyp gcp givp hp yyp cp ivp k1p yybackp cbackp ivbackp gbackp kp dp gp ap];

lv = [logvar logvarp];

%Make f a function of the logarithm of the state and control vector
f = subs(f, lv, exp(lv));

%Compute analytical derivatives of f
[fx,fxp,fy,fyp]=anal_deriv(f,x,y,xp,yp);