function [fx,fxp,fy,fyp,f] = financial_friction
%[fx,fxp,fy,fyp,f] = model computes a log-linear approximation to the  function f for a small open 
%economy with a debt-elastic interest-rate premium. The model is described in ``Real Business Cycles in 
%Emerging Countries?,'' by J. Garcia-Cicco, R. Pancrazi, and M. Uribe (AER forthcoming). . The 
%function f defines  the DSGE model (a p denotes next-period variables) : 
%  E_t f(yp,y,xp,x) =0. 
%
%Inputs: none
%
%Output: Numerical first derivative of f
%
%Calls: anal_deriv.m and anal_deriv_print2f.m. These files are taken from the suite of programs for 
%first- and second-order approximation of DSGE models prepared by S. Schmitt-Grohé and Martin 
%Uribe and available online at www.columbia.edu/~mu2166/2nd_order.htm
%
%(c) Martin Uribe
%
%Date September 2009.

%Define parameters	 
syms BETTA GAMA DELTA ALFA PSSI RHOA RHOG OMEGA PHI RSTAR DY G SIGMAG SIGMAA	

%Define variables 
syms c cp k kp k1 k1p a ap h hp d dp  yy yyp ivv ivvp tb tbp la lap tby tbyp gy gc giv yyback cback ivvback gyp gcp givp yybackp cbackp ivvbackp g gp gback gbackp

%Give functional form for utility, production, and interest-rate premium functions

%Argument in expressions related to capital adjustment costs
ac = kp/k*g-G;
acp = k1p/k1*gp-G;

syms s sp
%Trade balance
tb = yy - c - s - ivv - PHI/2 * (kp/k*g -G)^2*k;

%Interest Rate
syms mu mup
syms YY
%note r is the gross interest rate 
r = RSTAR + PSSI * (exp(dp/YY-DY) - 1) + exp(mu-1)-1; 

%Marginal utility of consumption
syms nu nup
la = nu * (c - 1/OMEGA * h^OMEGA)^(-GAMA);
lap = nup * (cp - 1/OMEGA * hp^OMEGA)^(-GAMA);

%Write equations (e1, e2,...en)
%Note: we take a linear, rather than log-linear, approximation with respect to tb, the trade balance)
e1 = d - tb - dp*g/r;

e2 = -yy + a * k^ALFA *  (g*h)^(1-ALFA);

e3 = -ivv + kp*g - (1-DELTA) *k;

e4 = - la  + BETTA / g^GAMA * r * lap;

syms ETA
e5 = - h^(OMEGA-1) + (1-ALFA) * a * g^(1-ALFA) * (k/h)^ALFA / (1+ETA*(r-1)/r);

e6 = -la * (1+ PHI * ac) + BETTA/ g^GAMA * lap * (1 - DELTA + ALFA * ap * (gp*hp/kp)^(1-ALFA) + PHI * k1p/k1*gp * acp - PHI/2 * acp^2 );

e7 = -k1 + kp;

e8 = -log(ap) + RHOA * log(a); 

e9 = -log(gp/G) + RHOG * log(g/G); 

e10 = -tby + tb / yy; 

e11 = -gy + yy/yyback*gback;

e12 = -gc + c/cback*gback;

e13 = -giv + ivv/ivvback*gback;

e14 = -yybackp + yy;

e15 = -cbackp + c;

e16 = -ivvbackp + ivv;

e17 = -gbackp + g;

%Preference shock
syms RHONU
e18 = -log(nup) + RHONU * log(nu);

%Interest rate shock 
syms RHOMU
e19 = -log(mup) + RHOMU * log(mu);

%Speinding shock
syms S RHOS
e20 = -log(sp/S) + RHOS * log(s/S);

%Create function f
f = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20];


% Define the vector of controls, y, and states, x
statevar_cu=[yyback cback ivvback gback k d nu mu s g a];
statevar_cup = [yybackp cbackp ivvbackp gbackp kp dp nup mup sp gp ap];

controlvar_cu = [gy gc giv tby h yy c ivv k1];
controlvar_cup = [gyp gcp givp tbyp hp yyp cp ivvp k1p];

log_linearize_cu = [statevar_cu gy gc giv h yy c ivv k1];
log_linearize_cup = [statevar_cup gyp gcp givp hp yyp cp ivvp k1p];

%variables to substitute from levels to logs
log_linearize = [log_linearize_cu log_linearize_cup];

f = subs(f, log_linearize, exp(log_linearize));

%Compute analytical derivatives of f
[fx,fxp,fy,fyp]=anal_deriv(f,statevar_cu,controlvar_cu,statevar_cup,controlvar_cup);

%Make f and its derivatives a function of the level of its arguments rather than the log
f = subs(f, log_linearize, log(log_linearize));
fx = subs(fx, log_linearize, log(log_linearize));
fy = subs(fy, log_linearize, log(log_linearize));
fxp = subs(fxp, log_linearize, log(log_linearize));
fyp = subs(fyp, log_linearize, log(log_linearize));

%Eliminate future variables
cu = [statevar_cu controlvar_cu];
cup = [statevar_cup controlvar_cup];

f = subs(f, cup,cu);
fx = subs(fx, cup,cu);
fy = subs(fy, cup,cu);
fxp = subs(fxp, cup,cu);
fyp = subs(fyp, cup,cu);

%Construct ETASHOCK matrix: the vector innov indicates which rows in the system statevar_cup = hx * statevar_cu + ETAMATRIX innovationp carry an innovaition

syms SIGMAG SIGMAA SIGMANU SIGMAMU SIGMAS
ETASHOCK(find(statevar_cu=='g'),1)=SIGMAG;
ETASHOCK(find(statevar_cu=='a'),2)=SIGMAA;
ETASHOCK(find(statevar_cu=='nu'),3)=SIGMANU;
ETASHOCK(find(statevar_cu=='mu'),4)=SIGMAMU;
ETASHOCK(find(statevar_cu=='s'),5)=SIGMAS;

ns = length(statevar_cu);

if size(ETASHOCK,1)<ns
    ETASHOCK(ns,1) = 0;
end

varshock = ETASHOCK*ETASHOCK';

syms STDmey STDmec STDmeiv STDmetby

varme = diag([STDmey; STDmec; STDmeiv; STDmetby].^2);

%Print derivatives to file `model_num_eval.m'  for model evaluation
filename = 'financial_friction';
anal_deriv_print2f(filename,fx,fxp,fy,fyp,f,ETASHOCK, varme);

save financial_friction.mat 