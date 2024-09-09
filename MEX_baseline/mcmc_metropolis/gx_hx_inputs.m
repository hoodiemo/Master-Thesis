%M: keep the Handbook debt SS, but take tby definition from paper
function [nfx, nfy, nfxp, nfyp, nvarshock, nvarme, nETASHOCK,tby,ALFA] = gx_hx_inputs(param_estimate, iteration);

%
%[nfx, nfy, nfxp, nfyp, nvarshock, nvarme, nETASHOCK,tby,ALFA] = gx_hx_inputs(param_estim) 
% evaluates at the deterministic steady state the equlibrium conditioins, nf, first derivatives, 
% nfx, nfy, nfxp, nfyp,  and variance matrices, nvarshock \equiv nETASHOCK*nETASHOCK' of innovations 
% associated with a DSGE model. The variance matrix of measurement errors varme is also evaluated. 
%
%This particular applicaiton is
%the version of the Garcia-Cicco, Pancrazi, and Uribe 
%(GPU) (AER, 2010) model presented in the book 
%``Open Economy Macroeconomics,'' in the 
%chapter ``Emerging-Market Business Cycles 
%Through the Lens of the SOE-RBC Model.''  The 
%
%input: param_estim is a vector containing parameter values
%
%(c)Martin Uribe, December 2008

if nargin == 1
    iteration = 0;
end

SIGMAG = param_estimate(1);
RHOG = param_estimate(2);
SIGMAA = param_estimate(3);
RHOA = param_estimate(4);
PHI = param_estimate(5);
STDmey = param_estimate(6);
STDmec = param_estimate(7);
STDmeiv = param_estimate(8);
STDmetby = param_estimate(9);


PSSI = 0.001; %debt elasticity 

GAMA = 2; %intertemporal elasticity of substitution

G =    1.016; %mean growth rate of output in Mexico for the period 1901 to 2023. 

DELTA = 0.082; %Depreciation rate

ALFA = 0.32; %Capital elasticity of the production function

OMEGA = 1.6; %exponent of labor in utility function

BETTA = 0.98^4; %discount factor

RSTAR = 1/BETTA * G^GAMA; %world interest rate

r = RSTAR;  %Country interest rate

k_over_gh = ((G^GAMA/BETTA - 1 + DELTA) / ALFA)^(1/(ALFA-1)); %K/(G*H)

h = ((1-ALFA)* G * k_over_gh^ALFA)^(1/(OMEGA-1)); %hours

k = k_over_gh * G * h; %capital

ivv = (G-1+DELTA) * k; %investment

yy = k^ALFA * (h*G)^(1-ALFA); %output

YY = yy; %parameter of the premium function

DY = 0.095; %net external debt detrended

d = DY * YY; %parameter of interest-rate function, equal to the steady state level of net foreign debt

c = (G/r-1) * d + yy - ivv; %Consumption

la = (c - 1/OMEGA*h^OMEGA)^(-GAMA); %marginal utility of wealth

tb = yy - c - ivv; %trade balance

tby = tb / yy; %trade balance to GDP ratio

k1 = k; %Auxiliary variable for capital

a = 1; %productivity shock 

g = G; %Growth rate of nonstationary productivity shock

gc = g;
giv = g;
yyback = yy;
cback = c;
ivvback = ivv;
gback = g;
gy = g;

nu = 1;
mu = 1;

model_num_eval; %this is a .M file produced by running model.m (c) M. Uribe

nvarshock = nETASHOCK*nETASHOCK';

if iteration > 999990
    disp("Evaluating yy ivv betta k/gh h k d DY c")
    eval(["save gx_hx_inputs.mat yy ivv BETTA k_over_gh h k d DY c tby"])
end