%RBC.M computes unconditional moments implied by the RBC model in ``Real Business Cycles 
%in Emerging Countries?'' by Javier Garcia-Cicco, Roberto Pancrazi, and Martin Uribe (AER, 
%forthcoming). The parameters defining the stochastic processes of the two technology shocks 
%are assigned values corresponding to the median of their estimated posterior distribution.  This 
%distribution was estimated using Bayesian methods and annual Argentine data on  output 
%growth, consumption growth, investment growth, and the trade-balance-to-output ratio  from 1900 to 
%2005. The other parameters of the model are calibrated. Their values are contained  in the 
%program rbc_ss.m,  which also computes the steady state of the model. The dynamic system and 
%its derivatives are computed in the program growoth.m. 
%
%Calls: rbc.m rbc_ss.m gx_hx.m mom.m num_eval.m 
%
%The auxiliary programs gx_hx.m, num_eval.m, and mom.m are taken from the package for computing first- and second-order 
%accurate dynamics of DSGE models produced by Stephanie Schmitt-Grohe and Martin Uribe, and 
%available online at www.columbia.edu/~mu2166/2nd_order.htm.  
%
%(c) Martin Uribe, September 2009. 

format compact
% recall that g is the non-stationary productivity shock, a the stationary
% and phi capital adjustment cost parameter
%Estimated Structural Parameters (sigma_g, rho_g, sigma_a, rho_a, phi) ->
%for consistency I take median values for all calculations

% b as per 3rd September

b = [
0.0234460287073452 %0
0.870256145108376 %0
0.0313004817075586
0.821737524728437
2.18350460757999
];

%The Linearized Theoretical model
[fx,fxp,fy,fyp,f] = rbc;


%Assign values to remaining parameters and compute steady-state
[GAMA, DELTA, ALFA, PSSI, OMEGA, G, SIGMAG, RHOG, SIGMAA, RHOA, PHI, RSTAR, DY, BETTA, c, cp, h, hp, k, kp, k1, k1p, d, dp, iv, ivp, tb, tbp, la, lap, a, ap, tby, tbyp, g, gp, gy, gyp, gc, gcp, giv, givp, yy, YY, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, gback, gbackp]=rbc_ss(b);
%Numerical evaluation of linearized model
num_eval

%Policy Functions
[gx,hx, deterflag] = gx_hx(nfy,nfx,nfyp,nfxp); 

%MOMENTS

%Construct Variance/Covariance Matrix of Innovations to State Process
varshock = zeros(size(hx,1));
varshock(end-1:end,end-1:end)=[SIGMAG^2 0;0 SIGMAA^2]; 

%Compute Variance/covariance matrix of controls (var_y) and states (var_x)
[var_y,var_x] = mom(gx,hx,varshock,0);


%SECOND MOMENTS OF of gY, gC, gI, and TB/Y
disp('From left to right, the figures correspond to output growth, consumption growth, investment growth and the trade-balance-to-output ratio')

disp('Standard Deviations')
STDEV = ((diag(var_y(1:4,1:4)).^(1/2))*100)'


disp('Correlations with Output Growth')
CGY = (var_y(1:4,1) ./ (diag(var_y(1:4,1:4)).^(1/2)) /var_y(1,1)^(1/2))'


disp('Correlations with  the Trade-Balance-to-GDP Ratio')
CTBY = (var_y(1:4,4) ./ (diag(var_y(1:4,1:4)).^(1/2)) /var_y(4,4)^(1/2))'


%First-order autocovariance
[var_y1,var_x1] = mom(gx,hx,varshock,1);
disp('First-Order Autocorrelations')
autocorre1 =  diag(var_y1)./ diag(var_y);
FOAC = autocorre1(1:4)'

% disp("Remainder shows 2nd until 10th autocorrelation of tby only")
% 
% %Second-order autocovariance
% [var_y2,var_x2] = mom(gx,hx,varshock,2);
% disp('Second-Order Autocorrelations')
% autocorre2 = diag(var_y2)./ diag(var_y);
% SOAC = autocorre2(4)'
% 
% %Third-order autocovariance
% [var_y3,var_x3] = mom(gx,hx,varshock,3);
% disp('Third-Order Autocorrelations')
% autocorre3 = diag(var_y3)./ diag(var_y);
% TOAC = autocorre3(4)'
% 
% %Fourth-order autocovariance
% [var_y4,var_x4] = mom(gx,hx,varshock,4);
% disp('Fourth-Order Autocorrelations')
% autocorre4 = diag(var_y4)./ diag(var_y);
% FFOAC = (autocorre4(4))'
% 
% %Fifth-order autocovariance
% [var_y5,var_x5] = mom(gx,hx,varshock,5);
% disp('Fifth-Order Autocorrelations')
% autocorre5 = diag(var_y5)./ diag(var_y);
% FOAC = autocorre5(4)'
% 
% %Sixth-order autocovariance
% [var_y6,var_x6] = mom(gx,hx,varshock,6);
% disp('Sixth-Order Autocorrelations')
% autocorre6 = diag(var_y6)./ diag(var_y);
% SOAC = autocorre6(4)'
% 
% %Seventh-order autocovariance
% [var_y7,var_x7] = mom(gx,hx,varshock,7);
% disp('Seventh-Order Autocorrelations')
% autocorre7 = diag(var_y7)./ diag(var_y);
% SOAC = (autocorre7(4))'
% 
% %Eighth-order autocovariance
% [var_y8,var_x8] = mom(gx,hx,varshock,8);
% disp('Eighth-Order Autocorrelations')
% autocorre8 = diag(var_y8)./ diag(var_y);
% EOAC = autocorre8(4)'
% 
% %Nineth-order autocovariance
% [var_y9,var_x9] = mom(gx,hx,varshock,9);
% disp('Nineth-Order Autocorrelations')
% autocorre9 = diag(var_y9)./ diag(var_y);
% NOAC = (autocorre9(4))'
% 
% %Tenth-order autocovariance
% [var_y10,var_x10] = mom(gx,hx,varshock,10);
% disp('Tenth-Order Autocorrelations')
% autocorre10 = diag(var_y10)./ diag(var_y);
% TOAC = autocorre10(4)'


