%FINANCIAL_FRICTION_RUN.M computes unconditional moments implied by the financial-
%friction model in the paper  ``Real Business Cycles 
%in Emerging Countries?'' by Javier Garcia-Cicco, Roberto Pancrazi, and Martin Uribe (AER, 
%forthcoming). The estimated parameters 
%are assigned values corresponding to the median of their estimated posterior distribution.  This 
%distribution was estimated using Bayesian methods and annual Argentine data on  output 
%growth, consumption growth, investment growth, and the trade-balance-to-output ratio  from 1900 to 
%2005. The other parameters of the model are calibrated. Their values are contained  in the 
%program financial_friction_ss.m,  which also computes the steady state of the model and 
%evaluates the first derivatives of the equilibrium conditions . The dynamic system and 
%its derivatives are computed analytically  in the program financial_friction.m. 
%
%The auxiliary programs gx_hx.m, num_eval_print2f.m, and mom.m are taken from the package for computing first- and second-order 
%accurate dynamics of DSGE models produced by Stephanie Schmitt-Grohe and Martin Uribe and 
%available online at www.columbia.edu/~mu2166/2nd_order.htm
%
%(c) Martin Uribe, September 2009. 
format compact
%Median of posterior distribution of estimated parameters

% b 3rd September
b =[
    0.00600573639305974
    0.165101144792673
    0.0310785914346065
    0.824777731758842
    0.447608554159049
    0.840818844490291
    0.0471780959100502
    0.537120427881169
    0.0968845301681774
    0.928956390643194
    5.24055474199999
    1.19272506845779
    0.469047836951961
    0.00431555399294389
    0.0101609513494971
    0.0447445105306922
    0.00294153428339901
    ];

financial_friction

% calibrate parameters, compute the steady state, and evaluate the first
% derivatives of the equilibrium conditions
[nfx, nfy, nfxp, nfyp, nvarshock, nETASHOCK] = financial_friction_ss(b);

%Policy Functions
[gx,hx, deterflag] = gx_hx(nfy,nfx,nfyp,nfxp); 


%MOMENTS


%Compute Variance/covariance matrix of controls (var_y) and states (var_x)
[var_y,var_x] = mom(gx,hx,nvarshock,0);

clc

%SECOND MOMENTS OF of gY, gC, gI, and TB/Y
disp('From left to right, the figures correspond to output growth, consumption growth, investment growth and the trade-balance-to-output ratio')

disp('Standard Deviations')
STDEV = ((diag(var_y(1:4,1:4)).^(1/2))*100)'


disp('Correlations with Output Growth')
CGY = (var_y(1:4,1) ./ (diag(var_y(1:4,1:4)).^(1/2)) /var_y(1,1)^(1/2))'


disp('Correlations with  the Trade-Balance-to-GDP Ratio')
CTBY = (var_y(1:4,4) ./ (diag(var_y(1:4,1:4)).^(1/2)) /var_y(4,4)^(1/2))'


%First-order autocovariance
[var_y1,var_x1] = mom(gx,hx,nvarshock,1);
disp('First-Order Autocorrelations')
autocorre1 =  diag(var_y1)./ diag(var_y);
FOAC = autocorre1(1:4)'

% disp("Remainder shows 2nd until 10th autocorrelation of tby only")
% 
% %Second-order autocovariance
% [var_y2,var_x2] = mom(gx,hx,nvarshock,2);
% disp('Second-Order Autocorrelations')
% autocorre2 = diag(var_y2)./ diag(var_y);
% SOAC = autocorre2(4)'
% 
% %Third-order autocovariance
% [var_y3,var_x3] = mom(gx,hx,nvarshock,3);
% disp('Third-Order Autocorrelations')
% autocorre3 = diag(var_y3)./ diag(var_y);
% TOAC = autocorre3(4)'
% 
% %Fourth-order autocovariance
% [var_y4,var_x4] = mom(gx,hx,nvarshock,4);
% disp('Fourth-Order Autocorrelations')
% autocorre4 = diag(var_y4)./ diag(var_y);
% FFOAC = (autocorre4(4))'
% 
% %Fifth-order autocovariance
% [var_y5,var_x5] = mom(gx,hx,nvarshock,5);
% disp('Fifth-Order Autocorrelations')
% autocorre5 = diag(var_y5)./ diag(var_y);
% FOAC = autocorre5(4)'
% 
% %Sixth-order autocovariance
% [var_y6,var_x6] = mom(gx,hx,nvarshock,6);
% disp('Sixth-Order Autocorrelations')
% autocorre6 = diag(var_y6)./ diag(var_y);
% SOAC = autocorre6(4)'
% 
% %Seventh-order autocovariance
% [var_y7,var_x7] = mom(gx,hx,nvarshock,7);
% disp('Seventh-Order Autocorrelations')
% autocorre7 = diag(var_y7)./ diag(var_y);
% SOAC = (autocorre7(4))'
% 
% %Eighth-order autocovariance
% [var_y8,var_x8] = mom(gx,hx,nvarshock,8);
% disp('Eighth-Order Autocorrelations')
% autocorre8 = diag(var_y8)./ diag(var_y);
% EOAC = autocorre8(4)'
% 
% %Nineth-order autocovariance
% [var_y9,var_x9] = mom(gx,hx,nvarshock,9);
% disp('Nineth-Order Autocorrelations')
% autocorre9 = diag(var_y9)./ diag(var_y);
% NOAC = (autocorre9(4))'
% 
% %Tenth-order autocovariance
% [var_y10,var_x10] = mom(gx,hx,nvarshock,10);
% disp('Tenth-Order Autocorrelations')
% autocorre10 = diag(var_y10)./ diag(var_y);
% TOAC = autocorre10(4)'
