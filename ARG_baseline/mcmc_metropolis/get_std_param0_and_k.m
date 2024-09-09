%get_std_param0_and_k.m
%produces the files param0.xls, std_param0.xls, and k.xls, which are used as inputs in bayes.m, these files contain, respectively:
%param0 = initial value of the vector of estimated parameters. 
%std_param0 = matrix such that std_param0*std_param0' is the variance/covariance matrix of the stand-in distribution. 
%k = scale factor of the stand_in distribution calibrated to deliver an acceptance rate of about 25 percent. 
%The algorithm involves running the program 
%bayes1.m twice and the program bayes.m twice as well. 
%Running time: On a top-of-the-line desktop computer vintage 2011, this program ran 4, when the lengths of the MCMC chains in bayes.m and bayes1.m were set at 1 million.  
%(c) M. Uribe and S. Schmitt-Grohe, January 2014. 

clear all, format compact

%Upper and lower bounds of supports of parameter distributions
D =[     0    0.2000 %sigmag
   -0.9900    0.9900 %rhog
         0    0.2000 %sigmaa
   -0.9900    0.9900 %rhoa
         0    8.0000 %phi
    0.0001    0.0135 %STDmey
    0.0001    0.0184 %STDmec
    0.0001    0.049  %STDmeiv
    0.0001    0.0121 %STDmetby
];
%last 4 elements were obtained as std(Y)'*0.25;  This implies that the upper bound of the prior distribution of the standard  deviation of measurement errors is 25 percent of the standard deviation of the observables.  Equivalently, measurement errors are allowed to explain at most 6.25 percent of the variance of the observable. 

LB = D(:,1); %lower bound
UB = D(:,2);%upper bound

param0 = (LB+UB)/2;

np = length(param0); %length of parameter vector

k = 0.03; %parameter scaling the variance of the stand-in distribution

%iniital standard deviation of stand-in  distribution
std_param0 =  eye(np); 
if exist('std_param0.xls')
    !del std_param0.xls
end
xlswrite('std_param0.xls',std_param0);

disp('running bayes1.m for the first time')
bayes1(std_param0,param0,k,LB,UB);
load bayes1 mean_k mean_param 
k = mean_k;
if exist('k.xls')
    !del k.xls
end
xlswrite('k.xls',k);
param0 = mean_param;
if exist('param0.xls')
    !del param0.xls
end
xlswrite('param0.xls',param0);

disp('running bayes.m for the first time')
bayes
std_param0 = std_param;
if exist('std_param0.xls')
    !del std_param0.xls
end
xlswrite('std_param0.xls',std_param0);
param0 = mean_param;

disp('running bayes1.m for the second time')
bayes1(std_param0,param0,k,LB,UB);
load bayes1 mean_k mean_param 
k = mean_k;
if exist('k.xls')
    !del k.xls
end
xlswrite('k.xls',k);
param0 = mean_param;
if exist('param0.xls')
    !del param0.xls
end
xlswrite('param0.xls',param0);

disp('running bayes.m for the second time.')
bayes
param0 = mean_param;
if exist('param0.xls')
    !del param0.xls
end
xlswrite('param0.xls',param0);