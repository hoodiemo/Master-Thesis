%bayes.m
%Estimates the model with financial frictions presented in ``Open Economy 
%Macroeconomics,'' chapter ``Emerging-Country Business Cycles Through the Lens 
%of the SOE-RBC Model.''
%
%output (save in the file bayes.mat):
%param_chain is a matrix containing 1 million draws from the posterior distribution of the vector of estimated parameters
%mean_param mean of param_chain
%logposterior_chain vector containing 1 million draws from the posterior log lilelihood function
%var_param variance/covariance matrix of the vector of estimated parameters, computed using the 1 million MCMC chain. 
%std_param  lower choleski decomposition of std_param
%std_param0 lower choleski decomposition of the  variance matrix of the stand-in distribution. 
%accept_rate fraction of draws accepted
%exitflag_chain is a 1 million vector of 0s, 1s, 2s, 3s, 4s, or 5s:
%0 the parameter draw delivered no equilibrium in the vicinity of the steady state
%1 the draw was rejected even though the equilibrium was unique (the key feature of the mcmc algorithm)
% 2 the draw was rejected because it implied multiple equilibria
%3 the draw was rejected because it resulted in the message ``invertibility violated.''
%4 the draw was rejected because it fell outside of the support of the prior distribution
%5 the draw was accepted.
%Running time: On a top-of-the-line desktop computer in vintage 2011, this program ran 55 minutes to produce an MCMC of 1 million draws. 
%(c)  M. Uribe and S. Schmitt-Grohe, January 2014. 

clear all, format compact

%reset random generator
s=sum(100*clock);
randn('state',s)

%Number of draws
ndraw = 1e+6;

%scaling factor of var-cov matrix of normal distribution from which parameter values are drawn
k = xlsread('k.xls'); %this xls file is created by running get_std_param0_and_k.m

%load data
Y1 = xlsread('mexico_1980q2_2003q2.xlsx');
%Y1 is a matrix with columns: 
%growth rate of output, 
%growth rate of consumption
%growth rate of investment
%trade-balance-to-GDP ratio

T = size(Y1,1); %number of observations
%Demean the data
Y = Y1 - ones(T,1)*mean(Y1);

%The vector of parameters to be estimated, param,  has the following order: 
%SIGMAG = param(1);
%RHOG = param(2);
%SIGMAA = param(3);
%RHOA = param(4);
%PHI = param(5);
%STDmey = param(6);
%STDmec = param(7);
%STDmeiv = param(8);
%STDmetby = param(9);

%Upper and lower bounds of supports of parameter distributions (Mexico)
% D =[    0    0.2000
%    -0.9900    0.9900
%          0    0.2000
%    -0.9900    0.9900
%          0    8.0000
%     0.0001    0.0147
%     0.0001    0.0206
%     0.0001    0.0655
%     0.0001    0.0141
% ];

D =[     0    0.2000
   -0.9900    0.9900
         0    0.2000
   -0.9900    0.9900
         0    8.0000
    0.0001    0.0038
    0.0001    0.0048
    0.0001    0.0142
    0.0001    0.0099
]; %D 0.25

%last 4 elements were obtained as std(Y)'*0.35;  This implies that the upper bound of the prior distribution of the standard deviation of measurement errors is 35 percent of the standard deviation of the observables.  Equivalently, measurement errors are allowed to explain at most 12.25 percent of the variance of the observable. 

LB = D(:,1); %lower bound
UB = D(:,2);%upper bound

%Initial guess for estimated parameter
param0 = xlsread('param0.xls'); %this xls file is created by running get_std_param0_and_k.m

%Number of estimated parameters
np = length(param0);

%Iniital variance of stand-in normal density
std_param0 = xlsread('std_param0.xls'); %this xls file is created by running get_std_param0_and_k.m

pick_variables = [1 2 3 4]; %picks the elements of the policy funciton (rows of gx) corresponding to the observables.

%Initial Prior likelihood
logprior0 = log(prod( duniform(param0,LB,UB)));
 
%Calculate the initial log likelihood of the model for param0
loglikelihood0 =  log_likelihood(param0,pick_variables,Y);
 
%calculate initial posterior density value
logposterior0 =  loglikelihood0 + logprior0;

%Initialize the vectors used in the chain: 

%Matrix collecting the history of posterior draws
param_chain = zeros(np,ndraw);

%Vector collecting history of log posterior distribution
logposterior_chain = zeros(ndraw,1);

exitflag_chain = zeros(ndraw,1); %Vector collecting record of reason for rejection of draw: 0=no equilibrium, 1=draw rejected; 2=multiple equilibriu; 3=invertibility violated; 4=bound violated; 5=draw accepted

%The following two variables initialize the computatiion of the posterior mean and variance of the vector of estimated parameters
mean_param = 0;
var_param = 0;
accept_rate = 0;

%start drawing from posterior distribution
for j=1:ndraw
    
    % counter to see status
    if rem(j, 10000) == 0
        fprintf('bayes percent done: %f \n', j/ndraw);
    end

    exitflag = 4;
    
    %Draw parameter value
    param = param0 + k * std_param0 * randn(np,1);
    
    %Make sure the draw is within bounds
    
    if (param>=LB) & (param <= UB) 
    
        %Compute derivatives of eq'm conditions
        [nfx, nfy, nfxp, nfyp, nvarshock, nvarme] = gx_hx_inputs(param, j); 
        
        %Compute policy functions
        [gx,hx,exitflag]=gx_hx(nfy,nfx,nfyp,nfxp);
        
        %Make sure eq'm is locally unique
        if exitflag==1
        
            %Pick only policy functions of variables with observable empirical counterparts
            GX = gx(pick_variables,:);
            
            %evaluate log of posterior density
            logposterior = log_like_kalman(Y, hx, nvarshock, GX', nvarme) + logprior0;
            
            %Determine whether to accept or reject the draw
            if rand <= exp(logposterior-logposterior0)
            
                exitflag = 5; 
                
                param0 = param;
                
                logposterior0 = logposterior;
    
            end %if rand
    
        end %if exitflag
    
    end %^if param>=LB ...
    
    param_chain(:,j) = param0;
    
    logposterior_chain(j) = logposterior0;
    
    exitflag_chain(j) = exitflag;
    
    accept_rate =  accept_rate  * (j-1)/j + (exitflag==5)/j;
    
    mean_param = mean_param * (j-1)/j + param0/j;
    
    var_param = var_param * (j-1)/j + param0*param0'/j;
    
end %for j=1:ndraw

var_param = var_param - mean_param*mean_param';

std_param = chol(var_param,'lower');

eval(['save bayes.mat param_chain logposterior_chain logprior0 exitflag_chain mean_param var_param std_param std_param0 accept_rate LB UB k Y1 Y'])