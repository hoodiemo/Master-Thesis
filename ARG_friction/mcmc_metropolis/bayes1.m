function bayes1(std_param0,param0,k,LB,UB)
%bayes1(std_param0,param0,k,LB,UB);
%computes the value of k that delivers an 
%acceptance rate of about 25 percent. The parameter k is  a scale factor of the variance of the stand in distribution.  
%(c) M. Uribe and S. Schmitt-Grohe, January 2014. 

format compact

%reset the seed of the random number generator
s=sum(100*clock);
randn('state',s)

%Number of draws
ndraw = 1e+6;

%Number of estimated parameters
np = length(param0);

%load data
Y1 = xlsread('argentina_data.xls');
%Note for the authors: argentina_data.xls was created by running mat2xls.m 

T = size(Y1,1); %number of observations

%Demean the data
Y = Y1 - ones(T,1)*mean(Y1);

pick_variables = [1 2 3 4]; %picks the elements of the policy function (rows of gx) corresponding to the observables.

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

mean_param = 0;
var_param = 0;
accept_rate = 0;
mean_k = k;

%start drawing from posterior distribution
for j=1:ndraw

    exitflag = 4;

    % counter to see status
    if rem(j, 10000) == 0
        fprintf('bayes1 percent done: %f \n', j/ndraw);
    end
    
    %Draw parameter value
    param = param0 + k * std_param0 * randn(np,1);
    
    %Make sure the draw is within bounds
    
    if (param>=LB) & (param <= UB) 
    
        %Compute derivatives of eq'm conditions
        [nfx, nfy, nfxp, nfyp, nvarshock, nvarme] = gx_hx_inputs(param); 
        
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
    
    mean_param = mean_param * (j-1)/j + param0/j;
    
    var_param = var_param * (j-1)/j + param0*param0'/j;
    
    accept_rate = accept_rate * (j-1)/j + (exitflag==5)/j;
    
    mean_k = mean_k*(j-1)/j + k/j;
    
    %update the scaling factor of the variance of the stand-in distribution
    k = k * (1 + 0.01 * ((exitflag==5)-0.25));

end %for j=1:ndraw

var_param = var_param - mean_param*mean_param';

std_param = chol(var_param,'lower');

eval(['save bayes1.mat  mean_param std_param0 var_param std_param accept_rate mean_k k'])