function [LL,exitflag] =  log_likelihood(param,pick_variables,Y)
%LL =  log_likelihood(param,pick_variables,Y) 
%computes the logarithm of the likelihood of a DSGE model for parameter vector 
%param,   the ata Y, and the vector pick_variable, which picks the  elements of 
%the policy function that correspond to the observables. 
%[LL,exitflag] =  log_likelihood(param,pick_variables,Y)  
%produces the indicator variable exitflag that takes the value 1 if the equilibrium is 
%locally unique and 0 otherwise. 
%(c) S. Schmitt-Grohe and M. Uribe, January 2013. 


%Compute derivatives of eq'm conditions

[nfx, nfy, nfxp, nfyp, nvarshock, nvarme] = gx_hx_inputs(param); 

%Compute policy functions
[gx,hx,exitflag]=gx_hx(nfy,nfx,nfyp,nfxp);

%Make sure eq'm is locally unique
if exitflag==1

    %Pick only policy functions of variables with observable empirical counterparts
    GX = gx(pick_variables,:);
    
    %evaluate log of posterior density
    LL = log_like_kalman(Y, hx, nvarshock, GX', nvarme);

else
    LL = -inf;

end %if exitflag