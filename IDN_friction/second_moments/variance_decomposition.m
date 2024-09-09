%M: adapted to only calculate the variance decomposition of the four
%M: control variables of interest
function [Vyr,Vxr,Vy,Vx]=variance_decomposition(gx,hx,ETA1,ETA2)
%[Vy,Vx]=variance_decomposition(gx,hx,ETA1) computes 
%the variance decomposition of y and x, 
%where y and x evolve according to
% x_t+1 = hx x_t + ETA1 epsilon_t+1
% y_t = gx x_t
% epsilon_t ~ iidN(0,I)
%[Vy,Vx]=variance_decomposition(gx,hx,ETA1,ETA2) 
%computes the variance decomposition of y and x, 
%where y and x evolve according to
% x_t+1 = hx x_t + ETA1 epsilon_t+1
% y_t = gx x_t + ETA2 mu_t
% epsilon_t ~ iidN(0,I) and mu_t ~ iidN(0,I) 
% E(epsilon_t mu_T) = 0 for all t,T
% (c) Martin Uribe, February 2009

Vy = [];


n1 = size(ETA1,2);

for j=1:n1
    % creates a n1 x n1 (5x5) matrix of zeros and fills the diagonal jth
    % element with 1
    I1 = zeros(n1);
    I1(j,j) = 1;
    
    V1 = ETA1*I1*ETA1';
    
    [sigy,~]=mom(gx,hx,V1);
    % append the variances
    Vy = [Vy;diag(sigy)'];

end % for j

% Only keep the first 4 columns of Vy
Vy = Vy(:, 1:4);
Vyr = Vy*diag(sum(Vy))^(-1);