%log_like = log_like_kalman(Y,F,Q,H,R)
%Evaluates the log likelihood function at the sample Y using the Kalman filter. The model is:
% x_t+1 = F x_t + v_t+1
% y_t = H'x_t + w_t 
%var(v_t) = Q
%var(w_t) = R
% x_t is nx by 1
% F is nx by nx
% v_t is nx by 1
% y_t is ny by 1
% H' is ny by nx
% w_t is ny by 1
% Q is nx by nx
% R is ny by ny
%Y is a T by ny matrix containing T observations of y_t
% See Hamilton, page 372
%(c) S. Schmitt_Grohe and M. Uribe, December 11, 2007

function log_like = log_like_kalman(Y,F,Q,H,R)

%Number of Observations
T = size(Y,1);

%Rows of transition matrix
r = size(F,1);
n = size(Y,2);

%Initial forecast
XI_10 = zeros(r,1);

%MSE of initial forecast
%P_10 = zeros(r,r);
%P_10(:) = (eye(r^2) - kron(F,F))\Q(:);
%Doubling algorithm
F_old=F;
Q_old=Q;
P_10_old=eye(size(F));
diferenz=0.1;
while diferenz>1e-25;
P_10 =F_old*P_10_old*F_old' + Q_old;

diferenz = max(max(abs(P_10-P_10_old)));
Q_old=F_old*Q_old*F_old' + Q_old;
F_old = F_old * F_old;
P_10_old=P_10;
end    %while diferenz

XIold = XI_10;
Pold = P_10;
log_like = 0;

tmp_const = n*0.5 * log(2*pi);

Fprime = F';
Hprime = H';
Yprime = Y';

for t=1:T
    y = Yprime(:,t);

    dev_y = y - Hprime*XIold;

    tmp1 = Hprime*Pold;
    var_y = tmp1*H+R;

    ivar_y = var_y \ eye(n);

    log_like = log_like - tmp_const -0.5 * log(det(var_y)) - 0.5 *dev_y'*ivar_y*dev_y;

    tmp_fp = F * Pold;
    tmp2 = tmp_fp * H * ivar_y;
    XIold = F*XIold  + tmp2 *dev_y;

    Pold = (tmp_fp - tmp2 * tmp1) * Fprime + Q;
end