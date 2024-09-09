%M: modified to only use 1st derivative and not determining the size of x
%M: etc.
function [fx,fxp,fy,fyp]=anal_deriv(f,x,y,xp,yp) 
% Computes analytical first order derivatives of the function f(yp,y,xp,x) 
% with respect to x, y, xp, and yp.  
% For documentation, see the paper ``Solving Dynamic General Equilibrium Models Using a Second-Order Approximation to the Policy Function,'' by Stephanie Schmitt-Grohe and Martin Uribe, 2001). 
%
%Inputs: f, x, y, xp, yp
%
%Output: Analytical first derivatives of f. 
%
%Note: This program requires MATLAB's Symbolic Math Toolbox
%
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001

%Compute the first derivatives of f
fx = jacobian(f,x);
fxp = jacobian(f,xp);
fy = jacobian(f,y);
fyp = jacobian(f,yp);

end 