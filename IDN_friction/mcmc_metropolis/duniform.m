function [y,m,s] = duniform(x,a,b);
%y = duniform(x,a,b) evaluates at x the 
%probability density function of the uniform 
%distribution with lower bound a and upper bound b. 
%[y,m,s] = duniform(x,a,b) produeces the mean and standard deviation. 

m = (a+b)/2;
s = 1/6*3^(1/2)*(b-a);
y = 1./(b-a) .* (0*x+1).* ((x>=a)&(x<=b));