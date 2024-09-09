%This program evaluates the analytical first derivatives of f numerically. The parameters and steady state values of the arguments of the function f are assumed to be in the workspace. Also, the order of approximation must be in the workspace.
%
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001
%Changed on September 25, 2001 to replace subs with eval and make it no longer a function.


nfx = zeros(size(fx));
nfx(:) = eval(fx(:));

nfxp = zeros(size(fxp));
nfxp(:)= eval(fxp(:));

nfy = zeros(size(fy));
nfy(:) = eval(fy(:));

nfyp = zeros(size(fyp));
nfyp(:)= eval(fyp(:));

nf = zeros(size(f));
nf(:)=eval(f(:));