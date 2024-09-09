% File name: financial_friction_num_eval.m 
% File generated by anal_deriv_print2f.m Date: 05-Sep-2024

nfx=[0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -(gback*yy)/yyback; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -(c*gback)/cback; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -(gback*ivv)/ivvback; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; (gback*yy)/yyback; (c*gback)/cback; (gback*ivv)/ivvback; 0; 0; 0; 0; 0; 0; 0; (PHI*k*(G - g)^2)/2 + PHI*g*k*(G - g); ALFA*a*k*k^(ALFA - 1)*(g*h)^(1 - ALFA); k*(DELTA - 1); 0; -(ALFA*a*g^(1 - ALFA)*k*(k/h)^(ALFA - 1)*(ALFA - 1))/(h*((ETA*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) + 1)); (PHI*g*nu)/(c - h^OMEGA/OMEGA)^GAMA; 0; 0; 0; -((PHI*k*(G - g)^2)/2 + PHI*g*k*(G - g))/yy; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; d; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -nu/(c - h^OMEGA/OMEGA)^GAMA; 0; (nu*(PHI*(G - g) - 1))/(c - h^OMEGA/OMEGA)^GAMA; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; RHONU; 0; 0; (d*g*mu*exp(mu - 1))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1)^2; 0; 0; (BETTA*mu*nu*exp(mu - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); (a*g^(1 - ALFA)*((ETA*mu*exp(mu - 1))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) - (ETA*mu*exp(mu - 1)*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1)^2)*(k/h)^ALFA*(ALFA - 1))/((ETA*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) + 1)^2; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; RHOMU; 0; s; 0; 0; 0; 0; 0; 0; 0; 0; -s/yy; 0; 0; 0; 0; 0; 0; 0; 0; 0; RHOS; - (d*g)/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) - PHI*g*k*(G - g); -(a*g*h*k^ALFA*(ALFA - 1))/(g*h)^ALFA; g*k; -(BETTA*GAMA*g*nu*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1))/(g^(GAMA + 1)*(c - h^OMEGA/OMEGA)^GAMA); (a*g*(k/h)^ALFA*(ALFA - 1)^2)/(g^ALFA*((ETA*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) + 1)); (BETTA*GAMA*g*nu*(DELTA + (PHI*(G - g)^2)/2 - ALFA*a*((g*h)/k)^(1 - ALFA) + PHI*g*(G - g) - 1))/(g^(GAMA + 1)*(c - h^OMEGA/OMEGA)^GAMA) - (PHI*g*nu)/(c - h^OMEGA/OMEGA)^GAMA; 0; 0; RHOG; (PHI*g*k*(G - g))/yy; 0; 0; 0; 0; 0; 0; g; 0; 0; 0; 0; a*k^ALFA*(g*h)^(1 - ALFA); 0; 0; -(a*g^(1 - ALFA)*(k/h)^ALFA*(ALFA - 1))/((ETA*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) + 1); 0; 0; RHOA; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
nfxp=[0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -yyback; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -cback; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -ivvback; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -gback; 0; 0; 0; -PHI*g*k*(G - g); 0; g*k; 0; 0; (ALFA*BETTA*a*g*h*nu*(ALFA - 1))/(g^GAMA*k*(c - h^OMEGA/OMEGA)^GAMA*((g*h)/k)^ALFA) - (PHI*g*nu)/(c - h^OMEGA/OMEGA)^GAMA; k; 0; 0; (PHI*g*k*(G - g))/yy; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; (PSSI*d^2*g*exp(d/YY - DY))/(YY*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1)^2) - (d*g)/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1); 0; 0; (BETTA*PSSI*d*nu*exp(d/YY - DY))/(YY*g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); (a*g^(1 - ALFA)*(k/h)^ALFA*((ETA*PSSI*d*exp(d/YY - DY))/(YY*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1)) - (ETA*PSSI*d*exp(d/YY - DY)*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(YY*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1)^2))*(ALFA - 1))/((ETA*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) + 1)^2; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; (BETTA*nu*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); 0; -(BETTA*nu*(DELTA + (PHI*(G - g)^2)/2 - ALFA*a*((g*h)/k)^(1 - ALFA) + PHI*g*(G - g) - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0; (BETTA*nu*(PHI*g^2 - (ALFA*a*g*h*(ALFA - 1))/(k*((g*h)/k)^ALFA)))/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); 0; 0; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; (ALFA*BETTA*a*nu*((g*h)/k)^(1 - ALFA))/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); 0; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
nfy=[0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -gy; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -gc; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -giv; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -(a*g*h*k^ALFA*(ALFA - 1))/(g*h)^ALFA; 0; -(GAMA*h*h^(OMEGA - 1)*nu)/(c - h^OMEGA/OMEGA)^(GAMA + 1); (ALFA*a*g^(1 - ALFA)*k*(k/h)^(ALFA - 1)*(ALFA - 1))/(h*((ETA*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) + 1)) - h*h^(OMEGA - 2)*(OMEGA - 1); (GAMA*h*h^(OMEGA - 1)*nu*(PHI*(G - g) - 1))/(c - h^OMEGA/OMEGA)^(GAMA + 1); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -yy; -yy; 0; 0; 0; 0; 0; 0; 0; (c + ivv + s - yy + (PHI*k*(G - g)^2)/2)/yy + 1; (gback*yy)/yyback; 0; 0; yy; 0; 0; 0; 0; 0; 0; c; 0; 0; (GAMA*c*nu)/(c - h^OMEGA/OMEGA)^(GAMA + 1); 0; -(GAMA*c*nu*(PHI*(G - g) - 1))/(c - h^OMEGA/OMEGA)^(GAMA + 1); 0; 0; 0; -c/yy; 0; (c*gback)/cback; 0; 0; c; 0; 0; 0; 0; 0; ivv; 0; -ivv; 0; 0; 0; 0; 0; 0; -ivv/yy; 0; 0; (gback*ivv)/ivvback; 0; 0; ivv; 0; 0; 0; 0; 0; 0; 0; 0; 0; -(BETTA*PHI*g^2*nu)/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); -k1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
nfyp=[0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; (BETTA*GAMA*h*h^(OMEGA - 1)*nu*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^(GAMA + 1)); 0; - (BETTA*GAMA*h*h^(OMEGA - 1)*nu*(DELTA + (PHI*(G - g)^2)/2 - ALFA*a*((g*h)/k)^(1 - ALFA) + PHI*g*(G - g) - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^(GAMA + 1)) - (ALFA*BETTA*a*g*h*nu*(ALFA - 1))/(g^GAMA*k*(c - h^OMEGA/OMEGA)^GAMA*((g*h)/k)^ALFA); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -(BETTA*GAMA*c*nu*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^(GAMA + 1)); 0; (BETTA*GAMA*c*nu*(DELTA + (PHI*(G - g)^2)/2 - ALFA*a*((g*h)/k)^(1 - ALFA) + PHI*g*(G - g) - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^(GAMA + 1)); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; (BETTA*PHI*g^2*nu)/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
nf=[c + d + ivv + s - yy - (d*g)/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) + (PHI*k*(G - g)^2)/2; a*k^ALFA*(g*h)^(1 - ALFA) - yy; g*k - ivv + k*(DELTA - 1); (BETTA*nu*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA) - nu/(c - h^OMEGA/OMEGA)^GAMA; - h^(OMEGA - 1) - (a*g^(1 - ALFA)*(k/h)^ALFA*(ALFA - 1))/((ETA*(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 2))/(RSTAR + exp(mu - 1) + PSSI*(exp(d/YY - DY) - 1) - 1) + 1); (nu*(PHI*(G - g) - 1))/(c - h^OMEGA/OMEGA)^GAMA - (BETTA*nu*(DELTA + (PHI*(G - g)^2)/2 - ALFA*a*((g*h)/k)^(1 - ALFA) + PHI*g*(G - g) - 1))/(g^GAMA*(c - h^OMEGA/OMEGA)^GAMA); k - k1; RHOA*log(a) - log(a); RHOG*log(g/G) - log(g/G); - tby - (c + ivv + s - yy + (PHI*k*(G - g)^2)/2)/yy; (gback*yy)/yyback - gy; (c*gback)/cback - gc; (gback*ivv)/ivvback - giv; yy - yyback; c - cback; ivv - ivvback; g - gback; RHONU*log(nu) - log(nu); RHOMU*log(mu) - log(mu); RHOS*log(s/S) - log(s/S)];
nETASHOCK=[0; 0; 0; 0; 0; 0; 0; 0; 0; SIGMAG; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; SIGMAA; 0; 0; 0; 0; 0; 0; SIGMANU; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; SIGMAMU; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; SIGMAS; 0; 0];
nvarme=[STDmey^2; 0; 0; 0; 0; STDmec^2; 0; 0; 0; 0; STDmeiv^2; 0; 0; 0; 0; STDmetby^2];

nfx=reshape(nfx,[20  11]);
nfxp=reshape(nfxp,[20  11]);
nfy=reshape(nfy,[20   9]);
nfyp=reshape(nfyp,[20   9]);
nf=reshape(nf,[20   1]);
nETASHOCK=reshape(nETASHOCK,[11   5]);
nvarme=reshape(nvarme,[4  4]);
