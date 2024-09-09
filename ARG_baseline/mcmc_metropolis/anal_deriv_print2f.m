%M: modified to only use 1st order derivative and fixed contents of
%M: resulting matrix -> I use the entirety of S, not S(8:end-1)
%
% anal_deriv_print2f.m
% Converting symbolic objects to an m-file for numeric avaluation.
% inputs: anal_deriv.m (symbolic) output
%         approx: 1 or 2 approximation order
%         filename: ' ...' name of the file to append to '_num_eval'
%         byrow: 1 or 0. If 1 the printed file is with many lines, if 0 viceversa
% output: m-file "filename_num_eval.m" 

% by Andrea Pescatori Jun 18 2008, modified Aug 08
%Modified by SSG on 10/15/2008
function anal_deriv_print2f(filename,fx,fxp,fy,fyp,f,ETASHOCK,varme);

fid=fopen([filename,'_num_eval.m'],'w');
fprintf(fid,'%% File name: %s_num_eval.m \n',filename);
fprintf(fid,'%% File generated by anal_deriv_print2f.m Date: %s\n\n',datetime("today"));

% write all derivatives of state/control vector to file
% turn to characters to manipulate later
S = char(fx(:));
S=['nfx=', S,';\n'];  
% replace all "," with return 
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

S = char(fxp(:));
S=['nfxp=', S,';\n']; 
S=regexprep(S, ',', '\r');
fprintf(fid,S);

S = char(fy(:));
S=['nfy=', S,';\n'];
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

S = char(fyp(:));
S=['nfyp=', S,';\n']; 
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

S = char(f(:));
S=['nf=', S,';\n']; 
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

S = char(ETASHOCK(:));
S=['nETASHOCK=', S,';\n'];  
S=regexprep(S, ',', '\r');
fprintf(fid,S);

S = char(varme(:));
S=['nvarme=', S,';\n\n']; 
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

% write the reshape dimensions for the resulting array to put it back to
% matrix format depending on the size of the derivative matrix

S2F = ['nfx=reshape(nfx,[', num2str(size(fx)),']);\n'];  
fprintf(fid,S2F);

S2F = ['nfxp=reshape(nfxp,[', num2str(size(fxp)),']);\n']; 
fprintf(fid,S2F);

S2F = ['nfy=reshape(nfy,[', num2str(size(fy)),']);\n'];
fprintf(fid,S2F);

S2F = ['nfyp=reshape(nfyp,[', num2str(size(fyp)),']);\n'];  
fprintf(fid,S2F);

S2F = ['nf=reshape(nf,[', num2str(size(f)),']);\n'];  
fprintf(fid,S2F);

S2F = ['nETASHOCK=reshape(nETASHOCK,[', num2str(size(ETASHOCK)),']);\n'];  
fprintf(fid,S2F);

S2F = ['nvarme=reshape(nvarme,[', num2str(size(varme)),']);\n']; 
fprintf(fid,S2F);

fclose(fid);
end