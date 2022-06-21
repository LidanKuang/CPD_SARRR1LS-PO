function [pS1] = phase_de_ambiguity(cS)
%cS: phase corrected SM
%pS1: phase de-noised SM
pS1=zeros(size(cS)); SS1 = zeros(size(cS));
out_p=angle(cS);
ind1=find(out_p<(pi/4)&out_p>(-pi/4));
SS1(ind1) = 1;
pS1(ind1)=cS(ind1);