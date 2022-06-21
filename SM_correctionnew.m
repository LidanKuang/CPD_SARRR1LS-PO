function [cS,cS1,pS1,pS11,Rs21,Rs22] = SM_correctionnew(S,Sm2)
inr = find(zscore(Sm2)>2.5);
for m=0:2*512
    cS=S.*exp(sqrt(-1)*(-m*pi/512));
    vn1(m+1) = length(find(abs(angle(cS(inr)))<pi/4));
    [pS1] = phase_de_ambiguity(cS);
    psm = pS1; in = abs(psm)<0; psm(in) = 0; p = corrcoef(abs(psm),Sm2); Rs(m+1) =(p(1,2));
end
[vnm,bbm1]=max(vn1); [rmm,bbm]=max(Rs); 
cS=S.*exp(sqrt(-1)*(-(bbm-1)*pi/512));
[pS1] = phase_de_ambiguity(cS);
psm = pS1; %psm = zscore(pS1);
in = abs(psm)<0; psm(in) = 0; p = corrcoef(abs(psm),Sm2); Rs21 = abs(p(1,2));
cS1=S.*exp(sqrt(-1)*(-(bbm1-1)*pi/512));
[pS11] = phase_de_ambiguity(cS1);
psm = pS11; %psm = zscore(pS1);
in = abs(psm)<0; psm(in) = 0; p = corrcoef(abs(psm),Sm2); Rs22 = abs(p(1,2));