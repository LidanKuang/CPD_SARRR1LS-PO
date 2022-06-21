function [Sf]=updateS(Sf,Xfp,A,T,f,krprtt,krkrt,Nf,Lambda)  
  
  Sfg=zeros(size(Sf));
  for m = 1:Nf(2)
      comp=[1:size(Sf,1)];
      Xft = Xfp(:,m,:);
      F = A(:,comp).*exp(T(:,comp)*f(m));
      krpr = krprod(krprtt(:,comp),F);
      krkr = krkrt(comp,comp).*(F'*F);
      Sf(comp,m) = (krkr+(diag(Lambda(comp)+j*Lambda(comp))))\krpr'*Xft(:);
  end 
  
   