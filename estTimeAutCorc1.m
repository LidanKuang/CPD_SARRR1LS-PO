function [T,A] = estTimeAutCorc11(Xf,A,B,C,krSf,krf,T,Nf,TauW)
noc = size(A,2); sSf = Nf(2); Sf = fft(B.',[],2);
t1 = randperm(size(A,1));
t2 = randperm(noc);
Br = real(B);  Bi = imag(B); Bfr=fft(Br',[],2);Bfi=fft(Bi',[],2);
Cr = real(C);  Ci = imag(C);
for k = t1   
    Resf = Xf(k,:)-A(k,:)*(krSf.*exp(T(k,:)'*krf));
    for d = t2     
       Resfud = Resf+A(k,d)*(krSf(d,:).*exp(T(k,d)*krf)); 
       Xft = squeeze(unmatricizing(Resfud,1,[1 Nf(2) prod(Nf(3:end))]));
      if Nf(3)>1, Resd = ifft(Xft.',[],2);  else Resd = ifft(Xft,[],2); end
       Rr = real(Resd);  Ri = imag(Resd); Rfr=fft(Rr,[],2);Rfi=fft(Ri,[],2);
       Ar = real(A(k,d));  Ai = imag(A(k,d)); 
       CRrr = conj(Cr(:,d).'*Rfr); CRii = conj(Ci(:,d).'*Rfi); CRri = conj(Cr(:,d).'*Rfi); CRir = conj(Ci(:,d).'*Rfr); 
       Xdf1 = Ar*(CRrr.*Bfr(d,:)-CRir.*Bfi(d,:)+CRii.*Bfr(d,:)+CRri.*Bfi(d,:));
       Xdf2 = Ai*(CRri.*Bfr(d,:)-CRrr.*Bfi(d,:)-CRir.*Bfr(d,:)-CRii.*Bfi(d,:));
       Xdf = Xdf1+Xdf2; 
       Xd = (ifft(Xdf,[],2));  Xd1 = (ifft(Xdf1,[],2)); Xd2 = (ifft(Xdf2,[],2));
       Xd  = Xd .*TauW(d,:);       
       [y,ind] = max(abs(Xd));
       T(k,d) = (ind-sSf)-1;
       A(k,d) = real(Xd1(ind)/(sum((krSf(d,:).*conj(krSf(d,:))))/(sSf)))+...
                 i*real(Xd2(ind)/(sum((krSf(d,:).*conj(krSf(d,:))))/(sSf)));
       if abs(T(k,d))>(sSf/2)
           if T(k,d)>0
               T(k,d) = T(k,d)-sSf;
           else
               T(k,d) = T(k,d)+sSf;
           end
       end
       Resf = Resfud-A(k,d)*(krSf(d,:).*exp(T(k,d)*krf));
   end
end