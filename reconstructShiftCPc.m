function [Rec,Stt] = reconstructShiftCPc(FACT,T)

% Reconstructs the data from FACT obtained from the shiftCP model
% Written by Morten M鴕up
noc=size(FACT{1},2);
nrmodes=length(FACT);
krprt = ones(1,noc);
krkrt = ones(noc);
for k = 3:nrmodes
    krprt = krprod(FACT{k},krprt);
end
for k=1:nrmodes
    N(k)=size(FACT{k},1);
end
if nrmodes<3
    N(3)=1;
end

Ns = N(2);
f = sqrt(-1)*2*pi*[0:Ns-1]/Ns;
Nf = length(f);
Sf = fft(FACT{2}.',[],2);

Rec=zeros(N(1:3));
Stt=cell(N(1),noc);

for i = 1:N(1)
   Sft = Sf.*exp(T(i,:)'*f); 
   St = ifft(Sft,[],2);
   Stt{i} = St;
   Rec(i,:,:) = St.'*krprod(FACT{1}(i,:),krprt).';
end
