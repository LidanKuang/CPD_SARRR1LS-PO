 %shiftmax
function [wt1,PP,kk1] = shiftmax(w,ref_new)
wf=fft(w);
Ns=length(w);
f = sqrt(-1)*2*pi*[0:Ns-1]/Ns;
Nf = floor(length(f)/2)+1;
wf = wf(1:Nf);
f = f(1:Nf);
for k=0:length(ref_new)-1
    wft=wf.*exp(k.*f); 
    if mod(length(w),2) == 0
        wft = [wft conj(wft(:,end-1:-1:2))];
    else
        wft = [wft conj(wft(end:-1:2))];
    end
    wt1 = real(ifft(wft,[],2));
    p=corrcoef(wt1,ref_new');
    P(k+1,:)=p(1,2);
end
[PP,kk]=max(abs(P)); PP = P(kk); kk1=kk-1;
wft=wf.*exp((kk-1).*f); 
wft = [wft conj(wft(:,end:-1:2))];
wt1 = real(ifft(wft,[],2));
     