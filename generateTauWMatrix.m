function TauWMatrix=generateTauWMatrix(TauW,N2)
    TauWMatrix=zeros(size(TauW,1),N2);
    for d=1:size(TauW,1)
       TauWMatrix(d,1:TauW(d,2))=1;
       TauWMatrix(d,end:-1:end+TauW(d,1))=1;
    end
    