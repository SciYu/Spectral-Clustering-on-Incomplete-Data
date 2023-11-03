% Given inputs with data points to be estimated denoted by -999
% this code will fill them in using the least squares estimate E(y|x). This
% can be used to perform function approximation by setting the output
% variable columns to -999.
%
% the input data should be in the variable    inputs.
% the output is returned in the variable      cinputs.

h=zeros(N,M);

for i=1:M
  Ri = R((i-1)*D+1:i*D,:);
  Cov=inv(Ri);
  for tp=1:N
    input=inputs(tp,:);
    z=~(input==-999); lz=sum(z);
    if(lz)
      t1=(input(z)-Mean(i,z))';
      Sz=Cov(z,z);
      dist=t1'*inv(Sz)*t1;
      h(tp,i) =(2*pi)^(-lz/2)/sqrt(det(Sz))*exp(-dist/2);
    else
      h(tp,i)=1.0;
    end;
  end;
end;

h=h+epsi;

for tp=1:N
    h(tp,:)=h(tp,:)/sum(h(tp,:));
end;

cinputs=zeros(inputs);

for tp=1:N
  input=inputs(tp,:);
  z=~(input==-999); 
  for i=1:M  
    Ri=R((i-1)*D+1:i*D,:);
    Cov=inv(Ri);   
    if(z) 
      Rxx=inv(Cov(z,z));  
      input(~z)=Mean(i,~z)+Cov(~z,z)*Rxx*(input(z)-Mean(i,z))';
    else
      input=Mean(i,:);
    end;
    cinputs(tp,:)=cinputs(tp,:)+h(tp,i)*input;
  end; 
end;

