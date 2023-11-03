% Given inputs with data points to be estimated denoted by -999
% this will fill them in using the least squares estimate E(y|x). This
% code can also be used to classify a complete or incomplete
% data set, simply by setting the last column of the data set to -999.
% For missing values in this class column the code will fill in the class
% with the highest posterior probability. 
%
% The conditional Mean of class k can be obtained by setting
% the input to be the vector [-999 -999 ... -999 k].
%
% the input data should be in the variable    inputs.
% the output is returned in the variable      cinputs.
 
h=zeros(N,M);

for i=1:M
  Ri = R((i-1)*(D-1)+1:i*(D-1),:);
  Cov=inv(Ri);
  for tp=1:N
    input=inputs(tp,1:D-1);
    z=~(input==-999); lz=sum(z);
    class=inputs(tp,D);
    zc=~(class==-999);
    if (z)
	    t1=(input(z)-Mean(i,z))';
	    Sz=Cov(z,z);
	    dist=t1'*inv(Sz)*t1;
            h(tp,i) = (2*pi)^(-lz/2)/sqrt(det(Sz))*exp(-dist/2);
    else
	    h(tp,i) = 1.0;
    end;
    if (zc) 
	h(tp,i) = mu(i,class)*h(tp,i);
    end;
  end;
end;

h=h+epsi;

for tp=1:N
    h(tp,:)=h(tp,:)/sum(h(tp,:));
end;

cinputs=inputs(:,1:D-1);
classes=zeros(N,1);

% classify the incomplete data

for tp=1:N
  clprob=h(tp,:)*mu; % class probabilities
  [dummy,class]=max(clprob);
  classes(tp)=class; % maximum likelihood class
end;

% complete the data using regression

for tp=1:N
  input=inputs(tp,1:D-1);
  z=~(input==-999); 
  for i=1:M  
    Ri=R((i-1)*(D-1)+1:i*(D-1),:);
    Cov=inv(Ri);   
    if(sum(z))
      Rxx=inv(Cov(z,z));  
    end;
    if(sum(~z))
      if(sum(z))    
	input(~z)=Mean(i,~z)+Cov(~z,z)*Rxx*(input(z)-Mean(i,z))';
      else
	input(~z)=Mean(i,~z)
      end;
    end;
    cinputs(tp,:)=cinputs(tp,:)+h(tp,i)*input;
  end; 
end;

cinputs=[cinputs classes];
