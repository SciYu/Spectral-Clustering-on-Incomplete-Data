% BATCH EM ALGORITHM
% MODIFIED TO DEAL WITH INCOMPLETE DATA
% Zoubin Ghahramani
% May 21, 1992

% EXPECTATION STEP -- COMPUTE E(Zij)
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

% MAXIMIZATION STEP

% MEANS AND VARIANCES

cinputs=inputs;

for i=1:M
  Cov=zeros(D);
  Ri=R((i-1)*D+1:i*D,:);
  oCov=inv(Ri);   

  % complete data for Gaussian i
  
  for tp=1:N
    input=inputs(tp,:);
    z=~(input==-999); 
    if(z)
      Rxx=inv(oCov(z,z));           
      input(~z)=Mean(i,~z)+(oCov(~z,z)*Rxx*(input(z)-Mean(i,z))')';
    else
      input=Mean(i,:);
    end;
    cinputs(tp,:)=input;
  end;  

  % estimate Mean
  
  Z=sum(h(:,i));
  Mean(i,:)=(h(:,i)/Z)'*cinputs;

  % estimate variance

  for tp=1:N
    cinput=cinputs(tp,:);
    z=~(inputs(tp,:)==-999); 
    Rxx=inv(oCov(z,z));           
    t2=(cinput-Mean(i,:));
    Cov = Cov + h(tp,i)*(t2'*t2); 
    if(sum(~z))
      Cov(~z,~z) = Cov(~z,~z)+h(tp,i)*(oCov(~z,~z) ...
	  - oCov(~z,z) * Rxx * oCov(z,~z));
    end;
  end;
  Cov=Cov/Z;
  
  R((i-1)*D+1:i*D,:) = inv(Cov+epsi*eye(D));
end;



