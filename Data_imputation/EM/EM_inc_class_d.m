% BATCH EM ALGORITHM FOR CLASSIFICATION OF REAL DATA
% Zoubin Ghahramani
% December 21, 1993

h=zeros(N,M);

% EXPECTATION STEP

for i=1:M
  Ri = R((i-1)*(D-1)+1:i*(D-1));
  for tp=1:N
    input=inputs(tp,1:D-1);
    z=~(input==-999); lz=sum(z);
    class=inputs(tp,D);
    zc=~(class==-999);
    if (z)
	    t1=(input(z)-Mean(i,z))';
	    Rz=Ri(z);
	    dist=t1'*(Rz.*t1);
            h(tp,i) = (2*pi)^(-lz/2) *sqrt(prod(Rz)) *exp(-dist/2);
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

% MAXIMIZATION STEP

% MEANS AND VARIANCES

cinputs=inputs(:,1:D-1);

for i=1:M
  Cov=zeros(1,D-1);
  nMean=zeros(1,D-1);
  Ri=R((i-1)*(D-1)+1:i*(D-1));
  oCov=(Ri).^(-1);
  
  % complete data for Gaussian i
  
  for tp=1:N
    input=inputs(tp,1:D-1);
    z=~(input==-999); 
    if(sum(~z))   
      input(~z)=Mean(i,~z);
    end;    
    cinputs(tp,:)=input;
  end;  

  % estimate Mean
  
  Z=sum(h(:,i));
  Mean(i,:)=(h(:,i)/Z)'*cinputs;
  
  % estimate variance
  
  for tp=1:N
    cinput=cinputs(tp,1:D-1);
    z=~(inputs(tp,1:D-1)==-999); 
    t2=cinput-Mean(i,:);
    Cov = Cov + h(tp,i)*((t2.*t2) + (~z).*oCov'); 
  end;
  Cov=Cov/Z;
  R((i-1)*(D-1)+1:i*(D-1)) = (Cov+epsi).^(-1);
end;

% CLASSES

for i=1:M
  Z=sum(h(:,i));
  for k=1:nclass
    mu(i,k)=(h(:,i)/Z)'* ...
	((inputs(:,D)==k)+(inputs(:,D)==-999)*mu(i,k));
  end;
end;

for i=1:M;
  mu(i,:)=mu(i,:)/sum(mu(i,:));
end;
