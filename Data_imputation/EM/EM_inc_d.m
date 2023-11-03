% BATCH EM ALGORITHM
% Zoubin Ghahramani
% DIAGONAL VERSION
% May 15, 1992

% EXPECTATION STEP -- COMPUTE E(Zij)
h=zeros(N,M);

for i=1:M
  Ri = R((i-1)*D+1:i*D);
  for tp=1:N
    input=inputs(tp,:);
    z=~(input==-999); lz=sum(z);
    t1=(input(z)-Mean(i,z))';
    Rz=Ri(z);
    dist=t1'*(Rz.*t1);
    h(tp,i) = (2*pi)^(-lz/2)*sqrt(prod(Rz))*exp(-dist/2);
  end;
end;

h=h+epsi;

for tp=1:N
    h(tp,:)=h(tp,:)/sum(h(tp,:));
end;

% MAXIMIZATION STEP

% MEANS

for i=1:M
  nMean=zeros(input);
  for tp=1:N;
    input=inputs(tp,:);
    z=~(input==-999);
    nMean=nMean+h(tp,i)*((input.*z)+Mean(i,:).*(~z));
  end;
  Z=sum(h(:,i));
  Mean(i,:)=nMean/Z;
end;

% VARIANCES

for i=1:M
  Cov=zeros(1,D);
  oldCov=(R((i-1)*D+1:i*D)).^(-1);
  for tp=1:N
    input=inputs(tp,:);
    z=~(input==-999); 
    t2=(input.*z-Mean(i,:).*z);
    Cov = Cov + h(tp,i)*(t2.*t2 + (~z).*oldCov');
  end;
  Z=sum(h(:,i));
  Cov=Cov/Z;
  R((i-1)*D+1:i*D) = (Cov+epsi).^(-1);
end;




