% BATCH EM ALGORITHM FOR CLASSIFICATION OF REAL DATA
% Zoubin Ghahramani
% December 21, 1993

h=zeros(N,M);

% EXPECTATION STEP

for i=1:M
  Ri = R((i-1)*(D-1)+1:i*(D-1),:);
  for tp=1:N
    input=inputs(tp,1:D-1);
    class=inputs(tp,D);
    t1=(input-Mean(i,:))';
    dist=t1'*Ri*t1;
    h(tp,i) = mu(i,class)*(2*pi)^(-(D-1)/2)*sqrt(det(Ri))*exp(-dist/2);
  end;
end;

h=h+epsi;

for tp=1:N
  h(tp,:)=h(tp,:)/sum(h(tp,:));
end;

% MAXIMIZATION STEP

% MEANS

for i=1:M
  Z=sum(h(:,i));
  Mean(i,:)=(h(:,i)/Z)'*inputs(:,1:D-1);
end;

% VARIANCES

for i=1:M
  Z=sum(h(:,i));
  Cov=zeros(D-1);
  for tp=1:N
    input=inputs(tp,1:D-1);
    Cov = Cov + (h(tp,i)/Z)* ...
	(input-Mean(i,:))'*(input-Mean(i,:));
  end;
  R((i-1)*(D-1)+1:i*(D-1),:) = inv(Cov+epsi*eye(D-1));
end;

% CLASSES

for i=1:M
  Z=sum(h(:,i));
  for k=1:nclass
    mu(i,k)=(h(:,i)/Z)'*(inputs(:,D)==k);
  end;
end;

for i=1:M;
   mu(i,:)=mu(i,:)/sum(mu(i,:));
end;
