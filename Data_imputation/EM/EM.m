% BATCH EM ALGORITHM
% Zoubin Ghahramani
% April 8, 1992

h=zeros(N,M);

% EXPECTATION STEP

for i=1:M
  Ri = R((i-1)*D+1:i*D,:);
  for tp=1:N
    input=inputs(tp,:);
    t1=(input-Mean(i,:))';
    dist=t1'*Ri*t1;
    h(tp,i) = (2*pi)^(-D/2)*sqrt(det(Ri))*exp(-dist/2);
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
  Mean(i,:)=(h(:,i)/Z)'*inputs;
end;

% VARIANCES

for i=1:M
  Z=sum(h(:,i));
  Cov=zeros(D);
  for tp=1:N
    input=inputs(tp,:);
    Cov = Cov + (h(tp,i)/Z)* ...
	(input-Mean(i,:))'*(input-Mean(i,:));
  end;
  R((i-1)*D+1:i*D,:) = inv(Cov+epsi*eye(D));
end;


