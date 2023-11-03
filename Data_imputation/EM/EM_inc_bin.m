% BATCH EM ALGORITHM
% Bernoulli Mixture
% Zoubin Ghahramani
% June 1, 1993
% note inputs have to be binary vectors

h=zeros(N,M);

% EXPECTATION STEP

for i=1:M
  for tp=1:N
    input=inputs(tp,:);    
    h(tp,i) = prod(Mean(i,input==1));
  end;
end;

h=h+epsi;

for tp=1:N
  h(tp,:)=h(tp,:)/sum(h(tp,:));
end;

% MAXIMIZATION STEP

% MEANS

cinputs=inputs;

for i=1:M
  Z=sum(h(:,i));
  % complete data
  cinputs=inputs.*(inputs~=-999)+(ones(N,1)*Mean(i,:)).*(inputs==-999);
  % re-estimate Mean
  Mean(i,:)=(h(:,i)/Z)'*cinputs;
end;

