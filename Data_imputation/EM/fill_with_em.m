function [cinputs] = fill_with_em(inputs)
% function [cinputs] = fill_with_em(inputs)
%
% Each row of the input is a data sample.
%
% Sample Script for running EM algorithm
%
% These variables must be set before running the program:
% full=1 (full Gaussians) or full=0 (diagonal Gaussians)
% reset=1 (randomize Mean and covariance) or 0 (keep old values)
% cycles (number of cycles to run)
%
% Also needs the files data.mat and tdata.mat containing the 
% training and test data.

M=50; 			% number of Gaussians
D=length(inputs(1,:)); 	% dimensionality of the input
N=length(inputs(:,1)); 	% number of training patterns
epsi=10e-8; 		% small number for preventing divide by zeros
cycles = 250;

Mean=rand(M,D);

for cycle=1:cycles
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

  if (mod(cycle, 50) == 0)  
%     fprintf('cycle %g \n',cycle);
  end;
end;
