% Sample Script for running EM algorithm
%
% These variables must be set before running the program:
% full=1 (full Gaussians) or full=0 (diagonal Gaussians)
% reset=1 (randomize Mean and covariance) or 0 (keep old values)
% cycles (number of cycles to run)
%
% Also needs the files data.mat and tdata.mat containing the 
% training and test data.

disp('loading training data');
load data;
inputs=data;

disp('loading test data');
load tdata;
testdata=tdata;

M=3; 			% number of Gaussians
D=length(inputs(1,:)); 	% dimensionality of the input
N=length(inputs(:,1)); 	% number of training patterns
epsi=10e-8; 		% small number for preventing divide by zeros
scale=0.1; 		% initial scale of the Gaussian variances

if(reset)		% randomize Means and variances
  Mean=rand(M,D);
  if (full) 		% full Gaussians
    R=zeros(M*D,D);
    for i=1:M
      R((i-1)*D+1:i*D,:)=eye(D)*scale;
    end
  else 			% diagonal Gaussians
    R=ones(M*D,1)*scale;                
  end;
end;

for cycle=1:cycles;
  if (full)		% full Gaussians
    EM_inc;    
  else			% diagonal Gaussians
    EM_inc_d;
  end;

  fprintf('cycle %g \n',cycle);
  Mean
  R
end;

