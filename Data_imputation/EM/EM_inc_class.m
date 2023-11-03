% BATCH EM ALGORITHM FOR CLASSIFICATION OF REAL DATA
% Zoubin Ghahramani
% December 21, 1993

h=zeros(N,M);

% EXPECTATION STEP

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

% MAXIMIZATION STEP

% MEANS AND VARIANCES

cinputs=inputs(:,1:D-1);

for i=1:M
  Cov=zeros(D-1);
  Ri=R((i-1)*(D-1)+1:i*(D-1),:);
  oCov=inv(Ri);   
  Z=sum(h(:,i));
 
  % complete data for Gaussian i
  
   for tp=1:N
    input=inputs(tp,1:D-1);
    z=~(input==-999); 
    if(sum(z))
	    Rxx=inv(oCov(z,z));  
    end;
    if(sum(~z))   
      if(sum(z))
	input(~z)=Mean(i,~z)+(oCov(~z,z)*Rxx*(input(z)-Mean(i,z))')';
      else 
	input(~z)=Mean(i,~z);
      end
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
    if(sum(z))
      Rxx=inv(oCov(z,z));  
    end;
    t2=cinput-Mean(i,:);
    Cov = Cov + h(tp,i)*(t2'*t2); 
    if(sum(~z))
      if(sum(z))
	Cov(~z,~z) = Cov(~z,~z)+h(tp,i)*(oCov(~z,~z) ...
	    - oCov(~z,z) * Rxx * oCov(z,~z));
      else
	Cov(~z,~z) = Cov(~z,~z)+h(tp,i)*(oCov(~z,~z));
      end;
    end;
  end;
  Cov=Cov/Z;
  R((i-1)*(D-1)+1:i*(D-1),:) = inv(Cov+epsi*eye(D-1));
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
	if(sum(mu(i,:)))  mu(i,:)=mu(i,:)/sum(mu(i,:));
	else mu(i,:)=ones(1,nclass)/nclass;
	end;
end;
