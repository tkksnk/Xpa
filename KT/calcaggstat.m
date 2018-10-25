function [Yvec Ivec Nvec Kvec Cvec Kpvec Zvec] = calcaggstat(json,simT,drop,GAMY,DELTA)

% load variables: not here?
% json = jsondecode(fileread('./output.json'));
Kvec = json.output.Kvec;
Zvec = json.output.Zvec;
Yvec = json.output.Yvec;
Ivec = json.output.Ivec;
Nvec = json.output.Nvec;
Cvec = json.output.Cvec;
Kpvec = json.output.Kpvec;
Xvec = json.output.Xvec;
% eval(['load ' dir 'Kvec.txt;']);
% eval(['load ' dir 'Zvec.txt;']);
% eval(['load ' dir 'Yvec.txt;']);
% eval(['load ' dir 'Ivec.txt;']);
% eval(['load ' dir 'Nvec.txt;']);
% eval(['load ' dir 'Cvec.txt;']);

% Xvec = Ivec./Kvec;
Kvec = Kvec(drop+1:simT+drop);
Zvec = Zvec(drop+1:simT+drop);
Yvec = Yvec(drop+1:simT+drop);
Ivec = Ivec(drop+1:simT+drop);
Nvec = Nvec(drop+1:simT+drop);
Cvec = Cvec(drop+1:simT+drop);
Kpvec = Kpvec(drop+1:simT+drop);
% Xvec = Xvec(drop+1:simT+drop);
% Kpvec = (Ivec+(1-DELTA)*Kvec)/GAMY;

so = log([Yvec Cvec Ivec Nvec Kvec Zvec]);
sd = hpfilter(so,100);
%sd = so-sf;
sd = sd(9:simT-8,:);
std0 = std(sd);

% TODO: export to csv?
disp('  Standard deviation');
disp('    Y         C         I         N         K         Z');
disp([std0(1)*100 std0(2:6)/std0(1)]);

corr0 = corr(sd);
disp('  Output correlation');
disp('    Y         C         I         N         K         Z');
disp([corr0(1,1:6)]);

disp('  Aggregate investment rate');
% aggregate investment rate nonlinearities
xmom(1) = sum(Xvec(drop+1:simT+drop),1)/simT;
xmom(2) = sum((Xvec(drop+1:simT+drop)-xmom(1)).^2,1)/simT;
xmom(3) = sum((Xvec(drop+1:simT+drop)-xmom(1)).^3,1)/simT;
xmom(4) = sum((Xvec(drop+1:simT+drop)-xmom(1)).^4,1)/simT;
% persistence
rho1 = sum((Xvec(drop+1:simT-1)-xmom(1))'*(Xvec(drop+2:simT)-xmom(1)),1)/sum((Xvec(drop+1:simT-1)-xmom(1)).^2,1);
% standard deviation
sig1 = xmom(2)^0.5d0;
% skewness
g1   = xmom(3)/(xmom(2)^1.5d0);
% excess kurtosis
g2   = xmom(4)/(xmom(2)^2) - 3.0d0;
disp('    Persist   Stddev    Skew      Exc.Kur');
disp([rho1 sig1 g1 g2]);