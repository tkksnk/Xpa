function [Yvec Ivec Nvec Kvec Cvec Kpvec Zvec] = calcaggstat(json,simT,drop)

% load variables
Kvec = json.output.Kvec;
Zvec = json.output.Zvec;
Yvec = json.output.Yvec;
Ivec = json.output.Ivec;
Nvec = json.output.Nvec;
Cvec = json.output.Cvec;
Kpvec = json.output.Kpvec;
Xvec = json.output.Xvec;

Kvec = Kvec(drop+1:simT+drop);
Zvec = Zvec(drop+1:simT+drop);
Yvec = Yvec(drop+1:simT+drop);
Ivec = Ivec(drop+1:simT+drop);
Nvec = Nvec(drop+1:simT+drop);
Cvec = Cvec(drop+1:simT+drop);
Kpvec = Kpvec(drop+1:simT+drop);

KYmean = mean(Kvec)/mean(Yvec);
corr1 = corr([Yvec Cvec]);
std1 = std(Ivec);
ac1 = corr([Yvec(1:end-3) Yvec(4:end)]);

% so = log([Yvec Cvec Ivec Nvec Kvec Zvec]);
% sd = hpfilter(so,100);
% %sd = so-sf;
% sd = sd(9:simT-8,:);
% std0 = std(sd);
% corr0 = corr(sd);

disp('    K/Y       Corr(Y,C) Std(I)    AC(Y)');
disp([KYmean corr1(1,2) std1 ac1(1,2)]);
% % TODO: export to csv?
% disp('  Standard deviation');
% disp('    Y         C         I         N         K         Z');
% disp([std0(1)*100 std0(2:6)/std0(1)]);
% 
% disp('  Output correlation');
% disp('    Y         C         I         N         K         Z');
% disp([corr0(1,1:6)]);