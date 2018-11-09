% aggresults.m: Display results about aggregate variables
clear all;

%% parameters
jsonXpa = jsondecode(fileread('./results_calibKMP_Xpa.json'));
jsonKS = jsondecode(fileread('./results_calibKMP_KS.json'));
% jsonXpa = jsondecode(fileread('./results_calibhety_Xpa.json'));
% jsonKS = jsondecode(fileread('./results_calibhety_KS.json'));
% jsonXpa = jsondecode(fileread('./results_calibKS_Xpa.json'));
% jsonKS = jsondecode(fileread('./results_calibKS_KS.json'));

drop = jsonXpa.input.drop;
simT = jsonXpa.input.simT;
irfdrop = jsonXpa.input.irfdrop;
nm = size(jsonXpa.input.knotsm,1);
nz = size(jsonXpa.input.Gz,1);

fitflag = 0; % =1 fit linear forecasting rule for static errors
priflag = 1; % =1 print eps files


%% business cycle statistics
disp(' ');
disp(' Xpa:business cycle statistics');
[YvecXpa IvecXpa NvecXpa KvecXpa CvecXpa KpvecXpa Zvec] = calcaggstat(jsonXpa,simT,drop);
disp(' ');
disp(' KS:business cycle statistics');
[YvecKS IvecKS NvecKS KvecKS CvecKS KpvecKS Zvec] = calcaggstat(jsonKS,simT,drop);

figure;
time = 1001:1200;

subplot(221);
plot(time,log(KvecXpa(time)),'b-o');
hold on;
plot(time,log(KvecKS(time)),'k-x');
title('Capital');
ylabel('Log');
xlim([time(1) time(end)]);
xticks([1 40 80 120 160 200]+1000);
xticklabels([1 40 80 120 160 200]);

subplot(222);
plot(time,log(CvecXpa(time)),'b-o');
hold on;
plot(time,log(CvecKS(time)),'k-x');
title('Consumption');
xlim([time(1) time(end)]);
xticks([1 40 80 120 160 200]+1000);
xticklabels([1 40 80 120 160 200]);

subplot(223);
plot(time,log(YvecXpa(time)),'b-o');
hold on;
plot(time,log(YvecKS(time)),'k-x');
title('Output');
xlabel('Quarter');
ylabel('Log');
legend('Xpa','KS');
xlim([time(1) time(end)]);
xticks([1 40 80 120 160 200]+1000);
xticklabels([1 40 80 120 160 200]);

subplot(224);
plot(time,log(IvecXpa(time)),'b-o');
hold on;
plot(time,log(IvecKS(time)),'k-x');
title('Investment');
xlabel('Quarter');
xlim([time(1) time(end)]);
xticks([1 40 80 120 160 200]+1000);
xticklabels([1 40 80 120 160 200]);

if (priflag); print -depsc2 simall_bench.eps; end;
%if (priflag); print -depsc2 simall_extend.eps; end;

izvec = jsonXpa.input.izvec;
izvec = izvec(drop+1:simT+drop);

for iz = 1:nz

    Diffmax(iz,1) = 100.0*max(abs(log(KpvecXpa(izvec==iz)./KpvecKS(izvec==iz))));
    Diffmean(iz,1) = 100.0*sum(abs(log(KpvecXpa(izvec==iz)./KpvecKS(izvec==iz))))/simT;
   
end

disp(' ');
disp(' Difference between Xpa and KS');
disp('    max         mean');
disp([Diffmax Diffmean]);


%% IRFs
disp(' ');
disp(' Xpa:Stochastic steady state');
[YirvecXpa IirvecXpa NirvecXpa KirvecXpa CirvecXpa Zirvec shareWirvecXpa] = calcstochss(jsonXpa,irfdrop);
disp(' ');
disp(' KS:Stochastic steady state');
[YirvecKS IirvecKS NirvecKS KirvecKS CirvecKS Zirvec shareWirvecKS] = calcstochss(jsonKS,irfdrop);

figure;
st = irfdrop+1;
ed = irfdrop+7;
time = st:ed;
subplot(221);
plot(time,100*log(Zirvec(st:ed)/Zirvec(st)),'r-');
hold on;
plot([st ed],[0 0],'k-');
title('TFP');
ylabel('Percent');
xlim([time(1) time(end)]);
xticks([0 1 2 3 4 5]+(st+1)); % st is set to -1
xticklabels([0 1 2 3 4 5]);
subplot(222);
plot(time,100*log(CirvecXpa(st:ed)/CirvecXpa(st)),'b-o');
hold on;
plot(time,100*log(CirvecKS(st:ed)/CirvecKS(st)),'k-x');
plot([st ed],[0 0],'k-');
title('Consumption');
xlim([time(1) time(end)]);
xticks([0 1 2 3 4 5]+(st+1)); % st is set to -1
xticklabels([0 1 2 3 4 5]);
subplot(223);
plot(time,100*log(YirvecXpa(st:ed)/YirvecXpa(st)),'b-o');
hold on;
plot(time,100*log(YirvecKS(st:ed)/YirvecKS(st)),'k-x');
plot([st ed],[0 0],'k-');
title('Output');
xlabel('Quarter');
ylabel('Percent');
xlim([time(1) time(end)]);
xticks([0 1 2 3 4 5]+(st+1)); % st is set to -1
xticklabels([0 1 2 3 4 5]);
subplot(224);
plot(time,100*log(IirvecXpa(st:ed)/IirvecXpa(st)),'b-o');
hold on;
plot(time,100*log(IirvecKS(st:ed)/IirvecKS(st)),'k-x');
plot([st ed],[0 0],'k-');
title('Investment');
xlabel('Quarter');
xlim([time(1) time(end)]);
xticks([0 1 2 3 4 5]+(st+1)); % st is set to -1
xticklabels([0 1 2 3 4 5]);

if (priflag); print -depsc2 irfall_bench.eps; end;


%% Accuracy Statistics
disp(' ');
disp(' Xpa:Accuracy statistics');
[mpvec0Xpa mpvec1Xpa] = calcaccurat(jsonXpa,simT,drop,nm,nz,fitflag);
disp(' ');
disp(' KS:Accuracy statistics');
[mpvec0KS mpvec1KS] = calcaccurat(jsonKS,simT,drop,nm,nz,fitflag);

figure;
time = 1001:1500;

subplot(211);
plot(time,mpvec0KS(time),'m-.');
hold on;
plot(time,mpvec1KS(time),'g--');
plot(time,KpvecKS(time),'k-');
title('KS:Capital');
ylabel('Log');
legend('Dynamic','Static','Actural','Location','NorthEast');
xlim([time(1) time(end)]);

subplot(212);
plot(time,mpvec0Xpa(time),'m-.');
hold on;
plot(time,mpvec1Xpa(time),'g--');
plot(time,KpvecXpa(time),'k-');
title('Xpa:Capital');
xlabel('Quarter');
ylabel('Log');
xlim([time(1) time(end)]);

if (priflag); print -depsc2 DHall_bench.eps; end;


disp([mean(jsonKS.output.eptimein) mean(jsonKS.output.eptimeout) size(jsonKS.output.eptimeout,1) sum(jsonKS.output.eptimein)+sum(jsonKS.output.eptimeout)]);
disp([mean(jsonXpa.output.eptimein) mean(jsonXpa.output.eptimeout) size(jsonXpa.output.eptimeout,1) sum(jsonXpa.output.eptimein)+sum(jsonXpa.output.eptimeout)]);

[11.06265669  0.94051792  0.04897656  0.88511369]
[11.06235819  0.9391175   0.04936251  0.88555691]
[11.06721553  1.03328721  0.34264785]
[ 0.11580618  0.88331998  3.63404863 12.65910948 82.70771572 19.2581517
 29.2371157  15.78997009]
[ 0.78278043 -2.31790857]
[11.06782983  1.03328721  0.34264443]
[ 0.11588526  0.88391026  3.63591532 12.6631013  82.70118786 19.25868397
 29.23370142 15.78344641]
[ 0.78272229 -2.27940205]
[[[0.34043067]
  [0.34454578]]
