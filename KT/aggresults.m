% Aggresults.m: Display results about aggregate variables from KT model

clear all;

% parameters
jsonKS = jsondecode(fileread('./results_extend_KS.json'));
jsonXpa = jsondecode(fileread('./results_extend_Xpa.json'));
% jsonKS = jsondecode(fileread('./results_traditional_KS.json'));
% jsonXpa = jsondecode(fileread('./results_traditional_Xpa.json'));
drop = 500;
simT = 2000;
simTT = drop+simT;
irfdrop = 200;
irfT = 50;
irfTT = irfdrop+irfT;
nm = size(jsonXpa.input.knotsm,1);
nz = size(jsonXpa.input.Gz,1);
GAMY = jsonXpa.input.param(1);
DELTA = jsonXpa.input.param(3);
fitflag = 0;
priflag = 1;

% business cycle statistics
disp(' ');
disp(' Xpa:business cycle statistics');
[YvecXpa IvecXpa NvecXpa KvecXpa CvecXpa KpvecXpa Zvec] = calcaggstat(jsonXpa,simT,drop,GAMY,DELTA);
disp(' ');
disp(' KS:business cycle statistics');
[YvecKS IvecKS NvecKS KvecKS CvecKS KpvecKS Zvec] = calcaggstat(jsonKS,simT,drop,GAMY,DELTA);

% simulated sequence of Y, I, N, and C for 50 years
figure;
time = 1001:1050;

subplot(231);
plot(time,log(YvecXpa(time)),'b-o');
hold on;
plot(time,log(YvecKS(time)),'k-x');
title('Output');
ylabel('Log');
legend('Xpa','KS');
xlim([time(1) time(end)]);
ylim([-0.8 -0.1]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

subplot(232);
plot(time,log(IvecXpa(time)),'b-o');
hold on;
plot(time,log(IvecKS(time)),'k-x');
title('Investment');
xlim([time(1) time(end)]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

subplot(233);
plot(time,log(NvecXpa(time)),'b-o');
hold on;
plot(time,log(NvecKS(time)),'k-x');
title('Labor');
xlim([time(1) time(end)]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

subplot(234);
plot(time,log(KvecXpa(time)),'b-o');
hold on;
plot(time,log(KvecKS(time)),'k-x');
title('Capital');
xlabel('Year');
ylabel('Log');
xlim([time(1) time(end)]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

subplot(235);
plot(time,log(CvecXpa(time)),'b-o');
hold on;
plot(time,log(CvecKS(time)),'k-x');
title('Consumption');
xlabel('Year');
xlim([time(1) time(end)]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

subplot(236);
plot(time,log(Zvec(time)),'r-');
title('TFP');
xlabel('Year');
xlim([time(1) time(end)]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

%if (priflag); print -depsc2 simall_extend.eps; end;
% if (priflag); print -depsc2 simall_trad.eps; end;

izvec = jsonXpa.input.izvec;
izvec = izvec(drop+1:simT+drop);

for iz = 1:nz

    Diffmax(iz,1) = 100.0*max(abs(log(KpvecXpa(izvec==iz)./KpvecKS(izvec==iz))));
    Diffmean(iz,1) = 100.0*sum(abs(log(KpvecXpa(izvec==iz)./KpvecKS(izvec==iz))))/simT;
    Diffmax(iz,2) = 100.0*max(abs(log(CvecXpa(izvec==iz)./CvecKS(izvec==iz))));
    Diffmean(iz,2) = 100.0*sum(abs(log(CvecXpa(izvec==iz)./CvecKS(izvec==iz))))/simT;
   
end

disp(' ');
disp(' Difference between Xpa and KS');
disp('    max                 mean');
disp('    K''        p         K''        p');
disp([Diffmax Diffmean]);


% IRFs
disp(' ');
disp(' Xpa:Stochastic steady state');
[YirvecXpa IirvecXpa NirvecXpa KirvecXpa CirvecXpa Zirvec ikirvecXpa] = calcstochss(jsonXpa,irfdrop);

disp(' ');
disp(' KS:Stochastic steady state');
[YirvecKS IirvecKS NirvecKS KirvecKS CirvecKS Zirvec ikirvecKS] = calcstochss(jsonKS,irfdrop);

figure;
st = irfdrop+1;
ed = irfdrop+22;
time = st:ed;
subplot(231);
% NOTE: denominator should be stoch SS???
plot(time,100*log(YirvecXpa(time)/YirvecXpa(end)),'b-o');
hold on;
plot(time,100*log(YirvecKS(time)/YirvecKS(end)),'k-x');
plot([st ed],[0 0],'k-');
title('Output');
ylabel('Percent');
legend('Xpa','KS');
xlim([time(1) time(end)]);
%ylim([-0.8 -0.1]);
xticks([0 5 10 15 20]+(st+1)); % st is set to -1
xticklabels([0 5 10 15 20]);

subplot(232);
plot(time,100*log(IirvecXpa(time)/IirvecXpa(end)),'b-o');
hold on;
plot(time,100*log(IirvecKS(time)/IirvecKS(end)),'k-x');
plot([st ed],[0 0],'k-');
title('Investment');
xlim([time(1) time(end)]);
xticks([0 5 10 15 20]+(st+1)); % st is set to -1
xticklabels([0 5 10 15 20]);

subplot(233);
plot(time,100*log(NirvecXpa(time)/NirvecXpa(end)),'b-o');
hold on;
plot(time,100*log(NirvecKS(time)/NirvecKS(end)),'k-x');
plot([st ed],[0 0],'k-');
title('Labor');
xlim([time(1) time(end)]);
xticks([0 5 10 15 20]+(st+1)); % st is set to -1
xticklabels([0 5 10 15 20]);

subplot(234);
plot(time,100*log(KirvecXpa(time)/KirvecXpa(end)),'b-o');
hold on;
plot(time,100*log(KirvecKS(time)/KirvecKS(end)),'k-x');
plot([st ed],[0 0],'k-');
title('Capital');
xlabel('Year');
ylabel('Percent');
xlim([time(1) time(end)]);
xticks([0 5 10 15 20]+(st+1)); % st is set to -1
xticklabels([0 5 10 15 20]);

subplot(235);
plot(time,100*log(CirvecXpa(time)/CirvecXpa(end)),'b-o');
hold on;
plot(time,100*log(CirvecKS(time)/CirvecKS(end)),'k-x');
plot([st ed],[0 0],'k-');
title('Consumption');
xlabel('Year');
xlim([time(1) time(end)]);
xticks([0 5 10 15 20]+(st+1)); % st is set to -1
xticklabels([0 5 10 15 20]);

subplot(236);
plot(time,100*log(Zirvec(time)/Zirvec(end)),'r-');
hold on;
plot([st ed],[0 0],'k-');
title('TFP');
xlabel('Year');
xlim([time(1) time(end)]);
xticks([0 5 10 15 20]+(st+1)); % st is set to -1
xticklabels([0 5 10 15 20]);

%if (priflag); print -depsc2 irfall_extend.eps; end;
% if (priflag); print -depsc2 irfall_trad.eps; end;


% Accuracy Statistics
disp(' ');
disp(' Xpa:Accuracy statistics');
[mpvec0Xpa pvec0Xpa mpvec1Xpa pvec1Xpa] = calcaccurat(jsonXpa,simT,drop,nm,nz,GAMY,DELTA,fitflag);
disp(' ');
disp(' KS:Accuracy statistics');
[mpvec0KS pvec0KS mpvec1KS pvec1KS] = calcaccurat(jsonKS,simT,drop,nm,nz,GAMY,DELTA,fitflag);

figure;
time = 1001:1050;

subplot(221);
plot(time,log(mpvec0KS(time)),'m-^');
hold on;
plot(time,log(mpvec1KS(time)),'c-o');
plot(time,log(KpvecKS(time)),'k-x');
title('KS:Capital');
ylabel('Log');
legend('Dynamic','Static','Actural','Location','SouthEast');
xlim([time(1) time(end)]);
ylim([-0.2 0.4]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

subplot(222);
plot(time,-log(pvec0KS(time)),'m-^');
hold on;
plot(time,-log(pvec1KS(time)),'c-o');
plot(time,log(CvecKS(time)),'k-x');
title('KS:Consumption');
xlim([time(1) time(end)]);
ylim([-0.95 -0.75]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

subplot(223);
plot(time,log(mpvec0Xpa(time)),'m-^');
hold on;
plot(time,log(mpvec1Xpa(time)),'c-o');
plot(time,log(KpvecXpa(time)),'k-x');
title('Xpa:Capital');
xlabel('Year');
ylabel('Log');
xlim([time(1) time(end)]);
ylim([-0.2 0.4]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

subplot(224);
plot(time,-log(pvec0Xpa(time)),'m-^');
hold on;
plot(time,-log(pvec1Xpa(time)),'c-o');
plot(time,log(CvecXpa(time)),'k-x');
title('Xpa:Consumption');
xlabel('Year');
xlim([time(1) time(end)]);
ylim([-0.95 -0.75]);
xticks([1 10 20 30 40 50]+1000);
xticklabels([1 10 20 30 40 50]);

%if (priflag); print -depsc2 DHall_extend.eps; end;
% if (priflag); print -depsc2 DHall_trad.eps; end;


disp([mean(jsonKS.output.eptimein) mean(jsonKS.output.eptimeout) size(jsonKS.output.eptimeout,1) sum(jsonKS.output.eptimein)+sum(jsonKS.output.eptimeout)]);
disp([mean(jsonXpa.output.eptimein) mean(jsonXpa.output.eptimeout) size(jsonXpa.output.eptimeout,1) sum(jsonXpa.output.eptimein)+sum(jsonXpa.output.eptimeout)]);
