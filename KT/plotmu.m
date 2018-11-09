clear all;

% distribution
ne = 5;
nb = 2001;
GAMY = 1.016;
DELTA = 0.069;
json = jsondecode(fileread('./results_extend_KS.json'));
mu0 = json.ss.muss;
mpmat = json.ss.mpmat;
ymat = json.ss.ymat;
nmat = json.ss.nmat;
knotsb = json.input.knotsb;
%dir = './';
% %dir = './results/traditional/SS/';
%eval(['load ' dir 'mu0.txt;']);
% eval(['load ' dir 'mpmat.txt;']);
% eval(['load ' dir 'ymat.txt;']);
% eval(['load ' dir 'nmat.txt;']);
% eval(['load ' dir 'knotsb.txt;']);

mnow = 0.0;
mp = 0.0;
ynow = 0.0;
nnow = 0.0;
for ie=1:ne
    
    for ib=1:nb
        
        mumat1(ib,ie) = mu0(nb*(ie-1)+ib);
        mpmat1(ib,ie) = mpmat(nb*(ie-1)+ib);
        ymat1(ib,ie)  = ymat(nb*(ie-1)+ib);
        nmat1(ib,ie)  = nmat(nb*(ie-1)+ib);
        ikmat1(ib,ie) = GAMY*mpmat1(ib,ie)/knotsb(ib) - (1-DELTA);
        
    end
    
    mnow = mnow + mumat1(:,ie)'*knotsb;
    mp   = mp   + mumat1(:,ie)'*mpmat1(:,ie);
    ynow = ynow + mumat1(:,ie)'*ymat1(:,ie);
    nnow = nnow + mumat1(:,ie)'*nmat1(:,ie);

end

inow = GAMY*mp - (1-DELTA)*mnow;
%ik = GAMY-1+DELTA
mp
mnow
cnow = ynow-inow;
disp('    K/Y       N         p');
disp([mnow/ynow nnow 1/cnow]);
% disp('    Y         C         I         N         K');
% disp([ynow cnow inow nnow mnow]);

ikvec = zeros(7,1);
% mean
ikvec(1) = sum(sum(mumat1.*ikmat1));

for ie = 1:ne

    for ib = 1:nb

        ikvec(2) = ikvec(2) + mumat1(ib,ie)*(ikmat1(ib,ie)-ikvec(1))^2;
        if (abs(ikmat1(ib,ie))<0.01d0); ikvec(3) = ikvec(3) + mumat1(ib,ie); end;
        if (ikmat1(ib,ie)>0.20d0);      ikvec(4) = ikvec(4) + mumat1(ib,ie); end;
        if (ikmat1(ib,ie)<-0.20d0);     ikvec(5) = ikvec(5) + mumat1(ib,ie); end;
        if (ikmat1(ib,ie)>=0.01d0);     ikvec(6) = ikvec(6) + mumat1(ib,ie); end;
        if (ikmat1(ib,ie)<=-0.01d0);    ikvec(7) = ikvec(7) + mumat1(ib,ie); end;

    end

end

disp('    Mean      Stddev    Inaction  Spike+    Spike-    Invest+   Invest-');
disp(ikvec');


% figure;
% subplot(211);
% index = 1:0.4*(nb-1)+1;
% plot(knotsb(index),mumat1(index,1),'LineWidth',2.0);
% hold on;
% plot(knotsb(index),mumat1(index,2),'LineWidth',2.0);
% xlim([knotsb(index(1)) knotsb(index(end))]);
% % xlabel('Capital');
% ylabel('Density');
% subplot(212);
% index = 1:0.8*(nb-1)+1;
% plot(knotsb(index),mumat1(index,3),'LineWidth',2.0);
% hold on;
% plot(knotsb(index),mumat1(index,4),'LineWidth',2.0);
% plot(knotsb(index),mumat1(index,5),'LineWidth',2.0);
% xlim([knotsb(index(1)) knotsb(index(end))]);
% xlabel('Capital');
% ylabel('Density');

% figure;
% index = 1:nb;
% plot(knotsb(index),sum(mumat1(index,:),2),'LineWidth',2.0);
% xlim([knotsb(index(1)) knotsb(index(end))]);
% xlabel('Capital');
% ylabel('Density');

for ie=1:ne
    if (ne>1); subplot(2,3,ie); end;
%    index = 1:(nb-1)*0.4+1;
    index = 1:nb;
    plot(knotsb(index),mumat1(index,ie),'LineWidth',2.0);
    xlim([knotsb(index(1)) knotsb(index(end))]);
    xlabel('Capital');
    ylabel('Density');
    title(sprintf('\\epsilon_%d',ie),'FontWeight','Normal');
end