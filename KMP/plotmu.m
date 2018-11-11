clear all;

% distribution
drop = 500;
simT = 2000;
ne = 2;
ny = 1;
nd = 1;
nx = ne*ny*nd;
nb = 2001;
nm = 5;
nz = 2;
%json = jsondecode(fileread('./results_calibKMPbenchnk2500_KS.json'));
json = jsondecode(fileread('./results_calibKS_KS.json'));
%json = jsondecode(fileread('./results_calibKMPhety_KS.json'));
%mu0 = json.ss.muss;
%mux = json.input.mux;
Ge = json.input.Ge;
mu0 = json.irf.mu0;
mpmat = json.ss.mpmat;
knotsb = json.input.knotsb;
shareWvec = json.ss.shareWvec;
gini = json.ss.gini;
%THETA = json.input.param(1);
% load mu0.txt;
% load knotsb.txt;

mnow = 0.0;
mp = 0.0;
for ix=1:nx
    
    for ib=1:nb
        
        mumat1(ib,ix) = mu0(nb*(ix-1)+ib);
        mpmat1(ib,ix) = mpmat(nb*(ix-1)+ib);
        
    end
    
    mnow = mnow + mumat1(:,ix)'*knotsb;
    mp   = mp   + mumat1(:,ix)'*mpmat1(:,ix); %*THETA;

end

if (nx>2)
    muu = sum(mumat1(:,[1:2:nx-1]),2);
    mue = sum(mumat1(:,[1:2:nx-1]+1),2);
else
    muu = mumat1(:,1);
    mue = mumat1(:,2);
end
lnow = sum(muu)*Ge(1) + sum(mue)*Ge(2);

disp([shareWvec']);
disp([gini (mnow/lnow)^(1-0.36)]);
figure;
plot(knotsb,muu,'LineWidth',2.0);
hold on;
plot(knotsb,mue,'LineWidth',2.0);
xlim([knotsb(1) knotsb(end)]);


% %for id=1:nd
%     figure;
%     for iy=1:ny
%          if (ny>1); subplot(3,3,iy); end;
% %         if (iy==1)
% %             index = 1:51;
% %         elseif (iy==2)
% %             index = 1:101;
% %         else
%             index = 1:nb;
% %         end
%         plot(knotsb(index),sum(mumat1(index,[1 15 29]+iy-1),2),'LineWidth',2.0);
%         xlim([knotsb(index(1)) knotsb(index(end))]);
% %        plot(knotsb(index),log10(mumat(index,ne*ny*(id-1)+iy)));
%         hold on;
% % figure;
%         plot(knotsb(index),sum(mumat1(index,[1 15 29]+iy),2),'LineWidth',2.0);
% %        plot(knotsb(index),log10(mumat(index,ne*ny*(id-1)+iy+1)));
%         xlim([knotsb(index(1)) knotsb(index(end))]);
%     end
% % end