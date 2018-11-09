clear all;

jsonKS = jsondecode(fileread('./results_traditional_KS.json'));
jsonXpa = jsondecode(fileread('./results_traditional_Xpa.json'));
knotsk = jsonKS.input.knotsk;
knotsb = jsonKS.input.knotsb;
knotsm = jsonKS.input.knotsm;
Gz = jsonKS.input.Gz;
nk = size(knotsk,1);
nm = size(knotsm,1);
nz = size(Gz,1);

gmat0Xpa = reshape(jsonXpa.output.gmat0,[nk,nm,nz]);
gmat0KS = reshape(jsonKS.output.gmat0,[nk,nm,nz]);
mpmat0Xpa = reshape(jsonXpa.output.mpmat0,[nm nz]);
mpmat0KS = reshape(jsonKS.output.mpmat0,[nm nz]);

% z = 1.0, k = K (linear interp)
iz = 3; % Gz(iz) = 1.0
nm1 = 51;
knotsm1 = linspace(knotsm(1),knotsm(end),nm1)';

for ik = 1:nk
    for im1 = 1:nm1
        gmat1Xpa(ik,im1) = interp2(knotsk,knotsm,gmat0Xpa(:,:,iz)',knotsk(ik),knotsm1(im1));
        gmat1KS(ik,im1) = interp2(knotsk,knotsm,gmat0KS(:,:,iz)',knotsk(ik),knotsm1(im1));
    end
end

for im1 = 1:nm1

    gvec1Xpa(im1,1) = interp1(knotsk,gmat1Xpa(:,im1),knotsm1(im1),'linear','extrap');
    gvec1KS(im1,1) = interp1(knotsk,gmat1KS(:,im1),knotsm1(im1),'linear','extrap');

end

gmatss = jsonKS.ss.gmatss;
muss = jsonKS.ss.muss;
mss = muss'*knotsb
zetak = mss - interp1(knotsk,gmatss,mss)

figure;
plot3(knotsm1,knotsm1,gvec1Xpa,'r-','LineWidth',2.0);
hold on;
mesh(knotsk,knotsm1,gmat1Xpa');
xlabel('k'); ylabel('K_t'); zlabel('k''');
title('g_k(k;z_t,K_t)','FontWeight','Normal');
%legend('g_k(k;z_t,K_t)','g_k(K_t;z_t,K_t)');
xlim([0.5,3.0]);
view(-20,35);
%axis square;
print -depsc2 pf1.eps

figure;
plot(knotsm1,gvec1Xpa,'r-','LineWidth',2.0);
hold on;
plot(knotsm1,interp1(knotsm,mpmat0Xpa(:,3),knotsm1),'b--','LineWidth',2.0)
%plot(knotsm1,gvec1KS,'r-','LineWidth',1.0);
%plot(knotsm1,interp1(knotsm,mpmat0KS(:,3),knotsm1),'b:','LineWidth',1.0)
plot([mss mss],[knotsm1(1) knotsm1(end)],'k-','LineWidth',2.0,'Color',[.5 .5 .5]);
plot(knotsm1,knotsm1,'k-');
xlabel('K_t'); ylabel('K_{t+1}');
legend('g_k(K_t;z_t,K_t)','\Gamma_K(z_t,K_t)',...
    'Location','SouthEast');
%    'g_k^{KS}(K_t;z_t,K_t)','\Gamma_K^{KS}(z_t,K_t)',...
xlim([knotsm1(1) knotsm1(end)]);
ylim([knotsm1(1) knotsm1(end)]);
axis square;

print -depsc2 pf2.eps

figure
[AX,H1,H2] = plotyy(knotsk,gmatss,knotsb,muss);
%plot(knotsk,gmatss);
hold on;
plot(mss,mss,'k*');
plot(mss,interp1(knotsk,gmatss,mss),'bo');
%xlim([0.5 3.0]);
xlim(AX(1),[0.5 2.5]);
xlim(AX(2),[0.5 2.5]);
legend('g(k)','K_{ss}','g(K_{ss})','Location','SouthEast');

print -depsc2 pf3.eps
