function [mpvec0 pvec0 mpvec1 pvec1] = calcaccurat(json,simT,drop,nm,nz,GAMY,DELTA,fitflag)

% load variables
% json = jsondecode(fileread('./output.json'));
knotsm = json.input.knotsm;
izvec = json.input.izvec;
mpmat0 = json.output.mpmat0;
pmat0 = json.output.pmat0;
Kvec = json.output.Kvec;
Kpvec = json.output.Kpvec;
Cvec = json.output.Cvec;
% dir = './';
% eval(['load ' dir 'knotsm.txt;']);
% eval(['load ' dir 'mpmat0.txt;']);
% eval(['load ' dir 'pmat0.txt;']);
% eval(['load ' dir 'izvec.txt;']);
% eval(['load ' dir 'Kvec.txt;']);
% eval(['load ' dir 'Ivec.txt;']);
% eval(['load ' dir 'Cvec.txt;']);

for iz=1:nz
    for im=1:nm
        
        mpmat(im,iz) = mpmat0(nm*(iz-1)+im);
        pmat(im,iz) = pmat0(nm*(iz-1)+im);
        
    end
end

[mpvec0 pvec0] = calcDH(knotsm,mpmat,pmat,izvec,Kvec,0);
[mpvec1 pvec1] = calcDH(knotsm,mpmat,pmat,izvec,Kvec,1);

izvec = izvec(drop+1:simT+drop);
%Kpvec = (Ivec(drop+1:simT+drop) + (1-DELTA)*Kvec(drop+1:simT+drop))/GAMY;
Kvec = Kvec(drop+1:simT+drop);
Kpvec = Kpvec(drop+1:simT+drop);
Cvec = Cvec(drop+1:simT+drop);
mpvec0 = mpvec0(drop+1:simT+drop);
pvec0 = pvec0(drop+1:simT+drop);
mpvec1 = mpvec1(drop+1:simT+drop);
pvec1 = pvec1(drop+1:simT+drop);

for iz = 1:nz

    DHmax(iz,1) = 100.0*max(abs(log(mpvec0(izvec==iz)./Kpvec(izvec==iz))));
    DHmean(iz,1) = 100.0*sum(abs(log(mpvec0(izvec==iz)./Kpvec(izvec==iz))))/simT;
    DHmax(iz,2) = 100.0*max(abs(log(pvec0(izvec==iz).*Cvec(izvec==iz))));
    DHmean(iz,2) = 100.0*sum(abs(log(pvec0(izvec==iz).*Cvec(izvec==iz))))/simT;

    y = [log(Kpvec(izvec==iz)) log(1./Cvec(izvec==iz))];
    ymean = mean(y,1);
    
    if (fitflag)
        % log-linear forecasts by fitting to minimize least sqaures 
        X = [ones(sum(izvec==iz),1) log(Kvec(izvec==iz))];
        Beta = inv(X'*X)*X'*y;
        yhat = zeros(sum(izvec==iz),2);
        yhat(:,1) = Beta(1,1) + Beta(2,1)*X(:,2);
        yhat(:,2) = Beta(1,2) + Beta(2,2)*X(:,2);
%        yhat = X*Beta;
        e1 = y(:,1)-yhat(:,1);
        e2 = y(:,2)-yhat(:,2);
    else
        % nonlinear forecasts ("static" forecast errors)
        e1 = log(Kpvec(izvec==iz)./mpvec1(izvec==iz));
        e2 = log((1./Cvec(izvec==iz))./pvec1(izvec==iz));
    end
    
    RMSE(iz,1) = 100*sqrt(mean(e1.^2));
    RMSE(iz,2) = 100*sqrt(mean(e2.^2));
    Rsq(iz,1) = 1-e1'*e1/((y(:,1)-ymean(1))'*(y(:,1)-ymean(1)));
    Rsq(iz,2) = 1-e2'*e2/((y(:,2)-ymean(2))'*(y(:,2)-ymean(2)));
    
end

% TODO: export to csv?
disp('    DH max              DH mean             RMSE                Rsq');
disp('    K''        p         K''        p         K''        p         K''       p');
disp([DHmax DHmean RMSE Rsq]);


function [mpvec pvec] = calcDH(knotsm,mpmat0,pmat0,izvec,Kvec,statflag)


tott = size(izvec,1);

mp = Kvec(1);
mpvec = zeros(tott,1);
pvec = zeros(tott,1);

for tt = 1:tott

    if (statflag==1)
        mnow = Kvec(tt); 
    else
        mnow = mp;
    end
    iz = izvec(tt);
%     % linear interpolation
%     im = gridlookup2(mnow,knotsm);
%     wm = (knotsm(im+1)-mnow)/(knotsm(im+1)-knotsm(im));
%     mp = wm*mpmat0(im,iz) + (1.0d0-wm)*mpmat0(im+1,iz);
%     cnow = 1.0/(wm*pmat0(im,iz) + (1.0-wm)*pmat0(im+1,iz));
    % log-linear interpolation
    im = gridlookup2(log(mnow),log(knotsm));
    wm = log(knotsm(im+1)/mnow)/log(knotsm(im+1)/knotsm(im));
    mp = exp(wm*log(mpmat0(im,iz)) + (1.0d0-wm)*log(mpmat0(im+1,iz)));
    cnow = 1.0/exp(wm*log(pmat0(im,iz)) + (1.0-wm)*log(pmat0(im+1,iz)));

    mpvec(tt) = mp;
    pvec(tt) = 1.0/cnow;

end


%! NOTE: 070718 from mod_utils.f90 in the JSY code
function ix = gridlookup2(x0,xgrid)

nx = size(xgrid,1);

ix = 0;
for jx = 1:nx
    if (x0<=xgrid(jx)); break; end;
    ix = ix+1;
end

ix = min(max(1,ix),nx-1);
