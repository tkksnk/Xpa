function [mpvec0 mpvec1] = calcaccurat(json,simT,drop,nm,nz,fitflag)

% load variables
knotsm = json.input.knotsm;
izvec = json.input.izvec;
mpmat0 = json.output.mpmat0;
Kvec = json.output.Kvec;
Kpvec = json.output.Kpvec;

for iz=1:nz
    for im=1:nm
        
        mpmat(im,iz) = mpmat0(nm*(iz-1)+im);
        
    end
end

% dynamic forecast
mpvec0 = calcDH(knotsm,mpmat,izvec,Kvec,0);
% static forecast
mpvec1 = calcDH(knotsm,mpmat,izvec,Kvec,1);

izvec = izvec(drop+1:simT+drop);
Kvec = Kvec(drop+1:simT+drop);
Kpvec = Kpvec(drop+1:simT+drop);
mpvec0 = mpvec0(drop+1:simT+drop);
mpvec1 = mpvec1(drop+1:simT+drop);

for iz = 1:nz

    DHmax(iz,1) = 100.0*max(abs(log(mpvec0(izvec==iz)./Kpvec(izvec==iz))));
    DHmean(iz,1) = 100.0*sum(abs(log(mpvec0(izvec==iz)./Kpvec(izvec==iz))))/simT;

    y = log(Kpvec(izvec==iz));
    ymean = mean(y,1);
    
    if (fitflag)
        % log-linear forecasts by fitting to minimize least sqaures 
        X = [ones(sum(izvec==iz),1) log(Kvec(izvec==iz))];
        Beta = inv(X'*X)*X'*y;
        yhat = zeros(sum(izvec==iz),2);
        yhat(:,1) = Beta(1,1) + Beta(2,1)*X(:,2);
%        yhat = X*Beta;
        e1 = y(:,1)-yhat(:,1);
    else
        % nonlinear forecasts ("static" forecast errors)
        e1 = log(Kpvec(izvec==iz)./mpvec1(izvec==iz));
    end
    
    RMSE(iz,1) = 100*sqrt(mean(e1.^2));
    Rsq(iz,1) = 1-e1'*e1/((y(:,1)-ymean(1))'*(y(:,1)-ymean(1)));
    
end

% TODO: export to csv?
disp('    DH max    DH mean   RMSE      Rsq');
%disp('    K''        K''        K''        K''');
disp([DHmax DHmean RMSE Rsq]);


function mpvec = calcDH(knotsm,mpmat0,izvec,Kvec,statflag)


tott = size(izvec,1); % error to be fixed in fortran code

mp = Kvec(1);
mpvec = zeros(tott,1);

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
    % log-linear interpolation
    im = gridlookup2(log(mnow),log(knotsm));
    wm = log(knotsm(im+1)/mnow)/log(knotsm(im+1)/knotsm(im));
    mp = exp(wm*log(mpmat0(im,iz)) + (1.0d0-wm)*log(mpmat0(im+1,iz)));

    mpvec(tt) = mp;

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
