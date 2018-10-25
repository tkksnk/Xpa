function [Yirvec Iirvec Nirvec Kirvec Cirvec Zirvec ikirvec] = calcstochss(json,irfdrop)

Kirvec = json.irf.Kvec;
Zirvec = json.irf.Zvec;
Yirvec = json.irf.Yvec;
Iirvec = json.irf.Ivec;
Nirvec = json.irf.Nvec;
Cirvec = json.irf.Cvec;
ikirvec(:,1) = json.irf.ikmean;
ikirvec(:,2) = json.irf.ikstddev;
ikirvec(:,3) = json.irf.ikinaction;
ikirvec(:,4) = json.irf.ikspikepos;
ikirvec(:,5) = json.irf.ikspikeneg;
ikirvec(:,6) = json.irf.ikpos;
ikirvec(:,7) = json.irf.ikneg;
% eval(['load ' dir 'Kirvec.txt;']);
% eval(['load ' dir 'Zirvec.txt;']);
% eval(['load ' dir 'Yirvec.txt;']);
% eval(['load ' dir 'Iirvec.txt;']);
% eval(['load ' dir 'Nirvec.txt;']);
% eval(['load ' dir 'Cirvec.txt;']);
mnow = Kirvec(irfdrop);
ynow = Yirvec(irfdrop);
nnow = Nirvec(irfdrop);
cnow = Cirvec(irfdrop);
disp('    K/Y       N         p');
disp([mnow/ynow nnow 1/cnow]);
disp('    Mean      Stddev    Inaction  Spike+    Spike-    Invest+   Invest-');
disp(ikirvec(irfdrop,:));
% disp('    Y         C         I         N         K');
% disp([Yirvec(irfdrop) Cirvec(irfdrop) Iirvec(irfdrop) Nirvec(irfdrop) Kirvec(irfdrop)]);