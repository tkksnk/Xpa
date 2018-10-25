function [Yirvec Iirvec Nirvec Kirvec Cirvec Zirvec shareWirvec] = calcstochss(json,irfdrop)

Kirvec = json.irf.Kvec;
Zirvec = json.irf.Zvec;
Yirvec = json.irf.Yvec;
Iirvec = json.irf.Ivec;
Nirvec = json.irf.Nvec;
Cirvec = json.irf.Cvec;
shareWirvec(:,1) = json.irf.shareW1;
shareWirvec(:,2) = json.irf.shareW2;
shareWirvec(:,3) = json.irf.shareW3;
shareWirvec(:,4) = json.irf.shareW4;
shareWirvec(:,5) = json.irf.shareW5;
shareWirvec(:,6) = json.irf.shareW9095;
shareWirvec(:,7) = json.irf.shareW9599;
shareWirvec(:,8) = json.irf.shareWT1;
gini = json.irf.gini;
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
disp('    K/Y');
disp([mnow/ynow]);
disp('    Q1        Q2        Q3        Q4        Q5        9095      9599      T1%');
disp(shareWirvec(irfdrop,:));
disp('    Gini      DeltaC');
disp([gini(irfdrop) 100*log(Cirvec(irfdrop+2)/Cirvec(irfdrop+1))]);

% disp('    Y         C         I         N         K');
% disp([Yirvec(irfdrop) Cirvec(irfdrop) Iirvec(irfdrop) Nirvec(irfdrop) Kirvec(irfdrop)]);
