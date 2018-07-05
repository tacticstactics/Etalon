% etalon_func.m
function [wlcol, PTetacol, PRetacol] = etalon_func(mm, wl0, stepwl);

if nargin ==0
 mm = 2^7; % points
 stepwl = 0.0001;	%um
 wl0 = 0.6500;	% um

end

global c
 c = 2.99792458 * 10^14; % um / s


%______________________________________________________________

% etalon
 Re1 = 0.8;  Re2 = 0.8;
 re1 = sqrt(Re1);  re2 = sqrt(Re2); 
 te1 = sqrt(1-Re1); te2 = sqrt(1-Re2);

% etalen = 1501.590; % um -> odd(red)
%etalen = 1501.9775; % um + 0.25 * 1.55 um -> even(green)
%etalen = 1502.340; % um + 0.5*1.55um -> odd(red)
%etalen = 1502.7525; % um + 0.75 * 1.55 um -> even(green)
%etalen = 1503.090; % um + 1* 1.55um -> odd(red)

 %etalen = 3003.18; % um ->
 %etalen = 6006.36;% + 0.25*wavel0; % um ->
%etalen = 9009.54 ;%+ 0.25*wavel0; % um ->
%etalen = 12012.66 + 0.25*1.545; % um ->


%etalen = 0.5;% + 0.25*wavel0; % um ->
etalen = 10.0;% + 0.25*wavel0; % um ->

%______________________________________________________________

 endwl = wl0 + stepwl*mm;
 wlcol = zeros(mm,1); omega1col = zeros(mm,1);


 ETeta_col=zeros(mm,1);EReta_col=zeros(mm,1);
 angleE1_col=zeros(mm,1);angleE2_col=zeros(mm,1);
 PTetacol=zeros(mm,1);angleTetacol=zeros(mm,1);
 PRetacol=zeros(mm,1);angleRetacol=[];

%__________________________________________________________________

for ii = 1:mm;
wavel1 = wl0 + ii .* stepwl;
 wlcol(ii,1) = wavel1;
omega1 = 2*pi*c/wavel1; % rad/ps
omega1col(ii,1) = omega1;

 % Trans!!!
 ETeta =(te1*te2*exp(-1j*4*pi*10^0*(1.00*etalen)./ wavel1)) ./ (1+re1*re2*exp(-1j*4*pi*10^0 *(1.00*etalen)./ wavel1));
 ETeta_col(ii,1) = ETeta;

 PTeta = (ETeta * conj(ETeta));
 PTetacol(ii,1) = PTeta; % Trans !!!
 
 angleTeta = (angle(ETeta));
 angleTetacol(ii,1) = angleTeta;
 
%__________________________________________________________________
 % Reflect!!!
 EReta =(re1+re2*exp(-1j*4*pi*10^0*(1.00*etalen)./ wavel1)) ./ (1+re1*re2*exp(-1j*4*pi*10^0 *(1.00*etalen)./ wavel1));
 EReta_col(ii,1) = EReta;

 PReta = (EReta * conj(EReta));
 PRetacol(ii,1) = PReta; % reflect !!!
 
angleReta = unwrap(angle(EReta));
angleRetacol = [angleRetacol;angleReta];

end 

 %Tau, delay _________________________________________________________________________  

%tau1_col=[]; tau2_col=[];tauTeta_col=[];tauReta_col=[];
%for qq = 1:mm-1
   
   % Reflect
 %  tauReta = 10^12 .*((angleRetacol(qq+1)) - (angleRetacol(qq))) ./...
 %     (wavel1col(qq+1) - wavel1col(qq)) .* (wavel1col(qq)).^2 ./ (2*pi*c); 

%   tauReta_col=[tauReta_col;tauReta];
   
   %_________________________________________________________________________  
   
   % Trans
 %  tauTeta = 10^12 .*((angleTetacol(qq+1)) - (angleTetacol(qq))) ./...
 %     (wavel1col(qq+1) - wavel1col(qq)) .* (wavel1col(qq)).^2 ./ (2*pi*c); 

%   tauTeta_col=[tauTeta_col;tauTeta];

% end
 
 % Chromatic Dispersion
% cd1_col = [];cd2_col=[];cdReta_col=[];cdTeta_col=[];

%for ww = 1:mm-2
% cdReta = (tauReta_col(ww+1) - tauReta_col(ww)) ./ (1000 .* (wavel1col(ww+1) - wavel1col(ww)));   
% cdReta_col = [cdReta_col;cdReta];
 
% cdTeta = (tauTeta_col(ww+1) - tauTeta_col(ww)) ./ (1000 .* (wavel1col(ww+1) - wavel1col(ww)));   
% cdTeta_col = [cdTeta_col;cdTeta];
 
% end
%_________________________________________________________________________________________________

 
%fig2=figure(2);
%plot(wavel1col,PRetacol,'r.-',wavel1col,PTetacol,'g.-');
% grid on;
% xlabel('wavelength'); ylabel('Power'); ylim([-0.1 1.1]);
 
%fig3=figure(3);
%plot(wavel1col,angleRetacol,'r.-',wavel1col,angleTetacol,'g.-');
%grid on;
%xlabel('wavelength'); ylabel('angle');
%  ylim([-5 0.5]);
  
% fig4=figure(4);
%plot(wavel1col([1:mm-1]), tauReta_col, 'r.', wavel1col([1:mm-1]), tauTeta_col, 'g.');
%ylim([-2 2]);
%grid on;
%xlabel('wavelength'); ylabel('delay');


% fig5=figure(5);
% plot(wavel1col([1:mm-2]), cdReta_col, 'r.-',    wavel1col([1:mm-2]), 1.05 .* cdTeta_col, 'g.-');
% ylim([-.5 .5]);
% grid on;
% xlabel('wavelength');
% ylabel('Chromatic Dispersion')
 
%_________________________________________________________________________________________________

%wlcol = wavel1col;
%[b,a] = itlvr_getba(wlcol,E1_col);

%tfba = tf(b,a)
%rootsb = roots(b)
%rootsa = roots(a)

% fig6=figure(6);
% subplot(211);stem(b);
%  subplot(212);stem(a);
  
%minrealtfba = minreal(tfba)

%[H,T] = impz(b,a);
%fig7=figure(7);
%impulse(minrealtfba,10)
%stem(T, H .* conj(H));
%xlim([0 30])

%fig8 = figure(8);
%set(fig8,'Position',[10 200 500 500])
%zplane(b,a);
%grid on

%fig9 = figure(9);
%set(fig9,'Position',[10 200 300 300])
%freqz(b,a);

%
%____________________________________________
%for ii = 1:mm;
%%   t = 1e12.* (1/Fs) *ii;
 %  tcol(ii,1) = t;
   
   % 1 or 244 is good to see what's goin'on
 %if  ii >= 1 & ii <= 1
   %x1 = sin(t*1*2*pi); %red
   % x2 = sin(t*0.1 *2*pi);	%green
  % x1 = 0.999;
  % x2 = 0.66;
   
   %Pass frequency : 0.35;	x3 = sin(t* 0.35 *2*pi);	%blue
   %x3 = 0.33;
   %x3 = sin(t * (1/mm).* 100 *2*pi);	%green

% else 
%    x1 = 0;
%    x2 = 0;
%    x3 = 0;
%    x4 = 0;    
    
 %   end
 %   x1col(ii,1) = x1;   
 %  x2col(ii,1) = x2;   
 %  x3col(ii,1) = x3;   

%end

%fig5 = figure(5);
%set(fig5,'Position',[50 200 900 700])

%subplot(211);
%stem(tcol,x1col,'r.-');hold on;
%stem(tcol,x2col,'g.-');
%stem(tcol,x3col,'b.-');hold off;

%xlim([0  100]);

%y1 = filter(b,a,x1col);
%y2 = filter(b,a,x2col);
%y3 = filter(b,a,x3col);
%y4 = filter(b,a,x4col);

% since inverse of 200THz is 5 fsec, you can take square !!!
%subplot(212);
%stem(tcol,(y1).^2,'r.-','filled');hold on;
%stem(tcol,(y2).^2,'g.-','filled');
%stem(tcol,(y3).^2,'b.-','filled');hold off;




%xlim([0 10])



