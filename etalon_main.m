% etalon_main.m

 clear all; close all; warning on
 
 mm = 2^9;
 wl0 = 0.6500;	% um
 stepwl = 0.0001;	%um

 [wlcol, PTetacol, PRetacol] = etalon_func(mm, wl0, stepwl);
 
 fig1=figure(1);
 
 plot(wlcol,PTetacol,'r-', wlcol, PRetacol,'g-');


