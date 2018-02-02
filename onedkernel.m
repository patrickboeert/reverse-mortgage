%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reproduces kernel density estimates for age / collateral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function onedkernel
  clear all
% generate hecm data
  load pattern
  age(age==-1912)=63;
% call the routine
 [fc, xic, uc] = ksdensity(collateral, 'support', [0,2000000], 'function', 'pdf');
 [fa, xia, ua] = ksdensity(age, 'support', [62,130], 'function', 'pdf'); 
% plot
 subplot(2,1,1)
 plot(xic,fc);
 subplot(2,1,2)
 plot(xia,fa);
 
end