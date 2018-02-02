%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reproduces 2d kernel density estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function twodkernel
  clear all
% generate hecm data
  load pattern
  data=[collateral, age];
% call the routine
   n=2^7;
   
   %MAX=max(data,[],1); MIN=min(data,[],1); Range=MAX-MIN;
   %MAX_XY=MAX+Range/4; MIN_XY=MIN-Range/4;   
   MAX_XY = [1200000, 126];
   MIN_XY = [0,62];
   
  [bandwidth,density,X,Y]=kde2d(data,n,MIN_XY,MAX_XY);
% plot the data and the density estimate
  contour3(X,Y,density,50), hold off
  %plot(data(:,1),data(:,2),'r.','MarkerSize',1)
    
  % surf(X,Y,density,'LineStyle','none'), view([0,60])
  % colormap hot, hold on, alpha(.8)
  % set(gca, 'color', 'blue');
  % plot(data(:,1),data(:,2),'w.','MarkerSize',5)
end