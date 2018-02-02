function plotsimu(sim)

% arrange data
examplepath = 35;

claim = zeros(size(sim.b_e));
for i=1:100
 claim(:,i) = sum(sim.claim{i,1},2);   % cumulate claims
end

data{1}.data = sim.fund;         
data{1}.des = 'nominal net receivables of insurance fund';  % net receivables fund (cum. premiums - cumm. claims)
data{2}.data = cumsum(claim,1);  
data{2}.des = 'cumulated claims';         % cum. claims
data{3}.data = double(sim.b_e);          
data{3}.des = 'outstanding loan balance'; % outstanding loan balance
data{4}.data = sim.advances;     
data{4}.des = 'advances to borrowers';    % advances
data{5}.data = cumsum(sim.premiums,1);     
data{5}.des = 'cumulated insurance premiums'; % cumm premiums
data{6}.data = cumsum(sim.interest,1);     
data{6}.des = 'interest payments';        % interest

% compute data
for i = 1:6  
  
  % mean
  data{i}.mean = mean(data{i}.data,2);
  
  % upper and lower bound
  sorted = sort(data{i}.data,2);
  data{i}.lower = sorted(:,6);
  data{i}.upper = sorted(:,95);
  
  % example path
  data{i}.ex = data{i}.data(:,examplepath);  
end               

% plot data
years = (0:1/12:60-1/12)'; i = 6;

h(1) = plot(years,data{i}.mean, '-k'); hold on   % mean
h(2) = plot(years,data{i}.lower, '--b'); hold on  % lower
h(3) = plot(years,data{i}.upper, '--b'); % upper
%h(4) = plot(years,data{i}.ex);             % example path

set(h(:), 'LineWidth', 1);
    
  legend(h, {'Mean of 100 trials',...
             'Lower and Upper 5% Bound of 100 Trials'},...
             'Location','NorthWest');  
  xlabel('Year of Simulation'); 

end
       
