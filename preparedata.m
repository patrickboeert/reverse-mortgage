function varinput = preparedata()

% PREPAREDATA  transforms original data for use in VAR(p) model

%% prepare data
  % load data
  load ir_bundeswertpapier
  load cpi
  load bulwien
  clear city westcities

  % create varinput
  varinput{1,1} = national;
  varinput{2,1} = bundeswertpapier(2);
  varinput{3,1} = bundeswertpapier(1);
  varinput{4,1} = cpi;

%% truncate & transform data
  
  % truncate for Bundeswertpapiere
  for i=2:3
    varinput{i}.month = {varinput{i}.month{41:(end-2),1}}';
    varinput{i}.index =  varinput{i}.index(41:(end-2),1);
  end

  % truncate for CPI
  varinput{4}.year = varinput{4}.year(28:end,1);
  varinput{4}.index = varinput{4}.index(28:end,1);
  
  % transform: compute yearly mean for Bundeswertpapier returns
  year  = (1:33)';
  start = 12*year-11;
  final = 12*year;

  for i=2:3   
    monthly = varinput{i}.index;

    % produce mean
    for y = 1:length(year)
      yearly(y,1) = mean(monthly(start(y):final(y)));
    end
    
    % reorder structure & reassign data
    varinput{i} = rmfield(varinput{i},'month');
    varinput{i}.name = {'year' 'index'};
    varinput{i}.index =[];  
    varinput{i}.index(:,1) = varinput{1}.index(2:end,1);    
    varinput{i}.index(:,2) = yearly;
    
  end
  
%% difference in logarithms of house prices
  varinput{1}.index = [varinput{1}.index(2:end,1), ...
                       100 * diff(log(varinput{1}.index(:,2:end)))];
                     
%% change cpi data
    varinput{4}.index = [varinput{4}.year(2:end), 100*diff(log(varinput{4}.index))];
    varinput{4}.description = varinput{4}.name;
    varinput{4} = rmfield(varinput{4},'name');
    varinput{4} = rmfield(varinput{4},'year');    
    varinput{4}.name = {'year' 'index'};
     
%% plot
  plot(datenum(varinput{1}.index(:,1),1,1), [varinput{1}.index(:,2), varinput{2}.index(:,2), varinput{3}.index(:,2), varinput{4}.index(:,2)]); 
  legend({'house' '1-year' '10-year' 'cpi'})
  grid on
  datetick('x','yyyy')	
                     
end

