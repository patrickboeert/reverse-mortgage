classdef simulation < handle
% SIMULATION  generates & simulates reverse mortgage portfolio based on the
% simulated paths of a previously estimated ARX(p) and MCMC time series 
% model

   properties
      
      % object instances
      model         = [];   % model
      
      % monte carlo markov chain - interest rates 
      r_prior       = 5.508789954; % initial rate of 1-year rate for MCMC simulation (mean of Bundesbank WZ9808 9/1972-2/2009)
      rl_prior      = 6.733630137; % initial rate of 10-year rate for MCMC simulation (mean of Bundesbank WZ9826 9/1972-2/2009)          
      margin        = 1.5;         % lender's profit margin
      r_sim         = [];          % simulated interest rates from MCMC, nRUN x (12xnYears)      
      state         = [];          % markov chain states (upper, mid, lower)      
      trans         = [];          % markov chain transition matrix
      
      % arx model - house prices  
      mu_sim      = [];     % simulated mu from ARX(p), nRUN x nYears
      arx         = [];     % estimated ARX(p) model of house returns 
      AIC         = [];     % Akaike information criterion
      reject      = [];     % results for LLR hypothesis testing of lag
      pvalue      = [];     % p-values for LLR hypothesis testing   
                 
      % demand data
      demand      = [];     % demand pattern (nContracts)
           
      % aggregated cashflows
      fund                = [];   % fund: cum(premiums)-claims; (12 * nYears x nRUN)
      claim               = {};   % claim: regular, claim 1, claim 2, {nRUN} cell of [12*nYears x nContracts]
      b_e                 = [];   % outstanding balance (end of period); (12 * nYears x nRUN)
      advances            = [];   % advances (12 * nYears x nRUN)    
      premiums            = [];   % premiums (12 * nYears x nRUN)
      interest            = [];   % interest (12 * nYears x nRUN)

      % other results
      period_origination  = [];   % period of origination
      period_terminated   = [];   % period of termination (due to death, moveout)
      reason_terminated   = {};   % reason for termination: 'death', 'moveout'
      type_terminated     = {};   % 'Payoff', 'Claim 1 Type', 'Claim 2 Type' (assignment)      
      period_assigned     = [];   % period of assignment to insurer       
      
      % simulation
      nRUN          = 10;       % # of simulations      
      nContracts    = 15;       % # of contracts in each single portfolio
      nYears        = 60;        % # simulated years in every simulation 
      
      % status information
      madedir          = [];  % has directory already been created?
      savedorigination = [];  % have originated contracts been stored?      
      savedrun         = [];  % has simulation run already been stored?
      savedsim         = [];  % has simulation already been stored?      

      % storage
      path   = ''; % storage path
      dir    = ''; % storage directory
      file   = ''; % file name
      vars   = {}; % var names for run
       
   end % end

   methods   
     
     %% CONSTRUCTOR
     function obj = simulation(model)  
       % SIMULATION  constructs a simulation object for a reverse mortgage
       % portfolio simulation
                  
         % model configuration
         obj.model = model;           % set model instance      
         
         % time series simulation
         obj = get_arx(obj);          % set ARX time series instance
         obj = mcmc(obj,6,6);       % simulate short-term interest rates
         obj = sim_arx(obj,...        % simulate house appreciation rates
                       obj.arx,...
                       obj.r_sim);                                
         
         % set statuses
         obj.madedir          = false;
         obj.savedsim         = false;  
         obj.savedorigination = false;         
         obj.savedrun         = false(obj.nRUN,1); 
          
         % set storage parameters         
         obj.path = strcat(obj.model.basedir);         
         obj.dir  = strcat('simulation_', datestr(clock, 30));         
         obj.file = strcat(obj.dir, '.mat');
         
         % set variable names for storage
         prefix   = 'run_';           
         suffix   = strtrim(cellstr(int2str((1:obj.nRUN)')));         
         obj.vars = strcat(prefix, suffix);  
                                                       
     end % end SIMULATION
     
     %% SIMULATION
     function obj = get_arx(obj) 
      % GET_ARX  defines and estimates ARX(p) model
      % for house price returns
      %
      % vgxset       - defines model specification
      % vgxdisp      - displays model information
      % vgxvarx      - estimate calibration of VAR or VARX models
      % vgxcount     - count number of estimated parameters
      % lratiotest   - do likelihood ratio test on model specifications
      % aicbic       - Akaike information criterion
      % vgxqual      - qualify multivariate time series process (stability)
      % vgxplot      - plot multivariate time series process
      % vgxsim       - simulate multivariate time series process

      %% parameter setup
      desiredp    = 4;           % selected lag p of AR(p)          
      maxP        = 12;          % max. p of AR(p) to estimate
      nSeries     = 1;           % # of time series
      nX          = 1;           % # of exg. time series
      alpha       = 0.01;        % alpha level for log-likelihood test

      %% load data
      p  = (1:maxP)';   % number of lags p         
      load vardata.mat vardata
      Y   = vardata{1}.index(:,12);       % real estate price index   
      X   = vardata{2}.index(:,2);        % 1-year risk free rate 

      %% estimate (V)ARX(p) model      
      % preallocate
      spec = vgxset.empty(length(p),0);
      
      % define & estimate
      for i = 1:maxP
        
        % define model
        spec(i) = vgxset('n',nSeries,...
                         'nAR',p(i),...
                         'nX',nX,...
                         'Constant', true,...                              
                         'Series',{vardata{1}.name{12}});
                                                                           

        % get log-likelihood statistic to determine number of lags        
        % with lratiotest, estimation with maximum likelihood
        % (for further estimation, see vgxvarx help)
        Ypresample = Y(1:p(i),:);
        Ysample    = Y((p(i)+1):end,:);
        Xsample    = num2cell(X((p(i)+1):end,:));
        
        [Estspec(i),...
          EstStdErrors(i),...
           LLF(i),...
            W{i}] = vgxvarx(spec(i),...
                            Ysample,...
                            Xsample,...
                            Ypresample);

        % get number of all and unrestricted parameters
        [NumParam(i),NumActive(i)] = vgxcount(Estspec(i));
        
      end

      %% test for best lag structure      
      % get test results of likelihood ratio test at level alpha
      reject = NaN(length(p),length(p));
      pvalue = NaN(length(p),length(p));

      for i = 1:length(p)
        if p(i)>0
          for j = 1:(i-1)
            % reject is vector of Boolean decisions the same size as 
            % NullLLF. A 0 indicates acceptance of the restricted model under the null
            % hypothesis. 1 indicates rejection of the restricted, 
            % null hypothesis model relative to the unrestricted
            % alternative associated with BaseLLF.
            %
            % [H,pValue,Ratio,CriticalValue] = ...
            % lratiotest(BaseLLF,NullLLF,DoF,Alpha)
            
            [reject(i,j), pvalue(i,j)] =...
              lratiotest(LLF(i),LLF(j),NumActive(i) - NumActive(j),alpha);
          end
        end
      end

      % get Akaike information criterion
      AIC = aicbic(LLF,NumActive);
     
      %% return results
      obj.arx     = Estspec(desiredp);
      obj.AIC     = AIC;
      obj.reject  = reject;
      obj.pvalue  = pvalue;
      
     end % end GET_ARX
     function obj = sim_arx(obj, arx, r_sim)
      % SIM_ARX  simulates ARX(p) model data of house price returns given
      % exogenous, simulated paths of interest rates
      %
      % INPUT:  arx - ARX(p) time series model stored in object
      %       r_sim - simulates short-term interest rate paths
      %
      % OUTPUT: object 'simulation' with updated mu_sim      
     
        % average interest rate simulation data to yearly values
        subs     = floor((1:1/12:(obj.nYears+1-1/12))');
        r_sim_y  = zeros(obj.nRUN, obj.nYears);
        for run = 1:obj.nRUN
          vals = r_sim(run,:)';
          r_sim_y(run,:) = accumarray(subs,vals)/12;
        end
        X = num2cell(r_sim_y');       
        
        % load presample data       
        load vardata.mat vardata
        Y   = vardata{1}.index(:,12);         % real estate price index   
        Ypresample = Y(end-arx.nAR:end,:);    % presample data
        
        % simulate ARX for desired p
        mu_sim = vgxsim(arx,...                  % model
                        obj.nYears,...           % # obs
                        X,...                    % exogenous series
                        Ypresample,...           % presample 
                        [],...                   % no innovations
                        obj.nRUN);               % # paths
        
        mu_sim = squeeze(mu_sim); mu_sim = mu_sim';
                      
        % return simulated mu
        obj.mu_sim = mu_sim;
                        
     end % end SIM_ARX
     function obj = mcmc(obj, ir_dim, ir_diff_dim)
       % MCMC  Monte Carlo Markov Chain model for short-term interest
       % rates; estimates transition matrix from German interest rate 
       % data and simulates interest rate paths 
       % 
       % INPUT: 
       %    ir_dim       - # of states for interest rate level for estimation of
       %                  transition matrix
       %    ir_diff_dim - # of states for interest rate change for 
       %                  estimation of transition matrix
       %
       % OUTPUT: 'simulation' object
       
       % load data
       load ir_bundeswertpapier
       ir = bundeswertpapier(2).index;    % 1-year risk free rate
       ir_diff = diff(ir);                % monthly change
       clear bundeswertpapier
       
       % construct ir grid
       ir_step   = (max(ir)-min(ir))/(ir_dim+1);   % steps
       ir_state  = ((min(ir)+ir_step):ir_step:(max(ir)-ir_step))';     % states
       
       ir_lower  = NaN(size(ir_state)); 
       ir_lower(2:end) = ir_state(2:end);        % lower bound of state (first is NaN)
       ir_upper  = NaN(size(ir_state)); 
       ir_upper(1:end-1) = ir_state(2:end);    % upper bound of state (last is NaN)                       
       
       ir_mid    = NaN(size(ir_state));                
       ir_mid(2:end-1) = ir_lower(2:end-1)+(ir_upper(2:end-1)-ir_lower(2:end-1))/2; % mid value of state (first & last NaN)
       step = unique(diff(ir_mid)); step = step(1);
       ir_mid(1) = ir_mid(2)-step; ir_mid(end) = ir_mid(end-1)+step;
           
       % construct ir_diff grid
       ir_diff_step   = 4.4/(ir_diff_dim-1);
       ir_diff_state  = (-2.2:ir_diff_step:2.2)';         % states
       
       ir_diff_lower  = NaN(size(ir_diff_state)); 
       ir_diff_lower(2:end) = ir_diff_state(2:end);        % lower bound of state (fir_diffst is NaN)
       ir_diff_upper  = NaN(size(ir_diff_state)); 
       ir_diff_upper(1:end-1) = ir_diff_state(2:end);    % upper bound of state (last is NaN)                                       
       
       ir_diff_mid    = NaN(size(ir_diff_state));                
       ir_diff_mid(2:end-1) = ir_diff_lower(2:end-1)+(ir_diff_upper(2:end-1)-ir_diff_lower(2:end-1))/2; % mid value of state (first & last NaN)
       step = unique(diff(ir_diff_mid)); step = step(1);
       ir_diff_mid(1) = ir_diff_mid(2)-step; ir_diff_mid(end) = ir_diff_mid(end-1)+step;
              
       % determine logical state of data
       % ir is in state i
       ir = ir(1:end-1);  % skip value to match periods       
       in_level    = cell(ir_dim,1);  % cell of logicals       
       
       in_level{1} = ir<ir_upper(1);
       for i = 2:ir_dim-1              
         in_level{i} = ir>=ir_lower(i) & ir<ir_upper(i);  % in level state
       end
       in_level{ir_dim} = ir>=ir_lower(ir_dim);
       
       % ir_diff is in state i
       in_diff = cell(ir_diff_dim,1);  % cell of logicals       
       
       in_diff{1} = ir_diff<ir_diff_upper(1);
       for i = 2:ir_diff_dim-1              
         in_diff{i} = ir_diff>=ir_diff_lower(i) & ir_diff<ir_diff_upper(i);  % in diff state
       end
       in_diff{ir_diff_dim} = ir_diff>=ir_diff_lower(ir_diff_dim);
       
       % compute transition matrix probs
       trans  = zeros(ir_dim, ir_diff_dim);   % probability of transition      
       Ntrans = zeros(size(trans));           % frequency of transition       
       
       for i=1:ir_dim
         for j=1:ir_diff_dim
           Ntrans(i,j) = sum(in_level{i} & in_diff{j});
           trans(i,j)  = Ntrans(i,j)/sum(in_level{i});
         end
       end
       
       % simulate data with functionality from Kevin Murphy's 'pmtk3'
       % toolbox; see Google Code page http://code.google.com/p/pmtk3/
       
       % parameters
       prior  = obj.r_prior;           % start value
       T      = 12*obj.nYears;         % # months to simulate
       trials = obj.nRUN;              % # of paths to simulate
       ir_lower(1)   = -inf;           % for logical comparison
       ir_upper(end) = inf;            % for logical comparison
       
       % simulate interest rate paths 
       S = zeros(trials,T);   % state matrix
       V = zeros(trials,T);   % value matrix
       C = zeros(trials,T);   % change state to current period matrix
              
       for trial = 1:trials                 
         % set initial state
         S(trial,1) = find(prior>=ir_lower & prior<ir_upper);
         V(trial,1) = prior;
          for t = 2:T            
            C(trial,t) = sampleDiscrete(trans(S(trial,t-1),:));
            V(trial,t) = V(trial,t-1) + ir_diff_mid(C(trial,t));
                if V(trial,t) < min(ir), V(trial,t) = min(ir);
                elseif V(trial,t) > max(ir), V(trial,t) = max(ir); end
            S(trial,t) = find(V(trial,t)>=ir_lower & V(trial,t)<ir_upper);
          end          
       end         
       
       % return outputs
       obj.trans       = trans;           % transition matrix
       obj.r_sim       = V;               % simulated interest rate paths
       obj.state.ir_upper = ir_upper;
       obj.state.ir_mid   = ir_mid;       
       obj.state.ir_lower = ir_lower;                                  
       obj.state.ir_diff_upper = ir_diff_upper;
       obj.state.ir_diff_mid   = ir_diff_mid;       
       obj.state.ir_diff_lower = ir_diff_lower; 
       
     end % end MCMC
     function obj = get_demand(obj, fraction)
       % GET_DEMAND  construct demand model for simulation
       %
       % INPUT
       %      fraction: 2-dim vector of [tenure,term] payment scheme
       %                share in contract characteristics, e.g. [0.3,0.7]

       % error checking
       if sum(fraction)~=1
         error('ReverseMortgage:Simulation:GetDemand:Fraction',...
               'Elements of Fraction need to equal 1!');         
       end
       
       % preallocate cells
       age           = cell(obj.nContracts,1);
       gender        = cell(obj.nContracts,1);
       house         = cell(obj.nContracts,1);
       ir_exp        = cell(obj.nContracts,1);
       util          = cell(obj.nContracts,1);
       
       % expected rate = 10-year rate (long-term mean) + margin
       exp_rate = obj.rl_prior/100 + obj.margin/100;
       
       % sample from HECM data age, collateral and gender       
       
       % load data
       load age; load collateral; load gender            
       n = length(hecm_age);
       
       % find valid gender values
       ok  = find(strcmp('male',hecm_gender)|strcmp('female',hecm_gender));
       ind = unidrnd(length(ok),obj.nContracts,1);       
       
       for i = 1:obj.nContracts
           age{i}     = hecm_age(ok(ind(i)));
           gender{i}  = hecm_gender{ok(ind(i))};
           house{i}   = hecm_collateral(ok(ind(i)));
           ir_exp{i}  = exp_rate;   
           util{i}    = 1;           
       end
       
       clear hecm_age hecm_gender hecm_collateral
                      
       % preallocate demand structure (except scheme)

       obj.demand = ...
          struct('age', age,...
                 'gender', gender,...
                 'house', house,...
                 'util', util,...
                 'ir_exp', ir_exp);              

       % payment schemes
       nTenure = floor(fraction(1)*obj.nContracts); 
       nTerm   = ceil(fraction(2)*obj.nContracts);
       
       years = unidrnd(15,obj.nContracts,1)+5; % uniformly random periods between 5 and 20 years
       
              
       for i = 1:obj.nContracts
         if i <= nTenure
           obj.demand(i).scheme.name    = 'tenure';
         elseif i > nTenure
           obj.demand(i).scheme.name    = 'term';
           if (12*age{i} + 12*years(i)) < 1464               % can receive full term
               obj.demand(i).scheme.periods = 12*years(i);
           else                                              % must reduce term
               obj.demand(i).scheme.periods = 1464-12*age{i};
           end
         end
       end       

     end % end GET_DEMAND      
     
     %% BEHAVIOUR
     function obj = originate(obj)
       % ORIGINATE originate and save originated contracts without history
       
       % error checking: no origination if simulation already done
       if obj.savedorigination || obj.savedsim
         error('ReverseMortgage:Simulation:Originate:HasBeenOriginated',...
               ['This portfolio already has been originated or simulated. Another'... 
               ' origination is not possible. Create a new simulation instance!']);
       end   
       
       % create directory
       if ~obj.madedir, successful = mkdir(strcat(obj.path,obj.dir,'\'));
         if successful, obj.madedir = true;
         else error('ReverseMortgage:Simulation:Originate:CouldNotCreateDir',...
                    'Did not succeed at creating directory for storage !');  
         end
       end     
                 
       % originate contracts & save on disk     
       if ~obj.savedorigination
         
         % preallocate
         origination = cell(obj.nContracts,1);
         for i = 1:obj.nContracts, 
           origination{i} = contract.empty(1,0); 
         end
                
         % originate
         for i = 1:obj.nContracts         
                     
             % create contract
             origination{i} =...
               contract(obj.demand(i).age,...
                        obj.demand(i).gender,...
                        obj.demand(i).house,...
                        obj.demand(i).ir_exp,...
                        obj.demand(i).scheme,...
                        obj.model,...
                        obj.demand(i).util); 
                      
             % print status
             fprintf(1, 'contract %d of %d originated\n',...
                     i, obj.nContracts);                                                             
         end
         
         % save on disk
         save(strcat(obj.path, obj.dir, '\', 'origination.mat'),...
              'origination');          
         obj.savedorigination = true; 
         clear origination
         
       % save simulation object
         saveaway(obj);           
         obj.savedsim = true;                
       end       
       
     end % end ORIGINATE     
     function obj = simulate(obj)
       % SIMULATE  simulate contracts and store on disk 
       
       %% error checking: no simulation if there has been one before
       if obj.madedir && obj.savedsim && any(obj.savedrun)
         error('ReverseMortgage:Simulation:Simulate:HasBeenSimulated',...
               ['This portfolio already has been simulated. Another'... 
               ' simulation run is not possible, data has already been',...
               ' generated. Create a new simulation instance!']);
       end   
       
       %% error checking: no simulation without origination
       if ~obj.savedorigination
         error('ReverseMortgage:Simulation:Simulate:NoOrigination',...
               ['Need originated contracts before I can start to'... 
               ' simulate those contracts! Use method originate before',...
               ' simulate!']);
       end                                         
                    
       %% simulate contracts 
       for run = 1:obj.nRUN       % loop over sim runs
         
         tic
         load(strcat(obj.path, obj.dir, '\', 'origination.mat'));  % load var origination
         
         for rm = 1:obj.nContracts

           % simulate
           origination{rm} = origination{rm}.lifetime(obj, run); % termination events
           origination{rm} = origination{rm}.get_house_simulated(obj, run); % house price path
           origination{rm} = origination{rm}.sim_series(obj, run); % simulated time series paths    
           origination{rm} = origination{rm}.assignment; % determine time of assignment event                                
           
           % print status
           fprintf(1, 'contract %d/%d of run %d/%d simulated\n',...
                  rm, obj.nContracts,...
                  run, obj.nRUN);   

         end
         
         % store away run
         simulated = origination;
         clear origination
         save(strcat(obj.path, obj.dir, '\', char(obj.vars{run})), 'simulated');            
         obj.savedrun(run) = true;         
         clear simulated        
         
        toc

       end              
       
       %% save simulation object
       saveaway(obj);           
       obj.savedsim = true;        

     end % end SIMULATE
     function obj = close(obj)
       % CLOSE  compute claim to insurance fund and save run
       % variables
       
       %% error checking: have portfolios been simulated?
       if ~all(obj.savedrun)
         error('ReverseMortgage:Simulation:Close:PortfoliosNotSimulated',...
              ['To evaluate the aggregated values for simulated portfolios', ...
               ' you need to first simulate portfolio data!']);
       end
         
       %% loop over runs
       for run = 1:obj.nRUN
         
         % load simulation data
         load(strcat(obj.path, obj.dir, '\', obj.vars{run}));  % load simulation data
         data = simulated;
         clear simulated
         
         % preallocate         
         period_origination     = zeros(obj.nContracts,1);
         period_terminated      = zeros(obj.nContracts,1);
         type_terminated        = cell(obj.nContracts,1);
         reason_terminated      = cell(obj.nContracts,1);
         age_for_calculation    = zeros(obj.nContracts,1);        
         gender_for_calculation = cell(obj.nContracts,1);   
         plf_exogenous          = zeros(obj.nContracts,1); 
         house_max_claim_amount = zeros(obj.nContracts,1);      
         house_initialvalue     = zeros(obj.nContracts,1);      
         util_pl                = zeros(obj.nContracts,1);           
         premiums_expected_sim  = zeros(obj.nContracts,1);
         loss_expected_sim      = zeros(obj.nContracts,1);
         utilization_sim        = zeros(obj.nContracts,1);
         premiums_expected_pay  = zeros(obj.nContracts,1);
         loss_expected_pay      = zeros(obj.nContracts,1);
         utilization_pay        = zeros(obj.nContracts,1);                              
         payment_monthly        = zeros(obj.nContracts,1);
         
         house_simulated        = sparse(12*obj.nYears,obj.nContracts);         
         b_s                    = sparse(12*obj.nYears,obj.nContracts);
         b_e                    = sparse(12*obj.nYears,obj.nContracts);
         advances               = sparse(12*obj.nYears,obj.nContracts);
         interest               = sparse(12*obj.nYears,obj.nContracts);
         premiums               = sparse(12*obj.nYears,obj.nContracts); 
         claim                  = sparse(12*obj.nYears,obj.nContracts);           
           
         disp(run);
         
           %% loop over contracts          
           for rm = 1:obj.nContracts
           
           % fill data                     
           period_origination(rm)     = data{rm}.period_origination;
           period_terminated(rm)      = data{rm}.period_terminated;
           reason_terminated{rm}      = data{rm}.reason_terminated;
           age_for_calculation(rm)    = data{rm}.age_for_calculation;      
           gender_for_calculation{rm} = data{rm}.gender_for_calculation;  
           plf_exogenous(rm)          = data{rm}.plf_exogenous;
           house_max_claim_amount(rm) = data{rm}.house_max_claim_amount;    
           house_initialvalue(rm)     = data{rm}.house_initialvalue;
           util_pl(rm)                = data{rm}.util_pl;          
           premiums_expected_sim(rm)  = data{rm}.acc_simulation.premiums_expected;
           loss_expected_sim(rm)      = data{rm}.acc_simulation.loss_expected;
           utilization_sim(rm)        = data{rm}.acc_simulation.utilization;
           premiums_expected_pay(rm)  = data{rm}.acc_payment.premiums_expected;
           loss_expected_pay(rm)      = data{rm}.acc_payment.loss_expected;
           utilization_pay(rm)        = data{rm}.acc_payment.utilization;                             
           payment_monthly(rm)        = data{rm}.acc_simulation.payment_monthly;
         
           st = data{rm}.period_origination + 1;       
           fi = data{rm}.period_terminated; 

           house_simulated(:,rm)         = data{rm}.house_simulated;  
           b_s(1:data{rm}.s(fi),rm)      = data{rm}.acc_simulation.b_s(st:fi);
           b_e(1:data{rm}.s(fi),rm)      = data{rm}.acc_simulation.b_e(st:fi);
           advances(1:data{rm}.s(fi),rm) = data{rm}.acc_simulation.advances(st:fi);
           interest(1:data{rm}.s(fi),rm) = data{rm}.acc_simulation.interest(st:fi);
           premiums(1:data{rm}.s(fi),rm) = data{rm}.acc_simulation.premiums(st:fi);

           % determine type and value of payoff
             % arrange data         
             demand     = data{rm}.acc_simulation.b_e(fi);
             collateral = data{rm}.house_simulated(data{rm}.s(fi))...
                              .* (1-data{rm}.model.salesexpense);

             % find payoff type & payoff value
             pay = min(demand, collateral);
             switch pay
               case demand                             
                 type_terminated{rm} = 'Payoff';
                 claim(data{rm}.s(fi),rm) = 0;
               case collateral
                 type_terminated{rm} = 'Claim';
                 claim(data{rm}.s(fi),rm) = demand - collateral;
             end
             
           end % end loop contracts
         
           %% reduce dimensionality
           b_e      = nansum(b_e,2);
           b_s      = nansum(b_s,2);
           advances = nansum(advances,2); 
           house_simulated = nansum(house_simulated,2);
           interest = nansum(interest,2);    
           premiums = nansum(premiums,2);                                   
           
           %% save data for run
           % save away data and clear
           save(strcat(obj.path, obj.dir, '\', 'vars_', obj.vars{run},'.mat'),...
              'period_terminated',...
              'period_origination',...
              'type_terminated',...
              'reason_terminated',...
              'age_for_calculation',...
              'gender_for_calculation',...
              'plf_exogenous',...
              'house_max_claim_amount',...
              'house_initialvalue',...
              'util_pl',...
              'premiums_expected_sim',...
              'loss_expected_sim',...
              'utilization_sim',...
              'premiums_expected_pay',...
              'loss_expected_pay',...
              'utilization_pay',...
              'payment_monthly',...
              'house_simulated',...
              'b_s',...
              'b_e',...
              'advances',...
              'interest',...
              'premiums',...
              'claim');   
            
           clear('period_terminated',...
              'period_origination',...
              'type_terminated',...
              'reason_terminated',...
              'age_for_calculation',...
              'gender_for_calculation',...
              'plf_exogenous',...
              'house_max_claim_amount',...
              'house_initialvalue',...
              'util_pl',...
              'premiums_expected_sim',...
              'loss_expected_sim',...
              'utilization_sim',...
              'premiums_expected_pay',...
              'loss_expected_pay',...
              'utilization_pay',...
              'payment_monthly',...
              'house_simulated',...
              'b_s',...
              'b_e',...
              'advances',...
              'interest',...
              'premiums',...
              'claim',...
              'data');              
           
         
       end % end loop run
     end % end CLOSE
     function obj = aggregate(obj)
       % AGGREGATE aggregate simulated time series variables into
       % simulation object
       
       % preallocate variables
       obj.claim      = cell(obj.nRUN,1);
       obj.b_e        = sparse(12*obj.nYears,obj.nRUN);
       obj.advances   = sparse(12*obj.nYears,obj.nRUN);    
       obj.premiums   = sparse(12*obj.nYears,obj.nRUN);
       obj.interest   = sparse(12*obj.nYears,obj.nRUN);
       obj.fund       = sparse(12*obj.nYears,obj.nRUN);              
       
       % loop over runs
       for run = 1:obj.nRUN
       
         % load data
         tic
         load(strcat(obj.path, obj.dir, '\', 'vars_', obj.vars{run}));  % load run variables

         % assign data
         obj.claim{run,1}      = claim;
         obj.b_e(:,run)        = b_e;
         obj.advances(:,run)   = advances;
         obj.premiums(:,run)   = premiums;
         obj.interest(:,run)   = interest;
         obj.fund(:,run)       = cumsum(premiums) - cumsum(nansum(claim,2));
         
         % clear data
         clear('period_terminated',...
              'period_origination',...
              'type_terminated',...
              'reason_terminated',...
              'age_for_calculation',...
              'gender_for_calculation',...
              'plf_exogenous',...
              'house_max_claim_amount',...
              'house_initialvalue',...
              'util_pl',...
              'premiums_expected_sim',...
              'loss_expected_sim',...
              'utilization_sim',...
              'premiums_expected_pay',...
              'loss_expected_pay',...
              'utilization_pay',...
              'payment_monthly',...
              'house_simulated',...
              'b_s',...
              'b_e',...
              'advances',...
              'interest',...
              'premiums',...
              'claim');                       
          toc
       end 
       
       % save simulation object
       saveaway(obj);           
       obj.savedsim = true;          
       
     end % end AGGREGATE
                              
     %% STORAGE
     function saveaway(obj)  
       % SAVEAWAY  stores simulation object on disk        
       sim = obj;
       save(strcat(obj.path, obj.dir, '\', 'simulation.mat'),...
            'sim');                 
     end % end SAVEAWAY
     
   end % methods  
end % end classdef