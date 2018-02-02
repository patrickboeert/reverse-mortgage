classdef probability
% PROBABILITY  computes customized probabilities needed by reverse mortgage
% pricing and simulation methods for a single reverse mortgage contract
%
%   PROBABILITY handles all probability data relevant for a specific
%   reverse mortgage contract. The probabilities are based on the specific
%   lifetable provided by the model settings.
%
%   USAGE:       OBJ = PROBABILITY(contract)
%
%   PROPERTIES:
%
%    p:              uncond. death probability (annual)
%    q:              uncond. survival probability (annual)
%
%    q_cond:         cond. death probability (annual)
%    p_cond:         cond. survival probability (annual)
%
%    t_cond_monthly: cond. contract termination probability (monthly, adj.)
%    s_cond_monthly: cond. contract survival probability (monthly, adj.)
%
%    t_death_cond_monthly:    cond. death termination probability (monthly)
%    s_death_cond_monthly:    cond. death survival probability (monthly)
%
%    t_moveout_cond_monthly:  cond. moveout termination probability (monthly)
%    s_moveout_cond_monthly:  cond. moveout survival probability (monthly)
%
%    loss:           analytical loss probability of HECM model (monthly)

   properties
       
       % associated reverse mortgage
       contract = [];           % reverse mortgage
       
       % lifetable data
       q = [];                
       p = [];                 
       q_cond = [];            
       p_cond = [];             
       
       % adj. contract lifetime probabilities
       t_cond_monthly = [];       
       s_cond_monthly = [];
       
       % termination probability by cause (death, moveout)
       t_death_cond_monthly = []; 
       s_death_cond_monthly = [];              
       t_moveout_cond_monthly = [];
       s_moveout_cond_monthly = [];
       
       % loss probability
       loss = [];
       
   end

   methods
       
       %% CONSTRUCTOR
       function obj = probability(contract)
           % PROBABILITY  constructor for probability instance
            
               % assign contract
               obj.contract = contract;    

               % import lifetable data
               obj = import(obj);

               % compute conditional probabilities
               obj = conditional(obj);

               % compute adj., monthly termination & survival probabilities
               obj = transform(obj);           
       
       end % end PROBABILITY
       
       %% IMPORT
       function obj = import(obj) 
           % IMPORT  Imports lifetable data
           % 
           % IMPORT computes uncond. death probabilities from lifetable 
           % and picks the right probability / age / gender combination
           % according borrower & coborrower profiles in contract
                
                % 1) import lifetable with uncond. death probs;
                %               age | men | women                    
                    importfile(obj.contract.model.lifetable);

                        q_male   = men;
                        q_female = women; 
                              
                 % 2) pick death probability 'probs' from lifetable given
                 %    contract's age, borrower, gender information
                    
                    % only one borrower
                    if isempty(obj.contract.coborrower_age)     % only one borrower

                        age    = obj.contract.borrower_age;
                        gender = obj.contract.borrower_gender;

                        switch obj.contract.borrower_gender
                            case 'male', probs = q_male;
                            case 'female', probs = q_female;
                        end

                    % two borrowers, pick the younger one according to his gender    
                    elseif ~isempty(obj.contract.coborrower_age) % two borrower, pick the younger
                        if obj.contract.borrower_age <= obj.contract.coborrower_age % borrower younger than coborrower

                            age    = obj.contract.borrower_age;
                            gender = obj.contract.borrower_gender;

                            switch obj.contract.borrower_gender
                                case 'male', probs = q_male;
                                case 'female', probs = q_female;
                            end
                        elseif obj.contract.coborrower_age < obj.contract.borrower_age % coborrower younger than borrower

                            age    = obj.contract.coborrower_age;
                            gender = obj.contract.coborrower_gender;

                            switch obj.contract.coborrower_age
                                case 'male', probs = q_male;
                                case 'female', probs = q_female;
                            end
                        end       
                    end
                    
                 % 3) assign values
                    obj.q = probs;
                    obj.contract.age_for_calculation = age;
                    obj.contract.gender_for_calculation = gender;
                    obj.contract.period_origination = 12*age;
                    
                    
       end % function IMPORT       
       
       %% COMPUTE PROPERTIES
       function obj = conditional(obj)
           % CONDITIONAL  Computes cond. death and survival probabilities
           %
           % CONDITIONAL computes cond. death and survival probabilities of
           % the contract, without move-out adjustments
           
                % 1) initialize vectors
                    p_cond    = zeros(length(obj.q),1);
                    q_cond    = zeros(length(obj.q),1);                     
                    
                    end_period_age = (1:(length(obj.q)))';  % age (end of period) 

                % 2) calculate pure (without move-out adjustment) 
                %    conditional death and survival probabilities q_cond,
                %    p_cond 
                    
                    % cond. survival probability
                    p_cond(end_period_age == obj.contract.age_for_calculation)...
                        = 1;

                    for year = (obj.contract.age_for_calculation + 1):(max(end_period_age)) 
                        p_cond(end_period_age == year) ...
                            = p_cond(end_period_age == (year-1))* obj.p(end_period_age == year);
                    end
                    
                    p_cond(end) = 0;
                    
                    % cond. death probability
                    for year = (obj.contract.age_for_calculation + 1):(max(end_period_age)); 
                       q_cond(end_period_age == year) ...
                           = p_cond(end_period_age == year-1) - p_cond(end_period_age == year);
                    end
                    
                % 3) assign values
                    obj.p_cond = p_cond;
                    obj.q_cond = q_cond;  
                    
       end % function CONDITIONAL       
       function obj = transform(obj) 
           % TRANSFORM  compute adjusted, monthly probabilities
           %
           % TRANSFORM computes adjusted, monthly termination and contract
           % survival probabilities based on the cond., yearly probabilities
           % provided by the lifetable. Adjustments for monthly values are
           % made by geometric interpolation.
           
              %% initialize vectors
                obj.s_cond_monthly = zeros(12 * length(obj.q),1);
                obj.t_cond_monthly = zeros(12 * length(obj.q),1);
                
                obj.s_death_cond_monthly   = zeros(12 * length(obj.q),1);                
                obj.t_death_cond_monthly   = zeros(12 * length(obj.q),1);
                
                obj.s_moveout_cond_monthly = zeros(12 * length(obj.q),1);               
                obj.t_moveout_cond_monthly = zeros(12 * length(obj.q),1);
                
              %% set transformation parameters
                m = ones(12 * length(obj.q),1) * obj.contract.model.moveout;    % moveout parameter
                                   
                month = (1 : (12 * length(obj.q)))';                            % lifetime months (at end of period)               
                start_period_month = 12 * obj.contract.age_for_calculation;     % month of start of contract
                running = (month >= start_period_month);                        % logical: contract in force, starts end of period
                
                j = (1/12:1/12:length(obj.q))';      
                j = floor(j);                    % age in years (at 'month')
                j(1:11) = 1;
                
                jplus = j+1;                     % age in years shifted one year (for indexing)
                jplus(end) = length(obj.q);
                    
                r_first_year  = (1:11)';         % months lived beyond current age count 'j' (end of period)
                r_other_years = (0:11)';                                                                               
                r = repmat(r_other_years,length(obj.q)-1,1);
                r = [r_first_year; r; 0];
                                             
              %% transform into monthly contract survival probabilities
              % through geometrics interpolation and adjustment for
              % moveouts
              
                obj.s_cond_monthly(month) = ...
                   running(month) .* ( obj.p_cond(j) .* ( ( obj.p_cond(jplus) ./ obj.p_cond(j) ) .^ (r(month)./12) ) ) .^ (1 + m(month));   
                 
                obj.s_death_cond_monthly =...
                   running(month) .* ( obj.p_cond(j) .* ( ( obj.p_cond(jplus) ./ obj.p_cond(j) ) .^ (r(month)./12) ) );
                 
                obj.s_moveout_cond_monthly =...
                   running(month) .* ( obj.p_cond(j) .* ( ( obj.p_cond(jplus) ./ obj.p_cond(j) ) .^ (r(month)./12) ) ) .^ (m(month));                  
                
              %% monthly contract termination probabilities                
                for i = (12 * obj.contract.age_for_calculation + 1):(12 * length(obj.q));
                  obj.t_cond_monthly(i)         = obj.s_cond_monthly(i-1)...
                                                    - obj.s_cond_monthly(i);
                  obj.t_death_cond_monthly(i)   = obj.s_death_cond_monthly(i-1)...
                                                    - obj.s_death_cond_monthly(i);
                  obj.t_moveout_cond_monthly(i) = obj.s_moveout_cond_monthly(i-1)...
                                                    - obj.s_moveout_cond_monthly(i);
                end 
                
                obj.t_moveout_cond_monthly(end-11) = 0;                                
                
       end % function TRANSFORM       
       function loss = get_loss(obj)
             % LOSS Compute loss probability of contract
             %
             % P(collateral value in t < outstanding balance in t), for
             % formula see 'price_collateral' of 'contract'
             
             [a, loss, b] = price_collateral(obj.contract.acc_basic.b_e, obj.contract);
             
         end % function GET_LOSS
                       
       %% GETS
       function p = get.p(obj)
           % GET.P  compute dependent p property if q has been set           
               if ~isempty(obj.q)
                   p = 1-obj.q;
               else
                   p = [];
               end              
       end % function GET P       
       function loss = get.loss(obj)
           % GET.LOSS  compute, if it has not been set yet           
               if isempty(obj.loss)
                   loss = get_loss(obj); return
               else
                   loss = obj.loss;
               end              
       end % function GET LOSS
                 
   end
end 
