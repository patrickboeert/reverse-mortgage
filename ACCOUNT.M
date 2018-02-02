classdef account
% ACCOUNT Payment account for a reverse mortgage contract
%
%   ACCOUNT = ACCOUNT(contract)
%
%   ACCOUNT manages the payment information for a reverse mortgage
%   contract. Depending on the payment scheme chosen by the borrower, a 
%   reverse mortgage contract may contain one or two payment accounts. 

   properties
               
       %% account information
       contract         = [];         % associated reverse mortgage contract
       type             = '';         % payment account type: 'tenure', 'term', 'loc'
       plf_mode         = '';         % mode of plf calculation: 'exogenous', 'endogenous'
 
       %% cashflow scenario
       scenario = '';                 % cashflow scenario: 'basic', 'payment', 'simulation'
       
       %% payment vectors
       principal_limit  = [];         % principal limit, T x 1       
       payment_monthly  = [];         % max. monthly payment at origination for type 'term', 'tenure'
       b_s       = [];                % start balance, T x 1
       b_e       = [];                % end balance, T x 1
       advances  = [];                % advances, T x 1
       interest  = [];                % interest, T x 1
       premiums  = [];                % premium fees, T x 1
       
       %% insurance
       premiums_expected = [];        % sum of expected, discounted premiums
       loss_expected     = [];        % sum of expected, dicounted losses; sum(E [B(t) - cond. H(t)])
       
   end
   
   properties (Dependent = true)
       
       %% insurance
       utilization = [];            % expected perc. utilization of the insurance := (exp.) loss / (exp.) premiums
   
   end

   methods
       
       %% CONSTRUCTOR
       function obj = account(contract, varargin)
           % ACCOUNT   creates new payment account for reverse mortgage

               %% get: number of accounts
               if any(strcmp(contract.payment_plan, {'tenure', 'term', 'loc'}))
                   dim = 1;
               elseif any(strcmp(contract.payment_plan, {'modifiedterm', 'modifiedtenure'}))
                   dim = 2;
               end
               
               %% handle 1x1 account
               if dim == 1
                   
                   % assign account type
                   obj(1).type = contract.payment_plan;
                   
                   % set plf_mode (because computed from plf_exogenous)
                   obj(1).plf_mode = 'exogenous';
                   
                   % assign contract
                   obj(1).contract = contract;                   
                   
                   % compute principal limit:
                   %  1) take from contract net pl (for 'term, tenure, loc')
                   %  2) but only fractional utilitzation level (for 'term,
                   %     tenure')
                   switch contract.payment_plan
                       case {'term','tenure'}
                          pl = obj.contract.util_pl * obj.contract.principal_limit_init_net;                   
                       case 'loc'
                          pl = obj.contract.principal_limit_init_net;
                   end
                               
                   obj(1).principal_limit = get_principal_limit(obj(1), pl); 
                           
                   % compute maximum monthly payment, (if not 'loc' scheme)
                   switch contract.payment_plan
                       case 'loc'
                               obj(1).payment_monthly = NaN;                                                    
                       case 'term'
                               obj(1).payment_monthly =...
                                   get_payment_monthly(obj, obj.principal_limit, 'term', obj.contract.payment_periods);
                       case 'tenure'
                               obj(1).payment_monthly =...
                                   get_payment_monthly(obj, obj.principal_limit, 'tenure');
                   end                   
               
              %% handle 2x1 account (first account)
              elseif dim == 2 && isempty(varargin)
                   
                   % assign account type
                   switch contract.payment_plan
                       case 'modifiedterm'
                           obj(1).type = 'term';                       
                       case 'modifiedtenure'
                           obj(1).type = 'tenure';
                   end

                   % set plf_mode (because computed from plf_exogenous)
                   obj(1).plf_mode = 'exogenous';
                   
                   % assign contract
                   obj(1).contract = contract;
                   
                   % compute principal limit:
                   %  1) take from contract net pl (for 'modifiedterm', 'modifiedtenure')
                   %  2) but only fractional utilitzation level for main account                                  
                   pl = obj.contract.util_pl * (obj.contract.principal_limit_init_net - obj.contract.loc_amount);
                   obj(1).principal_limit = get_principal_limit(obj(1), pl); 

                   % compute maximum monthly payment, (if not 'loc' scheme)
                   switch obj(1).type
                           case 'term'
                               obj(1).payment_monthly =...
                                   get_payment_monthly(obj(1), obj(1).principal_limit, 'term', obj(1).contract.payment_periods);
                           case 'tenure'
                               obj(1).payment_monthly =...
                                   get_payment_monthly(obj(1), obj(1).principal_limit, 'tenure');
                   end
                   
                   % call loop constructor
                   obj(2) = account(contract,'beenhere');

               %% handle 2x1 account (second account)                   
               elseif dim == 2 && ~isempty(varargin)
                   
                   % assign account type
                   obj.type = 'loc';

                   % set plf_mode (because computed from plf_exogenous)
                   obj.plf_mode = 'exogenous';
                   
                   % assign contract
                   obj.contract = contract;                   
                   
                   % compute principal limit:
                   %  1) all the rest for LOC account
                   %  2) and even more if fractional use
                   if obj.contract.util_pl == 1             % util factor standard: 1
                       pl = obj.contract.loc_amount;
                   elseif obj.contract.util_pl ~= 1         % other util factor
                       pl = obj.contract.principal_limit_net...
                            - obj.contract.util_pl...
                            *(obj.contract.principal_limit_init_net - obj.contract.loc_amount);                    
                   end
                       
                   obj.principal_limit = get_principal_limit(obj, pl); 
                   
                   % no monthly payment, this is LOC
                   obj.payment_monthly = NaN;
                   
               end % if

       end % function CONSTRUCTOR

       %% PAYMENT methods
       function obj = basic(obj)
           % BASIC   computes cashflows + insurance for the basic HECM model 
           %
           % BASIC computes the T x 1 cashflow vectors assuming a one-time
           % upfront lump-sum payment to the borrower in the amount of the
           % net principal limit. The borrower is assumed to outlive the
           % model's horizon.
           %
           % BASIC updates the expected premiums and expected losses using
           % the basic model and the basic model's plf.
           
               % update properties
               if length(obj) == 2, obj(2) = []; end
               obj.scenario        = 'basic';
               obj.plf_mode        = 'exogenous'; 
               obj.payment_monthly = NaN;           % not relevant anymore               
                                     
               % get cashflows for basic model
               cashflows = price_payments(obj.contract.plf_exogenous, obj.contract, 'exogenous');
           
               % assign cashflow properties
               obj.b_s      = cashflows.b_s;
               obj.b_e      = cashflows.b_e;
               obj.advances = cashflows.advances;
               obj.interest = cashflows.interest;
               obj.premiums = cashflows.premiums;               
               
               % compute insurance
               obj.premiums_expected = get_premiums_expected(obj); 
               obj.loss_expected     = get_loss_expected(obj);

       end % function BASIC
       function obj = payment(obj, plf_mode)
           % PAYMENT Computes cashflows + insurance for the HECM
           % payment model
           %
           % PAYMENT computes the cashflow vectors T x 1 for a cashflow
           % pattern consistent with the actual cashflows of a term or 
           % tenure loan if the borrower outlives the model's horizon.
           % PAYMENT also updates expected premiums and expected losses.
           %
           % OPTIONS:
           %
           %    PLF_MODE:
           %
           %        'exogenous': PAYMENT uses the plf factor given by the 
           %        equations of the basic pricing model. Assumed cashflow 
           %        pattern for the plf iteration is the maximal one-time
           %        upfront advance. Expected premiums and expected losses are
           %        computed using the basic model's plf as well. This is the
           %        way actual HECM contracts are computed.
           %
           %        'endogenous': PAYMENT uses a plf factor exlusively computed
           %        for this purpose. The plf is computed using the equations
           %        of the basic pricing model, but assuming the actual cashflow
           %        pattern demanded by the borrower. Expected premiums and
           %        expected losses are computed in the same way.
           % 
           % NOTE:
           %
           %        Outliving the horizon: In case of a tenure account, this means
           %        living to the end of the horizon. In case of a term account, 
           %        it means receiving payments until the end of the term with 
           %        the contract remaining in force up to the end of the
           %        horizon.

                %% error checking
                if strcmp(obj(1).scenario,'basic')
                    error('ReverseMortgage:Account:Payment:CannotChangeBasic',...
                        ['Cannot change previous payment plan BASIC to a payment plan '...
                        'of the HECM payment model. Create a new account instead !']);
                elseif nargin < 2
                    error('ReverseMortgage:Account:Payment:NumberOfInputs',...
                          'Too few input arguments !');                                               
                end
                
                %% not for LOC, because it doesn't make sense
                if any(strcmpi(obj(1).type,{'tenure', 'term'})) 
                
                    %% distinction between plf modes                
                    switch plf_mode

                        case 'exogenous' % we can fall back on the principal limit handed over to the account

                            % create monthly payment vectors
                            obj(1).b_s   = zeros(obj(1).contract.horizon_months,1);  
                            obj(1).premiums = zeros(obj(1).contract.horizon_months,1); 
                            obj(1).advances = zeros(obj(1).contract.horizon_months,1);
                            obj(1).interest = zeros(obj(1).contract.horizon_months,1);
                            obj(1).b_e      = zeros(obj(1).contract.horizon_months,1);

                            % set starting values: advances, premiums, end
                            % balance

                                % advances
                                obj(1).advances(12 * obj(1).contract.age_for_calculation + 1) =...
                                           obj(1).payment_monthly + obj(1).contract.closing_cost;

                                switch obj(1).type      
                                    case 'tenure'
                                      obj(1).advances((12*obj(1).contract.age_for_calculation+2):end) =...
                                          obj(1).payment_monthly;  

                                    case 'term'
                                      obj(1).advances((12*obj(1).contract.age_for_calculation+2):...
                                          (12*obj(1).contract.age_for_calculation+obj(1).contract.payment_periods)) =...
                                          obj(1).payment_monthly;                      
                                end

                                % premiums
                                obj(1).premiums(12 * obj(1).contract.age_for_calculation + 1) =...
                                    obj(1).contract.model.upfront * obj(1).contract.house_max_claim_amount; 

                                % end balance
                                obj(1).b_e(12 * obj(1).contract.age_for_calculation + 1) = ...
                                    obj(1).advances(12 * obj(1).contract.age_for_calculation + 1)...
                                    + obj(1).premiums(12 * obj(1).contract.age_for_calculation + 1);       

                            % recursive payment computation                    
                            for m = (12*obj(1).contract.age_for_calculation+2):obj(1).contract.horizon_months
                              obj(1).b_s(m)           = obj(1).b_e(m-1);
                              obj(1).interest(m)      = obj(1).b_s(m) .* (obj(1).contract.ir_expected_rate/12);             
                              obj(1).premiums(m)      = obj(1).contract.premiumstructure(m) .* obj(1).b_s(m);
                              obj(1).b_e(m)           = obj(1).b_s(m) + obj(1).advances(m)...
                                                     + obj(1).interest(m) + obj(1).premiums(m);
                            end
                            
                            % update insurance
                            obj(1).premiums_expected = get_premiums_expected(obj(1)); 
                            obj(1).loss_expected     = get_loss_expected(obj(1));

                        case 'endogenous' % we first need to compute and update plf

                          % get new plf first
                          
                                % get new plf
                                [obj(1).contract.plf_endogenous, newcontract] = price(obj(1).contract,'endogenous');

                                % get principal limit scalar (with new plf)
                                princ_lim_init_net =...
                                    obj(1).contract.plf_endogenous...
                                    .* obj(1).contract.house_max_claim_amount...
                                    - obj(1).contract.setaside(12 * obj(1).contract.age_for_calculation + 1);

                                % get principal limit vector (with new plf)
                                obj(1).principal_limit = get_principal_limit(obj(1), princ_lim_init_net);

                                % update monthly payment (with new plf)
                                switch obj(1).type
                                    case 'tenure'
                                        obj(1).payment_monthly = get_payment_monthly(obj(1), obj(1).principal_limit, 'tenure');
                                    case 'term'
                                        obj(1).payment_monthly = get_payment_monthly(obj(1), obj(1).principal_limit, 'term', obj(1).contract.payment_periods);
                                end

                          % now turn to computation of cashflows
                                % create monthly payment vectors
                                obj(1).b_s      = zeros(obj(1).contract.horizon_months,1);  
                                obj(1).premiums = zeros(obj(1).contract.horizon_months,1); 
                                obj(1).advances = zeros(obj(1).contract.horizon_months,1);
                                obj(1).interest = zeros(obj(1).contract.horizon_months,1);
                                obj(1).b_e      = zeros(obj(1).contract.horizon_months,1);

                                % set starting values: advances, premiums, end
                                % balance

                                    % advances
                                    obj(1).advances(12 * obj(1).contract.age_for_calculation + 1) =...
                                               obj(1).payment_monthly + obj(1).contract.closing_cost;

                                    switch obj(1).type      
                                        case 'tenure'
                                          obj(1).advances((12*obj(1).contract.age_for_calculation+2):end) =...
                                              obj(1).payment_monthly;  

                                        case 'term'
                                          obj(1).advances((12*obj(1).contract.age_for_calculation+2):...
                                              (12*obj(1).contract.age_for_calculation+obj(1).contract.payment_periods)) =...
                                              obj(1).payment_monthly;                      
                                    end

                                    % premiums
                                    obj(1).premiums(12 * obj(1).contract.age_for_calculation + 1) =...
                                        obj(1).contract.model.upfront * obj(1).contract.house_max_claim_amount; 

                                    % end balance
                                    obj(1).b_e(12 * obj(1).contract.age_for_calculation + 1) = ...
                                        obj(1).advances(12 * obj(1).contract.age_for_calculation + 1)...
                                        + obj(1).premiums(12 * obj(1).contract.age_for_calculation + 1);       

                                % recursive payment computation                    
                                for m = (12*obj(1).contract.age_for_calculation+2):obj(1).contract.horizon_months
                                  obj(1).b_s(m)           = obj(1).b_e(m-1);
                                  obj(1).interest(m)      = obj(1).b_s(m) .* (obj(1).contract.ir_expected_rate/12);             
                                  obj(1).premiums(m)      = obj(1).contract.premiumstructure(m) .* obj(1).b_s(m);
                                  obj(1).b_e(m)           = obj(1).b_s(m) + obj(1).advances(m)...
                                                           + obj(1).interest(m) + obj(1).premiums(m);
                                end
                                
                                % update insurance
                                obj(1).premiums_expected = get_premiums_expected(obj(1)); 
                                obj(1).loss_expected     = get_loss_expected(obj(1));                   

                    end % switch
                
                elseif strcmpi(obj(1).type,'loc')
                    
                    disp('Warning: Cannot use payment function on LOC Account! Nothing happened here!');
                    
                end % if
                
                %% update properties
                for i=1:length(obj)
                    obj(i).scenario        = 'payment';
                    obj(i).plf_mode        = plf_mode;                  
                end                
             
       end % function PAYMENT
       function obj = simulation(obj, simulation, run)
         % SIMULATE  computes cashflows for the HECM simulation of the
         % contract
         %
         % SIMULATE computes the cashflow vectors T x 1 resulting from the
         % projected payments, the stochastic termination event (death,
         % moveout) and the stochastic interest rate path. SIMULATE also 
         % updates the premiums and losses, which now have to be
         % interpreted as realized premiums and losses.

           %% error checking
            if strcmp(obj(1).scenario,'basic')
               error('ReverseMortgage:Account:Simulation:CannotChangeBasic',...
                    ['Cannot change previous payment plan BASIC to a payment plan '...
                     'of the HECM simulation model. Create a new account instead !']);
                   
            elseif ~any(strcmpi(obj(1).type,{'tenure', 'term'}))
               error('ReverseMortgage:Account:Simulation:CannotSimulateLOC',...
                     ['Cannot simulate LOC account! Need to implement ',...
                      'pattern for LOC drawing first']);
                    
            elseif isempty(simulation) || isempty(obj.contract.period_terminated)
              error('ReverseMortgage:Account:Simulation:NoSimulationData',...
                    ['For a simulation of an account, a simulation instance'...
                     ' needs to be set! Please create simulation directly'...
                     ' through simulation class.']);
            end              
 
           %% set loan interest rate
           ir_y = simulation.r_sim(run,:)/100 + simulation.margin/100;  % ass.: monthly adjustment
           ir_m = ir_y/12;
           ir   = ir_m;
                                  
           %% cashflow computation 
           
              % initialize vectors
              obj(1).b_s      = NaN(obj(1).contract.horizon_months,1);  
              obj(1).premiums = NaN(obj(1).contract.horizon_months,1); 
              obj(1).advances = NaN(obj(1).contract.horizon_months,1);
              obj(1).interest = NaN(obj(1).contract.horizon_months,1);
              obj(1).b_e      = NaN(obj(1).contract.horizon_months,1);

              % construct logical runtime period of contract
              month = (1:obj(1).contract.horizon_months)';
                     
              alive = (month > obj(1).contract.period_origination)...
                            & (month <= obj(1).contract.period_terminated);
              
              switch obj(1).type
                case 'tenure'
                  pay       = (month > obj(1).contract.period_origination)...
                            & (month <= obj(1).contract.period_terminated);
                case 'term'
                  pay       = (month > obj(1).contract.period_origination)...
                            & (month <= obj(1).contract.period_terminated)...
                            & (month <= obj(1).contract.period_origination + obj(1).contract.payment_periods);
              end
                                      
              ind = month(alive);  % index simulation month into lifetime month
              p   = month(pay);    % months of cash advances to borrower
              
              % set advances 
              obj(1).advances(p(1)) = obj(1).payment_monthly...
                                           + obj(1).contract.closing_cost;
              
              if any(pay ~= alive)
                obj(1).advances(p(2:end)) = obj(1).payment_monthly;
                obj(1).advances(p(end)+1:ind(end)) = 0;
              else
                obj(1).advances(p(2:end)) = obj(1).payment_monthly;
              end
                                
              % set start values
               % premiums
               obj(1).premiums(ind(1)) = obj(1).contract.model.upfront...
                                         * obj(1).contract.house_max_claim_amount; 
               % end balance
               obj(1).b_e(ind(1)) = obj(1).advances(ind(1))...
                                       + obj(1).premiums(ind(1));

           %% recursive payment computation                               
           for i = 2:length(ind)

             obj(1).b_s(ind(i))       = obj(1).b_e(ind(i-1));

             obj(1).interest(ind(i))  = obj(1).b_s(ind(i))...
                                              .* ir(i-1);

             obj(1).premiums(ind(i))  = obj(1).contract.premiumstructure(ind(i))...
                                              .* obj(1).b_s(ind(i));

             obj(1).b_e(ind(i))       = obj(1).b_s(ind(i))...
                                              + obj(1).advances(ind(i))...
                                              + obj(1).interest(ind(i))...
                                              + obj(1).premiums(ind(i));
           end             
                    
           %% update insurance
            obj(1).premiums_expected = get_premiums_expected(obj(1)); 
            obj(1).loss_expected     = get_loss_expected(obj(1));            
                
           %% update properties
            for i=1:length(obj)
               obj(i).scenario        = 'simulation';
               obj(i).plf_mode        = 'exogenous';                  
            end
            
       end % end SIMULATE
       
       %% PROPERTY methods
       function principal_limit = get_principal_limit(obj, princ_lim_init)
               % GET_PRINCIPAL_LIMIT   computes principal limit vector T x 1
               % (where T is model horizon in months), given initial
               % principal limit amount
               %
               % INPUT: 
               %   - initial amount of principal limit
               %
               % OUTPUT:
               %   - T x 1 vector of future principal limit values

                    % l: periods during which contract is in force
                    inforce = (1:obj.contract.horizon_months)'; 
                    inforce = (inforce > (12 * obj.contract.age_for_calculation));   

                    % t is array: number of months that contract is in
                    % force (at start of month when payments are made)
                    t = (1:obj.contract.horizon_months)' - (12 * obj.contract.age_for_calculation + 1); 
                    t =  t .* inforce;

                    % compute principal limit vector
                    principal_limit = inforce .* princ_lim_init .* ((1 + (obj.contract.ir_compounding/12)).^t);
                    principal_limit = single(principal_limit);
                
       end  % function GET_PRINCIPAL_LIMIT       
       function payment_monthly = get_payment_monthly(obj, principal_limit_vec, type, varargin)
              % GET_PAYMENT_MONTHLY Computes maximum monthly payments for term
              % contracts (with a given # of months) and tenure contracts 
              % at the time of the origination of the loan. Computation is 
              % based on a principal limit vector (e.g., that of the
              % account)
                     
                    % check inputs
                    if ~any(strcmp(type, {'term','tenure'}))
                        error('For calculation of maximum monthly payment, payment scheme needs to be term or tenure!');
                    end
                    if strcmp(type,'term') && isempty(varargin)
                        error('For maximum monthly payment with term scheme, please enter number of periods!');
                    end
                    
                    % get the whole vector of past and future max. monthly
                    % payments and assign value
                    if strcmp(type, 'term')
                        payments        = get_payment_monthly_vec(obj, principal_limit_vec, 'term', varargin{:});
                        payment_monthly = payments((12*obj.contract.age_for_calculation)+1);
                        return
                    elseif strcmp(type, 'tenure')
                        payments        = get_payment_monthly_vec(obj, principal_limit_vec, 'tenure');
                        payment_monthly = payments((12*obj.contract.age_for_calculation)+1);
                    end
                    
       end % function GET_PAYMENT_MONTHLY       
       function payment_monthly_vec = get_payment_monthly_vec(obj, principal_limit_vec, type, varargin)
              % GET_PAYMENT_MONTHLY_VEC Computes (T x 1) vector of (future)
              % maximum monthly payments for term contracts (with a given # of 
              % months) and tenure contracts. 
              %
              %  advance = GET_PAYMENT_MONTHLY(obj, 'tenure')
              %     returns maximum monthly payment vector for a tenure 
              %     payment scheme
              %
              %  advance = GET_PAYMENT_MONTHLY(obj, 'term', varargin)
              %     returns maximum monthly payment vector for a term payment
              %     scheme of a duration of month 'varargin'
              %
              %  for more details and explanation, see Szymanoski (1990), p.
              %     29
              
                    % check inputs
                    if ~any(strcmp(type, {'term','tenure'}))
                        error('For calculation of maximum monthly payment, payment scheme needs to be term or tenure!');
                    end
                    if strcmp(type,'term') && isempty(varargin)
                        error('For maximum monthly payment with term scheme, please enter number of periods!');
                    end                        
              
                    % inforce: periods during which contract is in force
                    inforce = (1:obj.contract.horizon_months)'; 
                    inforce = (inforce > (12 * obj.contract.age_for_calculation));   

                    % t: number of months that contract is in
                    % force (at start of month when payments are made)
                    t = (1:obj.contract.horizon_months)' - (12 * obj.contract.age_for_calculation + 1); 
                    t =  t .* inforce;
                    
                    % lastperiod: last period of monthly payments
                    if strcmpi(type, 'tenure') == 1
                        lastperiod = t(end);
                    elseif strcmpi(type, 'term') == 1
                        lastperiod = t(t==(varargin{:}-1));
                    end
                                        
                    % m: see Szymanoski (1990), p. 29
                    m = inforce .* max((lastperiod - t + 1),0);
                    
                    % compute monthly advance
                    payment_monthly_vec = inforce .* principal_limit_vec...
                        .* ((1+obj.contract.ir_compounding/12).^m)...
                        .* (obj.contract.ir_compounding/12)...
                        ./(((1+obj.contract.ir_compounding/12).^(m+1)) - (1+obj.contract.ir_compounding/12));

       end % function GET_PAYMENT_MONTHLY_VEC              
       function premiums_expected = get_premiums_expected(obj)
             % PREMIUMS_EXPECTED   Expected, discounted premiums given the 
             % casflow model
             %
             % Sum of expected, discounted premiums for insurer, given the
             % premiums of the cashflow account
           
                 % get them
                 premiums_expected = obj.contract.probs.s_cond_monthly ...
                                       .* obj.premiums;        

                 % discount them
                 premiums_expected_discounted = ...
                         obj.get_pv(premiums_expected,...
                         (obj.contract.model.ir_discount_basic/12),...
                         12 * obj.contract.age_for_calculation + 1);           

                 % sum them up
                 premiums_expected = nansum(premiums_expected_discounted);
             
       end % function GET_PREMIUMS_EXPECTED
       function loss_expected = get_loss_expected(obj)
             % LOSS_EXPECTED   Expected, discounted loss given the cashflow
             % model
             %
             % Sum of expected, discounted losses for insurer, given the
             % underlying end balance of the payment account; computation
             % is based on the conditional value of collateral, termination
             % probability and loss probability
             %
             % for the exact formula, see 'price_collateral' in 'contract'
             
                 % get expected, nominal loss
                 [loss_expected, a, b] = price_collateral(obj.b_e, obj.contract);             
             
                 % discount loss
                 losses_exp_discounted =...
                    obj.get_pv(loss_expected,...
                    (obj.contract.model.ir_discount_basic/12),...
                    12 * obj.contract.age_for_calculation + 1);
                 
                 % sum them up
                 loss_expected = nansum(losses_exp_discounted);             
             
       end % function GET_LOSS_EXPECTED          
   
   end % methods
   
   methods
      
       %% GETS
       function utilization = get.utilization(obj)
          % Compute dependent property of insurance utilization
          utilization = (obj.loss_expected./obj.premiums_expected).*100;
       end
           
   end % methods GETS
       
   methods (Static)
        
       %% STATIC methods
       
       % Compute present value
       function presentvalue = get_pv(nominalvalue,i,firstperiod)
             % GET_PV calculates present value cashflow vector for given 
             % nominal cashflow vector, given a periodic interest rate and
             % the start period t (in t: nv = pv)

                % 1) set interest rate (vector/scalar cases)
                if size(i,1) > 1
                    if size(i,1) == size(nominalvalue,1)
                        interestrate = i;
                    else
                        disp('ERROR: interest rate vector needs to have the')
                        disp('exact same length as the nominal value vector!')
                    end
                elseif size(i,1) == 1
                    interestrate = i;
                end
                
                % 2) set parameters of transformation                                  
                l = (1 : length(nominalvalue))';   
                l = (l >= firstperiod);               % l is logical: periods that get discounted
                
                t = (1 : length(nominalvalue))' - firstperiod;  % array t: # periods to be discounted
                t =  t.*l;

                % 3) create discount vector
                discountfactor = (1 ./ ((1 + interestrate).^t)) .* l;
                
                % 4) transform nominal into present value
                presentvalue = nominalvalue .* discountfactor;   
                        
         end % function GET_PV
         
    end % methods (Static)
    
end % classdef
