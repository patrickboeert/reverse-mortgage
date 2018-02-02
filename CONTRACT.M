classdef contract < handle
% CONTRACT  models a single reverse mortgage contract
%
%   CONTRACT models a reverse mortgage contract with the HECM pricing
%   procedure, cashflow modeling, plotting functions and the simulation 
%   mechanism for the stochastic simulation model
%
%   INPUTS:
%
%    name             description                         dim [class] 
%  ________________________________________________________________________
%    age              age at origination                  string or 2x1 [double]
%    gender           gender                              string or 2x1 [cellstring]
%    house            appraised value of collateral       1x1 [double]
%    expected_rate    expected rate at origination        1x1 [double]
%    scheme           payment scheme choice               3x1 [structure]
%      .name              payment scheme name             string
%      .loc               line of credit amount choice    1x1 double
%      .periods           payment period choice (term)    1x1 double
%    model            HECM model settings                 1x1 [model]
%    varargin         utilization % of principal limit    1x1 [double]
%
%
%    name             potential values 
%   _______________________________________________________________________             
%    gender           'female', 'male'   
%    scheme           
%      .name          'tenure', 'term', 'loc', 'modifiedtenure', 'modifiedterm'
%
%   EXAMPLE: 
%       
%    rm = contract(75, 'female', 100000, 0.1, scheme, hecm) 
%    returns a reverse mortgage contract for a 75 year old female with a
%    house valued at 100.000 in an interest rate environment of 10%. The
%    contract is computed for the payment scheme choice 'scheme', where
%    'scheme' is a structure with the fields 'name', 'loc' and 'periods' and
%    needs to be set before calling CONTRACT. CONTRACT returns a reverse
%    mortgage based on the model instance 'hecm' which has to be
%    initialized before calling CONTRACT.
%    
%   USAGE:
%
%    hecm = model();             % initialize mode configuration
%    scheme.name = 'term';       % set form of payment scheme
%    scheme.periods = 120;       % 10-year term payment
%    rm = contract(75, 'female', 100000, 0.1, scheme, hecm); % price contract

    properties        
             
        %% simulation: relations
        id         = [];           % simulation trial id            
      
        %% payments: accounts
        acc_simulation  = [];      % payment account: simulated HECM payment model        
        acc_payment     = [];      % payment account: actual HECM payment model
        acc_basic       = [];      % payment account: basic HECM payment model        
       
        %% payments: parameters
        payment_plan    = '';      % payment plan: 'tenure', 'term', 'loc', 'modifiedtenure', 'modifiedterm'
        payment_periods = [];      % payment periods (only 'term', 'modifiedterm')
        loc_amount      = [];      % loc_amt (only 'modifiedtenure', 'modifiedterm')
        util_pl         = 1;       % fractional utilization of principal limit, HECM payment account
               
        %% basic model: collateral
        house_initialvalue     = [];    % initial appraisal value of collateral         
        house_expected         = [];    % expected house price; E[H(t)]
        house_expected_cond    = [];    % expected house price, conditional on loss; E [H(t)|H(t)<B(t)]
        house_max_claim_amount = [];    % maximum claim amount
        house_simulated        = [];    % house price path simulated by GBM
        
        %% basic model: insurance
        premiumstructure     = [];      % time structure for premium fees
        
        %%  basic model: interest rates
        ir_expected_rate  = [];         % expected rate (fixed, yearly)
        ir_compounding    = [];         % compounding rate (fixed, yearly)
        ir_discount_basic = [];         % discount rate (fixed, yearly)
        
        %% basic model: pricing
        plf_exogenous             = [];  % principal limit factor, exogenous
        plf_endogenous            = [];  % principal limit factor, endogenous
        principal_limit_init      = [];  % initial principal limit
        principal_limit_init_net  = [];  % initial net principal limit (after setaside)
        
        %% basic model: iteration
        iter_npv      = [];             % net present value of insurance after iteration (plf_exogenous)
        iter_exitflag = [];             % exitflag of iteration (see MATLAB DOC 'fzero')
        iter_output   = [];             % output structure of iteration (see MATLAB DOC 'fzero')        
        
        %% basic model: fees & costs
        setaside        = [];           % set-aside of principal limit for servicing fees
        closing_cost    = 0;            % closing costs        
        
        %% properties: age
        borrower_age        = [];       % age of borrower
        coborrower_age      = [];       % age of coborrower
        age_for_calculation = [];       % age relevant for plf calculation
        
        %% properties: gender
        borrower_gender        = '';    % gender of borrower
        coborrower_gender      = '';    % gender of coborrower
        gender_for_calculation = '';    % gender relevant for plf calculation 
        
        %% properties: timing
        period_origination = [];        % period of origination
        period_terminated  = [];        % period of termination (due to death, moveout)
        period_assigned    = [];        % period of assignment to insurer   
        
        %% properties: status
        reason_terminated  = '';        % reason for termination: 'death', 'moveout'        
        
        %% objects: probability 
        probs = [];                     % probability data
                
        %% objects: model
        model               = [];       % model settings     
     
    end
    
    properties (Dependent = true)
        
        %% properties: parameters
        horizon_months = [];         % pricing horizon / months
        horizon_years  = [];         % pricing horizon / years 
        
    end
    
    methods
         
         %% CONSTRUCTOR
         function obj = contract(age, gender, house, expected_rate, scheme, model, varargin)             
          % CONTRACT  constructs a reverse mortgage contract
            
            % parse inputs
            if nargin >= 6
               
               %% set: age
               switch length(age)    % handle age for 1-dim, 2-dim input
                   case 1
                       obj.borrower_age = age;
                   case 2
                       obj.borrower_age   = age(1);
                       obj.coborrower_age = age(2);
                   otherwise
                       error('ReverseMortgage:Contract:Contract:AgeNotSpecifiedWell',...
                       ['Age variable at instantiation of contract must either be integer'...
                       '(one borrower) or 2x1,1x2 double (for borrower / coborrower!'])
                end
               
               %% set: gender
               if ischar(gender)     % handle gender for 1-dim, 2-dim input
                       obj.borrower_gender = gender;
               elseif iscell(gender) && max(size(gender)) == 2
                       obj.borrower_gender   = gender{1};
                       obj.coborrower_gender = gender{2};
               else
                       error('ReverseMortgage:Contract:Contract:GenderMisspecified',...
                       ['Gender at origination of contract must be either be string'...
                       '(female/male) or 2x1 / 1x2 cellstring for borrower/coborrower (female/male)!']);
               end
               
               %% set: house
               obj.house_initialvalue = house;
               
               %% set: expected & discount rate (Szymanoski (1990), p. 26)
               obj.ir_expected_rate  = expected_rate;                             
               
               %% set: model object
               if length(model) == 1
                   obj.model = model;
               else
                   error('ReverseMortgage:Contract:Contract:ModelInputTooBig',...
                   ['Model input at instantiation of contract must be 1-dim'...
                    'instance of class model !']);
               end
               
               %% set: utilization principal limit
               if ~isempty(varargin) && (varargin{1} >= 0 && varargin{1}<=1)
                   obj.util_pl = varargin{:};
               elseif ~isempty(varargin) && ~(varargin{1} >= 0 && varargin{1}<=1)
                   error('ReverseMortgage:Contract:Contract:WrongSpecificationOfUtilLevel',...
                         ['When creating a contract the utilization level'...
                          ' must be between 0 and 1']);
               end
               
               %% instantiate: probability object
               obj.probs = probability(obj);
               
               %% compute: premium structure
               obj.premiumstructure = get_premiumstructure(obj);
               
               %% compute: set aside
               obj.setaside = get_setaside(obj);
               
               %% compute: principal limit
                   % initial principal limit
                   obj.principal_limit_init =...
                       obj.plf_exogenous .* obj.house_max_claim_amount;

                   % initial net principal limit                   
                   obj.principal_limit_init_net =...
                       obj.plf_exogenous .* obj.house_max_claim_amount...
                             - obj.setaside(12 * obj.age_for_calculation + 1);               
               
               %% compute: payment accounts   
          
                   % set contract's payment scheme properties
                   
                       % payment scheme: name
                       if isfield(scheme,'name') && ~isempty(scheme.name)
                           obj.payment_plan = scheme.name;
                       else
                           error('ReverseMortgage:Contract:Contract:MustProvidePaymentSchemeName',...
                           ['For instantiation of a reverse mortgage contract,'...
                           'you must provide a payment scheme name for the constructor !']);
                       end

                       % payment scheme: periods
                       if ~isempty(findstr('term',scheme.name)) && (isfield(scheme,'periods') && ~isempty(scheme.periods));
                          obj.payment_periods = scheme.periods;
                       elseif ~isempty(findstr('tenure',scheme.name)) || ~isempty(findstr('loc',scheme.name))
                          obj.payment_periods = NaN;
                       else
                           error('ReverseMortgage:Contract:Contract:MustProvidePeriodForTermMortgage',...
                           ['For instantiation of a term reverse mortgage contract,'...
                           'you must provide a number of payment periods for the term mortgage !']);
                       end

                       % payment scheme: loc amount
                       if ~isempty(findstr('modified',scheme.name)) && (isfield(scheme,'loc') && ~isempty(scheme.loc));
                          obj.loc_amount = scheme.loc;
                       elseif any(strcmpi(scheme.name, {'tenure', 'term'}))
                           obj.loc_amount = NaN;                                                    
                       else
                           error('ReverseMortgage:Contract:Contract:MustProvideLOCAmount',...
                           ['For instantiation of a modified term or tenure reverse mortgage contract,'...
                            'you must provide an amount for the line of credit account !']);
                       end
                       
                  % check for sufficiency of principal limit for LOC
                       if any(strcmpi(scheme.name, {'modifiedtenure', 'modifiedterm'}))
                           if obj.principal_limit_init_net < obj.loc_amount
                               error('ReverseMortgage:Contract:Contract:InconsistentLOCAmount',...
                              ['For instantiation of a modified term or tenure reverse mortgage contract'...
                               'the principal limit must be bigger than the requested line of credit !']);
                           end
                       end

                  % open accounts: basic & payment (exogenous)
                  obj.acc_basic   = open_acc(obj,'basic');                  
                  obj.acc_payment = open_acc(obj,'payment');
                  

            elseif nargin < 6
                
                error('ReverseMortgage:Contract:Contract:TooFewInputs',...
                      'Too few input arguments to constructor of class contract!');
                
            end % nargin
         end % function CONSTRUCTOR
            
         %% PRICING           
         function [plf,obj] = price(obj,plf_mode)
             % PRICE calculates principal limit factor for reverse
             % mortgage contract through iteration
             %
             % principal limit factor is iterated until net present value of
             % expected losses equals net present value of expected premiums
             % based on the specifics of the model; in each plf iteration
             % 'price_payments' and 'price_collateral' are called
             %
             % OPTIONS:
             %
             %  PLF_MODE:
             %
             %      'exogenous': PRICE assumes cashflow pattern for the plf
             %      iteration is the maximal one-time upfront advance.
             %
             %      'endogenous': PRICE assumes a different cashflow pattern
             %      for the plf iteration. For the actual cashflow pattern,
             %      PRICE falls back on the properties of the payment 
             %      account that evaluates the endogenous plf.

                 % set: calling payment account (if mode endogenous)
                 if strcmpi(plf_mode,'endogenous') && ~isempty(evalin('caller','obj'))
                     acc_calling = evalin('caller','obj');
                 end
             
                 % set: iteration options
                 options  = optimset('Display',obj.model.iter_output, 'TolX',obj.model.iter_tolerance);            

                 % iterate: do net present value (npv) calculation in every iteration,
                 % using the probability model the and collateral model
                 [plf_iteration, npv, exitflag, output] = fzero(@(plf) get_npv(plf,obj,plf_mode),[0 1],options);

                 % store: iteration information locally
                 plf = plf_iteration;
                 obj.iter_npv      = npv;
                 obj.iter_exitflag = exitflag;
                 obj.iter_output   = output;  

             end % function PRICE
         function npv       = get_npv(plf, obj, plf_mode) 
             % GET_NPV calculates net present value of contract-specific
             % insurance
             % 
             % net present value (npv) of contract-specific insurance is
             % calculated for a random guess of a principal limit factor (plf); 
             % npv in this context means npv of expected discounted premiums 
             % - npv of expected discounted losses
             %
             % OPTIONS:
             %
             %  PLF_MODE:
             %
             %      'exogenous': GET_NPV assumes cashflow pattern for the plf
             %      iteration is the maximal one-time upfront advance.
             %
             %      'endogenous': GET_NPV assumes a different cashflow pattern
             %      for the plf iteration. For the actual cashflow pattern,
             %      PRICE falls back on the properties of the  payment 
             %      account that evaluates the endogenous plf
             
                 % set: calling payment account (if mode endogenous)
                 if strcmpi(plf_mode,'endogenous') && ~isempty(evalin('caller','acc_calling'))
                     acc_calling = evalin('caller','acc_calling');
                 end                

                 % 1) calculate cashflows based on plf guess 
                 cashflows = price_payments(plf, obj, plf_mode);
                 
                 % 2) calculate house price model based on plf guess and
                 %    resulting cashflows
                 [loss_expected, notneeded, notneeded2] =...
                     price_collateral(cashflows.b_e, obj);

                 % 2) premiums: get pv, exp_pv and sum_exp_pv of premiums
                 premiums_expected = ...
                     obj.probs.s_cond_monthly .* cashflows.premiums;             
                 
                 premiums_expected_discounted = ...
                     obj.get_pv(premiums_expected,(obj.model.ir_discount_basic/12),12*obj.age_for_calculation+1);           
                 
                 sum_premiums_expected_discounted =...
                     nansum(premiums_expected_discounted);
                 
                 % 3) losses: get exp_pv and sum_exp_pv of losses
                 losses_exp_discounted = ...
                     obj.get_pv(loss_expected,(obj.model.ir_discount_basic/12),12*obj.age_for_calculation+1);
                
                 sum_losses_exp_discounted =...
                     nansum(losses_exp_discounted);

                 % 4) return net present value of single insurance contract
                 npv = sum_premiums_expected_discounted - sum_losses_exp_discounted;

         end % function GET_NPV
         function cashflows = price_payments(plf, obj, plf_mode)
             % PRICE_PAYMENTS calculates payment vectors for a one-time upfront
             % payment for the purpose of calculating the contract's
             % principal limit factor
             % 
             % calculates nominal monthly cashflows (loan balance at start and
             % end of period, advances to borrower, premiums fees paid to insurer,
             % interest payments to the lender) of contract 
             %
             % OPTIONS:
             %
             %  PLF_MODE:
             %
             %      'exogenous': PRICE_PAYMENTS assumes cashflow pattern for the plf
             %      iteration is the maximal one-time upfront advance.
             %
             %      'endogenous': PRICE_PAYMENTS assumes a different cashflow pattern
             %      for the plf iteration. For the actual cashflow pattern,
             %      PRICE falls back on the properties of the payment account
             %      that evaluates the endogenous plf
             
                %% distinction between plf modes
                switch plf_mode
                    
                    case 'exogenous'
             
                        % create monthly payment vectors
                        cashflows.b_s           = zeros(obj.horizon_months,1);  
                        cashflows.premiums      = zeros(obj.horizon_months,1); 
                        cashflows.advances      = zeros(obj.horizon_months,1);
                        cashflows.interest      = zeros(obj.horizon_months,1);   
                        cashflows.b_e           = zeros(obj.horizon_months,1); 

                        % set starting values for arrays
                        cashflows.b_s(12 * obj.age_for_calculation + 1)=...
                            plf * obj.house_initialvalue;                  % first month's start balance is 0

                        cashflows.premiums(12 * obj.age_for_calculation + 1) = ...        % first month's upfront premium
                            obj.model.upfront * obj.house_initialvalue; 

                        cashflows.advances(12 * obj.age_for_calculation + 1) = ...        % first month's advances, upfront premium is paid implicitly in loan balance t=1
                            plf * obj.house_initialvalue - cashflows.premiums(12 * obj.age_for_calculation + 1);            

                        cashflows.b_e(12 * obj.age_for_calculation + 1) = ...     % first month's end balance
                            cashflows.advances(12 * obj.age_for_calculation + 1) + cashflows.premiums(12 * obj.age_for_calculation + 1);

                        % recursively calculate monthly nominal payments                
                        for period = (12 * obj.age_for_calculation + 2):obj.horizon_months
                          cashflows.b_s(period) = cashflows.b_e(period-1);
                          cashflows.interest(period)      = cashflows.b_s(period) .* (obj.ir_expected_rate/12);
                          cashflows.premiums(period)      = obj.premiumstructure(period) .* cashflows.b_s(period);
                          cashflows.b_e(period)   = cashflows.b_s(period) + cashflows.interest(period) + cashflows.premiums(period);
                        end;
                    
                    case 'endogenous'
                        
                        % 0) set: calling payment account (if in mode endogenous)
                        if strcmpi(plf_mode,'endogenous') && ~isempty(evalin('caller','acc_calling'))
                            acc_calling = evalin('caller','acc_calling');
                        end                          
                        
                        % 1) create monthly payments from principal limit
                        % factor guess

                            % get principal limit, principal limit vector
                            principal_limit_init_net =...
                                plf .* obj.house_initialvalue...
                                - obj.setaside(12 * obj.age_for_calculation + 1);
                            
                            principal_limit_vec =...
                                get_principal_limit(acc_calling(1), principal_limit_init_net);
                            
                            % get maximum monthly payment
                            % (with plf guess)
                            switch acc_calling(1).type
                                case 'tenure'
                                    payment_monthly =...
                                        get_payment_monthly(acc_calling(1), principal_limit_vec, 'tenure');
                                case 'term'
                                    payment_monthly =...
                                        get_payment_monthly(acc_calling(1), principal_limit_vec, 'term', obj.payment_periods);
                            end 
                        
                        % create monthly payment vectors
                        cashflows.b_s           = zeros(obj.horizon_months,1);  
                        cashflows.premiums      = zeros(obj.horizon_months,1); 
                        cashflows.advances      = zeros(obj.horizon_months,1);
                        cashflows.interest      = zeros(obj.horizon_months,1);   
                        cashflows.b_e           = zeros(obj.horizon_months,1);

                        % set starting values: advances, premiums, end
                        % balance

                                % advances
                                switch acc_calling(1).type 
                                    case 'tenure'
                                      cashflows.advances((12*obj.age_for_calculation+1):end) =...
                                          payment_monthly;  

                                    case 'term'
                                      cashflows.advances((12*obj.age_for_calculation+1):...
                                          (12*obj.age_for_calculation+obj.payment_periods)) =...
                                          payment_monthly;                      
                                end

                                % premiums
                                cashflows.premiums(12 * obj.age_for_calculation + 1) =...
                                    obj.model.upfront * obj.house_initialvalue; 

                                % end balance
                                cashflows.b_e(12 * obj.age_for_calculation + 1) = ...
                                    cashflows.advances(12 * obj.age_for_calculation + 1)...
                                    + cashflows.premiums(12 * obj.age_for_calculation + 1);       

                         % recursive payment computation                    
                         for m = (12*obj.age_for_calculation+2):obj.horizon_months
                           cashflows.b_s(m)           = cashflows.b_e(m-1);
                           cashflows.interest(m)      = cashflows.b_s(m) .* (obj.ir_expected_rate/12);             
                           cashflows.premiums(m)      = obj.premiumstructure(m) .* cashflows.b_s(m);
                           cashflows.b_e(m)           = cashflows.b_s(m) + cashflows.advances(m)...
                                                     + cashflows.interest(m) + cashflows.premiums(m);
                         end                                               
                end

         end % function PRICE_PAYMENTS 
         function [loss_expected, lossprobs, house_expected_cond] = price_collateral(balance_end, obj)
             % PRICE_COLLATERAL projects value of collateral, loss
             % probabilities and expected nominal loss given payment vectors
             % 
             % house price value over time is modelled as a geometric
             % Brownian motion process with constant drift and difussion; the
             % expected house value can hence be calculated analytically (note:
             % error in Szymanoski, 1994 on p. 354); returns the expected losses 
             % in present value terms for each year of the contract period.
             %
             % For a more extensive review of the HECM modeling technique,
             % see 'Risk and the Home Equity Conversion Mortgage', E.J. Szymanoski, Jr,
             % Real Estate Economics, 1994, Vol. 22, Iss. 2, pp. 347 - 366
             %
             % CONTENT  Code Overview
             %    1) path of expected (unconditional) house value
             %    2) loss probability for all t
             %    3) conditional expected collateral value
             %    4) nominal expected loss

             % (1) calculate path of expected collateral value 'house_expected'            
                 house_expected = get_house_expected(obj);

             % (2) loss probability: P(collateral(t) <= loanbalance(t))
             % loss probability lossprob is calculated analytically, where
             % H(t)/H(0) is distributed lognormal; p.14, equation (5)       
             
                 % l is logical: periods during which contract is in force
                 l = (1:obj.horizon_months)'; 
                 l = (l > (12 * obj.age_for_calculation));   

                 % t is array: number of months that contract is in force 
                 t = (1:obj.horizon_months)' - (12 * obj.age_for_calculation + 1); 
                 t =  t .* l;
             
                 %  upper integration bound
                 b = l .* balance_end ./ obj.house_initialvalue;   % outstanding balance in terms of house_initialvalue value
                 lnb = l .* log(b);                              % H(t)/H(0) is distributed lognormal

                 % transform mu, sigma into monthly values
                 mu    = obj.model.mu/12;
                 sigma = obj.model.sigma/sqrt(12);
                 
                 % (mu,sigma) of cumulated appreciation distribution (see 'help normcdf')
                 mu_cumulated_apprec = l .* mu .* t;
                 sigma_cumulated_apprec = l .* sigma .* sqrt(t);

                 % A(t) := P(collateral value in t < outstanding balance in t) using 
                 % cumulative normal distribution function
                 lossprobs = l .* normcdf(lnb,mu_cumulated_apprec,sigma_cumulated_apprec);
                 lossprobs(12 * obj.age_for_calculation + 1) = 0;

              % (3) conditional expected collateral value 'house_expected_cond'
              %     see p.15, equation (7)

                 %  upper integration bound U(t)
                 U = (lnb - mu_cumulated_apprec) ./ sigma_cumulated_apprec;

                 % (mu,sigma) of conditional expected value (see 'help normcdf')
                 mu_cond_expected_value    = sigma_cumulated_apprec;
                 sigma_cond_expected_value = l .* 1;                 

                 % calculate beta
                 beta = l .* (1./lossprobs) .* normcdf(U,mu_cond_expected_value,sigma_cond_expected_value); 

                 % conditional expected collateral value
                 house_expected_cond = l .* beta .* house_expected;

              %  (4) nominal expected loss
                 loss_expected = l .* obj.probs.t_cond_monthly .* lossprobs .* (balance_end - house_expected_cond);          

         end % function PRICE_COLLATERAL             
          
         %% CASHFLOW ANALYSIS
         function compare(obj)
             % COMPARE   Compare maximum monthly payment with two different
             % plf factors (exogen and endogen)
             %
             % COMPARE computes the maximum monthly payment for the
             % contract, given that the principal limit is used without a
             % line of credit. The maximum monthly payment is computed for
             % two different principal limit factors,
             % 
             %     1) the exogenous plf of the HECM model, where
             %        the advance is assumed to be a maximum 
             %        one-time upfront draw
             %
             %     2) the endogenous plf of the HECM Model, where the
             %        cashflow is modeled after the actual draws requested
             %        by the borrower (term or tenure)
             %
             % COMPARE displays the expected premiums and expected losses
             % associated with both principal limit factors.
             
                 %% check inputs
                    if strcmpi(obj.payment_plan,'loc')
                         error('ReverseMortgage:Contract:Compare:NoTermOrTenureAccount',...
                               ['Comparing payments for HECM contracts by exogenous / endogenous'...
                                'plf only makes sense for term and tenure contracts. This contract'...
                                'is a LOC account !']);
                    end
             
                 %% create payment accounts                          
                     % payment account with exogenous plf
                       exogenous = account(obj);                    % open
                       exogenous = payment(exogenous, 'exogenous'); % compute

                     % payment account with endogenous plf
                       endogenous = account(obj);                       % open
                       endogenous = payment(endogenous, 'endogenous');  % compute

                  %% display results
                     format bank

                     disp('      ');
                     disp('                 Monthly Payment    Utilization   E(premiums)    E(losses)'); 
                     disp('      ');
                     disp(strcat('HECM Payment:',['    ',num2str(exogenous(1).payment_monthly),'            ',num2str(round(exogenous(1).utilization)),'          ',num2str(exogenous(1).premiums_expected),'     ',num2str(exogenous(1).loss_expected)]));
                     disp(strcat('HECM Basic:',['      ',num2str(endogenous(1).payment_monthly),'            ',num2str(round(endogenous(1).utilization)),'          ',num2str(endogenous(1).premiums_expected),'     ',num2str(endogenous(1).loss_expected)]));                

         end % function COMPARE         
         function acc = open_acc(obj,kind,varargin)
            % OPEN_ACC  compute payment account for contract
            %
            % OPEN_ACC computes payment account for a given reverse mortgage
            % contract, either of kind 'basic' (HECM bascic model) or of 
            % kind 'payment' (HECM payment model).
            %
            % If account is of type 'payment' and called without another 
            % option, exogenous plf is used. 'exogenous' and 'endogenous'
            % can be passed as an option to open a HECM payment model
            % account with the respective plf factor.
                  
                  %% check inputs
                  if all(~strcmpi(kind,{'basic' 'payment' 'simulation'}))
                      error('ReverseMortgage:Contract:Open_Acc:NoSuchAccountExists',...
                            ['When opening an account, account needs to be either a'...
                             'HECM basic model, a HECM payment model or HECM simulation account!']);

                  elseif ~isempty(varargin) && (strcmpi(kind,'basic') || ~any(strcmpi(varargin{:},{'endogenous' 'exogenous'})))
                      error('ReverseMortgage:Contract:Open_Acc:WrongPaymentAccTypeSpecification',...
                            ['When opening HECM payment account, option needs to be'...
                             'one of endogenous / exogenous !']);                                                    
                  end
                  
                  %% set default option when no explicit call
                  if isempty(varargin) && strcmpi(kind,'payment')
                      varargin{1,1} = 'exogenous';
                  end
                  
                  %% open accounts
                  switch kind
                      case 'basic'          
                         
                         % open basic account
                         acc = account(obj); % open
                         acc = basic(acc);   % compute

                      case 'payment'
                          
                         % open payment account
                         acc = account(obj);  % open
                         acc = payment(acc, varargin{:});  % compute 
                         
                      case 'simulation'

                         % open simulation account
                         acc = account(obj); % open
                         acc = simulation(acc);   % compute                                            
                  end

         end % function OPEN_ACC          
         
         %% STOCHASTIC SIMULATION
         function obj = lifetime(obj, simulation, run)
           % LIFETIME  simulate stochastic reverse mortgage history
           %
           % LIFETIME generates lifetime events for reverse mortgages (e.g.
           % termination events, stochastic house price path) and handles
           % reverse mortgage assignment to insurer accordingly.
           
           %% error checking
            if isempty(simulation)
              error('ReverseMortgage:Contract:Lifetime:NoSimulationData',...
                    ['For a simulation of the contract history, a simulation instance'...
                     ' needs to be set! ']);
            end
            
           %% determine contract's termination time & type           
            % set contract run id
            obj.id = run;
           
            % get monthly probability data
            start       = obj.age_for_calculation*12 + 1;
            final       = length(obj.probs.s_cond_monthly);            
            PROBdeath   = obj.probs.t_death_cond_monthly(start:final);
            PROBmoveout = obj.probs.t_moveout_cond_monthly(start:final);
                        
            % determine time of event(s)
            drawnewsample = true;            
            while drawnewsample    
              
              % draw pseudo numbers              
              PSEUDOdeath   = rand(size(PROBdeath));
              PSEUDOmoveout = rand(size(PROBmoveout));              

              % determine realization of event
              EVENTdeath    = any(PSEUDOdeath <= PROBdeath);
              EVENTmoveout  = any(PSEUDOmoveout <= PROBmoveout);                                          
              
              % find respective month of event & reason
              if EVENTdeath || EVENTmoveout                
                drawnewsample = false;                
                if EVENTdeath && ~EVENTmoveout                  
                   month  = find((PSEUDOdeath <= PROBdeath), 1, 'first');
                   reason = 'death';                   
                elseif ~EVENTdeath && EVENTmoveout
                   month  = find((PSEUDOmoveout <= PROBmoveout), 1, 'first');
                   reason = 'moveout';
                elseif EVENTdeath && EVENTmoveout
                   month_death   = find((PSEUDOdeath <= PROBdeath), 1, 'first');
                   month_moveout = find((PSEUDOmoveout <= PROBmoveout), 1, 'first');
                   month = min(month_death, month_moveout);
                   if month_death == month
                     reason = 'death';
                   elseif month_moveout == month
                     reason = 'moveout';
                   end                   
                end                   
              end
            end % while
           
            %% set object's termination parameters 
            obj.period_terminated = 12*obj.age_for_calculation + month;
            obj.reason_terminated = reason;                        
         
         end % end LIFETIME  
         function obj = sim_series(obj, simulation, run)
           % SIM_SERIES compute simulated time series given termination
           % events and exogenous interest rate path
           
           acc = account(obj);                     % open account
           acc = acc.simulation(simulation, run);  % simulate account
           obj.acc_simulation = acc;
           
         end % end SIM_SERIES
         function obj = assignment(obj)
           % ASSIGNMENT   determines time of assignment event
              
            % create logical vector for life of contract         
            month = (1:obj.horizon_months)';            
            switch obj.acc_simulation(1).type
              case 'tenure'              
               notterminated = (month > obj.period_origination)...
                             & (month <= obj.period_terminated);                            
              case 'term'
               notterminated = (month > obj.period_origination)...
                            & (month <= obj.period_terminated)...
                            & (month <= obj.period_origination...
                                        + obj.payment_periods);
            end
            
            % get ob and mca across lifetime
            ob                 = obj.acc_simulation.b_e;
            mca                = NaN(size(obj.acc_simulation.b_e));
            mca(notterminated) = obj.model.threshold .* obj.house_max_claim_amount;
              
            % determine assignment period
            period_assigned    = (ob > mca);                            
            
            % this is Claim 2 Type / Claim 1 Type or Payoff not covered
            % here
            if any(period_assigned)
              obj.period_assigned = find(period_assigned,1,'first');            
            end              
            
         end % end ASSIGNMENT
         function obj = get_house_simulated(obj, simulation, run)
           % GET_HOUSE_SIMULATED  simulates house price path with Geometric
           % Brownian Motion process
           %
           % GET_HOUSE_SIMULATED simulates a single house price path for 
           % the contract based on an approximate solution to a geometric 
           % Brownian Motion process. The drift of the process is a time-
           % variant function based on the realized returns of the national
           % house price index simulated through the VAR(p) multivariate
           % time series model.
           %
           %    PARAMETERS:  nSteps            
           
            %% error checking
            if isempty(simulation)
              error('ReverseMortgage:Contract:GetHouseSimulated:NoSimulationData',...
                    ['For a simulation of the house price path, a simulation instance'...
                     ' needs to be set! ']);
            end
           
            %% get realized returns mu
            yields = simulation.mu_sim(run,:)/100;
            
            %% construct GBM process 
            r        = ts2func(yields,...                     % r is continous yield function
                               'Times', (1:length(yields))'); 
                             
            process  = gbm(r,...                                   % mu
                           simulation.model.sigma,...              % sigma
                           'StartTime', 1/12,...                   % start time of observation
                           'StartState', obj.house_initialvalue);  % start value
  
            %% simulate GBM process
            nPeriods = (12*length(yields))-1;                 % # of simulated observations at dt
            nSteps   = 10;                                    % # of intermediate time steps within each time dt increment
            dt       = 1/12;                                  % observation for each month            
                        
            obj.house_simulated = process.simulate(nPeriods,...
                                               'DeltaTime', dt,...
                                               'nSteps', nSteps);            
           
         end % end GET_HOUSE_SIMULATED
         function sub_sim = s(obj, sub_life)
           % S  transform lifetime index into simulation index           
           sub_sim = sub_life - 12 * obj.age_for_calculation;
         end % end S       
         
         %% PLOTTING
         function plot(obj, varargin)
             % PLOT  Plots contract's variables as time series
             %
             % If varargin is not used, PLOT will plot the contract's
             % time series advances, end balance, principal limit & expected
             % house value by default
             %
             % If the respective variable names are supplied as strings
             % in varargin, PLOT will plot at request: 
             %         
             %          advances ('advances')
             %          interest ('interest')
             %          premiums ('premiums')
             %          end balance ('b_e')
             %          principal limit ('principal_limit')
             %          initial house value ('house_initialvalue')
             %          expected house value ('house_expected')
             %          conditional house value ('house_expected_cond')
             %
             % EXAMPLES:
             %
             %    plot(contract) plots time series advances, end balance,
             %                   principal limit & expected house value
             %                   of contract
             %
             %    plot(contract, 'advances', 'premiums', 'interest',...
             %         'house_expected', 'b_e'} 
             %
             %                   plots advances, premiums, interest,
             %                   expected house value and end balance
             %                   of contract
             
                 %% define default time lines
                    % if no specific time series have been requested, plot
                    % defaults: advances, b_e, principal_limit,
                    % house_expected                    
                    if nargin == 1
                        varargin{1} = 'advances';
                        varargin{2} = 'b_e';
                        varargin{3} = 'principal_limit';
                        varargin{4} = 'house_expected';
                    end              
                    
                    % number of lines to plot               
                    nlines  = length(varargin);
                    
                 %% detect inputs             
                    var_name = cell(8,1);
                        var_name{1} = 'advances'; 
                        var_name{2} = 'interest';
                        var_name{3} = 'premiums';
                        var_name{4} = 'b_e';
                        var_name{5} = 'principal_limit';
                        var_name{6} = 'house_initialvalue';
                        var_name{7} = 'house_expected';
                        var_name{8} = 'house_expected_cond';
                        
                    for i = 1:nlines
                        var_ind(i) = strmatch(varargin{i}, var_name, 'exact');
                    end

                 %% legend description for plot
                    legend_name = cell(8,1);
                        legend_name{1} = 'Cumulated Advances';
                        legend_name{2} = 'Cumulated Interest';
                        legend_name{3} = 'Cumulated Premiums';
                        legend_name{4} = 'Outstanding Loan Balance';
                        legend_name{5} = 'Principal Limit';
                        legend_name{6} = 'Inital House Value';
                        legend_name{7} = 'Expected Value of Collateral';
                        legend_name{8} = 'Conditional Expected Value of Collateral (given a loss for the insurer)';                        

                 %% color for lines of plot
                    var_color = cell(8,1);
                        var_color{1} = [0.5 0.5 0.5];
                        var_color{2} = 'green'; 
                        var_color{3} = 'blue';
                        var_color{4} = 'cyan';
                        var_color{5} = 'magenta';
                        var_color{6} = 'yellow';
                        var_color{7} = 'black';
                        var_color{8} = [0.5 0.5 0.5];                   
             
                 %% arrange data
                     % index
                     start   = 12 * obj.age_for_calculation + 1;
                     year    = (1:(1/12):((obj.horizon_months-start+12)/12))';

                     % prepare data
                     if any(strcmpi(varargin,'advances')),
                        advances = cumsum(obj.acc_payment(1).advances(start:end)); end
                     if any(strcmpi(varargin,'interest'))
                        interest = cumsum(obj.acc_payment(1).interest(start:end)); end
                     if any(strcmpi(varargin,'premiums'))
                        premiums = cumsum(obj.acc_payment(1).premiums(start:end)); end
                     if any(strcmpi(varargin,'b_e'))
                        b_e = obj.acc_payment(1).b_e(start:end); end
                     if any(strcmpi(varargin,'principal_limit'))
                        principal_limit = obj.acc_payment(1).principal_limit(start:end); end
                     if any(strcmpi(varargin,'house_initialvalue'))
                        house_initialvalue = ones(obj.horizon_months-start+1,1) .* obj.house_initialvalue; end
                     if any(strcmpi(varargin,'house_expected'))
                        house_expected = obj.house_expected(start:end); end
                     if any(strcmpi(varargin,'house_expected_cond'))
                        house_expected_cond = obj.house_expected_cond(start:end); end
                    
                 %% plot                  
                    % create handle
                    h = zeros(nlines,1);
                    
                    % make plot
                    for i = 1:nlines
                       h(i) = plot(year,eval(var_name{var_ind(i)}));              % add line
                            set(h(i), 'LineWidth', 2);                            % set line width
                            set(h(i), 'Color', var_color{var_ind(i)});            % set line color                            
                            hold on                                               % wait to add more
                    end
                    hold off
                    
                    % title & label
                    legend(h, legend_name{var_ind},'Location','NorthWest'); % add legend to line
                    title('Stylized Cashflow for a Reverse Mortgage', 'FontWeight', 'bold');
                    xlabel('Years after Origination')  
         
         end % function PLOT
         function plotpaymentsensitivity(obj,deviation)
            % PLOTPAYMENTSENSITITY  visualize sensitivity of max. monthly 
            % payment to the assumption of the model's exogenous parameter
            %
            % PLOTPAYMENTSENSITITY plots the max. monthly payment for a tenure
            % contract for different parameters as a fraction of the 
            % baseline model's max. monthly payment. The baseline model is 
            % defined as the model underlying the current contract. The
            % domain of alternative parameters is defined as a percentage
            % deviation from the baseline model's parameter.

              % get data
              moveout = plotpaymentby(obj,'moveout',deviation);
              fee     = plotpaymentby(obj,'fee',deviation);
              mu      = plotpaymentby(obj,'mu',deviation);
              sigma   = plotpaymentby(obj,'sigma',deviation);
                                        
              % plot result
              h = plot(moveout.parameter,...
                      [moveout.tenurepayment, fee.tenurepayment,...
                       mu.tenurepayment, sigma.tenurepayment]);
                   set(h, 'LineWidth', 2);                             
 
              legend(h, {'by moveout rate',...
                         'by yearly insurance fee',...
                         'by drift parameter mu',...
                         'by diffusion parameter sigma'},...
                         'Location','NorthWest');                    
              title('sensitivity of max. monthly payment',...
                     'FontWeight', 'bold');
              ylabel('max. monthly payment as percentage of baseline model');
              xlabel('parameter percentage deviation from baseline model assumption');  
              
         end % function PLOTPAYMENTSENSITITY
         function plotinsurancesensitivity(obj,deviation)
            % PLOTINSURANCESENSITITY  visualize sensitivity of utilization
            % of net present value to deviation of parameters from the
            % model's assumption
            %
            % PLOTINSURANCESENSITITY plots the utilization level of the net
            % present value by different parameters, where parameters are 
            % used to compute expected premiums and expected losses using 
            % the basic model. In other words an answer to the question: How 
            % sensitive is the insurance to a wrong assumption AND which
            % parameter is more sensitive than others?            

              % get data
              discount = plotinsuranceby(obj,'discount rate',deviation);
              moveout  = plotinsuranceby(obj,'moveout',deviation);
              mu       = plotinsuranceby(obj,'mu',deviation);
              sigma    = plotinsuranceby(obj,'sigma',deviation);
                                        
              % plot result
              h = plot(moveout.parameter,...
                      [discount.eutilization, moveout.eutilization,...
                       mu.eutilization, sigma.eutilization]);
              set(h, 'LineWidth', 2);                             
 
              legend(h, {'by discount rate',...
                         'by moveout rate',...
                         'by drift parameter mu',...
                         'by diffusion parameter sigma'},...
                         'Location','NorthWest');                    
              title('sensitivity of utilization level of insurance net present value',...
                     'FontWeight', 'bold');
              ylabel('utilization level');
              xlabel('parameter percentage deviation from model assumption');  
              
         end % function PLOTINSURANCESENSITITY         
         function varargout = plotpaymentby(obj, factor, varargin)
            % PLOTPAYMENTBY  plot maximum monthly payment by exogenous
            % parameter
            %
            % PAYMENTBY performs comparative statics on the maximum monthly
            % payment for a tenure and 10-year term payment scheme that can
            % be supported through the insurance scheme. It plots the 
            % maximum monthly payment against one of the model's
            % parameters. Parameters can be age at origination, expected
            % rate at origination, the move-out factor, the monthly
            % insurance fee fraction, mu or sigma.
            %
            %    USAGE:
            %
            %       output = plotpaymentby(contract, factor, [option])
            %
            %    FACTOR is one string of:
            %
            %       age at origination           | 'age'
            %       expected interest rate       | 'interest rate'
            %       move-out factor              | 'moveout'
            %       insurance fee                | 'fee'
            %       drift of house price GBM     | 'mu'
            %       diffusion of house price GMB | 'sigma'
            %
            %   OPTIONS:
            %
            %       Option is a double input (e.g. 0.1 for 10%). It is the maximum
            %       % derivation of the parameter from the parameter
            %       of the underlying model.

              %% check inputs
              callforelasticities     = ~isempty(varargin);
              noelasticitiesavailable = (strcmp(factor,'age') || strcmp(factor,'interest rate'));
              
              if (callforelasticities && noelasticitiesavailable)
                warning('ReverseMortgage:Contract:PlotPaymentBy',...
                        'No elasticities can be computed for %s. Will switch to regular output without options.',...
                        factor);
              end                           
            
              %% get parameters
              switch factor
                  case 'age'
                      parameter = (obj.model.age_min:obj.model.age_max)';
                      label     = 'age';                                                  
                  case 'interest rate'
                      parameter  = obj.model.ir_expected_rate(30:end);
                      label      = 'expected rate at origination';                
                  case 'moveout'
                      if callforelasticities
                          base  = obj.model.moveout;
                          step  = base * 0.01; % 1% steps
                          start = base - (varargin{1}*100)*step;
                          final = base + (varargin{1}*100)*step;                          
                          parameter = (start:step:final)'; 
                          label      = 'moveout factor of probability model';                          
                      else
                          parameter  = (0:0.05:2)';
                          label      = 'moveout factor of probability model';
                      end                      
                  case 'fee'
                      if callforelasticities                          
                          base  = obj.model.monthly;
                          step  = base * 0.01; % 1% steps
                          start = base - (varargin{1}*100)*step;
                          final = base + (varargin{1}*100)*step;                          
                          parameter = (start:step:final)'; 
                          label     = 'yearly premium fee as % of the outstanding balance';                          
                      else                      
                          parameter = (0.001:0.001:0.02)';
                          label     = 'yearly premium fee as % of the outstanding balance';
                      end
                  case 'mu'
                      if callforelasticities
                          base  = obj.model.mu;
                          step  = base * 0.01; % 1% steps
                          start = base - (varargin{1}*100)*step;
                          final = base + (varargin{1}*100)*step;                          
                          parameter = (start:step:final)'; 
                          label     = 'mu: yearly drift parameter of house price model';
                      else                        
                          parameter = (-0.02:0.005:0.06)';
                          label     = 'mu: yearly drift parameter of house price model';
                      end
                  case 'sigma'
                      if callforelasticities
                          base  = obj.model.sigma;
                          step  = base * 0.01; % 1% steps
                          start = base - (varargin{1}*100)*step;
                          final = base + (varargin{1}*100)*step;                          
                          parameter = (start:step:final)'; 
                          label     = 'sigma: yearly diffusion parameter of house price model';
                      else                       
                          parameter = (0.01:0.01:0.2)';
                          label     = 'sigma: yearly diffusion parameter of house price model';
                      end
                  otherwise
                      error('ReverseMortgage:Contract:PlotPaymentBy:WrongInput',...
                            'The input argument %s is not a recognized parameter of the model. Use one of: age, interest rate, move out, fee, mu, sigma.',...
                            factor);                      
              end
              
              %% preallocate output
              termpayment   = zeros(length(parameter),1);
              tenurepayment = zeros(length(parameter),1);
              
              %% payment schemes
              termscheme.name    = 'term';
              termscheme.periods = 120;
              tenurescheme.name  = 'tenure';      
              
              %% get data
              switch factor
                  case 'age'
                    for i = 1:length(parameter)
                      % get contract
                      term = contract(parameter(i),...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    termscheme, obj.model);
                                
                      tenure = contract(parameter(i),...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    tenurescheme, obj.model);
                                
                      % get monthly payment
                      termpayment(i) = term.acc_payment.payment_monthly;
                      tenurepayment(i) = tenure.acc_payment.payment_monthly;
                    end  
                  case 'interest rate'
                    for i = 1:length(parameter)
                      % get contract
                      term = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    parameter(i),...
                                    termscheme, obj.model);
                                
                      tenure = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    parameter(i),...
                                    tenurescheme, obj.model);
                                
                      % get monthly payment
                      termpayment(i) = term.acc_payment.payment_monthly;
                      tenurepayment(i) = tenure.acc_payment.payment_monthly;
                    end                   
                  case 'moveout'
                    for i = 1:length(parameter)
                        
                      % set model
                      model         = obj.model;
                      model.moveout = parameter(i);
                      
                      switch obj.gender_for_calculation
                          case 'male', gender = 1;
                          case 'female', gender = 2;
                      end
                      
                      model.plf(obj.age_for_calculation, gender,...
                                (model.ir_expected_rate==obj.ir_expected_rate))...
                                    = NaN;
                      
                      % make contract
                      term = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    termscheme, model);
                                
                      tenure = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    tenurescheme, model);
                                
                      % get monthly payment
                      termpayment(i) = term.acc_payment.payment_monthly;
                      tenurepayment(i) = tenure.acc_payment.payment_monthly;
                    end 
                  case 'fee'
                    for i = 1:length(parameter)                        
                      % set model
                      model         = obj.model;
                      model.monthly = parameter(i)/12;
                      
                      switch obj.gender_for_calculation
                          case 'male', gender = 1;
                          case 'female', gender = 2;
                      end
                      
                      model.plf(obj.age_for_calculation, gender,...
                                (model.ir_expected_rate==obj.ir_expected_rate))...
                                    = NaN;
                      
                      % make contract
                      term = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    termscheme, model);
                                
                      tenure = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    tenurescheme, model);
                                
                      % get monthly payment
                      termpayment(i) = term.acc_payment.payment_monthly;
                      tenurepayment(i) = tenure.acc_payment.payment_monthly;
                    end                    
                  case 'mu'
                    for i = 1:length(parameter)                        
                      % set model
                      model         = obj.model;
                      model.mu      = parameter(i);
                      
                      switch obj.gender_for_calculation
                          case 'male', gender = 1;
                          case 'female', gender = 2;
                      end
                      
                      model.plf(obj.age_for_calculation, gender,...
                                (model.ir_expected_rate==obj.ir_expected_rate))...
                                    = NaN;
                      
                      % make contract
                      term = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    termscheme, model);
                                
                      tenure = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    tenurescheme, model);
                                
                      % get monthly payment
                      termpayment(i) = term.acc_payment.payment_monthly;
                      tenurepayment(i) = tenure.acc_payment.payment_monthly;
                    end 
                  case 'sigma'
                    for i = 1:length(parameter)                        
                      % set model
                      model         = obj.model;
                      model.sigma   = parameter(i);
                      
                      switch obj.gender_for_calculation
                          case 'male', gender = 1;
                          case 'female', gender = 2;
                      end
                      
                      model.plf(obj.age_for_calculation, gender,...
                                (model.ir_expected_rate==obj.ir_expected_rate))...
                                    = NaN;
                      
                      % make contract
                      term = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    termscheme, model);
                                
                      tenure = contract(obj.age_for_calculation,...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    tenurescheme, model);
                                
                      % get monthly payment
                      termpayment(i) = term.acc_payment.payment_monthly;
                      tenurepayment(i) = tenure.acc_payment.payment_monthly;
                    end                       
              end
                    
              %% adjust output according to option
              if callforelasticities
                  % x values
                  xvalues = ((-varargin{1}*100):1:(varargin{1}*100))';                  
                  
                  % rel. percentage deviation
                  termpayment = (termpayment/termpayment(parameter==base))*100;
                  tenurepayment = (tenurepayment/tenurepayment(parameter==base))*100;
                  
              elseif ~callforelasticities
                  % x values
                  xvalues = parameter;
              end
              
              %% plot result
              h(1) = plot(xvalues, termpayment);
                   set(h(1), 'LineWidth', 2);
                   set(h(1), 'Color', 'magenta');                               
                   hold on
              h(2) = plot(xvalues, tenurepayment); 
                   set(h(2), 'LineWidth', 2);
                   set(h(2), 'Color', 'black');                                                                                
                   hold off    
              legend(h, {'maximum monthly term payment',...
                         'maximum monthly tenure payment'},...
                         'Location','NorthWest');                    
              title('comparative statics of max. monthly payment',...
                     'FontWeight', 'bold');
              ylabel('payment');
              xlabel(label);
                                      
              %% output
              varargout{1}.parameter     = xvalues;
              varargout{1}.termpayment   = termpayment; 
              varargout{1}.tenurepayment = tenurepayment;                    
    
         end % function PLOTPAYMENTBY
         function varargout = plotinsuranceby(obj, factor, varargin)
            % PLOTINSURANCEBY  plots expected utilization rate u by the
            % value of an exogenous parameter that deviates from the true
            % parameter assumption in a baseline model
            %
            % The assumed cashflow is the cashflow derived from the principal
            % limit method, the assumed model is the 'model' instance
            % in the input instance 'contract'. 
            %
            % Expected utilization rate u can also be plotted for a 
            % percentage deviation of the parameter from
            % the true parameter (elasticity) with option 'varargin',
            % where 'varargin' is maximum percentage deviation from
            % baseline model's parameter value (e.g. 0.1 = 10% range). 
            %
            % Note: For parameter 'age' the plot shows the expected utilization
            % rate u for different ages at origination. This is not the
            % same logic as for the other parameters.
            %
            %   input 'factor' is one of:
            %
            %       'discount rate' | discount interest rate
            %       'moveout'       | move-out factor
            %       'mu'            | drift of house price GBM
            %       'sigma'         | diffusion of house price GMB
            %       'age'           | age at origination            
            %
            %   USAGE:
            %
            %    data = plotinsuranceby(contract, 'moveout', 0.1)
            %    data = plotinsuranceby(contract, 'moveout')            


              %% called in elasticity mode?
              callforelasticities = (~isempty(varargin) &&...
                                      isnumeric(varargin{1}));
            
              %% check inputs
              noparameter = (~strcmp(factor,'age') &&...
                             ~strcmp(factor,'discount rate') &&...
                             ~strcmp(factor,'moveout') &&...
                             ~strcmp(factor,'mu') &&...
                             ~strcmp(factor,'sigma'));              
              nonumber    = (~isempty(varargin) && ~isnumeric(varargin{1}));              
              nodeviation = (strcmp(factor,'age') && ~isempty(varargin));              
                            
              if noparameter
                error('ReverseMortgage:Contract:PlotInsuranceBy:WrongInput',...
                      'Input %s is not a recognized parameter. Use one of: age, discount rate, move out, fee, mu, sigma.',...
                      factor);
              elseif nonumber
                error('ReverseMortgage:Contract:PlotInsuranceBy:WrongInput',...
                      'Input for deviation must be a number.');                                   
              elseif nodeviation
                warning('ReverseMortgage:Contract:PlotInsuranceBy:DeviationNotPossible',...
                        'No elasticities can be computed for %s. Will switch to regular output.',...
                        factor);                
              end                           
            
              %% get parameters
              switch factor
                  case 'age'
                      parameter = (obj.model.age_min:obj.model.age_max)';
                      label     = 'age';                                                  
                  case 'discount rate'
                      if callforelasticities
                          base  = obj.model.ir_discount_basic;
                          step  = base * 0.01; % 1% steps
                          start = base - (varargin{1}*100)*step;
                          final = base + (varargin{1}*100)*step;                          
                          parameter = (start:step:final)'; 
                          label      = 'discount rate deviation';                          
                      else                    
                        parameter  = ((obj.model.ir_discount_basic-0.03):0.002:(obj.model.ir_discount_basic+0.05))';
                        label      = 'discount rate';
                      end
                  case 'moveout'
                      if callforelasticities
                          base  = obj.model.moveout;
                          step  = base * 0.01; % 1% steps
                          start = base - (varargin{1}*100)*step;
                          final = base + (varargin{1}*100)*step;                          
                          parameter = (start:step:final)'; 
                          label      = 'moveout factor deviation';                          
                      else
                          parameter  = (0:0.05:2)';                          
                          label      = 'moveout factor';
                      end                      
                  case 'mu'
                      if callforelasticities
                          base  = obj.model.mu;
                          step  = base * 0.01; % 1% steps
                          start = base - (varargin{1}*100)*step;
                          final = base + (varargin{1}*100)*step;                          
                          parameter = (start:step:final)'; 
                          label     = 'mu: yearly drift parameter deviation';
                      else                        
                          parameter = (-0.02:0.0005:0.06)';                          
                          label     = 'mu: yearly drift parameter of house price model';
                      end
                  case 'sigma'
                      if callforelasticities
                          base  = obj.model.sigma;
                          step  = base * 0.01; % 1% steps
                          start = base - (varargin{1}*100)*step;
                          final = base + (varargin{1}*100)*step;                          
                          parameter = (start:step:final)'; 
                          label     = 'sigma: yearly diffusion parameter deviation';
                      else                       
                          parameter = (0.01:0.005:0.2)';                          
                          label     = 'sigma: yearly diffusion parameter of house price model';
                      end                    
              end
              
              %% preallocate output
              epremiums    = zeros(length(parameter),1);
              elosses      = zeros(length(parameter),1);
              eutilization = zeros(length(parameter),1);
              
              %% construct payment scheme              
              if strcmpi(obj.payment_plan, 'loc')
                error('ReverseMortgage:Contract:PlotInsuranceBy:NoLOCAccountComputable',...
                      'Cannot plot this analysis for a LOC account. Create different contract!');
              end
              if any(strcmpi(obj.payment_plan, {'modifiedtenure','modifiedterm'}))
                scheme.loc = obj.loc_amount;
              end 
              if any(strcmpi(obj.payment_plan, {'modifiedterm','term'}))
                scheme.periods = obj.payment_periods;
              end              
              scheme.name = obj.payment_plan;             
                       
              %% get data
              switch factor
                  case 'age'
                    for i = 1:length(parameter)
                      % create contracts
                      mortgage = contract(parameter(i),...
                                    obj.gender_for_calculation,...
                                    obj.house_initialvalue,...
                                    obj.ir_expected_rate,...
                                    scheme, obj.model);
                                
                      % assign data
                      epremiums(i)    = mortgage.acc_payment.premiums_expected;
                      elosses(i)      = mortgage.acc_payment.loss_expected;
                      eutilization(i) = mortgage.acc_payment.utilization;
                    end  
                  case 'discount rate'
                    
                    % save current parameter
                    old = obj.model.ir_discount_basic;
                    
                    for i = 1:length(parameter)               
                      % assign new temporary parameter
                      obj.model.ir_discount_basic = parameter(i);
                      
                      % get data
                      epremiums(i)    = get_premiums_expected(obj.acc_payment);
                      elosses(i)      = get_loss_expected(obj.acc_payment);
                      eutilization(i) = elosses(i)/epremiums(i).*100;
                    end  
                    
                    % restore current parameter
                    obj.model.ir_discount_basic = old;                    
                  case 'moveout'
                    % save current probability model + parameter
                    oldmoveout = obj.model.moveout;
                    oldprobs   = obj.probs;
                                        
                    for i = 1:length(parameter)               
                      % recompute
                      obj.model.moveout = parameter(i);
                      obj.probs = transform(obj.probs);
                      
                      % get data
                      epremiums(i) = get_premiums_expected(obj.acc_payment);
                      elosses(i)   = get_loss_expected(obj.acc_payment);
                      eutilization(i) = elosses(i)/epremiums(i).*100;
                    end  
                    
                    % restore current parameter
                    obj.model.moveout = oldmoveout;                    
                    obj.probs = oldprobs;                                       
                  case 'mu'
                    % save current parameters
                    oldmu = obj.model.mu;
                                        
                    for i = 1:length(parameter)               
                      % recompute
                      obj.model.mu = parameter(i);
                      
                      % get data
                      epremiums(i) = get_premiums_expected(obj.acc_payment);
                      elosses(i)   = get_loss_expected(obj.acc_payment);
                      eutilization(i) = elosses(i)/epremiums(i).*100;
                    end  
                    
                    % restore current parameter
                    obj.model.mu = oldmu;                    
                  case 'sigma'
                    % save current parameters
                    oldsigma = obj.model.sigma;
                                        
                    for i = 1:length(parameter)               
                      % recompute
                      obj.model.sigma = parameter(i);
                      
                      % get data
                      epremiums(i) = get_premiums_expected(obj.acc_payment);
                      elosses(i)   = get_loss_expected(obj.acc_payment);
                      eutilization(i) = elosses(i)/epremiums(i).*100;
                    end  
                    
                    % restore current parameter
                    obj.model.sigma = oldsigma;                       
              end
              
              %% adjust output according to option
              if callforelasticities
                  xvalues = ((-varargin{1}*100):1:(varargin{1}*100))';                                    
              elseif ~callforelasticities
                  xvalues = parameter;
              end
              
              %% plot result
              h = plot(xvalues, eutilization);
              set(h, 'LineWidth', 2);                              
              legend(h, 'utilization level in % of net present value');                    
              title('sensitivity of utilization',...
                     'FontWeight', 'bold');
              ylabel('utilization level');
              xlabel(label);
                                      
              %% output
              varargout{1}.parameter    = xvalues;
              varargout{1}.epremiums    = epremiums;
              varargout{1}.elosses      = elosses;              
              varargout{1}.eutilization = eutilization;                                 
              
         end % function PLOTINSURANCEBY         
         
         %% PROPERTIES                       
         function house_expected      = get_house_expected(obj)
             % HOUSE_EXPECTED calculate path of expected collateral value
             % see Szymanoski, 1990, p.14, equation (6) for further detail
                 
                 % l is logical: periods during which contract is in force
                 l = (1 : obj.horizon_months)'; 
                 l = (l > (12 * obj.age_for_calculation));   

                 % t is array: number of months that contract is in force 
                 t = (1:obj.horizon_months)' - (12 * obj.age_for_calculation + 1); 
                 t =  t .* l;

                 % transform mu, sigma into monthly values
                 mu    = obj.model.mu/12;
                 sigma = obj.model.sigma/sqrt(12);
                 
                 % calculate end-of-month expected house price E [ H(t) ]
                 house_expected  = ...
                     l .* obj.house_initialvalue .* exp(mu * t + 1/2 * (sigma.^2) * t);
             
         end % function HOUSE_EXPECTED
         function house_expected_cond = get_house_expected_cond(obj)
             % HOUSE_EXPECTED_COND calculate path of expected collateral
             % value conditioned on a loss (balance > collateral value)
             % assuming a geometric Brownian motion process and basic
             % payment model
             %
             % for more information and formula, see method
             % 'price_collateral'
             
             [loss_expected, lossprobs, house_expected_cond] =...
                 price_collateral(obj.acc_basic.b_e, obj);
         
         end % function GET_HOUSE_EXPECTED_COND        
         function premiumstructure    = get_premiumstructure(obj)
             % GET_PREMIUMSTRUCTURE Builds / Updates time stucture of
             % premiums
             %
             % returns array which sets the payment structure of the 
             % insurance premium fees each month as a percentage of initial
             % appraised value and the monthly outstanding balance

                premiumstructure = zeros(obj.horizon_months,1);
                premiumstructure((12 * obj.age_for_calculation + 2):end) = obj.model.monthly;

         end % function GET_PREMIUMSTRUCTURE           
         function setaside            = get_setaside(obj)
             % GET_SETASIDE present value of future fees
             
                    % l: periods during which contract is in force
                    policymonth = (1:obj.horizon_months)'; 
                    inforce = (policymonth > (12 * obj.age_for_calculation));   

                    % t is array: number of months that contract is in
                    % force (at start of month when payments are made)
                    t = (1:obj.horizon_months)' - (12 * obj.age_for_calculation + 1); 
                    t =  t .* inforce;
                    
                    % m is the number of remaining months that the servicing fee
                    % could be collected
                    remaining = obj.horizon_months - 12*obj.age_for_calculation;                    
                    m         = zeros(obj.horizon_months,1);
                    m((12*obj.age_for_calculation+1):end) = (remaining:-1:1)';
                    
                    % compute set aside: HUD Handbook 4325.1, Appendix 22,
                    % p.1
                    setaside = (m>0) .* obj.model.setaside .* ...
                             ( ((1+obj.ir_compounding/12).^(m+1)) - (1+obj.ir_compounding/12) ) ...
                             ./ ((obj.ir_compounding/12) .* ((1+obj.ir_compounding/12).^m));
                           
        end
                  
    end % methods
    
    methods   
        
         %% GET methods 
             
         % If property is empty, then compute
         function plf_exogenous       = get.plf_exogenous(obj)
             % GET.PLF_EXOGENOUS If empty, then iterate
             %
             % PLF_EXOGENOUS is determined by characteristics of contract and model's
             % assumptions; it will be constant over lifetime of the model
             %
             % GET.PLF_EXOGENOUS checks for existence of a previous iteration in the 
             % model object. If no such value exists, then the principal
             % limit factor is iterated and stored in the model for future
             % use
             
                 % set: indexes for indexing into plf matrix
                 if strcmpi(obj.gender_for_calculation, 'male')==1
                   gender = 1;
                 elseif strcmpi(obj.gender_for_calculation, 'female')==1
                   gender = 2; 
                 end
                 age     = obj.age_for_calculation;
                 diffabs = -abs(obj.model.ir_expected_rate-obj.ir_expected_rate);
                 maximum = max(diffabs);
                 ind     = find(diffabs==maximum);
                                  
                 % has plf already been set in the 'model' instance?                                     
                     % yes! --> get and save value 
                     if ~isnan(obj.model.plf(age,gender,ind));                    
                     plf_exogenous = obj.model.plf(age,gender,ind); 
                     return
                     % no!  plf has not been set: price this contract
                     % through iteration
                     else                    
                     [plf_exogenous, obj] = price(obj,'exogenous'); 
                     obj.model.plf(age,gender,ind) = plf_exogenous;   
                     end            
             
         end % function GET.PLF                               
         function house_expected      = get.house_expected(obj)
                 % GET.HOUSE_EXPECTED If empty, then fill with current value
                 %
                 % HOUSE_EXPECTED will only change if assumptions of model
                 % change; hence, value is not expected to change

                    if isempty(obj.house_expected)
                        house_expected = get_house_expected(obj); return
                    else
                        house_expected = obj.house_expected;
                    end
         end % function GET.HOUSE_EXPECTED
         function house_expected_cond = get.house_expected_cond(obj)
                 % GET.HOUSE_EXPECTED_COND If empty, then fill with current value
                 %
                 % HOUSE_EXPECTED_COND will only change if assumptions of model
                 % change; hence, value is not expected to change

                    if isempty(obj.house_expected_cond)
                        house_expected_cond = get_house_expected_cond(obj); return
                    else
                        house_expected_cond = obj.house_expected_cond;
                    end
         end % function GET.HOUSE_EXPECTED_COND                                                                           
         function setaside            = get.setaside(obj)   
             % GET.SETASIDE If empty, then compute
             
                    if isempty(obj.setaside)
                        setaside = get_setaside(obj); return
                    else
                        setaside = obj.setaside;
                    end                        
         end % function GET.SETASIDE
         
         % Dependent properties
         function house_max_claim_amount = get.house_max_claim_amount(obj)
             % GET.HOUSE_MAX_CLAIM_AMOUNT
             
                house_max_claim_amount = obj.model.util_collateral * obj.house_initialvalue;
                
         end % function GET.HOUSE_MAX_CLAIM_AMOUNT
         function ir_compounding         = get.ir_compounding(obj)
             
                % compounding interest rate;
                % see Szymanoski (1990), p. 27
                ir_compounding = obj.ir_expected_rate + 12*obj.model.monthly;
                
         end % function GET.IR_COMPOUNDING         
         function horizon_months         = get.horizon_months(obj)
                 
                  horizon_months = obj.horizon_years*12;
                  
         end % function GET.HORIZON_MONTHS
         function horizon_years          = get.horizon_years(obj)
                
                horizon_years = length(obj.probs.p);              
                
         end  % function GET.HORIZON_YEARS
           
    end % methods GET
  
    methods (Static)        
         function presentvalue = get_pv(nominalvalue,i,firstperiod)
             % GET_PV  computes present value cashflow vector for given 
             % nominal cashflow vector, given a periodic interest rate i and
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