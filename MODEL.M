classdef model
% MODEL  defines the exogenous parameters and other configuration options
%        of the HECM Program's risk management model
%
%   NOTE: before using class model, set property 'basedir' to directory
%         of source code
%
%   INPUT: MODEL is initialized without input with default parameters. If
%          you want to change the parameters, change them directly through
%          assignment of values after creating the instance
%           
%   OUTPUT: HECM Program risk management model parameters and configuration
%
%   FUNCTIONALITY: MODEL separates the parameters and configuration of the 
%                  HECM Program from the logic of the model. MODEL is
%                  required input for all other classes, especially class
%                  CONTRACT

    properties        
        
        %% config
        basedir = 'J:\simulation\';    % storage path
    
        %% lifetable
        lifetable       = '';  % lifetable

        % POSSIBLE ALTERNATIVES:
        %
        % US Decennial Lifetable 1979-1981 | 'survival_prob_US_7981.csv'           
        %
        % US Decennial Lifetable 1999-2001 | 'survival_prob_US_9901.csv'     
        %
        % DAV 2004 R Sterbetafel, Basistafel 1. Ordung, Aggregattafel (mit
        % Sicherheitsabschlaegen) | 'survival_prob_2004R_firstorder.csv'
        %    
        % DAV 2004 R Sterbetafel, Basistafel 2. Ordung, Aggregattafel (ohne
        % Sicherheitsabschlaege) | 'survival_prob_2004R_secondorder.csv' 
        %    
        % Sterbetafel Deutschland 2005-2007, Statistisches Bundesamt |
        % 'survival_prob_destatis_2007.csv'
           
        %% probability model           
        moveout  = 0.3;    % paper: move-out factor (as fraction of death probability)
          
        %% costs
        setaside = 0;             % monthly HECM set aside fee vector
        cost     = [];            % cost for origination of contract
        
        %% premium structure
        upfront = 0.02;           % paper: \alpha; premium rate t=1, as fraction of maximum claim amount 
        monthly = 0.005/12;       % paper: \beta/12; monthly premium rate t>1, as fraction of outstanding balance(t)
          
        %% collateral model
        util_collateral = 0.85;   % fractional utilization of collateral: max_claim_amount policy
        threshold       = 0.98;   % fraction of max claim amount at which assignment to insurer occurs
        salesexpense    = 0.00;   % sale expense for collateral as fraction of house price        
        mu              = 0.00;   % paper: \mu; drift of geometrics Brownian motion process (yearly)
        sigma           = 0.10;   % paper: \sigma^2; diffusion of geometric Brownian motion (yearly)
              
        %% principal limit factor - matrix           
        % plf matrix
        plf = [];                 % 3-dim plf matrix (age,gender,interestrate)
        
        % plf matrix - age
        age_min = 60;             % minimum age for plf matrix
        age_max = 130;            % maximum age for plf matrix 
            
        % plf matrix - interest rates
        ir_expected_rate       = [];        % expected interest rate (vector)
        ir_expected_rate_min   = 0.00;      % minimum expected rate for plf matrix
        ir_expected_rate_max   = 0.16;      % maximum expected rate for plf matrix
        ir_expected_rate_steps = 0.00125;   % intermediate steps for expected rate of plf matrix, 0.125 (1/8) is HUD measure        
        ir_discount_basic      = 0.06733630137;  % discount rate, yearly
        
        % plf  - iteration options
        iter_output    = 'notify';   % output options for plf iteration step:
                                        % 'off'    displays no output for iteration output;
                                        % 'iter'   displays output at each iteration;
                                        % 'final'  displays just the final output;
                                        % 'notify' displays output only if the function does not converge                           
        iter_tolerance = 1e-6;       % termination tolerance of principal limit factor iteration       
                
    end % properties
    
    methods
        
        %% CONSTRUCTOR
        function obj = model()
            % MODEL  create model instance from input parameters
        
                 %% set expected rate vector 
                 obj.ir_expected_rate = ...
                     (obj.ir_expected_rate_min : obj.ir_expected_rate_steps : obj.ir_expected_rate_max)';                  

                 %% preallocate plf matrix               
                 obj.plf = NaN(obj.age_max, 2, length(obj.ir_expected_rate));
                 
                 %% assign lifetable
                 obj.lifetable = 'survival_prob_2004R_firstorder.csv';

        end % function MODEL
        
        %% PRICING
        function obj = price(obj)
             % PRICE  computes all principal limit factors for the model's
             % plf matrix
                                           
                 % set parameters
                 gender = {'male' 'female'}';
                 scheme.name = 'tenure';                 
                 
                 % fill plf matrix                            
                 for age = obj.age_min : obj.age_max
                     for gender_ind = 1 : length(gender)
                        for ir_index = 1 : length(obj.ir_expected_rate)                                                                                                     
                              
                            % create reverse mortgage with right parameters
                            rm = contract(age, gender{gender_ind}, 100000, obj.ir_expected_rate(ir_index),...
                                scheme, obj);                                    
                            
                            % set output: plf at 'model' level
                            obj.plf(age,gender_ind,ir_index) = rm.plf_exogenous;                                                                          
                        end
                     end
                 end

        end % function PRICE
        
        %% ANALYSIS
        function binary = ismonotonic(obj)
            % ISMONOTONIC  Is the plf matrix strictly monotonic in age and
            % the expected rate?
            
                % check for monotonicity
                age = ismonotonic_age(obj);
                ir  = ismonotonic_interestrate(obj);
                
                % output
                binary = all([age, ir]);
        end % function ISMONOTONIC              
        function binary = ismonotonic_gender(obj)
            % ISMONOTONIC_GENDER  Is the plf strictly monotonic by gender 
            % (male > female)?
                
                % get data
                plf_diff = get_plf_diff(obj);
                
                % compute
                binary = all(all(plf_diff > 0));
                
        end % function ISMONOTONIC_GENDER
        function binary = ismonotonic_interestrate(obj)
            % ISMONOTONIC_INTERESTRATE  Is the plf strictly monotonically 
            % decreasing in the expected rate?
            
                % arrange data
                plf{1} = squeeze(obj.plf(obj.age_min:obj.age_max,1,:));
                plf{2} = squeeze(obj.plf(obj.age_min:obj.age_max,2,:));
                
                % check for male, female
                for i=1:2
                    monotonic(i) = all(ismonotonic(plf{i}, 1, 'DECREASING', 2));
                end
                
                binary = all(monotonic);
     
        end % function ISMONOTONIC_INTERESTRATE
        function binary = ismonotonic_age(obj)
            % ISMONOTONIC_AGE  Is the plf strictly monotonically 
            % increasing in the age of origination?
            
                % arrange data
                plf{1} = squeeze(obj.plf(obj.age_min:obj.age_max,1,:));
                plf{2} = squeeze(obj.plf(obj.age_min:obj.age_max,2,:));
                
                % check for male, female
                for i=1:2
                    monotonic(i) = all(ismonotonic(plf{i}, 1, 'INCREASING'));
                end
                
                binary = all(monotonic);
                
        end % function ISMONOTONIC_AGE             
        
        %% VARIABLES
        function plf_diff = get_plf_diff(obj)
            % PLF_DIFF  computes difference of plf between males and
            % females, by age and expected interest rate
            
                % compute difference
                male   = squeeze(obj.plf(obj.age_min:obj.age_max,1,:));
                female = squeeze(obj.plf(obj.age_min:obj.age_max,2,:));
                plf_diff = male - female;
                
        end % function PLF_DIFF
                
        %% PLOTTING
        function plot(obj)
            % MESH  creates surf plot of principal limit factors
            %
            % MESH creates surf plot of principal limit factors for female
            % and male borrowers, as well as for the difference of
            % principal limit factors between male and female borrowers
            
                % get data by gender, squeeze
                male   = squeeze(obj.plf(obj.age_min:obj.age_max,1,:));
                female = squeeze(obj.plf(obj.age_min:obj.age_max,2,:));
                diff   = male - female;
                
                % prepare grid axes
                [age_axes, ir_axes] = meshgrid(obj.age_min:obj.age_max ,...
                                               (obj.ir_expected_rate.*100)');
                                           
                % plot surf                
                subplot(3,1,1);                            
                    h(1) = surf(age_axes, ir_axes, male');
                    title('Principal Limit Factor, Male',...
                              'FontWeight', 'bold');

                    % labels
                    zlabel('Principal Limit Factor');                        
                    xlabel('age at origination','fontsize',10);
                    ylabel('expected interest rate','fontsize',10);                        

                    % colorbar
                    colormap Jet
                    %colorbar  
                    
                    % others
                    set(h(1),'EdgeColor','none'); % no grid on data
                    axis tight square % data close to axes
                    view(-60,10) % view direction
                        
                subplot(3,1,2);                            
                    h(2) = surf(age_axes, ir_axes, female');
                    title('Principal Limit Factor, Female',...
                              'FontWeight', 'bold');

                    % labels
                    zlabel('Principal Limit Factor');
                    xlabel('age at origination','fontsize',10);
                    ylabel('expected interest rate','fontsize',10);

                    % colorbar
                    colormap Jet
                    %colorbar

                    % others
                    set(h(2),'EdgeColor','none'); % no grid on data
                    axis tight square % data close to axes
                    view(-60,10) % view direction                    
                      
                subplot(3,1,3);
                    h(3) = surf(age_axes, ir_axes, diff');
                    title('Principal Limit Factor Difference',...
                              'FontWeight', 'bold');

                    % labels
                    zlabel('Principal Limit Factor Difference');
                    xlabel('age at origination','fontsize',10);
                    ylabel('expected interest rate','fontsize',10);                        

                    % colorbar
                    colormap Jet
                    %colorbar
                    
                    % others
                    set(h(3),'EdgeColor','none'); % no grid on data
                    axis tight square % data close to axes
                    view(-60,10) % view direction                    
        
        end % function SURF
             
    end % methods
end % classdef