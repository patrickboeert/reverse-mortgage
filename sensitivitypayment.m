%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reproduces max. monthly payment sensitivity analysis in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model 1st German
german            = model();
german.lifetable  = 'survival_prob_2004R_firstorder.csv'; 
german.mu         = 0.024;

% contract, tenure for 65-year old female
age         = 65;
gender      = 'female';
i           = 0.07;
c           = 200000;
scheme.name = 'tenure';
vertrag     = contract(age, gender, c, i, scheme, german);

% get insurance sensitivity data

% age
subplot(3,2,1);
byage         = plotpaymentby(vertrag, 'age');

h(1) = plot(byage.parameter, byage.termpayment);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'magenta');                               
         hold on
h(2) = plot(byage.parameter, byage.tenurepayment); 
         set(h(2), 'LineWidth', 2);
         set(h(2), 'Color', 'black');                                                                                
         hold off    
 legend(h, {'10-year term payment plan',...
        'tenure payment plan'},...
        'Location','NorthWest');                    
 ylabel('maximum monthly cash advance');
 xlabel('initial age of borrower');

% expected rate
subplot(3,2,2);
interestrate  = plotpaymentby(vertrag, 'interest rate');

h(1) = plot(interestrate.parameter, interestrate.termpayment);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'magenta');                               
         hold on
h(2) = plot(interestrate.parameter, interestrate.tenurepayment); 
         set(h(2), 'LineWidth', 2);
         set(h(2), 'Color', 'black');                                                                                
         hold off    
 legend(h, {'10-year term payment plan',...
        'tenure payment plan'},...
        'Location','NorthWest');                    
 ylabel('maximum monthly cash advance');
 xlabel('expected rate');
 
% move out rate
subplot(3,2,3);
moveout = plotpaymentby(vertrag, 'moveout');

h(1) = plot(moveout.parameter, moveout.termpayment);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'magenta');                               
         hold on
h(2) = plot(moveout.parameter, moveout.tenurepayment); 
         set(h(2), 'LineWidth', 2);
         set(h(2), 'Color', 'black');                                                                                
         hold off    
 legend(h, {'10-year term payment plan',...
        'tenure payment plan'},...
        'Location','NorthWest');                    
 ylabel('maximum monthly cash advance');
 xlabel('move-out rate');
 
% monthly premium fee
subplot(3,2,4);
fee  = plotpaymentby(vertrag, 'fee');

h(1) = plot(fee.parameter, fee.termpayment);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'magenta');                               
         hold on
h(2) = plot(fee.parameter, fee.tenurepayment); 
         set(h(2), 'LineWidth', 2);
         set(h(2), 'Color', 'black');                                                                                
         hold off    
 legend(h, {'10-year term payment plan',...
        'tenure payment plan'},...
        'Location','NorthWest');                    
 ylabel('maximum monthly cash advance');
 xlabel('monthly premium fee');
 
 % mu
subplot(3,2,5);
mu  = plotpaymentby(vertrag, 'mu');

h(1) = plot(mu.parameter, mu.termpayment);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'magenta');                               
         hold on
h(2) = plot(mu.parameter, mu.tenurepayment); 
         set(h(2), 'LineWidth', 2);
         set(h(2), 'Color', 'black');                                                                                
         hold off    
 legend(h, {'10-year term payment plan',...
        'tenure payment plan'},...
        'Location','NorthWest');                    
 ylabel('maximum monthly cash advance');
 xlabel('drift parameter \mu');
 
 % sigma
subplot(3,2,6);
sigma  = plotpaymentby(vertrag, 'sigma');

h(1) = plot(sigma.parameter, sigma.termpayment);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'magenta');                               
         hold on
h(2) = plot(sigma.parameter, sigma.tenurepayment); 
         set(h(2), 'LineWidth', 2);
         set(h(2), 'Color', 'black');                                                                                
         hold off    
 legend(h, {'10-year term payment plan',...
        'tenure payment plan'},...
        'Location','NorthWest');                    
 ylabel('maximum monthly cash advance');
 xlabel('diffusion parameter \sigma');
 
% plot


