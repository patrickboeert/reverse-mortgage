%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reproduces insurance sensitivity analysis in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model
german           = model();
german.lifetable = 'survival_prob_2004R_firstorder.csv'; 
german.mu        = 0.024;

% contract
age         = 65;
gender      = 'female';
i           = 0.07;
c           = 200000;
scheme.name = 'tenure';
vertrag     = contract(age, gender, c, i, scheme, german);

%% plot data
% mu
subplot(2,2,1);
mu = plotinsuranceby(vertrag, 'mu');

h(1) = plot(mu.parameter, mu.eutilization);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'black');                               
    legend(h(1), {'expected utilization rate u (%)'},...
        'Location','NorthWest');                    
    ylabel('expected utilization rate u (%)');
    xlabel('assumed drift \mu');
    
% sigma
subplot(2,2,2);
sigma = plotinsuranceby(vertrag, 'sigma');

h(1) = plot(sigma.parameter, sigma.eutilization);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'black');                               
    legend(h(1), {'expected utilization rate u (%)'},...
        'Location','NorthWest');                    
    ylabel('expected utilization rate u (%)');
    xlabel('assumed diffusion \sigma^2');
   
% moveout
subplot(2,2,3);
moveout  = plotinsuranceby(vertrag, 'moveout');

h(1) = plot(moveout.parameter, moveout.eutilization);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'black');                               
    legend(h(1), {'expected utilization rate u (%)'},...
        'Location','NorthWest');                    
    ylabel('expected utilization rate u (%)');
    xlabel('assumed move-out factor m');    

% discount interest rate
subplot(2,2,4);
discount = plotinsuranceby(vertrag, 'discount rate');

h(1) = plot(discount.parameter, discount.eutilization);
         set(h(1), 'LineWidth', 2);
         set(h(1), 'Color', 'black');                               
    legend(h(1), {'expected utilization rate u (%)'},...
        'Location','NorthWest');                    
    ylabel('expected utilization rate u (%)');
    xlabel('assumed discount interest rate i^d');  
    

