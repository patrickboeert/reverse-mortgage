classdef insurer < handle
  % INSURER  insurance object
  
  properties
    
    %% simulation data
    id         = [];            % simulation id    
    simulation = [];            % simulation data
    lender     = [];            % lender instance    
    
    %% assignment data
    assignments = {};           % assigned contracts
        
    %% business indicators
    % +
    insurancefund       = [];   % insurance premiums (12 * nSimHorizon x nSimContracts)
    payoff              = [];   % claim 2 payoffs (assignment), (12 * nSimHorizon x nSimContracts)
    
    % -
    acquisitioncost     = [];   % acquisition cost (12 * nSimHorizon x nSimContracts)  
    advancestoborrower  = [];   % advances (12 * nSimHorizon x nSimContracts)    
    
    % others
    loantoborrower      = [];   % outstanding balance (end of period); (12 * nSimHorizon x nSimContracts)    
    interesttoinsurer   = [];   % interest (12 * nSimHorizon x nSimContracts)
    nAssigned           = [];   % # of assigned contracts after simulation    
    
    %% cashflow results
    inflow     = [];   % premiums + payoff
    outflow    = [];   % acquisition cost for loan + advances
    cashflow   = [];   % inflow - outflow
    
    inflow_y   = [];   % premiums + payoff
    outflow_y  = [];   % acquisition cost for loan + advances
    cashflow_y = [];   % inflow - outflow    
  
  end
  
  methods
    
     %% CONSTRUCTOR
     function obj = insurer(sim, id)
       % INSURER  constructor
       
       % set and preallocate
       obj.simulation  = sim;
       obj.id          = id;
       obj.assignments = cell(obj.simulation.nSimContracts,1);

       % preallocate business time series
       obj.loantoborrower     = NaN(12*obj.simulation.nSimHorizon,...
                                    obj.simulation.nSimContracts);
       obj.advancestoborrower = NaN(size(obj.loantoborrower));
       obj.interesttoinsurer  = NaN(size(obj.loantoborrower));       
       obj.payoff             = NaN(size(obj.loantoborrower));
       obj.acquisitioncost    = NaN(size(obj.loantoborrower));     
       
     end % end INSURER
     
     %% BEHAVIOUR
     function obj = takeassignment(obj, contract)
       % TAKEASSIGNMENT  take assignment from lender and process contract
       
       %% assign contract
       start = contract.period_assigned + 1;
       final = contract.period_terminated;       
       
       obj.assignments{contract.id}...                    % assign
         = contract; 
       
       obj.acquisitioncost(...                            % acqisition cost
               contract.s(start),...
               contract.id)...
                    = contract.acc_simulation.b_e(start);      
       
       %% set cashflows      
       % loan balance after assignment
       obj.loantoborrower(contract.s(start):contract.s(final))...
         = contract.acc_simulation.b_e(start:final);
       
       % advances after assignment
       obj.advancestoborrower(contract.s(start):contract.s(final))...
         = contract.acc_simulation.advances(start:final); 
       
       % interest after assignment
       obj.interesttoinsurer(contract.s(start):contract.s(final))...
         = contract.acc_simulation.interest(start:final);            
       
       %% compute payoff at Claim 2 Type
       house       = contract.house_simulated(contract.s(final));
       
       outstanding = contract.acc_simulation.b_e(final)...
                                .* (1 - contract.model.salesexpense);
       
       % compute payoff at termination
       obj.payoff(contract.s(final)) = min(house, outstanding);                      

     end % end TAKEASSIGNMENT
     
     %% TRANSFORMATION
     function yearlydata = yearly(obj, monthlydata)
       % YEARLYDATA  computes outflow property
       
       % prepare logical indexing
       continous  = (1/12:1/12:obj.simulation.nSimHorizon);
       year       = ceil(continous);
       yearlydata = NaN(obj.simulation.nSimHorizon, 1);
       
       % sum over contracts for each year
       for y = 1:numel(unique(year))
         yearlydata(y) =...
           nansum(nansum(monthlydata(year==y,:),2)); % sum over contracts
       end
       
     end % end YEARLYDATA      
     
     %% GET
     function nAssigned = get.nAssigned(obj)
       % NASSIGNED  compute number of assigned contracts after simulation
                          
       % compute
       assigned = zeros(size(obj.assignments));
       
       for i = 1:length(obj.assignments)
         assigned(i) = isa(obj.assignments{i}, 'contract');
       end  
       
       nAssigned = sum(assigned);
     
     end % end NASSIGNED       
     
     function insurancefund = get.insurancefund(obj)
      % GET.INSURANCEFUND  computes dependent insurancefund property
         insurancefund = obj.lender.premiumstoinsurer;
                    
     end % function GET.INSURANCEFUND
     function inflow = get.inflow(obj)
       % GET.INFLOW  computes inflow property
       inflow = nansum(obj.insurancefund,2)...
                  + nansum(obj.payoff, 2);
       
     end % end GET.INFLOW
     function outflow = get.outflow(obj)
       % GET.OUTFLOW  computes outflow property
       outflow = nansum(obj.acquisitioncost,2)...
                  + nansum(obj.advancestoborrower,2);
                
     end % end GET.OUTFLOW
     function cashflow = get.cashflow(obj)
       % GET.CASHFLOW  computes cashflow property
       cashflow = obj.inflow - obj.outflow;
       
     end % end CASHFLOW 
     
     function outflow_y = get.outflow_y(obj)
       % GET.OUTFLOW_Y  computes outflow_y property
       outflow_y = obj.yearly(obj.outflow);
       
     end % end GET.OUTFLOW     
     function inflow_y = get.inflow_y(obj)
       % GET.INFLOW_Y  computes inflow property
       inflow_y = obj.yearly(obj.inflow);
       
     end % end GET.INFLOW_Y
     function cashflow_y = get.cashflow_y(obj)
       % GET.CASHFLOW_Y  computes cashflow property
       cashflow_y = obj.yearly(obj.cashflow);
       
     end % end CASHFLOW_Y      
  
  end
end

