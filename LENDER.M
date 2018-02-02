classdef lender < handle
%LENDER   bank object

   properties
         
    %% simulation: relations
    id         = [];  % portfolio id    
    simulation = [];  % simulation instance
       
    %% business indicators 
    payoff              = [];   % payoff: regular, claim 1, claim 2, (12 * nYears x nContracts)    
    loantoborrower      = [];   % outstanding balance (end of period); (12 * nYears x nContracts)
    advancestoborrower  = [];   % advances (12 * nYears x nContracts)    
    premiumstoinsurer   = [];   % premiums (12 * nYears x nContracts)
    interestforlender   = [];   % interest (12 * nYears x nContracts)
    type_terminated     = '';   % 'Payoff', 'Claim 1 Type', 'Claim 2 Type' (assignment)    
    
    %% cashflows
    outflow     = [];   % advances + premiums [+ interest for refinancing (not implemented)]
    inflow      = [];   % payoff
    cashflow    = [];   % inflow - outflow
    
    outflow_y   = [];   % advances + premiums [+ interest for refinancing (not implemented)]
    inflow_y    = [];   % payoff
    cashflow_y  = [];   % inflow - outflow    

   end

   methods
      
     %% CONSTRUCTOR
     function obj = lender(sim, id)
       % LENDER  constructs lender instance
       
         % assign simulation data
         obj.id         = id;
         obj.simulation = sim;

         % preallocate business indicators
         obj.loantoborrower     = NaN(12*obj.simulation.nYears,...
                                      obj.simulation.nContracts);                           
         obj.advancestoborrower = NaN(size(obj.loantoborrower));
         obj.premiumstoinsurer  = NaN(size(obj.loantoborrower));
         obj.interestforlender  = NaN(size(obj.loantoborrower));
         obj.payoff             = NaN(size(obj.loantoborrower));
                         
     end % end LENDER 
     
     %% BEHAVIOUR
     function obj = assign(obj, contract)
       % ASSIGN  assigns contract to insurer, computes lender's cashflow
       % and calls on insurer to close the contract
       
           %% error checking: has contract been simulated?
            if isempty(contract.period_terminated)
              error('ReverseMortgage:Lender:Assign:ContractNotSimulated',...
                    ['To assign a contract to the insurer, the contract'...
                     ' first has to be simulated to generate a history.'...
                     ' Call contract-->simulate first!']);
            end
            
           %% assign cashflows for lender
            origination = contract.period_origination;
            assigned    = contract.period_assigned;
           
            obj.loantoborrower(...
              contract.s(origination+1):contract.s(assigned),...
              contract.id)...
              = contract.acc_simulation.b_e((origination+1):assigned);
            
            obj.advancestoborrower(...
              contract.s(origination+1):contract.s(assigned),...
              contract.id)...
              = contract.acc_simulation.advances((origination+1):assigned);
            
            obj.premiumstoinsurer(...
              contract.s(origination+1):contract.s(assigned),...
              contract.id)...
              = contract.acc_simulation.premiums((origination+1):assigned);
            
            obj.interestforlender(...
              contract.s(origination+1):contract.s(assigned),...
              contract.id)...
              = contract.acc_simulation.interest((origination+1):assigned);                       
            
           %% handle payoff at assignment
            contract.type_terminated = 'Claim 2 Type';           
           
            obj.payoff(contract.s(assigned+1), contract.id)...
              = contract.acc_simulation.b_e(assigned);
            
           %% clear contract from portfolio & force insurer to assign
           %% contract 
            obj.insurer = takeassignment(obj.insurer, contract);            
  
     end % end ASSIGN
     function obj = close(obj, contract)
       % CLOSE  close contract through 1) Claim 1 Type or Payoff (if not
       % assigned) or 2) Claim 2 Type (if assigned) and compute lender´s
       % properties
       
       %% error checking: has contract been simulated?
         if isempty(contract.period_terminated)
           error('ReverseMortgage:Lender:Close:ContractNotSimulated',...
                ['To assign a contract to the insurer, the contract'...
                 ' first has to be simulated to generate a history.'...
                 ' Call contract-->simulate first!']);
         end
         
       %% assign cashflows for lender
         start       = contract.period_origination + 1;       
         final       = contract.period_terminated;          
           
         obj.loantoborrower(contract.s(start):contract.s(final), contract.id)...
           = contract.acc_simulation.b_e(start:final);
           
         obj.advancestoborrower(contract.s(start):contract.s(final), contract.id)...
           = contract.acc_simulation.advances(start:final);
            
         obj.premiumstoinsurer(contract.s(start):contract.s(final), contract.id)...
           = contract.acc_simulation.premiums(start:final);
            
         obj.interestforlender(contract.s(start):contract.s(final), contract.id)...
           = contract.acc_simulation.interest(start:final);           
            
       %% determine type and value of payoff
         % arrange data         
         outstanding = contract.acc_simulation.b_e(final);
         housesale   = contract.house_simulated(contract.s(final))...
                          .* (1-contract.model.salesexpense);
       
         % fix minimum and find payoff type
         minimum = min(outstanding, housesale);
               
         % determine payoff value and payoff type
         switch minimum
           case outstanding           
             contract.type_terminated                  = 'Payoff';
             obj.payoff(contract.s(final),contract.id) = minimum;
           case housesale
             contract.type_terminated                  = 'Claim 1 Type';
             obj.payoff(contract.s(final),contract.id) = minimum;                      
         end
   
     end % end CLOSE
     
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
     function outflow = get.outflow(obj)
       % GET.OUTFLOW  computes outflow property
       outflow = nansum(obj.advancestoborrower,2)...
                   + nansum(obj.premiumstoinsurer,2);
       
     end % end GET.OUTFLOW
     function inflow = get.inflow(obj)
       % GET.INFLOW  computes inflow property
       inflow = nansum(obj.payoff,2);
       
     end % end GET.INFLOW
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
