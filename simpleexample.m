%% BASIC USE CASES

%% VALUE A SINGLE CONTRACT
%
% EXAMPLE and USAGE: 
%       
%  rm = contract(75, 'female', 100000, 0.1, scheme, hecm)
%
%  returns a reverse mortgage contract for a 75 year old female with a
%  house valued at 100.000 in an interest rate environment (expected rate)
%  of 10%. The contract is computed for the payment scheme choice 'scheme',
%  where 'scheme' is a structure with the fields 'name', 'loc' and 'periods'
%  and needs to be set before calling CONTRACT. The model configuration
%  'hecm' has to be initialized before calling CONTRACT. CONTRACT returns 
%  a reverse mortgage based on the model instance 'hecm' which has to be
%  initialized before calling CONTRACT.

hecm = model();             % initialize mode configuration
scheme.name = 'term';       % set form of payment scheme
scheme.periods = 120;       % 10-year term payment
rm = contract(75, 'female', 100000, 0.1, scheme, hecm); % price contract

%% RUN A SIMULATION
%
%  runs a simulation for model 'hecm'; change properties in 'model' and save
%  class file to run simulation with different:
%    nRUN          = 10;  % # of simulation trials
%    nContracts    = 15;  % # of contracts in each single simulation trial
%    nYears        = 60;  % # simulated years in every simulation 
%    r_prior       = 5.508789954; % initial rate of 1-year rate for MCMC simulation
%    rl_prior      = 6.733630137; % initial rate of 10-year rate for MCMC simulation
%    margin        = 1.5;         % lender's profit margin

hecm = model();                     % initialize mode configuration
sim  = simulation(hecm);            % initialize simulation
sim  = sim.get_demand([0.5,0.5]);   % sample characteristics of contracts for simulation
sim  = sim.originate;               % originate contracts
sim  = sim.simulate;                % simulate contract lifetime & time series
sim  = sim.close;                   % extract data from portfolios
sim  = sim.aggregate;               % aggregate contract data across portfolios

