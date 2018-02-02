%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reproduces table of utilization results in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate models
one = model; one.lifetable='survival_prob_US_9901.csv';
two = model;
three = model; three.lifetable='survival_prob_2004R_firstorder.csv'; three.mu=0.024;
four = model; four.lifetable='survival_prob_2004R_secondorder.csv'; four.mu=0.024; 
five = model; five.lifetable='survival_prob_destatis_2007.csv'; five.mu=0.024;

% set parameters
age    = 75;
gender = 'male';
i      = 0.07;
c      = 200000;
level  = 1;                   % utilization level of principal limit

% how much of the collateral is accounted for by bank?
use = 0.8;  % bank takes x% of market collateral value into account
one.util_collateral = use;
two.util_collateral = use;
three.util_collateral = use; 
four.util_collateral = use;
five.util_collateral = use;

% set payment scheme
scheme.name='tenure';

% originate contracts
c1 = contract(age, gender, c, i, scheme, one, level);
c2 = contract(age, gender, c, i, scheme, two, level);
c3 = contract(age, gender, c, i, scheme, three, level);
c4 = contract(age, gender, c, i, scheme, four, level);
c5 = contract(age, gender, c, i, scheme, five, level);

% show maximum payment
disp(gender)
disp(age)
disp(scheme.name)
if strcmp(scheme.name,'term'), disp(scheme.periods), end;
[c1.acc_payment.utilization, c2.acc_payment.utilization,...
c3.acc_payment.utilization,...
c4.acc_payment.utilization,...
c5.acc_payment.utilization]

clear scheme c1 c2 c3 c4 c5

% set payment scheme
scheme.name='term';
scheme.periods=120;

% originate contracts
c1 = contract(age, gender, c, i, scheme, one, level);
c2 = contract(age, gender, c, i, scheme, two, level);
c3 = contract(age, gender, c, i, scheme, three, level);
c4 = contract(age, gender, c, i, scheme, four, level);
c5 = contract(age, gender, c, i, scheme, five, level);

% show maximum payment
disp(gender)
disp(age)
disp(scheme.name)
if strcmp(scheme.name,'term'), disp(scheme.periods), end;
[c1.acc_payment.utilization, c2.acc_payment.utilization,...
c3.acc_payment.utilization,...
c4.acc_payment.utilization,...
c5.acc_payment.utilization]

clear scheme c1 c2 c3 c4 c5

% set payment scheme
scheme.name='term';
scheme.periods=240;

% originate contracts
c1 = contract(age, gender, c, i, scheme, one, level);
c2 = contract(age, gender, c, i, scheme, two, level);
c3 = contract(age, gender, c, i, scheme, three, level);
c4 = contract(age, gender, c, i, scheme, four, level);
c5 = contract(age, gender, c, i, scheme, five, level);

% show maximum payment
disp(gender)
disp(age)
disp(scheme.name)
if strcmp(scheme.name,'term'), disp(scheme.periods), end;
[c1.acc_payment.utilization, c2.acc_payment.utilization,...
c3.acc_payment.utilization,...
c4.acc_payment.utilization,...
c5.acc_payment.utilization]

clear scheme c1 c2 c3 c4 c5




