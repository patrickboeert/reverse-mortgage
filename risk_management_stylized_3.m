function risk_management_stylized_3

%Figure 2.1 shows the relevant, stylized time series for a reverse mortgage. The
%example was computed for a 75-year old female (at the time of the origination
%of the loan) with a property initially valued at 100,000 in a 10% interest rate
%environment and with a tenure payment plan.

hecm         = model();
scheme.name  = 'tenure';
rm           = contract(75, 'female', 100000, 0.1, scheme, hecm);
st           = rm.period_origination+1;
fi           = rm.period_origination+25*12;

y = (0:1/12:1463/12);
years = y(st:fi)'- rm.period_origination/12;

b_e            = rm.acc_payment.b_e(st:fi);
advances       = cumsum(rm.acc_payment.advances(st:fi));
interest       = cumsum(rm.acc_payment.interest(st:fi));
premiums       = cumsum(rm.acc_payment.premiums(st:fi));
house_expected = rm.house_expected(st:fi);

h(1) = plot(years,b_e);
    hold on
h(2) = plot(years,advances);
    hold on
h(3) = plot(years,interest);
    hold on
h(4) = plot(years,premiums);
    hold on
h(5) = plot(years,house_expected);

set(h(:), 'LineWidth', 2);
    
legend(h, {'Outstanding Loan Balance',...
                'Cumulated Cash Advances',...
                'Cumulated Interest',...
                'Cumulated Premiums',...                
                'Expected Value of Collateral'},...                
                'Location','NorthWest');                    
     xlabel('Years after Origination'); 
       

end