
rm = cases();
rm = rm.compute;

n = 258051;

% collateral values
collateral = zeros(n,1);
for i=1:n, collateral(i)=rm.data{i,31}; end
save collateral.mat
clear collateral

% age
age = zeros(n,1);
for i=1:n, age(i)=rm.data{i,44}; end
save age.mat
clear age

% gender
gender = cell(n,1);
  for i=1:n, gender{i} = char(cellstr(rm.data{i,45})); end
save gender.mat
clear gender




