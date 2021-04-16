load 'result.mat'
data = readtable('data_Performance.csv');
data_fixed = readtable('data_fixed.csv');

% holdgamma = zeros(size(gamma, 1) * size(gamma, 2), size(gamma, 3));
% count = 1;
% for i = 1:size(gamma, 1)
%     for j = 1:size(gamma, 2)
%         holdgamma(count, :) = gamma(i, j, :);
%         count = count + 1;
%     end
% end

x_age = (20:40)';
x_agestar = x_age - 24.8417;
data_pred = [ones(length(x_age), 1) x_agestar x_agestar.^2 x_agestar.^3 x_agestar.^4 zeros(length(x_age), 5)];

f = data_pred * gamma';

figure(1)
plot(x_age, median(f, 2) - max(median(f, 2)), 'b-')
hold on
x = quantile(f, [0.025 0.975], 2);
low = x(:, 1);
upp = x(:, 2);
plot(x_age, upp - max(median(f, 2)), 'b--')
plot(x_age, low - max(median(f, 2)), 'b--')

data_pred = [zeros(length(x_age), 5) ones(length(x_age), 1) x_agestar x_agestar.^2 x_agestar.^3 x_agestar.^4];

f = data_pred * gamma';

plot(x_age, median(f, 2) - max(median(f, 2)),'r-')
x = quantile(f, [0.025 0.975], 2);
low = x(:, 1);
upp = x(:, 2);
plot(x_age, upp - max(median(f, 2)), 'r--')
plot(x_age, low - max(median(f, 2)), 'r--')
hold off

xlabel 'age'
ylabel 'performance'

print -dpdf 'SP_male_female.pdf'

data_pred = [ones(length(x_age), 1) x_agestar x_agestar.^2 x_agestar.^3 x_agestar.^4 zeros(length(x_age), 5)];
f = data_pred * gamma';

%figure(1)
%plot(data.Age(data_fixed.Var1 == 1), data.Performance(data_fixed.Var1 == 1), 'b.');
%hold on
plot(x_age, median(f, 2), 'r-')
hold on
x = quantile(f, [0.025 0.975], 2);
low = x(:, 1);
upp = x(:, 2);
plot(x_age, upp, 'r--')
plot(x_age, low, 'r--')
hold off
print -dpdf 'SP_male.pdf'


data_pred = [zeros(length(x_age), 5) ones(length(x_age), 1) x_agestar x_agestar.^2 x_agestar.^3 x_agestar.^4];
f = data_pred * gamma';

%plot(data.Age(data_fixed.Var1 == 0), data.Performance(data_fixed.Var1 == 0), 'b.');
%hold on
plot(x_age, median(f, 2),'r-')
hold on
x = quantile(f, [0.025 0.975], 2);
low = x(:, 1);
upp = x(:, 2);
plot(x_age, upp, 'r--')
plot(x_age, low, 'r--')
hold off

xlabel 'age'
ylabel 'performance'

print -dpdf 'SP_female.pdf'