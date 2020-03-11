clc
clear all

load MeasuredData/Test.txt;
data = Test;
clear Test;

term1 = data(:,1) .* data(:,1) .* data(:,5);
term2 = data(:,2) .* data(:,2) .* data(:,5);

sum = +term1 + term2;
ax = data(:,5);

figure(1)
plot(1:rows(data),ax,'k--', 1:rows(data), term1, 'b', 1:rows(data), term2, 'g', 2:rows(data), sum(2:end), 'r--');
legend("a" , "x²a", "y²a", "x²a+y²a");
xlabel("number of state sample")

