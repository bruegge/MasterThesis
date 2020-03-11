clc
clear all

load MeasuredData/pendulumSpring_h0_00001_T1000_k10.txt;
k10 = pendulumSpring_h0_00001_T1000_k10;
clear pendulumSpring_h0_00001_T1000_k10;

load MeasuredData/pendulumSpring_h0_00001_T1000_k100.txt;
k100 = pendulumSpring_h0_00001_T1000_k100;
clear pendulumSpring_h0_00001_T1000_k100;

load MeasuredData/pendulumSpring_h0_00001_T1000_k1000.txt;
k1000 = pendulumSpring_h0_00001_T1000_k1000;
clear pendulumSpring_h0_00001_T1000_k1000;

load MeasuredData/pendulumSpring_h0_00001_T1000_k10000.txt;
k10000 = pendulumSpring_h0_00001_T1000_k10000;
clear pendulumSpring_h0_00001_T1000_k10000;

load MeasuredData/pendulumSpring_h0_00001_T1000_k100000.txt;
k100000 = pendulumSpring_h0_00001_T1000_k100000;
clear pendulumSpring_h0_00001_T1000_k100000;

load MeasuredData/pendulumSpring_h0_00001_T1000_k1000000.txt;
k1000000 = pendulumSpring_h0_00001_T1000_k1000000;
clear pendulumSpring_h0_00001_T1000_k1000000;
  
r = rows(k10);  

#cut
k10 = k10(1:580,:);
k100 = k100(1:580,:);
k1000 = k1000(1:580,:);
k10000 = k10000(1:580,:);

unitCircle = zeros(360,2);
for i = 1:360
  unitCircle(i,1) = cos(i/180*3.1415926);
  unitCircle(i,2) = sin(i/180*3.1415926);
  
endfor

figure(1)
plot(unitCircle(:,1), unitCircle(:,2),'-', k10(:,1),k10(:,2),'-',k100(:,1),k100(:,2),'-',k1000(:,1),k1000(:,2),'-',k10000(:,1),k10000(:,2),'-');
legend("unit circle" , "k = 10", "k = 100", "k = 1000", "k = 10000");
xlabel("x")
ylabel("y")

