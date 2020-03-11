clc
clear all

load MeasuredData/LorenzAttractor_T10_h0_001.txt;
T100 = LorenzAttractor_T10_h0_001;
clear LorenzAttractor_T10_h0_001;

load MeasuredData/LorenzAttractor_T100_h0_0001.txt;
T10 = LorenzAttractor_T100_h0_0001;
clear LorenzAttractor_T100_h0_0001;
 
r = rows(T10);  

#cut
T10 = T10(1:3000,:);
T100 = T100(1:3000,:);


figure(1)
plot3(T10(:,1),T10(:,2), T10(:,3),'--',T100(:,1),T100(:,2), T100(:,3),'--');
legend("x = 1, y = 1, z = 1", "x = 1.01, y = 1, z = 1");
xlabel("x")
ylabel("y")
zlabel("z")

