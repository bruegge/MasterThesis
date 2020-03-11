clc
clear all

load MeasuredData/Pendulum_Spring_Diff_L2_K100000_h0_000001.txt;
dataSpring = Pendulum_Spring_Diff_L2_K100000_h0_000001;
clear Pendulum_Spring_Diff_L2_K100000_h0_000001;

load MeasuredData/TestL2.txt;
dataAngular = TestL2;
clear TestL2;
  

figure(1)
plot(1:rows(dataAngular), dataAngular(:,5), 1:rows(dataSpring),dataSpring(:,5));
legend("unit circle" , "k = 10", "k = 100", "k = 1000", "k = 10000");
xlabel("x")
ylabel("y")

