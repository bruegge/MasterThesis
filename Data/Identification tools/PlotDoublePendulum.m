clc
clear all

function dataOut = ModAngle(data)
  dataOut = zeros(rows(data),columns(data));

  Pi = 3.1415926;
  
  for j = 1:columns(data)
    for i = 1:rows(data)
      dataOut(i,j) = data(i,j);      
      while dataOut(i,j) > Pi
        dataOut(i,j) = dataOut(i,j) - 2*Pi;   
      endwhile
      while dataOut(i,j) < -Pi
        dataOut(i,j) = dataOut(i,j) + 2*Pi;
      endwhile
    endfor
  endfor
endfunction


load MeasuredData/DoublePendulum_Plot.txt;
data = DoublePendulum_Plot;
clear DoublePendulum_Plot;

load MeasuredData/DoublePendulum_Plot2.txt;
data2 = DoublePendulum_Plot2;
clear DoublePendulum_Plot2;



dataCarthesian = zeros(rows(data),4);

dataCarthesian(:,1) = sin(data(:,1));
dataCarthesian(:,2) = -cos(data(:,1));
dataCarthesian(:,3) = sin(data(:,2)) + dataCarthesian(:,1);
dataCarthesian(:,4) = -cos(data(:,2)) + dataCarthesian(:,2);

dataCarthesian2 = zeros(rows(data),4);

dataCarthesian2(:,1) = sin(data2(:,1));
dataCarthesian2(:,2) = -cos(data2(:,1));
dataCarthesian2(:,3) = sin(data2(:,2)) + dataCarthesian2(:,1);
dataCarthesian2(:,4) = -cos(data2(:,2)) + dataCarthesian2(:,2);


dataOut1 = ModAngle(data(:,1:2));
dataOut2 = ModAngle(data2(:,1:2));

dataCarthesian = dataCarthesian(1:end/4,:);
dataCarthesian2 = dataCarthesian2(1:end/4,:);



figure(1)
plot(dataCarthesian(:,1),dataCarthesian(:,2),'g-',dataCarthesian2(:,1),dataCarthesian2(:,2),'k--',dataCarthesian(:,3),dataCarthesian(:,4),'b',dataCarthesian2(:,3),dataCarthesian2(:,4),'r')
#plot(1:rows(dataOut1),dataOut1(:,1),1:rows(dataOut2), dataOut2(:,1));
legend("inner bob" , "double pendulum 1", "double pendulum 2");
xlabel("x")
ylabel("y")

figure(2)
plot(dataCarthesian(:,1),dataCarthesian(:,2),'g',dataCarthesian(:,3),dataCarthesian(:,4),'b',dataCarthesian2(:,3),dataCarthesian2(:,4),'r')
#plot(1:rows(dataOut1),dataOut1(:,2),1:rows(dataOut2), dataOut2(:,2));
legend("inner bob" , "double pendulum 1", "double pendulum 2");
xlabel("x")
ylabel("y")

