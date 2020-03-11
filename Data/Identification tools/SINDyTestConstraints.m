clc
clear all

#change the file to load here
function dataOut = LoadPreparedData()
  load MeasuredData/Pendulum_Spring_Diff_L2_K10_h0_0001.txt;
  dataOut = Pendulum_Spring_Diff_L2_K10_h0_0001;
  clear Pendulum_Spring_Diff_L2_K10_h0_0001;
endfunction



function Xi = SparseRegression(Theta, dXdt, lambda)
  n = columns(dXdt);
  %% compute Sparse regression: sequential least squares
  Xi = Theta\dXdt; % initial guess: Least-squares
  % lambda is our sparsification knob.
  for k=1:10
    smallinds = (abs(Xi)<lambda); % find small coefficients
    Xi(smallinds)=0; % and threshold
    for ind = 1:n % n is state dimension
      biginds = smallinds(:,ind);
      % Regress dynamics onto remaining terms to find sparse Xi
      Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
    endfor
  endfor
endfunction

function matrix = DeleteSmaller(A, minimalNumber)
  matrix = A;
  for n = 1:rows(A)
    for m = 1:columns(A)
      if abs(matrix(n,m)) < minimalNumber
        matrix(n,m) = 0;
      endif
    endfor
  endfor
endfunction

function [dataOut, titlesOut] = Monomes3(data, titles)
  dataOut = data;
  titlesOut = [0; 0; 0];
  c = columns(data);
  
  for i =1:c
    titlesOut(:,i) = [titles(1,i); 0; 0];
  endfor
  column = c+1;
  
  %degree 2
  for i = 1:c
    for j = i:c
      dataOut(:,column) = data(:,i) .* data(:,j);
      titlesOut(:,column) = [0; 0; 0];
      titlesOut(1,column) = [titles(1,j)];
      titlesOut(2,column) = [titles(1,i)];
      column = column + 1;  
    endfor
  endfor
  
  cOut = columns(dataOut);
  %degree 3
  for i = 1:c
    ende = 1;
    for k = 1:cOut
      if titlesOut(1,k) == i && titlesOut(2,k) == i
        ende = k;
        break
      endif
    endfor
    for j = ende:cOut
      dataOut(:,column) = data(:,i) .* dataOut(:,j);
      titlesOut(:,column) = [0; 0; 0];
      titlesOut(1,column) = titlesOut(1,j);
      titlesOut(2,column) = titlesOut(2,j);
      titlesOut(3,column) = titles(1,i);
      
      column = column + 1;  
    endfor
  endfor
  
endfunction

function error = GetError(theta, Xi, responseVector)
  error = 100000;
  Y = theta * Xi;
  Y = responseVector - Y;
  Y = Y.*Y;
  r = rows(Y);
  c= columns(Y);
  summe = 0;
  for i = 1:r
    summe = summe + Y(i);
  endfor    
  error = sqrt(summe) * 1 / sqrt(rows(theta));
endfunction

function dataOut = OnlyDataOfInterest(data)
  %the only needed values are 5, 6, 7, 8, 13, 39, 43, 54, 58
  dataOut = zeros(rows(data),4);
  #dataOut(:,1) = data(:,5);
  #dataOut(:,2) = data(:,6);
  dataOut(:,1) = data(:,7);
  dataOut(:,2) = data(:,8);
  #dataOut(:,5) = data(:,13);
  dataOut(:,3) = data(:,39);
  dataOut(:,4) = data(:,43);
  #dataOut(:,8) = data(:,54);
  #dataOut(:,9) = data(:,58);
  
  
endfunction

function [dataOut, l] = NormalizeColumns(data, lNorm)
  dataOut = zeros(rows(data), columns(data));
  l = zeros(columns(data),1);
  if lNorm == 2
    for i = 1:columns(data)
      l(i,1) = sqrt(sum(data(:,i) .* data(:,i)));
      dataOut(:,i) = data(:,i) / l(i,1);
    endfor
  endif
  if lNorm == 1
    for i = 1:columns(data)
      l(i,1) = sum(abs(data(:,i)));
      dataOut(:,i) = data(:,i) / l(i,1);
    endfor
  endif
  if lNorm == 1000
    for i = 1:columns(data)
      l(i,1) = max(abs(data(:,i)));
      dataOut(:,i) = data(:,i) / l(i,1);
    endfor
  endif
endfunction

dataOut = LoadPreparedData();

[monomes, titlesOut] = Monomes3(dataOut, [1, 2, 3, 4, 5, 6]);

responseVector = -dataOut(:,5);
theta= dataOut = OnlyDataOfInterest(monomes); 
lambda = 0;  
  
Xi = SparseRegression(theta, responseVector, lambda);
