clc
clear all

function matrix = LoadSINDyDoublePendulum()
  load MeasuredData/DoublePendulum_908770.txt;
  matrix = DoublePendulum_908770;
  clear DoublePendulum_908770;
endfunction

function value = factorial(x)
  value = 1;
  for k = 1:x
      value = k * value;
  end
endfunction

function value = BinomialCoefficient(n, k)
  value = factorial(n) / (factorial(k) * factorial(n-k));
endfunction

function matrix = AddOnes(X)
  X = X';
  r = rows(X);
  c = columns(X);
  matrix = zeros(r+1,c);
  matrix(1:r,:) = X;
  one = ones(1, c);
  matrix(r+1,:) = one;
  matrix = matrix';
endfunction

function matrix = AddMonome(A, degree) % degree 1 till 3
  A = A';
  r = rows(A);
  c = columns(A);
  if degree == 1
    matrix = A;
  endif
  
  if degree == 2
    matrix = zeros(BinomialCoefficient(r+degree-1,degree) + r,c);
    matrix(1:r, :) = A; %copy matrix A to the beginning
    matrix(r+1:r+r, :) = A .* A; %copy A² (element wise) to the matrix
    
    count = r * 2 + 1; 
    for n = 1:r
      for m = n + 1:r
        matrix(count, :) = A(n, :) .* A(m, :);
        count = count +1;
      endfor
    endfor
  endif
  
  if degree == 3
    matrix = zeros(BinomialCoefficient(r+1,2) + BinomialCoefficient(r+2,3) + r,c);
    matrix(1:r, :) = A; %copy matrix A to the beginning (degree 1)
    matrix(r+1:r+r, :) = A .* A; %copy A² (element wise) to the matrix (degree 2)
    matrix(2*r+1:2*r+r, :) = A .* A .* A; %copy A³ (element wise) to the matrix (degree 3)
    
    %monome degree 2  Other * Others
    count = 3 * r + 1; 
    for n = 1:r
      for m = n + 1:r
        matrix(count, :) = A(n, :) .* A(m, :);
        count = count +1;
      endfor
    endfor
    
    %monome degree 3  
    %first use the A² * others to build degree 3
    for n = 1:r
      for m = 1:r
        if n != m
          matrix(count, :) = matrix(r + n, :) .* A(m, :);
          count = count +1;
        endif
      endfor
    endfor
    
    startOthersOthersOthers = count;
   
    %second use the rest combinations
    for n = 1:r
      for m = n + 1:r
        for o = m + 1:r
          matrix(count, :) = A(n, :) .* A(m, :) .* A(o, :);
          count = count +1;
        endfor
      endfor
    endfor
  endif  
  matrix = matrix';
endfunction


function matrix = AddLags(A, lags)
  A = A';
  r = rows(A);
  c = columns(A);
  cNew = c - lags;
  matrix = zeros(r * (lags + 1), cNew);
  
  for n = 1:lags+1
    matRowStart = (n-1)*r+1;
    matRowEnd = (n-1)*r+r;
    AColumnStart = lags + 2 - n;
    AColumnEnd = lags + 1 - n + cNew;
    matrix(matRowStart:matRowEnd,:) = A(1:r, AColumnStart:AColumnEnd);
  endfor
  matrix = matrix';
endfunction 

function matrix = AddDifference(A, h, row)
  A = A';
  r = rows(A);
  c = columns(A);
  matrix = zeros(r+1,c-1);
  matrix(1:r,:) = A(:,1:c-1);
  for n = 1: columns(row)
    matrix(r+n,1:c-1) = (A(row(n),2:c) - A(row(n),1:c-1)) * (1 / h);  
  endfor  
  matrix = matrix';
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


function matrix = DeleteRows(X, del)
  matrix = X;
  for m = 1:columns(del)
    index = 0;
    max = 0;    
    for n = 1:columns(del)
      if del(1,n) > max
        max = del(1,n);
        index = n;
      endif
    endfor
    if index > 0 && max > 0
      matrix(max,:) = [];
      del(1,index) = 0;
    endif
  endfor
endfunction

function X = AddMonomeLagsOnes(X, monomeDegree, lags, activeOnes)
  X = AddMonome(X, monomeDegree);
  X = AddLags(X, lags);
  if activeOnes != 0
    X = AddOnes(X);
  endif
endfunction

function [matrix] = IteratedSindyReconstruction(observables, monomeDegree, lags, activeOnes, Xi, iterations, h)
  matrix = zeros(iterations,columns(Xi));
  
  for n = 1:lags+2
    matrix(n,:) = observables(n,1:columns(matrix));
  endfor
  
  for n = 2+lags:iterations
    theta = AddMonomeLagsOnes(matrix(n-lags-1:n-1,:), monomeDegree, lags,activeOnes);
    derivative = theta * Xi;
    matrix(n,:) = matrix(n-1,:) + derivative * h;
  endfor
endfunction

function positions = AngleToPosition(values)
  positions = zeros(rows(values),4);

  positions(:,1) = sin(values(:,1));
  positions(:,2) = -cos(values(:,1));

  positions(:,3) = sin(values(:,2))+ positions(:,1);
  positions(:,4) = -cos(values(:,2))+ positions(:,2);

endfunction

lags = 0;
monomeDegree = 3;
activeOnes = 1;
 
values = LoadSINDyDoublePendulum();

figure(1);
positions = AngleToPosition(values);

plot(positions(1:2100,1), positions(1:2100,2), '-', positions(1:2100,3), positions(1:2100,4));

figure(2);
plot(values(:,1), values(:,3));
xlabel("Theta1");
ylabel("d Theta1/dt");
figure(3);
plot(values(:,2), values(:,4));
xlabel("Theta2");
ylabel("d Theta2/dt");

Theta = AddMonomeLagsOnes(values, monomeDegree, lags, activeOnes);
Theta = DeleteRows(Theta, rows(Theta));
dXdt = AddDifference(values,0.01,1:4);
dXdt = dXdt(:,5:8);
dXdt = DeleteRows(dXdt, rows(dXdt) - lags+1:rows(dXdt));

Xi = SparseRegression(Theta, dXdt, 100);

[matrix] = IteratedSindyReconstruction(values, monomeDegree, lags, activeOnes, Xi, 2100, 0.01);
positions2 = AngleToPosition(matrix);
figure(4);
plot(positions2(:,1), positions2(:,2), '-', positions2(:,3), positions2(:,4));

distanceM1toM2 = sqrt((positions2(:,1) - positions2(:,3)) .* (positions2(:,1) - positions2(:,3)) + (positions2(:,2) - positions2(:,4)) .* (positions2(:,2) - positions2(:,4)));

figure(5);
plot(1:rows(distanceM1toM2), distanceM1toM2);
