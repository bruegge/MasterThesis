clc
clear all

function matrix = LoadSINDyRigidBody1()
  load MeasuredData/RigidBody3Points_C_S1_P1.txt;
  matrix = RigidBody3Points_C_S1_P1;
  clear RigidBody3Points_C_S1_P1;
endfunction

function matrix = LoadSINDyRigidBody2()
  load MeasuredData/RigidBody3Points_C_S3_P3.txt;
  matrix = RigidBody3Points_C_S3_P3;
  clear RigidBody3Points_C_S3_P3;
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

function ShowTrajectorySindy(X, dimension)
  if dimension == 1
    plot(X(:,1));
  endif
  if dimension == 2
    plot(X(:,1),X(:,2));
  endif
  if dimension == 3
    plot3(X(:,1),X(:,2), X(:,3));  
  endif
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

function matrix = DeleteColumns(X, del)
  X = X';
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
  matrix = matrix';
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

function [observables, matrix] = SindyRigidBody(particlenumber)
  numberLags = 1;

  observables = LoadSINDyRigidBody1();
  if particlenumber == 1
    observables = observables(:,1:3);
  endif
  if particlenumber == 2
    observables = observables(:,4:6);
  endif
  if particlenumber == 3
    observables = observables(:,7:9);
  endif

  theta = observables;
  dXdt = AddDifference(theta, 0.01, 1:3);
  dXdt = DeleteColumns(dXdt, 1:3);
  
  dXdt = DeleteRows(dXdt, 1:numberLags);

  theta = AddMonome(theta, 3);
  theta = AddLags(theta,numberLags);
  theta = AddOnes(theta);
  theta = DeleteRows(theta, rows(theta));

  Xi = SparseRegression(theta, dXdt, 0);
  %Xi = DeleteSmaller(Xi, 0.0001);

  matrix = IteratedSindyReconstruction(observables, 3, numberLags, 1, Xi, 1000, 0.01);
endfunction

[observables1, matrix1] = SindyRigidBody(1);
[observables2, matrix2] = SindyRigidBody(2);
[observables3, matrix3] = SindyRigidBody(3);

figure(1)
plot3(observables1(:,1),observables1(:,2), observables1(:,3),'-');

figure(2)
plot3(matrix1(:,1),matrix1(:,2), matrix1(:,3),'--');

diff = observables1(1:1000,:) - matrix1;

figure(3)
plot(1:rows(diff), diff(:,1),'-');

figure(4)
plot(1:rows(diff), diff(:,2),'-');

figure(5)
plot(1:rows(diff), diff(:,3),'-');
