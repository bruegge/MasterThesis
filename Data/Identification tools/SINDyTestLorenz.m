clc
clear all

function matrix = LoadSINDyLorenz()
  load MeasuredData/LorenzAttractor_251924.txt;
  matrix = LorenzAttractor_251924;
  clear LorenzAttractor_251924;
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
    matrix(r+1:r+r, :) = A .* A; %copy A� (element wise) to the matrix
    
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
    matrix(r+1:r+r, :) = A .* A; %copy A� (element wise) to the matrix (degree 2)
    matrix(2*r+1:2*r+r, :) = A .* A .* A; %copy A� (element wise) to the matrix (degree 3)
    
    %monome degree 2  Other * Others
    count = 3 * r + 1; 
    for n = 1:r
      for m = n + 1:r
        matrix(count, :) = A(n, :) .* A(m, :);
        count = count +1;
      endfor
    endfor
    
    %monome degree 3  
    %first use the A� * others to build degree 3
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

function X = AddMonomeLags(X, monomeDegree, lags, activateOnes)
  X = AddLags(X, lags);
  X = AddMonome(X, monomeDegree);
  if activateOnes == 1
    X = AddOnes(X);
  endif
endfunction

function [matrix] = IteratedSindyReconstruction(observables, monomeDegree, lags, activateOnes, Xi, iterations, h)
  matrix = zeros(iterations,columns(Xi));
  
  for n = 1:lags+2
    matrix(n,:) = observables(n,1:columns(matrix));
  endfor
  
  for n = 2+lags:iterations
    theta = AddMonomeLags(matrix(n-lags-1:n-1,:), monomeDegree, lags, activateOnes);
    derivative = theta * Xi;
    matrix(n,:) = matrix(n-1,:) + derivative * h;
  endfor
endfunction

function [observables, matrix, Xi] = SindyLorenz1Observablen(observableNumber)
  numberLags = 4;
  degreeMonomes = 3;
  activateOnes = 1;
  
  observables = LoadSINDyLorenz();
  theta = observables;
  dXdt = AddDifference(theta, 0.01, 1:3);

  dXdt = DeleteColumns(dXdt, 1:3);
  if observableNumber == 1
    dXdt = DeleteColumns(dXdt, 2:3); 
  endif
  if observableNumber == 2
    dXdt = DeleteColumns(dXdt, 3);
    dXdt = DeleteColumns(dXdt, 1); 
  endif
  if observableNumber == 3
    dXdt = DeleteColumns(dXdt, 1:2);
  endif
  
  dXdt = DeleteRows(dXdt, 1:numberLags);

  if observableNumber == 1
    theta = DeleteColumns(theta,2:3);
  endif
  if observableNumber == 2
    theta = DeleteColumns(theta, 3);
    theta = DeleteColumns(theta, 1);
  endif
  if observableNumber == 3
    theta = DeleteColumns(theta,1:2);
  endif
  
  theta = AddLags(theta,numberLags);
  theta = AddMonome(theta, degreeMonomes);
  if activateOnes == 1
    theta = AddOnes(theta);
  endif
  theta = DeleteRows(theta, rows(theta));

  Xi = SparseRegression(theta, dXdt, 1);
  %Xi = DeleteSmaller(Xi, 0.0001);
  matrix = 1;
  %[matrix]  = IteratedSindyReconstruction(observables, degreeMonomes, numberLags, activateOnes, Xi, 1000, 0.01);
endfunction

function [observables, Xi, theta, dXdt] = SindyLorenz3Observablen()  
  numberLags = 0;
  degreeMonomes = 2;
  activateOnes = 1;
  
  observables = LoadSINDyLorenz();
  theta = observables;
  dXdt = AddDifference(theta, 0.01, 1:3);

  dXdt = DeleteColumns(dXdt, 1:3);
  %dXdt = DeleteColumns(dXdt, 2:3);
  dXdt = DeleteRows(dXdt, 1:numberLags);
  dXdt = dXdt(:,3);
  
  %theta = DeleteColumns(theta,2:3);
  theta = AddLags(theta,numberLags);
  theta = AddMonome(theta, degreeMonomes);
  if activateOnes == 1
    theta = AddOnes(theta);
  endif  
  
  theta = DeleteRows(theta, rows(theta));
  theta = theta(:,[3 7 10]);
  
  Xi = SparseRegression(theta, dXdt, 1000);
  %Xi = DeleteSmaller(Xi, 0.0001);

 # [matrix]  = IteratedSindyReconstruction(observables, degreeMonomes, numberLags, activateOnes, Xi, 20000, 0.01);

 # figure(1)
#  plot(observables(1:2000,1),'-',matrix(1:2000,1),'--');
#  figure(2)
 # plot3(observables(1:2000,1),observables(1:2000,2), observables(1:2000,3), '-',matrix(1:2000,1),matrix(1:2000,2), matrix(1:2000,3),'--'); 
#  figure(3)
#  plot3(matrix(:,1),matrix(:,2), matrix(:,3),'-');   
endfunction

%SindyLorenz3Observablen();
[observables1, Xi, theta, dXdt] = SindyLorenz3Observablen();
Xis = DeleteSmaller(Xi, 0.0001);
%[observables2, matrix2] = SindyLorenz1Observablen(2);
%[observables3, matrix3] = SindyLorenz1Observablen(3);

