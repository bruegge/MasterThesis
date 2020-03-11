clc
clear all

function matrix = LoadSINDySingleParticle()
  load MeasuredData/SinglePoint_h0_01_ax5_ay20.txt;
  matrix = SinglePoint_h0_01_ax5_ay20;
  clear SinglePoint_h0_01_ax5_ay20;
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

function X = AddMonomeLagsOnes(X, monomeDegree, lags, activeOnes, derivative,h)
  X = AddMonome(X, monomeDegree);
  X = AddLags(X, lags);
  if activeOnes != 0
    X = AddOnes(X);
  endif
  if derivative != 0
    X = AddDifference(X,h,derivative);
  endif
endfunction

function [matrix] = IteratedSindyReconstruction(observables, monomeDegree, lags, activeOnes, derivatives, Xi, iterations, h)
  matrix = zeros(iterations,columns(Xi));
  
  for n = 1:lags+2
    matrix(n,:) = observables(n,1:columns(matrix));
  endfor
  for n = 2+lags:iterations
    theta = AddMonomeLagsOnes(matrix(n-lags-1:n-1,:), monomeDegree, lags,activeOnes,derivatives, h);
    derivative = theta * Xi;
    matrix(n,:) = matrix(n-1,:) + derivative * h;
  endfor
endfunction

function [Xi,matrix, values] = TestBothValuesTogether()
  lags = 0;
  monomeDegree = 2;
  activeOnes = 1;
  
  values = LoadSINDySingleParticle();

  Theta = AddMonomeLagsOnes(values, monomeDegree, lags, activeOnes,0,0);
  Theta = DeleteRows(Theta, rows(Theta));
  dXdt = AddDifference(values,0.01,1:4);
  dXdt = dXdt(:,5:8);
  dXdt = DeleteRows(dXdt, rows(dXdt) - lags+1:rows(dXdt));

  Xi = SparseRegression(Theta, dXdt, 100);

  [matrix] = IteratedSindyReconstruction(values, monomeDegree, lags, activeOnes,0, Xi, 1000, 0.01);

  figure(1)
  plot(1:1000,matrix(:,1),'-', 1:1000, values(1:1000,1), '--');
 
  figure(2)
  plot(1:1000,matrix(:,2),'-', 1:1000, values(1:1000,2), '--');
  
endfunction

function [Xi,matrix, values] = TestfirstValue()
  lags = 0;
  monomeDegree = 2;
  activeOnes = 1;
  values = LoadSINDySingleParticle();

  values = values(:,1);
  Theta = AddMonomeLagsOnes(values, monomeDegree, lags, activeOnes,0,0);
  Theta = DeleteRows(Theta, rows(Theta));
  dXdt = AddDifference(values,0.01,1);
  dXdt = dXdt(:,2);
  dXdt = DeleteRows(dXdt, rows(dXdt) - lags+1:rows(dXdt));

  Xi = SparseRegression(Theta, dXdt, 100);

  [matrix] = IteratedSindyReconstruction(values, monomeDegree, lags, activeOnes,0, Xi, 1000, 0.01);

  figure(3)
  plot(1:1000,matrix(:,1),'-', 1:1000, values(1:1000,1), '--');
  
endfunction

function [Xi,matrix, values,dXdt] = TestSecondValue()
  lags = 0;
  monomeDegree = 2;
  activeOnes = 1;
 
  values = LoadSINDySingleParticle();
  v = values;
  values = values(:,[2,4]);
  Theta = AddMonomeLagsOnes(values, monomeDegree, lags, activeOnes,0,0.01);
  Theta = DeleteRows(Theta, rows(Theta));
  dXdt = AddDifference(values,0.01,1:2);
  dXdt = dXdt(:,3:4);
  dXdt = DeleteRows(dXdt, rows(dXdt) - lags+1:rows(dXdt));

  Xi = SparseRegression(Theta, dXdt, 0.001);

  [matrix] = IteratedSindyReconstruction(values, monomeDegree, lags, activeOnes, 0, Xi, 1000, 0.01);

  figure(4)
  plot(1:1000,matrix(:,1),'-', 1:1000, values(1:1000,1), '--'); 
endfunction

#[Xi,matrix,values] = TestBothValuesTogether();
#[Xi2,matrix2, values2] = TestfirstValue();
#[Xi3,matrix3, values3,dXdt] = TestSecondValue();


values = LoadSINDySingleParticle();
figure(1)
plot(values(:,1),values(:,2));

dXdt = AddDifference(values,0.01,1:4);
dXdt(:,1:6) = [];

#Theta = values(:,3:4);
Theta(:,1) = ones(rows(values),1);

Theta = Theta(1:end-1,:);
Xi = SparseRegression(Theta, dXdt, 0.001);
