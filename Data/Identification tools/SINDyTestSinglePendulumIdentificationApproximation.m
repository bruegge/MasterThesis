clc
clear all

function dataOut = Smooth(data, windowSize)
  dataOut = zeros(rows(data)-windowSize+1,1);
  for i = 1:rows(dataOut)
    
    result = 0;
    for j = 1:windowSize
      result = result + data(i+j-1);
    endfor
    dataOut(i) = result / windowSize;
    
  endfor
endfunction



function u = TVRegDiff( data, iter, alph, u0, scale, ep, dx, plotflag, diagflag )
  % u = TVRegDiff( data, iter, alph, u0, scale, ep, dx, plotflag, diagflag );
  % u = tvdiff( e, dx, iter, ep, alph );
  % Rick Chartrand (rickc@lanl.gov), Apr. 10, 2011
  % Please cite Rick Chartrand, "Numerical differentiation of noisy,
  % nonsmooth data," ISRN Applied Mathematics, Vol. 2011, Article ID 164564, 
  % 2011. 
  %
  % Inputs:  (First three required; omitting the final N parameters for N < 7
  %           or passing in [] results in default values being used.) 
  %       data        Vector of data to be differentiated.
  %
  %       iter        Number of iterations to run the main loop.  A stopping
  %                   condition based on the norm of the gradient vector g
  %                   below would be an easy modification.  No default value.
  %
  %       alph        Regularization parameter.  This is the main parameter
  %                   to fiddle with.  Start by varying by orders of
  %                   magnitude until reasonable results are obtained.  A
  %                   value to the nearest power of 10 is usally adequate.
  %                   No default value.  Higher values increase
  %                   regularization strenght and improve conditioning.
  %
  %       u0          Initialization of the iteration.  Default value is the
  %                   naive derivative (without scaling), of appropriate
  %                   length (this being different for the two methods).
  %                   Although the solution is theoretically independent of
  %                   the intialization, a poor choice can exacerbate
  %                   conditioning issues when the linear system is solved.
  %
  %       scale       'large' or 'small' (case insensitive).  Default is
  %                   'small'.  'small' has somewhat better boundary
  %                   behavior, but becomes unwieldly for data larger than
  %                   1000 entries or so.  'large' has simpler numerics but
  %                   is more efficient for large-scale problems.  'large' is
  %                   more readily modified for higher-order derivatives,
  %                   since the implicit differentiation matrix is square.
  %
  %       ep          Parameter for avoiding division by zero.  Default value
  %                   is 1e-6.  Results should not be very sensitive to the
  %                   value.  Larger values improve conditioning and
  %                   therefore speed, while smaller values give more
  %                   accurate results with sharper jumps.
  %
  %       dx          Grid spacing, used in the definition of the derivative
  %                   operators.  Default is the reciprocal of the data size.
  %
  %       plotflag    Flag whether to display plot at each iteration.
  %                   Default is 1 (yes).  Useful, but adds significant
  %                   running time.
  %
  %       diagflag    Flag whether to display diagnostics at each
  %                   iteration.  Default is 1 (yes).  Useful for diagnosing
  %                   preconditioning problems.  When tolerance is not met,
  %                   an early iterate being best is more worrying than a
  %                   large relative residual.
  %                   
  % Output:
  %
  %       u           Estimate of the regularized derivative of data.  Due to
  %                   different grid assumptions, length( u ) = 
  %                   length( data ) + 1 if scale = 'small', otherwise
  %                   length( u ) = length( data ).

  %% Copyright notice:
  % Copyright 2010. Los Alamos National Security, LLC. This material
  % was produced under U.S. Government contract DE-AC52-06NA25396 for
  % Los Alamos National Laboratory, which is operated by Los Alamos
  % National Security, LLC, for the U.S. Department of Energy. The
  % Government is granted for, itself and others acting on its
  % behalf, a paid-up, nonexclusive, irrevocable worldwide license in
  % this material to reproduce, prepare derivative works, and perform
  % publicly and display publicly. Beginning five (5) years after
  % (March 31, 2011) permission to assert copyright was obtained,
  % subject to additional five-year worldwide renewals, the
  % Government is granted for itself and others acting on its behalf
  % a paid-up, nonexclusive, irrevocable worldwide license in this
  % material to reproduce, prepare derivative works, distribute
  % copies to the public, perform publicly and display publicly, and
  % to permit others to do so. NEITHER THE UNITED STATES NOR THE
  % UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL
  % SECURITY, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY,
  % EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
  % RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF
  % ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR
  % REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED
  % RIGHTS. 

  %% BSD License notice:
  % Redistribution and use in source and binary forms, with or without
  % modification, are permitted provided that the following conditions
  % are met: 
  % 
  %      Redistributions of source code must retain the above
  %      copyright notice, this list of conditions and the following
  %      disclaimer.  
  %      Redistributions in binary form must reproduce the above
  %      copyright notice, this list of conditions and the following
  %      disclaimer in the documentation and/or other materials
  %      provided with the distribution. 
  %      Neither the name of Los Alamos National Security nor the names of its
  %      contributors may be used to endorse or promote products
  %      derived from this software without specific prior written
  %      permission. 
  %  
  % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
  % CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
  % INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  % MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
  % CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  % SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  % LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
  % USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
  % AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  % LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  % ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  % POSSIBILITY OF SUCH DAMAGE. 

  %% code starts here
  % Make sure we have a column vector.
  data = data( : );
  % Get the data size.
  n = length( data )

  % Default checking. (u0 is done separately within each method.)
  if nargin < 9 || isempty( diagflag )
      diagflag = 1;
  end
  if nargin < 8 || isempty( plotflag )
      plotflag = 1;
  end
  if nargin < 7 || isempty( dx )
      dx = 1 / n;
  end
  if nargin < 6 || isempty( ep )
      ep = 1e-6;
  end
  if nargin < 5 || isempty( scale )
      scale = 'small';
  end

  
  % Different methods for small- and large-scale problems.
  switch lower( scale )
      
    case 'small'
        % Construct differentiation matrix.
        c = ones( n + 1, 1 ) / dx;
        D = spdiags( [ -c, c ], [ 0, 1 ], n, n + 1 );
        clear c
        DT = D';
        % Construct antidifferentiation operator and its adjoint.
        A = @(x) chop( cumsum( x ) - 0.5 * ( x + x( 1 ) ),10 ) * dx;
        
        AT = @(w,n) ( sum( w ) * ones( n + 1, 1 ) - [ sum( w ) / 2; cumsum( w ) - w / 2 ] ) * dx;
        % Default initialization is naive derivative.
        if nargin < 4 || isempty( u0 )
            u0 = [ 0; diff( data ); 0 ];
        end
        u = u0;
        % Since Au( 0 ) = 0, we need to adjust.
        ofst = data( 1 );
        % Precompute.
        ATb = AT( ofst - data,n );
        
        % Main loop.
        for ii = 1 : iter
            % Diagonal matrix of weights, for linearizing E-L equation.
            Q = spdiags( 1 ./ ( sqrt( ( D * u ).^2 + ep ) ), 0, n, n );
            % Linearized diffusion matrix, also approximation of Hessian.
            L = dx * DT * Q * D;
            % Gradient of functional.
        
            rATB = rows(ATb)
            cATB = columns(ATb)
            
            g = AT( A( u ) , n+1);
            g(rows(g)) = [];
            g = g + ATb + alph * L * u;
            % Prepare to solve linear equation.
            tol = 1e-4;
            maxit = 100;
            % Simple preconditioner.
            P = alph * spdiags( spdiags( L, 0 ) + 1, 0, n + 1, n + 1 );
            if diagflag
                s = pcg( @(v) ( alph * L * v + AT( A( v ),n+1 )(1:end-1) ), g, tol, maxit, P );
                fprintf( 'iteration %4d: relative change = %.3e, gradient norm = %.3e\n', ii, norm( s ) / norm( u ), norm( g ) );
            else
                [ s, ~ ] = pcg( @(v) ( alph * L * v + AT( A( v ) ) ), g, tol, maxit, P );
            end
            % Update solution.
            u = u - s;
            % Display plot.
            if plotflag
                plot( u, 'ok' ), drawnow;
            end
        end
        
    case 'large'
        % Construct antidifferentiation operator and its adjoint.
        A = @(v) cumsum(v);
        AT = @(w) ( sum(w) * ones( length( w ), 1 ) - [ 0; cumsum( w( 1 : end - 1 ) ) ] );
        % Construct differentiation matrix.
        c = ones( n, 1 );
        D = spdiags( [ -c c ], [ 0 1 ], n, n ) / dx;
        D( n, n ) = 0;
        clear c
        DT = D';
        % Since Au( 0 ) = 0, we need to adjust.
        data = data - data( 1 );
        % Default initialization is naive derivative.
        if nargin < 4 || isempty( u0 )
            u0 = [ 0; diff( data ) ];
        end
        u = u0;
        % Precompute.
        ATd = AT( data );
        
        % Main loop.
        for ii = 1 : iter
            % Diagonal matrix of weights, for linearizing E-L equation.
            Q = spdiags( 1./ sqrt( ( D * u ).^2 +  ep ), 0, n, n );
            % Linearized diffusion matrix, also approximation of Hessian.
            L = DT * Q * D;
            % Gradient of functional.
            g = AT( A( u ) ) - ATd;
            g= g + alph * L * u;
            % Build preconditioner.
            c = cumsum( n : -1 : 1 ).';
            B = alph * L + spdiags( c( end : -1 : 1 ), 0, n, n );
            droptol = 1.0e-2;
            R = cholinc( B, droptol );
            % Prepare to solve linear equation.
            tol = 1.0e-4;
            maxit = 100;
            if diagflag
                s = pcg( @(x) ( alph * L * x + AT( A( x ) ) ), -g, tol, maxit, R', R );
                fprintf( 'iteration %2d: relative change = %.3e, gradient norm = %.3e\n', ii, norm( s ) / norm( u ), norm( g ) );
            else
                [ s, ~ ] = pcg( @(x) ( alph * L * x + AT( A( x ) ) ), -g, tol, maxit, R', R );
            end
            % Update current solution
            u = u + s;
            % Display plot.
            if plotflag
                plot( u, 'ok' ), drawnow;
            end
        end
   end
end

function matrixOut = KeepEveryNthRow(matrix, n)
  matrixOut = zeros(1,columns(matrix));
  counter = 1;
  for i = 1:n:rows(matrix)
    matrixOut(counter,:) = matrix(i,:);
    counter = counter + 1;
  endfor
  
endfunction

function matrix = LoadSINDySinglePendulum()
  load MeasuredData/pendulumSpring_n_h0_0001_T1_k100000.txt;
  matrix = pendulumSpring_n_h0_0001_T1_k100000;
  clear pendulumSpring_n_h0_0001_T1_k100000;
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

function [observables, Xi, theta, dXdt] = SindyPendulumObservablen()  
  numberLags = 0;
  degreeMonomes = 1;
  activateOnes = 1;
  
  observables = LoadSINDySinglePendulum();
  #observables(:,7) = [];
  theta = observables;
  dXdt = AddDifference(theta, 0.01, 1:columns(theta));

  dXdt = DeleteColumns(dXdt, 1:columns(theta));
  %dXdt = DeleteColumns(dXdt, 2:3);
  dXdt = DeleteRows(dXdt, 1:numberLags);
  %dXdt = dXdt(:,3);
  
  %theta = DeleteColumns(theta,2:3);
  theta = AddMonome(theta, degreeMonomes);
  theta = AddLags(theta,numberLags);
  if activateOnes == 1
    theta = AddOnes(theta);
  endif  
  
  theta = DeleteRows(theta, rows(theta));
  %theta = theta(:,[3 7 10]);
  
  Xi = SparseRegression(theta, dXdt, 100);
  %Xi = DeleteSmaller(Xi, 0.0001);

 # [matrix]  = IteratedSindyReconstruction(observables, degreeMonomes, numberLags, activateOnes, Xi, 20000, 0.01);

 # figure(1)
#  plot(observables(1:2000,1),'-',matrix(1:2000,1),'--');
#  figure(2)
 # plot3(observables(1:2000,1),observables(1:2000,2), observables(1:2000,3), '-',matrix(1:2000,1),matrix(1:2000,2), matrix(1:2000,3),'--'); 
#  figure(3)
#  plot3(matrix(:,1),matrix(:,2), matrix(:,3),'-');   
endfunction

function [observables, Xi, theta, dXdt] = SindyTest1(keepEveryNthElement)  
  degreeMonomes = 1;
  activateOnes = 1;
  
  observables = LoadSINDySinglePendulum();
  observables = KeepEveryNthRow(observables,keepEveryNthElement);
  theta = observables;
  c = columns(theta);
  r = rows(theta);
  
  dXdt = AddDifference(theta, 0.0001*keepEveryNthElement, 1:columns(theta));

  dXdt = DeleteColumns(dXdt, 1:columns(theta));
  
  
  theta = AddMonome(theta, degreeMonomes);
  c = columns(theta);
  
  theta(:,c+1) = theta(:,1) ./ sqrt(theta(:,1) .* theta(:,1) + theta(:,2) .* theta(:,2));
  theta(:,c+2) = theta(:,2) ./ sqrt(theta(:,1) .* theta(:,1) + theta(:,2) .* theta(:,2));
  
  observables = theta;
  if activateOnes == 1
    theta = AddOnes(theta);
  endif  
  
  theta = DeleteRows(theta, rows(theta));
  %theta = theta(:,[3 7 10]);
  
  Xi = SparseRegression(theta, dXdt, 0);
  %Xi = DeleteSmaller(Xi, 0.0001);

 # [matrix]  = IteratedSindyReconstruction(observables, degreeMonomes, numberLags, activateOnes, Xi, 20000, 0.01);

 # figure(1)
#  plot(observables(1:2000,1),'-',matrix(1:2000,1),'--');
#  figure(2)
 # plot3(observables(1:2000,1),observables(1:2000,2), observables(1:2000,3), '-',matrix(1:2000,1),matrix(1:2000,2), matrix(1:2000,3),'--'); 
#  figure(3)
#  plot3(matrix(:,1),matrix(:,2), matrix(:,3),'-');   
endfunction

function [Xi] = SindyTest2(theta, responseVector, lambda)  
  
  Xi = SparseRegression(theta, responseVector, lambda);
  
  
endfunction


%SindyLorenz3Observablen();
%[observables, Xi, theta, dXdt] = SindyPendulumObservablen();
%Xis = DeleteSmaller(Xi, 0.0001);
%[observables2, matrix2] = SindyLorenz1Observablen(2);
%[observables3, matrix3] = SindyLorenz1Observablen(3);



%g = 9.807;

%Test = zeros(10000,2);

%Test(:,1) = g * observables(:,2) .* observables(:,1);
%Test(:,2) = -g * observables(:,1) .* observables(:,1);

function angle = GetAngle(x,y)
  pi = 3.1415926;
  angle = 0;
  if x == 0 
    if y < 0
      angle = 0;
    else
      angle = 2* pi;
    endif
    return
  else   
    if y == 0
      if x < 0
        angle = -pi/2;
      else
        angle = pi/2;
      endif
      return
    endif
  endif
  if x < 0 && y < 0
    angle = -atan(-x/-y);
  endif
  if x < 0 && y > 0
    angle = -atan(y/-x)- pi/2;
  endif
  if x > 0 && y < 0
    angle = atan(x/-y);
  endif
  if x > 0 && y > 0
    angle = atan(y/x)+pi/2;
  endif
endfunction


function angle = clampAngle(angle)
  pi = 3.1415926;
  while angle > pi
    angle = angle - pi * 2;
  endwhile
  while angle < -pi
    angle = angle + pi * 2;
  endwhile
  
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
  dataOut = zeros(rows(data),9);
  dataOut(:,1) = data(:,5);
  dataOut(:,2) = data(:,6);
  dataOut(:,3) = data(:,7);
  dataOut(:,4) = data(:,8);
  dataOut(:,5) = data(:,13);
  dataOut(:,6) = data(:,39);
  dataOut(:,7) = data(:,43);
  dataOut(:,8) = data(:,54);
  dataOut(:,9) = data(:,58);
  #dataOut(:,10) = data(:,32);
  #dataOut(:,11) = data(:,52);
  
  
endfunction

function dataOut = LoadPreparedData()
  
  load MeasuredData/Test.txt;
  dataOut = Test;
  clear TestL2;
  
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


load MeasuredData/AngularPendulum_Dif_L1.txt;
data = AngularPendulum_Dif_L1;
clear AngularPendulum_Dif_L1;

data = data(1:end,:);

theta = zeros(rows(data),3);
theta(:,1) = sin(data(:,1));
theta(:,2) = data(:,1);
theta(:,3) = data(:,1).*data(:,1);
theta(:,4) = data(:,1).*data(:,1).*data(:,1);
theta(:,5) = data(:,1).*data(:,1).*data(:,1).*data(:,1);
theta(:,6) = data(:,1).*data(:,1).*data(:,1).*data(:,1).*data(:,1);

responseVector = data(:,3);

Xi1 = SindyTest2(theta(1:end/5,1), responseVector(1:end/5,:), 0);  
Xi2 = SindyTest2(theta(1:end/5,1:3), responseVector(1:end/5,:), 0);  #r
Xi3 = SindyTest2(theta(1:end/5,1:4), responseVector(1:end/5,:), 0);  #b
Xi4 = SindyTest2(theta(1:end/5,1:5), responseVector(1:end/5,:), 0);  

Result1 = zeros(rows(data),1);
Result2 = zeros(rows(data),1);
Result3 = zeros(rows(data),1);
Result4 = zeros(rows(data),1);
for i = 1:rows(Xi1)
  Result1 = Result1 + Xi1(i) .* theta(:,i);
endfor
for i = 1:rows(Xi2)
  Result2 = Result2 + Xi2(i) .* theta(:,i);
endfor
for i = 1:rows(Xi3)
  Result3 = Result3 + Xi3(i) .* theta(:,i);
endfor
for i = 1:rows(Xi4)
  Result4 = Result4 + Xi4(i) .* theta(:,i);
endfor


figure(1);
plot(1:rows(data)/5,responseVector(1:end/5,:), '.-k',
     rows(data)/5:rows(data) ,responseVector(end/5:end,:),'k--',
     1:rows(data), Result1, 'r--',
     #1:rows(data), Result2, 'r--',
     #1:rows(data), Result3, 'b--',
     1:rows(data), Result4, '--'
     
     );
legend("Training trajectory" , 
       "Continuing trajectory", 
       "Perfect identified trajectory",
       "Wrong identified trajectory"
       );
xlabel("number of state sample")
ylabel("angular acceleration in radiant")





