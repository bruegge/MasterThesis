function [fitness, xi, timingCreateTheta, timingSindy, timingCalcFitness] = fitnessEvaluation(obj, DataMatrix, dXdt, dXdtTest, bPrint, h, ea)

  function error = MeanGradientError(theta, dXdt, xi)
    #{
    theta = theta * xi;
    r = rows(theta);
    error = sum(abs(dXdt - theta)) / r;
    #}
    
    
    error = 100000;
    Y = theta * xi;
    Y = dXdt - Y;
    Y = Y.*Y;
    r = rows(Y);
    c= columns(Y);
    summe = sum(Y);
    
    error = sqrt(summe) * 1 / sqrt(rows(theta));
    #}
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
  
  col = columns(DataMatrix);
  tic 
  theta = calculate2(obj, DataMatrix, ea);   
  theta = theta(1:end-1, :);
  timingCreateTheta = toc;
 
  tic
  xi = SparseRegression(theta(1:end/2,:), dXdt, 0); 
  timingSindy = toc;
  tic
  fitness = -MeanGradientError(theta, dXdtTest, xi, h);
  timingCalcFitness = toc;
end





