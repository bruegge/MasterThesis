function obj = loop(obj, iterations, threshold, printTimings) 
  
  fitnessBest = displayFitness(obj,1) 
  eaBest = obj;
  ende = 0;
  i = 0;
  
  timings = zeros(10,1);
  while i < iterations && ende == 0 
    i = i + 1;
    ["Generation" num2str(i)]
    tic
    obj = deleteDuplicates(obj);
    timings(1,1) = toc;
    tic
    obj = externalSelection(obj);
    timings(2) = toc;
    count = size(obj.genomes,2);
    tic
    parents = parentSelection(obj);
    timings(3) = toc;
    tic
    obj = inheritance(obj, parents);
    timings(4) = toc;
    tic
    obj = mutation(obj, count, 10); 
    timings(5) = toc;
    
    [obj, timings(6), timings(7), timings(8)] = fitnessEvaluation(obj, obj.dataMatrix, 0); 
    
    
    tic
    fitness = displayFitness(obj,1);
    timings(9) = toc;
    if printTimings == 1
      ["deleteDuplicates " num2str(timings(1))]
      ["externalSelection " num2str(timings(2))]
      ["parentSelection " num2str(timings(3))]
      ["inheritance " num2str(timings(4))]
      ["mutation " num2str(timings(5))]
      ["CreateTheta " num2str(timings(6))]
      ["Sindy " num2str(timings(7))]
      ["CalcFitness " num2str(timings(8))]
      ["displayFitness " num2str(timings(9))]
      ["All together " num2str(sum(timings,1))]
    endif
    
    if isnan(fitness) == 0
      if fitness > fitnessBest
        eaBest = obj;
        fitness
      endif
      if fitness > threshold
        ende = 1;
      endif 
    endif
    ["---------------------------"]

  endwhile
  display(eaBest,1);
  displayEquation(eaBest,1);
endfunction