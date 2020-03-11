function fitness = displayFitness(obj, returnBest) 
  if returnBest == 0
    fitness = zeros(size(obj.genomes,2),1);
    
    for i = 1:size(obj.genomes,2)
      fitness(i,1) = obj.fitness(i);
      ["Fitness genome " num2str(i) " = " num2str(fitness(i,1))] 
    endfor
    ["----------------------"]   
  else
    fitnessIndex = 1;
    while isnan(obj.fitness(fitnessIndex)) == 1
      fitnessIndex = fitnessIndex + 1;
    endwhile
    for i = 2:size(obj.genomes,2)
      if isnan(obj.fitness(i)) == 0
        if obj.fitness(i) > obj.fitness(fitnessIndex)
          fitnessIndex = i;
        endif
      endif
    endfor
    fitness = obj.fitness(fitnessIndex);
  endif
end