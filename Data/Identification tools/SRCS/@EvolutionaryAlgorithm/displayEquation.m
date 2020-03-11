function fitness = displayEquation(obj, returnBest) 
  if returnBest == 0
      
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
    displayEquation(obj.genomes(fitnessIndex), obj.xi(fitnessIndex,:));
  endif
end