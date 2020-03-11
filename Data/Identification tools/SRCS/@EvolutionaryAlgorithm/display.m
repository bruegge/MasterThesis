function display(obj, returnBest) 
  if returnBest == 0
    for i = 1:size(obj.genomes,2) 
      i
      display(obj.genomes(i));
      ["k = " num2str(obj.xi(i))]
    endfor
    ["Fitness: " num2str(obj.fitness)]
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
    display(obj.genomes(fitnessIndex));
    ["Fitness: " num2str(obj.fitness(fitnessIndex))]
    ["k = " num2str(obj.xi(fitnessIndex,1:obj.xiCount))]
  endif
end