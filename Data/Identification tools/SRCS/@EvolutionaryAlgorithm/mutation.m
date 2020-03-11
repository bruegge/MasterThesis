function obj = mutation(obj, countExtSelection, mutationProbability) 
  for i =countExtSelection+1:obj.countGenomes
    obj.genomes(i) = mutate(obj.genomes(i), mutationProbability, obj.onlyMonomes);
  endfor
end