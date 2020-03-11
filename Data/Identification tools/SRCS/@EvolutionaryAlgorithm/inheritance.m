function obj = inheritance(obj, indexParents) 
  countInheritance = rows(indexParents);
  countGenomes = size(obj.genomes,2);
  for i=1:countInheritance
    out = crossover(obj.genomes(indexParents(i,1)), obj.genomes(indexParents(i,2)));
    countGenomes = countGenomes + 1;    
    obj.xi(countGenomes,:) = 0;
    obj.fitness(countGenomes) = 0;
    obj.genomes(countGenomes) = out; 
    
  endfor
end