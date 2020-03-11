function obj = externalSelection(obj) 
  
  while size(obj.genomes,2) > obj.countExtSelection
    smallestIndex = 1;  
    
    for i = 1:size(obj.genomes,2)
      if obj.fitness(i) < obj.fitness(smallestIndex) || isnan(obj.fitness(i)) == 1
        smallestIndex = i;
      endif
    endfor 
    deleteIndex = smallestIndex;
    #delete individual
    obj.xi(smallestIndex,:) = zeros(1,100);
    obj.fitness(smallestIndex) = [];
    obj.genomes(smallestIndex) = [];
  endwhile 
end