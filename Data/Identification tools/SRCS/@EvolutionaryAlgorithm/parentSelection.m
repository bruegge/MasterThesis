function genomeIDs = parentSelection(obj) 
  numGenomes = size(obj.genomes,2);
  numNewGenomes = obj.countGenomes - numGenomes;
  genomeIDs = zeros(numNewGenomes,2);
  for i = 1:numNewGenomes
    indexP1 = randi(numGenomes);
    indexP2 = randi(numGenomes);
    while indexP1 == indexP2
      indexP2 = randi(numGenomes);
    endwhile
    genomeIDs(i,:) = [indexP1 indexP2];
  endfor
end