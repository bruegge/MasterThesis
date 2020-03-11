function obj = EvolutionaryAlgorithm(countVariables, variableNumber, countGenomes, countExtSelection, dataMatrix, dXdt, h, onlyMonomes)

  #member variables
  mem.countVariables = countVariables;
  mem.variableNumber = variableNumber;
  mem.countGenomes = countGenomes;
  mem.countExtSelection = countExtSelection;
  mem.dataMatrix = dataMatrix;
  mem.onlyMonomes = onlyMonomes; 
  mem.h = h;
  mem.dXdt = dXdt(1:end/2,variableNumber);
  mem.dXdtTest = dXdt(:,variableNumber);
  mem.fitness = [0];
  mem.genomes = [""];
  mem.xi = zeros(countGenomes,100);
  mem.xiCount = 0;
  mem.FitnessResults = containers.Map('KeyType','char','ValueType','double');
  for i = 1:countGenomes
    genome = Genome(countVariables,variableNumber);
    genome = randomGenome(genome, onlyMonomes);
    mem.genomes(i) = genome;  
  endfor
  
  obj = class (mem, "EvolutionaryAlgorithm");

endfunction
  
