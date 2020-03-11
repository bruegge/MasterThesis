function [obj, timingCreateTheta, timingSindy, timingCalcFitness] = fitnessEvaluation(obj, DataMatrix, bPrint) 
  numGenomes = size(obj.genomes,2);
  timingCreateTheta = 0;
  timingSindy = 0;
  timingCalcFitness = 0;
  for i = size(obj.fitness,1):numGenomes
    
    [fitness, xi, timingCreateTheta1, timingSindy1, timingCalcFitness1] = fitnessEvaluation( obj.genomes(i) , DataMatrix, obj.dXdt, obj.dXdtTest, bPrint, obj.h, obj); 
    xi = xi';
    timingCalcFitness = timingCalcFitness + timingCalcFitness1;
    timingCreateTheta = timingCreateTheta + timingCreateTheta1;
    timingSindy = timingSindy + timingSindy1;
     
    obj.xi(i,:) = zeros(1,100);
    obj.xi(i,1:size(xi,2)) = xi;
    obj.xiCount = size(xi,2);
    obj.fitness(i) = fitness;
    
  endfor
end