function [out] = crossover(obj1, obj2)
  size1 = size(obj1.command,1);
  size2 = size(obj2.command,1);
  countVariables = obj1.countVariables;
 
  crossoverPoint = randi(min(size1-countVariables,size2-countVariables))+countVariables;
 
  out = obj1;
  
  out.command = out.command(1:crossoverPoint,:);
  out.commandNumber = out.commandNumber(1:crossoverPoint);
  out.parameters = out.parameters(1:crossoverPoint*2+1);
 
  out.command(crossoverPoint:size2,:) = obj2.command(crossoverPoint:size2,:);
  out.commandNumber(crossoverPoint:size2) = obj2.commandNumber(crossoverPoint:size2);
  out.parameters(crossoverPoint*2:size2*2+1) = obj2.parameters(crossoverPoint*2:size2*2+1); 
end