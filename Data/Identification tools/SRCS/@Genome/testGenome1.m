function obj = testGenome1(obj)
  listOfBuildingBlocks = ["+"; "-"; "*"; "sin"; "cos"; "const"; "variable"];
  obj.countVariables = 3;
  obj.variableNumber = 1;
  
  obj.command(1,:) = listOfBuildingBlocks(7,:);
  obj.parameters(2) = 1;   
  obj.parameters(3) = 0;  
  obj.command(2,:) = listOfBuildingBlocks(7,:);
  obj.parameters(4) = 2;   
  obj.parameters(5) = 0;  
  obj.command(3,:) = listOfBuildingBlocks(7,:);
  obj.parameters(6) = 3;   
  obj.parameters(7) = 0;  
 
  obj.command(4,:) = listOfBuildingBlocks(3,:);
  obj.parameters(8) = 1;   
  obj.parameters(9) = 2;  
 
  obj.command(5,:) = listOfBuildingBlocks(7,:);
  obj.parameters(10) = 3;   
  obj.parameters(11) = 0;  
 
 

  obj.command(6:rows(obj.command),:) = [];
  obj.parameters(12:rows(obj.parameters),:) = [];
end