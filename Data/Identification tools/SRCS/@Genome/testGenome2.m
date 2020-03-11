function obj = testGenome2(obj)
  listOfBuildingBlocks = ["+"; "-"; "*"; "sin"; "cos"; "const"; "variable"];
  
  obj.command(1,:) = listOfBuildingBlocks(7,:);
  obj.command(1,:) = listOfBuildingBlocks(7,:);
  obj.command(2:rows(obj.command),:) = [];
  obj.parameters(4:rows(obj.parameters),:) = [];
  obj.parameters(2) = 2;    
  obj.parameters(3) = 0;  
end