function obj = randomGenome(obj, onlyMonomes)
  
  countVariables = obj.countVariables;
  n = randi(20)+ countVariables;
  listOfBuildingBlocks = ["+"; "-"; "*"; "/"; "sin"; "cos"; "sqrt"; "const"; "variable"];
  
  for i = 1:countVariables      
    obj.command(i,:) = listOfBuildingBlocks(9,:);
    obj.commandNumber(i) = 9;
    obj.parameters(i*2) = i;
    obj.parameters(i*2+1) = 0; 
  endfor
  
  for i = countVariables+1:n  
    bb = randi(9);
    if onlyMonomes == 1
      bb = 3;
    endif
    obj.command(i,:) = listOfBuildingBlocks(bb,:);
    obj.commandNumber(i) = bb;
    if bb < 5
      obj.parameters(i*2) = randi(i-1);
      obj.parameters(i*2+1) = randi(i-1);      
    endif
    
    if bb >= 5 && bb <= 6
      obj.parameters(i*2) = randi(i-1);       
      obj.parameters(i*2+1) = 0;       
    endif
  
    if bb == 7
      obj.parameters(i*2) =randi(i-1);
      obj.parameters(i*2+1) = 0;
    endif
    
    if bb == 8
      obj.parameters(i*2) = normrnd(0,1)* 100; 
      obj.parameters(i*2+1) = 0;       
    endif
    
    if bb == 9
      obj.parameters(i*2) = randi(countVariables); 
      obj.parameters(i*2+1) = 0;       
    endif
  endfor  
  
end