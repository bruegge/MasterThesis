function obj = mutate(obj, probabilityInPercent, onlyMonomes)
  
  listOfBuildingBlocks = ["+"; "-"; "*"; "/"; "sin"; "cos"; "sqrt"; "const"; "variable"];
  countVariables = obj.countVariables;
 
  for i = countVariables+1:size(obj.command,1) 
    ran = randi(100);
    if ran <= probabilityInPercent
      command = obj.command(i,:);
    
      if command == "const   "
        obj.parameters(i*2) = obj.parameters(i*2) + normrnd(0,1);
      else
        ran = randi(9);
        if onlyMonomes == 1
          ran = 3;
        endif
        obj.command(i,:) = listOfBuildingBlocks(ran,:);
        obj.commandNumber(i) = ran;
        if ran < 5 
          obj.parameters(i*2) = randi(i-1);
          obj.parameters(i*2+1) = randi(i-1); 
        endif
        
        if ran >= 5 && ran <= 6
          obj.parameters(i*2) = randi(i-1);       
          obj.parameters(i*2+1) = 0;       
        endif
      
        if ran == 7
          obj.parameters(i*2) = randi(i-1);       
          obj.parameters(i*2+1) = 0;       
        endif
      
        if ran == 8
          obj.parameters(i*2) = normrnd(0,1)* 100; 
          obj.parameters(i*2+1) = 0;       
        endif
        
        if ran == 9
          obj.parameters(i*2) = randi(countVariables); 
          obj.parameters(i*2+1) = 0;       
        endif
      endif
  
    endif
  endfor
end