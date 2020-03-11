function values = calculate2(obj, variables, ea)
  values = zeros(size(variables,1),size(obj.command,1));
  treeInfo = zeros(size(obj.command,1),2);
  
  for i = 1:size(obj.command,1)
    command = obj.command(i,:);
    param1 = obj.parameters(i*2);
    param2 = obj.parameters(i*2 + 1);
    
    if command == "variable"
      values(:,i) = variables(:,param1);
      treeInfo(i,1) = 1;
      if i > obj.countVariables
        treeInfo(i,2) = 1;
      else
        treeInfo(i,2) = 0;
      endif
    endif
    
    if command == "const   "
      values(:,i) = param1;
      treeInfo(i,1) = 0;
      treeInfo(i,2) = 0;
    endif
    
    if command == "+       "
      values(:,i) = values(:,param1) + values(:,param2);
      
      treeInfo(i,1) = max(treeInfo(obj.parameters(i*2),1), treeInfo(obj.parameters(i*2+1),1));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || treeInfo(obj.parameters(i*2),1)> 0 || treeInfo(obj.parameters(i*2+1),1) >0
        treeInfo(i,2) = 1;
      endif
      treeInfo(obj.parameters(i*2),2) = 0;
      treeInfo(obj.parameters(i*2+1),2) = 0;
      if treeInfo(i,1) > 0
        treeInfo(i,1) = treeInfo(i,1) + 1;
      endif    
    endif
    
    if command == "-       "
      values(:,i) = values(:,param1) - values(:,param2);
      treeInfo(i,1) = max(treeInfo(obj.parameters(i*2),1), treeInfo(obj.parameters(i*2+1),1));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || treeInfo(obj.parameters(i*2),1)> 0 || treeInfo(obj.parameters(i*2+1),1) >0
        treeInfo(i,2) = 1;
      endif
      treeInfo(obj.parameters(i*2),2) = 0;
      treeInfo(obj.parameters(i*2+1),2) = 0;
      if treeInfo(i,1) > 0
        treeInfo(i,1) = treeInfo(i,1) + 1;
      endif    
    endif
    
    if command == "*       "
      values(:,i) = values(:,param1) .* values(:,param2);
      treeInfo(i,1) = max(treeInfo(obj.parameters(i*2),1), treeInfo(obj.parameters(i*2+1),1));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || treeInfo(obj.parameters(i*2),1)> 0 || treeInfo(obj.parameters(i*2+1),1) >0
        treeInfo(i,2) = 1;
      endif
      treeInfo(obj.parameters(i*2),2) = 0;
      treeInfo(obj.parameters(i*2+1),2) = 0;
      if treeInfo(i,1) > 0
        treeInfo(i,1) = treeInfo(i,1) + 1;
      endif    
    endif
    
    if command == "/       "
      values(:,i) = values(:,param1) ./ values(:,param2);
      treeInfo(i,1) = max(treeInfo(obj.parameters(i*2),1), treeInfo(obj.parameters(i*2+1),1));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || treeInfo(obj.parameters(i*2),1)> 0 || treeInfo(obj.parameters(i*2+1),1) >0
        treeInfo(i,2) = 1;
      endif
      treeInfo(obj.parameters(i*2),2) = 0;
      treeInfo(obj.parameters(i*2+1),2) = 0;
      if treeInfo(i,1) > 0
        treeInfo(i,1) = treeInfo(i,1) + 1;
      endif    
    endif
    
    if command == "sin     "
      values(:,i) = sin(values(:,param1));
      treeInfo(i,1) = treeInfo(obj.parameters(i*2),1);
      if obj.parameters(i*2) <= obj.countVariables || treeInfo(obj.parameters(i*2),1)>0
        treeInfo(i,2) = 1;
      endif
      treeInfo(obj.parameters(i*2),2) = 0;
      if treeInfo(i,1) > 0
        treeInfo(i,1) = values(i,1) + 1;
      endif 
    endif
    
    if command == "cos     "
      values(:,i) = cos(values(:,param1));
      treeInfo(i,1) = treeInfo(obj.parameters(i*2),1);
      if obj.parameters(i*2) <= obj.countVariables || treeInfo(obj.parameters(i*2),1)>0
        treeInfo(i,2) = 1;
      endif
      treeInfo(obj.parameters(i*2),2) = 0;
      if treeInfo(i,1) > 0
        treeInfo(i,1) = values(i,1) + 1;
      endif 
    endif
    
    if command == "sqrt    "
      values(:,i) = sqrt(abs(values(:,param1)));
      treeInfo(i,1) = treeInfo(obj.parameters(i*2),1);
      if obj.parameters(i*2) <= obj.countVariables || treeInfo(obj.parameters(i*2),1)>0
        treeInfo(i,2) = 1;
      endif
      treeInfo(obj.parameters(i*2),2) = 0;
      if treeInfo(i,1) > 0
        treeInfo(i,1) = values(i,1) + 1;
      endif 
    endif
    
    
    
  endfor
  
  i = 1;
  while i <= size(treeInfo,1) 
    
    if treeInfo(i,2) != 1
      values(:,i) = [];
      treeInfo(i,:) = [];
    else
      i = i + 1;
    endif  
  endwhile
  values(:, size(values,2)+1) = ones(size(values,1),1);
   
endfunction