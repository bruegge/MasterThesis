function obj = detectThermsOfInterest(obj)
  obj.terms = [0];
  
  values = zeros(size(obj.command,1),2);
  
  for i = 1:size(obj.command,1) 
    command = obj.command(i,:);
     
    if command == "variable"
      values(i,1) = 1;
      if i > obj.countVariables
        values(i,2) = 1;
      else
        values(i,2) = 0;
      endif
      
    endif
    if command == "const   "
      values(i,1) = 0;
      values(i,2) = 0;
    endif
    if command == "+       "
      values(i,1) = max(values(obj.parameters(i*2),1), values(obj.parameters(i*2+1),1));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || values(obj.parameters(i*2),1)> 0 || values(obj.parameters(i*2+1),1) >0
        values(i,2) = 1;
      endif
      values(obj.parameters(i*2),2) = 0;
      values(obj.parameters(i*2+1),2) = 0;
      if values(i,1) > 0
        values(i,1) = values(i,1) + 1;
      endif    
    endif
    if command == "-       "
      values(i,1) = max(values(obj.parameters(i*2),1), values(obj.parameters(i*2+1),1));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || values(obj.parameters(i*2),1)>0 || values(obj.parameters(i*2+1),1) >0
        values(i,2) = 1;
      endif
      values(obj.parameters(i*2),2) = 0;
      values(obj.parameters(i*2+1),2) = 0;
      if values(i,1) > 0
        values(i,1) = values(i,1) + 1;
      endif
    endif
    if command == "*       "
      values(i,1) = max(values(obj.parameters(i*2),1), values(obj.parameters(i*2+1),1));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || values(obj.parameters(i*2),1)>0 || values(obj.parameters(i*2+1),1) >0
        values(i,2) = 1;
      endif
      values(obj.parameters(i*2),2) = 0;
      values(obj.parameters(i*2+1),2) = 0;
      if values(i,1) > 0
        values(i,1) = values(i,1) + 1;
      endif
    endif
    if command == "/       "
      values(i,1) = max(values(obj.parameters(i*2),1), values(obj.parameters(i*2+1),1));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || values(obj.parameters(i*2),1)>0 || values(obj.parameters(i*2+1),1) >0
        values(i,2) = 1;
      endif
      values(obj.parameters(i*2),2) = 0;
      values(obj.parameters(i*2+1),2) = 0;
      if values(i,1) > 0
        values(i,1) = values(i,1) + 1;
      endif
    endif
    if command == "sin     "
      values(i,1) = values(obj.parameters(i*2),1);
      if obj.parameters(i*2) <= obj.countVariables || values(obj.parameters(i*2),1)>0
        values(i,2) = 1;
      endif
      values(obj.parameters(i*2),2) = 0;
      if values(i,1) > 0
        values(i,1) = values(i,1) + 1;
      endif 
    endif
    if command == "cos     "
      values(i,1) = values(obj.parameters(i*2),1);
      if obj.parameters(i*2) <= obj.countVariables || values(obj.parameters(i*2),1)>0
        values(i,2) = 1;
      endif
      values(obj.parameters(i*2),2) = 0;
      if values(i,1) > 0
        values(i,1) = values(i,1) + 1;
      endif
    endif
  endfor
 
  obj.terms = values;
endfunction