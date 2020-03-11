function result = calculate (obj, variables)
  values = zeros(size(obj.command,1),3);
  #----------------------------------------here is the error!!!!!!!!---------------#
  for i = 1:size(obj.command,1) 
    command = obj.command(i,:);
     
    if command == "variable"
      values(i,1) = variables(obj.parameters(i*2));
      values(i,2) = 1;
      if i > obj.countVariables
        values(i,3) = 1;
      else
        values(i,3) = 0;
      endif
      
    endif
    if command == "const   "
      values(i,1) = obj.parameters(i*2);
      values(i,2) = 0;
      values(i,3) = 0;
    endif
    if command == "+       "
      values(i,1) = values(obj.parameters(i*2),1) + values(obj.parameters(i*2+1),1);
      values(i,2) = max(values(obj.parameters(i*2),2), values(obj.parameters(i*2+1),2));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || values(obj.parameters(i*2),2)> 0 || values(obj.parameters(i*2+1),2) >0
        values(i,3) = 1;
      endif
      values(obj.parameters(i*2),3) = 0;
      values(obj.parameters(i*2+1),3) = 0;
      if values(i,2) > 0
        values(i,2) = values(i,2) + 1;
      endif    
    endif
    if command == "-       "
      values(i,1) = values(obj.parameters(i*2),1) - values(obj.parameters(i*2+1),1);
      values(i,2) = max(values(obj.parameters(i*2),2), values(obj.parameters(i*2+1),2));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || values(obj.parameters(i*2),2)>0 || values(obj.parameters(i*2+1),2) >0
        values(i,3) = 1;
      endif
      values(obj.parameters(i*2),3) = 0;
      values(obj.parameters(i*2+1),3) = 0;
      if values(i,2) > 0
        values(i,2) = values(i,2) + 1;
      endif
    endif
    if command == "*       "
      values(i,1) = values(obj.parameters(i*2),1) * values(obj.parameters(i*2+1),1);
      values(i,2) = max(values(obj.parameters(i*2),2), values(obj.parameters(i*2+1),2));
      if obj.parameters(i*2) <= obj.countVariables || obj.parameters(i*2+1) <= obj.countVariables || values(obj.parameters(i*2),2)>0 || values(obj.parameters(i*2+1),2) >0
        values(i,3) = 1;
      endif
      values(obj.parameters(i*2),3) = 0;
      values(obj.parameters(i*2+1),3) = 0;
      if values(i,2) > 0
        values(i,2) = values(i,2) + 1;
      endif
    endif
    if command == "sin     "
      values(i,1) = sin(values(obj.parameters(i*2),1));
      values(i,2) = values(obj.parameters(i*2),2);
      if obj.parameters(i*2) <= obj.countVariables || values(obj.parameters(i*2),2)>0
        values(i,3) = 1;
      endif
      values(obj.parameters(i*2),3) = 0;
      if values(i,2) > 0
        values(i,2) = values(i,2) + 1;
      endif 
    endif
    if command == "cos     "
      values(i,1) = cos(values(obj.parameters(i*2),1));
      values(i,2) = values(obj.parameters(i*2),2);
      if obj.parameters(i*2) <= obj.countVariables || values(obj.parameters(i*2),2)>0
        values(i,3) = 1;
      endif
      values(obj.parameters(i*2),3) = 0;
      if values(i,2) > 0
        values(i,2) = values(i,2) + 1;
      endif
    endif
  endfor
 
  result = 0;
  resultCount = 1;
  for i = 1:size(values,1) 
    if values(i,3) == 1
      result(1,resultCount) = values(i,1);
      resultCount = resultCount + 1;
    endif  
  endfor
  result(1,columns(result)+1) = 1;
endfunction