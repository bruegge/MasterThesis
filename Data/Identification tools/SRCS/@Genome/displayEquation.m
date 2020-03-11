function result = displayEquation(obj, xi)
  values = zeros(size(obj.command,1),2);
  
  for i = 1:size(obj.command,1) 
    command = obj.command(i,:);
     
    if command == "variable"
      #values(i,1) = variables(obj.parameters(i*2));
      values(i,1) = 1;
      if i > obj.countVariables
        values(i,2) = 1;
      else
        values(i,2) = 0;
      endif
      
    endif
    if command == "const   "
      #values(i,1) = obj.parameters(i*2);
      values(i,1) = 0;
      values(i,2) = 0;
    endif
    if command == "+       "
      #values(i,1) = values(obj.parameters(i*2),1) + values(obj.parameters(i*2+1),1);
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
      #values(i,1) = values(obj.parameters(i*2),1) - values(obj.parameters(i*2+1),1);
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
      #values(i,1) = values(obj.parameters(i*2),1) * values(obj.parameters(i*2+1),1);
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
      #values(i,1) = values(obj.parameters(i*2),1) / values(obj.parameters(i*2+1),1);
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
      #values(i,1) = sin(values(obj.parameters(i*2),1));
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
      #values(i,1) = cos(values(obj.parameters(i*2),1));
      values(i,1) = values(obj.parameters(i*2),1);
      if obj.parameters(i*2) <= obj.countVariables || values(obj.parameters(i*2),1)>0
        values(i,2) = 1;
      endif
      values(obj.parameters(i*2),2) = 0;
      if values(i,1) > 0
        values(i,1) = values(i,1) + 1;
      endif
    endif
    if command == "sqrt    "
      #values(i,1) = cos(values(obj.parameters(i*2),1));
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
  
  function string = rec(obj, commandNumber)
    command = obj.command(commandNumber,:);
    param1 = obj.parameters(commandNumber*2);
    param2 = obj.parameters(commandNumber*2+1);
    
    string = "";    
    if command == "variable"
      string = ["V" num2str(param1)];
    endif
    if command == "const   "
      string = ["(" num2str(param1) ")"];
    endif
    if command == "+       "
      string = ["(" rec(obj, param1) "+" rec(obj, param2) ")"];
    endif
    if command == "-       "
      string = ["(" rec(obj, param1) "-" rec(obj, param2) ")"];
    endif
    if command == "*       "
      string = ["(" rec(obj, param1) "*" rec(obj, param2) ")"];
    endif
    if command == "/       "
      string = ["(" rec(obj, param1) "/" rec(obj, param2) ")"];
    endif
    if command == "sin     "
      string = ["sin(" rec(obj, param1) ")"];
    endif
    if command == "cos     "
      string = ["cos(" rec(obj, param1) ")"];
    endif 
    if command == "sqrt    "
      string = ["sqrt(abs(" rec(obj, param1) "))"];
    endif 
  endfunction
  
  xiCounter = 1;
  for i= 1:size(values,1)
    if values(i,2) == 1
      [num2str(i) ":   " num2str(xi(1,xiCounter)) " * " rec(obj, i)]
      xiCounter = xiCounter + 1;
    endif
  endfor
  ["     " num2str(xi(1,xiCounter)) " * 1.0"]
endfunction