function display(obj) 
  ["-------------------------"]
  countVariables = obj.countVariables 
  variableNumber = obj.variableNumber
  
  for i = 1:size(obj.command,1) 
    command = obj.command(i,:);
    if command == "variable"
      [num2str(i) ": Variable " num2str(obj.parameters(i*2))]
    endif
    if command == "const   "
      [num2str(i) ": const " num2str(obj.parameters(i*2))]
    endif
    if command == "+       "
      [num2str(i) ": + " num2str(obj.parameters(i*2)) " , " num2str(obj.parameters(i*2+1))]
    endif
    if command == "-       "
      [num2str(i) ": - " num2str(obj.parameters(i*2)) " , " num2str(obj.parameters(i*2+1))]
    endif
    if command == "*       "
      [num2str(i) ": * " num2str(obj.parameters(i*2)) " , " num2str(obj.parameters(i*2+1))]
    endif
    if command == "/       "
      [num2str(i) ": / " num2str(obj.parameters(i*2)) " , " num2str(obj.parameters(i*2+1))]
    endif
 
    if command == "sin     "
      [num2str(i) ": sin " num2str(obj.parameters(i*2))]
    endif
    if command == "cos     "
      [num2str(i) ": cos " num2str(obj.parameters(i*2))]
    endif
    if command == "sqrt    "
      [num2str(i) ": sqrt(abs( " num2str(obj.parameters(i*2)) " ))"]
    endif
  endfor
 
end