function string = genomeTherm2String(obj, thermNumber)
  string = "";
  
  function string = rec(obj, commandNumber)
    command = obj.command(commandNumber,:);
    param1 = obj.parameters(commandNumber*2);
    param2 = obj.parameters(commandNumber*2+1);
    
    string = "";    
    if command == "variable"
      string = ["V" num2str(param1)];
    endif
    if command == "const   "
      string = [num2str(param1)];
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
  endfunction
  
  string = rec(obj, thermNumber);
endfunction