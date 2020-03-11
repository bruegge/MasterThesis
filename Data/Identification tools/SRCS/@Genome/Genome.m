function obj = Genome(countVariables, variableNumber)

  #member variables
  mem.countVariables = countVariables;
  mem.variableNumber = variableNumber;
  mem.parameters = [0];
  mem.command = ["",""];
  mem.commandNumber = [0];
  mem.terms = [0];
  obj = class (mem, "Genome");

endfunction
  
