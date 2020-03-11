function obj = deleteDuplicates(obj) 
  hashtable = containers.Map('KeyType','char','ValueType','double');

  i = 0;
  while i < size(obj.genomes,2) 
    i = i + 1;
    
    [string, string2] = genome2String(obj.genomes(i));
    #string = "test";
    
    if string
      if isKey(hashtable,string) == 1
        obj.genomes(i) = [];
        obj.fitness(i) = [];
        obj.xi(i) = 0;
        i = i - 1;
      else
        subsasgn(hashtable, substruct( '()', {string} ), 1.0 );
      endif 
    endif
  endwhile


endfunction



