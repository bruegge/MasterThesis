function [string2, string] = genome2String(obj)
  string = "";
  
 
  function [ o ] = myfastint2str( x )
    maxvalue=max(x(:));
    %maxvalue=intmax(class(x));%Alternative implementation based on class
    required_digits=ceil(log(double(maxvalue+1))/log(10));
    o=repmat(x(1)*0,size(x,1),required_digits);%initialize array of required size
    for c=size(o,2):-1:1
       o(:,c)=mod(x,10);
       x=(x-o(:,c))/10;
    end
    o=char(o+'0');
  endfunction
  
  count = size(obj.commandNumber,2);
  array = zeros(count*3,1);
  array(1:count,1) = obj.commandNumber(:);
  array(count+1:count*2,1) = obj.parameters(1:count);
  array(count*2+1:count*3,1) = obj.parameters(count+1:count*2);
  
  string2 = myfastint2str(array*1000)';
  stringRows = rows(string2);
  stringcols = columns(string2);
  
  stringComplete = "";
  col = 1;
  for i = 1: stringRows   
    stringComplete(1,col:col+stringcols-1) = string2(i,:);
    col = col + stringcols;
  endfor
  string2 = stringComplete;
  
  #{
  for i = 1:size(obj.command,1)   
    string = [string num2str(obj.command(i)) num2str(obj.parameters(i*2)) num2str(obj.parameters(i*2+1))];
  endfor
  #}
  string = string2;
endfunction