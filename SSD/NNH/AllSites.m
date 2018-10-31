function RCell = AllSites(z, N)
  RCell{1}{1} = [];
  
  for s = 1:N
    count = 1;
    for j = 1:NShell(z, s-1)
      for k = 1:NBranch(z, s-1) 
	RCell{s+1}{count} = [RCell{s}{j} k];
	count += 1;
      endfor
    endfor
  endfor
    
  RCell = [RCell{:}];
endfunction
