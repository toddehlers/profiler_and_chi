function rows=closeto_allvalues(value,col)

% closeto_allvalues.m finds ALL rows that contains a value 
%     in a specified column that is nearest to designated ‘value’.  Returns 
%     ALL matching indices (into vector col), not just one.
% 
% Usage:
% 	rows = closeto_allvalues(value,col)
% OUTPUT:
% 	rows - rows that has value in col closest to 'value'
% Input:  
% 	col - column to look at
% 	value - value want to find or be closest to

%set minormax=0 for low value if match, minormax=1 if highvalue to match.
minormax=1;

if isempty(col) | isempty(value),
%    disp('Empty variable passed to closeto_allvalues; returning no rows')
    rows=[];
    return
end

rows=find(col==value);
if length(rows)==0,
  gtval=min(col(find(col>value)));
  ltval=max(col(find(col<value)));
  if find(gtval) & isempty(ltval),
    rows=find(col==gtval);
  elseif find(ltval) & isempty(gtval),
    rows=find(col==ltval);
  elseif gtval-value<=value-ltval,
    rows=find(col==gtval);
  elseif gtval-value>=value-ltval,
    rows=find(col==ltval);
  else
    warning('error!')
    keyboard
  end
end

