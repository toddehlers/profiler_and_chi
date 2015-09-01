function row=closeto(value,col)

% closeto.m finds the row that contains a value in a 
%     specified column that is nearest to designated ‘value’. Only returns 
%     one value, even if there are multiple matching values.  Returns 
%     either the lowest or highest value if there are multiple matching 
%     values.
% 
% Usage:
% 	row = closeto(value,col)
% OUTPUT:
% 	row - row that has value in col closest to 'value'
% Input:  
% 	col - column to look at
% 	value - value want to find or be closest to

%set minormax=0 for low value if match, minormax=1 if highvalue to match.
minormax=1;

row=find(col==value);
if length(row)==0,
  gtval=min(col(find(col>value)));
  ltval=max(col(find(col<value)));
  if find(gtval) & isempty(ltval),  %case where value is below all values of col
    row=find(col==gtval);
    %  elseif find(ltval) & size(find(gtval))==0,
  elseif find(ltval) & isempty(gtval),  %case where value is above all values of col
    row=find(col==ltval);
  elseif gtval-value<=value-ltval,  %closer to gtval than ltval
    row=find(col==gtval);
  elseif gtval-value>=value-ltval,  %closer to ltval than gtval
    row=find(col==ltval);
  else
    warning('error in closeto.m')
    keyboard
  end
end

if length(row)>1,
    if minormax==0,
        row=row(1);
    else
        row=row(length(row));
    end
end