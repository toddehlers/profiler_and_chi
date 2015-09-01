function answ=answer_yn(questionstring)

% answer_yn is a simple function to return a binary flag in response 
% to a yes/no question 
% 
% Usage:
% 	answ=answer_yn(questionstring)
% Input
% 	A y/n answer to the input string,
% Output
% 	1 IF YES, 0 IF NO

%
% fcn requires one of these inputs: y/yes/n/no
loop=0;
while loop==0,
    yn=input([questionstring ' (y/n)....'],'s');
    if strcmp(yn,'y') | strcmp(yn,'yes'),
        answ=1;
        loop=1;
    elseif strcmp(yn,'n') | strcmp(yn,'no'),
        answ=0;
        loop=1;
    else
        disp('You must answer either y, n, yes, or no');
    end
end
