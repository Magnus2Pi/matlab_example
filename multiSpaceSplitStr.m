function resultList = multiSpaceSplitStr(s)
%multiSpaceSplitStr splits a string into a cell list based on multiple spaces.
%
% resultList = multiSpaceSplitStr(s) splits the string {s} into a cell list
% {resultList} using multiple connected spaces
%
% s                 - input string
%
%_______________________________________________________________________
% Examples
% s = 'This will be the first elemet.  This will be the second    Nr 3';
% resultList = multiSpaceSplitStr(s);
%_______________________________________________________________________
% Erling Hugo Jensen, 29/01/08
error(nargchk(1,1,nargin,'struct'));

%% Init
lngth = length(s);
lst = {};
index = 1;
lastSplitIndex = index;
tmpS = s;

%% Split string
while (index <= lngth)
    i = findSpace(tmpS, index);
    if (i < (lngth-1))
        j = findNonSpace(tmpS, i+1);
        index = j;
    else
        i = i + 1;   % adjustments made to include last character
        j = lngth;
        index = lngth+1;
    end
    if (((j-i) > 1) || (i >= (lngth)))
        lst = horzcat(lst,  ddewhite(tmpS(lastSplitIndex:(i-1))));
        lastSplitIndex = j;
    end
end

resultList = lst;

%% Finds index to first space in string starting from startIndex
function index = findSpace(s, startIndex)
    i = startIndex;
    lngth = length(s);
    while ((s(i:i) ~= ' ') && (i < lngth))
        i = i+1;
    end
    index = i;
    return

%% Finds index to first none-space in string starting from startIndex
function index = findNonSpace(s, startIndex)
    i = startIndex;
    lngth = length(s);
    while ((s(i:i) == ' ') && (i < lngth))
        i = i+1;
    end
    index = i;
    return
    
    
    
            
        