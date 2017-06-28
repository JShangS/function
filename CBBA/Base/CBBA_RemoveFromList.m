% Copyright 2010
% Massachusetts Institute of Technology
% All rights reserved
% Developed by the Aerospace Controls Lab, MIT

%---------------------------------------------------------------------%
% Remove item from list at location specified by index
%---------------------------------------------------------------------%

function newList = CBBA_RemoveFromList(oldList, index)

newList = -ones(1,length(oldList));

newList(1,1:index-1)   = oldList(1,1:index-1);

newList(1,index:end-1) = oldList(1,index+1:end);

end