function [S,S1,S2]=bsx_size(S1,S2)
% bsx_size gets two arrays S1,S2 and calculates the size of the result
% array after a bsxfun operation is performed. 
% The function verifies that the sizes are either identical, or for
% dimensions where thtey are differs, the dimension of one of the array
% must be scalar.
% The function also padding the size arrays with ones, in case one of the
% size arrays is shorter (i.e. its last dimensions are singleton.
    if (sum(S1(S1~=S2)~=1 & S2(S1~=S2)~=1)~=0)
        error('Matrix dimension must agree.');
    end

    D1 = length(S1);
    D2 = length(S2);
    if (D1>D2)
        S2(D1+1:D1)=1;
    elseif (D2>D1)
        S1(D2+1:D2)=1;
    end
    S = max(S1,S2);
end