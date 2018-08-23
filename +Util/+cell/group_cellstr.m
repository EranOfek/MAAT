function [GroupStr,GroupInd,GroupFlag]=group_cellstr(CellStr,CaseSens)
% Find and group distinct values in a cell array.
% Package: Util.cell
% Description: Given a cell array of strings create a cell array of
%              indices of distinct strings (see example).
% Input  : - Cell vector of strings.
%          - Case sensitive {'y' | 'n'}, default is 'y'.
% Output : - Cell vector of distinctive groups.
%          - Cell vector of indices of distinct groups (see example).
%          - Cell vector of flags indicating if a member of the input
%            list belong to this group (1) or not (0).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                   May 2006
%    URL : http://weizmann.ac.il/home/eran/matlab/
% Example: A={'try1';'try2';'hhh';'try2'}
%          [GroupStr,GroupInd]=Util.cell.group_cellstr(A);
%          % Group contains:
%          % GroupStr{1}='try1';  GroupInd{1}=[1];
%          % GroupStr{2}='try2';  GroupInd{2}=[2 4];
%          % GroupStr{3}='hhh';   GroupInd{3}=[3];
% Reliable: 1
%------------------------------------------------------------------------------
if (nargin==1),
   CaseSens = 'y';
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end

N    = length(CellStr);
VecN = [1:1:N].';
GroupStr = {};
GroupInd = {};
FI = 0;
for I=1:1:N,
   CurrentStr = CellStr{I};

   if (isempty(GroupStr)==0),
      % check if string already found
      switch CaseSens
       case 'y'
          Found      = sum(strcmp(CurrentStr,GroupStr));
       case 'n'
          Found      = sum(strcmpi(CurrentStr,GroupStr));
       otherwise
          error('Unknown CaseSens Option');
      end
   else
      Found = 0;
   end	
   if (Found>0),
      % string already found - go to next
   else
      switch CaseSens
       case 'y'
          Ind      = strcmp(CurrentStr,CellStr);
       case 'n'
          Ind      = strcmpi(CurrentStr,CellStr);
       otherwise
          error('Unknown CaseSens Option');
      end
      FI = FI + 1;
      GroupFlag{FI} = Ind==1 & VecN>=I;
      FoundInd = find(GroupFlag{FI});
      GroupStr{FI} = CurrentStr;
      GroupInd{FI} = FoundInd;
   end
end
