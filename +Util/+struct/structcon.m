function S=structcon(S1,S2,Dim)
% Concatenate two structures into one.
% Package: Util.struct
% Description: Concatenate two structures into one.
%              Example: S1.A S1.B, and S2.A, S2.B:
%                       S=structcon(S1,S2,1); will return a structure S,
%                       with fields A and B which contains the content
%                       of [S1.A;S2.A] and [S1.B;S2.B], respectively.
%              The concatantion can be done along the 1st or 2nd dimensions.
% Input  : - First structure.
%          - Second structure, containing the same fields as the
%            first structure.
%          - Dimension along to concatenate the fields.
%            1 - column concatenation;
%            2 - raw concatenation.
% Output : - New structure
% Tested : Matlab 7.3
%     By : Eran O. Ofek                  December 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: A.a.a = 1; A.a.b=[2 3]; A.c='ab'; A.d{1}='a'; A.d{2}=1;
%          S=Util.struct.structcon(A,A,1);
% Reliable: 2
%------------------------------------------------------------------------------

FN = fieldnames(S1);

NotExist = 1;
for I=1:1:length(FN)
   %if (eval(sprintf('isstruct(S1.%s)',FN{I}))==1),
   if (isstruct(S1.(FN{I}))==1)
      % field is a substructure - call structcon.m in recursive
      %Content = eval(sprintf('structcon(S1.%s,S2.%s,%d)',FN{I},FN{I},Dim));
      Content = structcon(S1.(FN{I}),S2.(FN{I}),Dim);

      if (NotExist==1)
         S = struct(FN{I},Content);
         NotExist = 0;
      else
          S.(FN{I}) = Content;
         %S=setfield(S,FN{I},Content);
      end
   else
      switch Dim
       case 1
          Content = [S1.(FN{I}); S2.(FN{I})];
       case 2
          Content = [S1.(FN{I}), S2.(FN{I})];
       otherwise
          error('Unknown dimension option');
      end

      if (NotExist==1),
         S = struct(FN{I},Content);
         NotExist = 0;
      else
          S.(FN{I}) = Content;
         %S=setfield(S,FN{I},Content);
      end
   end
end
