function catfile(OutputFile,varargin)
% Concatenate files into a single file.
% Package: Util.files
% Description: Concatenate files into a single file. Add a carriage
%              return at the end of each Concatenated file.
% Input  : - Output file name.
%          * Arbitrary number of strings containing files name.
%            Alternativel, a cell array containing files name.
% Output : null
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

N = length(varargin);

if (N==1),
   Files = varargin{1};
else
   Files = varargin(:);
end

FID = fopen(OutputFile,'w');

Nf = length(Files);
for If=1:1:Nf,
   FID_r = fopen(Files{If},'r');
   while (feof(FID_r)==0),
      Line = fgets(FID_r);
      fprintf(FID,'%s',Line);
   end
   fprintf(FID,'\n');
   fclose(FID_r);
end
fclose(FID);

