function Rand=rand_range(Size,Range,Type)
%--------------------------------------------------------------------------
% rand_range function                                            AstroStat
% Description: Generate uniformly random number in a given range.
%              The numbers can be uniform either in linear space or
%              log10 space.
% Input  : - Size of random number matrix [M by N].
%          - range of random values [Min Max].
%          - Uniform numbers in {'lin' | 'log'} space, default is 'lin'.
% Output : - Matrix of M by N random numbers.
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Oct 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Rand=rand_range([100000 1],[1 5]);
%          Rand=rand_range([100000 1],[10 1000],'log');
% Reliable: 2
%--------------------------------------------------------------------------

Def.Type = 'lin';
if (nargin==2),
   Type = Def.Type;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end


Rand = rand(Size);
switch lower(Type)
 case {'lin','linear'}
    Rand = Rand.*range(Range);
    Rand = Rand + min(Range);
 case 'log'
    Rand = Rand.*range(log10(Range));
    Rand = Rand + log10(min(Range));
    Rand = 10.^Rand;
 otherwise
    error('Unknown Type option');
end    
