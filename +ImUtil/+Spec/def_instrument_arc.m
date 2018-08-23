function [ArcKey,ArcKeyValue,ArcName]=def_instrument_arc(Instrum)
%--------------------------------------------------------------------------
% def_instrument_arc function                                      ImBasic
% Description: Definitions of how to identify spectroscopic arc accoring
%              to the image header information. Given an instrument,
%              this function return a header keyword name that contains
%              the arc information, as well as the values of this keyword
%              if the image is an arc image, and the translation of these
%              various values to standard arc name.
% Input  : - Instrument name. The following instruments are supported.
%            'P200-DBSP'  - Palomar 200" DBSP.
% Output : - The image header keyword to use in order to identify arcs.
%          - A cell array of the values of this keyword if the image is
%            an arc image.
%          - A cell array of translation between the values of the keyword
%            and standard arc name.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: is_arc_image.m
% Example: [ArcKey,ArcKeyValue,ArcName]=def_instrument_arc('P200-DBSP');
% Reliable: 2
%--------------------------------------------------------------------------

switch lower(Instrum)
    case 'p200-dbsp'
        ArcKey      = 'LAMPS';
        ArcKeyValue = {'1000000','0100000','0010000','0001000','0000100','0000010','0000001','0001100','0001110'};
        ArcName     = {'D',      'Fe+Ar',  'Hg',     'Ar',     'Ne',     'He',     'InCand' ,'Ne+Ar'  ,'Ne+Ar'};
    otherwise
        error('unknown Instrum option');
end

