function AstS=convert_flux(AstS,OutUnits)
% Convert flux units of AstSpec class object.
% Package: @AstSpec
% Description: Convert flux units of AstSpec class object.
% Input  : - AstSpec class object.
%          - Output units. See convert_flux.m for options.
% Outout : - AstSpec class object in which the .Int, .Err
%            and .Back fields units are converted to the
%            new system.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstS=load_pickles; AstS(1).IntUnits='cgs/A';
%          AstS=convert_flux(AstS,'mJy')
% Reliable: 
%--------------------------------------------------------------------------

Ns = numel(AstS);
for Is=1:1:Ns,
    if (~isempty(AstS(Is).Wave) && ~isempty(AstS(Is).Int)),
        ConvFactor   = convert.flux(1,AstS(Is).IntUnits,OutUnits,AstS(Is).Wave,AstS(Is).WaveUnits);
        AstS(Is).Int = AstS(Is).Int .* ConvFactor;
        AstS(Is).IntUnits  = OutUnits;
        if (~isempty(AstS(Is).Err)),
            AstS(Is).Err = AstS(Is).Err .* ConvFactor;
        end
        if (~isempty(AstS(Is).Back)),
            AstS(Is).Back = AstS(Is).Back .* ConvFactor;
        end
    end
end
