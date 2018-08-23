function Bit=def_bitmask_specpipeline(Field)
%--------------------------------------------------------------------------
% def_bitmask_specpipeline function                                 ImSpec
% Description: The spectroscopic pipeline bit mask definition.
%              Given the Bit mask name return the bit mask index.
% Input  : - Bit mask name.
% Output : - Bit index.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: maskflag_check.m, maskflag_set.m
% Example: Bit=def_bitmask_specpipeline('Bit_ImSaturated');
% Reliable: 2
%--------------------------------------------------------------------------

Ind = 0;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_ImSaturated';   % in raw image (1 - 1)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_BiasNoisy';     % in bias (2 - 2)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_BiasNoise0';    % std in bias is zero (bad pixel?) (3 - 4)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_BadPixel';      % flat ratio (4 - 8) 
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_FlatNaN';       % flat (5 - 16)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_FlatLowNim';    % small number of flat input images (6 - 32)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_MaxRelErr';     % flat large std (7-64)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_Divide0';       % need a new flag for low counts flat (28-128)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_CR';            % CR (9-256)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_TraceDev';      % large deviation between trace and fitted trace (10-512)
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_SpecHighBack';  % high back in 2D spec
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_WaveCalibDev';  % 
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_ResponseOutRange';     % 
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_Telluric';     % 
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_SkyEmissionLine';     % 
Map(Ind).Ind  = Ind;
Ind = Ind + 1;
Map(Ind).Name = 'Bit_StitchRegion';     % 
Map(Ind).Ind  = Ind;


Ind = find(strcmp({Map.Name},Field));
Bit = Map(Ind).Ind;
