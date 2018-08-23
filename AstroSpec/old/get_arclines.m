function List=get_arclines(ArcType);
%------------------------------------------------------------------------------
% get_arclines function                                              AstroSpec
% Description: Get list of arc/sky emission lines
%              Replaced by spec_get_arc.m
% Input  : - Arc type:
%            'Ar' - Ar arc [6965..10470]
%            'Cd  - Cd arc [3080..7345]
%            'Hg' - Hg arc [3125..5790]
%            'Ne' - Ne arc [5852..9665]
%            'Zn' - Zn arc [3018..7799]
%            'SA' - All sky lines [3142..10418]
%            'SS' - Selected sky lines [5577..10013]
% Output : - List of sky lines [Wavelength],
%            or [Wavelength, Intensity, FWHM, Flux].
% Tested : Matlab 7.0
% Needed : /data/AstroSpec/ArcLines directory.
%     By : Eran O. Ofek                     March 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: List=get_arclines('sa'); % Get all sky lines
% Reliable: 1
%------------------------------------------------------------------------------


switch lower(ArcType)
 case 'ar'
    [List] = textread('ArcLines_Ar.dat','%f\n','commentstyle','matlab');
 case 'cd'
    [List] = textread('ArcLines_Cd.dat','%f\n','commentstyle','matlab');
 case 'hg'
    [List] = textread('ArcLines_Hg.dat','%f\n','commentstyle','matlab');
 case 'ne'
    [List] = textread('ArcLines_Ne.dat','%f\n','commentstyle','matlab');
 case 'zn'
    [List] = textread('ArcLines_Zn.dat','%f\n','commentstyle','matlab');
 case 'sa'
    [W,Int,FWHM,F] = textread('ArcLines_SkyAll.dat','%f %f %f %f\n','commentstyle','matlab');
    List = [W,Int,FWHM, F];
 case 'ss'
    [List] = textread('ArcLines_SkySel.dat','%f\n','commentstyle','matlab');
 otherwise
    error('Unknown ArcType option');
end
