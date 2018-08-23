function []=spec_write_products(Spec,varargin)
%--------------------------------------------------------------------------
% spec_write_products function                                      ImSpec
% Description: 
% Input  : - Structure array containing all the spectroscopic products to
%            write.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.WaveCalibType = {'arc','sky'};  % preference order
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

Nct   = numel(InPar.WaveCalibType);


Nsp = numel(Spec);
for Isp=1:1:Nsp,
    for Ict=1:1:Nct,
        if (isfield(Spec(Isp).Info.FitWave,InPar.WaveCalibType{Ict})),
           WavePrim = Spec(Isp).Info.FitWave.(InPar.WaveCalibType{Ict}).SpecWave;
        
    Spec(Isp).Info.FitWave.arc.SpecWave
    
    