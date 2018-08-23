function Data=get_scattering_cross_section(Z)
% Get the attenuation and scattering cross-sections from the NIST website.
% Package: VO.NIST
% Description: Get the attenuation and scattering cross-sections from the
%              NIST website. Calculate opacity and cross section for
%              bound-free absorption.
% Input  : - Atomic number (Z).
% Output : - Structure with the following vectors.
%            .E - Energy [keV] (in the range ~0.01 to 430 keV).
%            .Kappa - Opacity [cm^2 g^-1]
%            .Sigma - Cross section per atom [cm^-2].
% Reference: https://physics.nist.gov/PhysRefData/FFast/html/form.html
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=VO.NIST.get_scattering_cross_section(1)
% Reliable: 2
%--------------------------------------------------------------------------


URL_low  = sprintf('https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z=%d&Formula=&gtype=2&range=L&lower=&upper=&density=&frames=no&htmltable=1',...
    Z);
URL_up   = sprintf('https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z=%d&Formula=&gtype=2&range=U&lower=&upper=&density=&frames=no&htmltable=1',...
    Z);

Str = webread(URL_low);
SearchStr = 'barns/atom) = [<i>&mu;</i>/<i>&rho;</i>](cm<sup>2</sup> g<sup>-1</sup>) &#160;&#215;&#160; ';
%D=regexp(Str,'barns/atom) = [<i>&mu;</i>/<i>&rho;</i>](cm<sup>2</sup> g<sup>-1</sup>) &#160;&#215;&#160; (?<val>\w+) ','names');
I = strfind(Str,SearchStr);
SubStr = Str(I+numel(SearchStr):I+numel(SearchStr)+20);
Iend = strfind(SubStr,'<dd>');
ConvLow = str2double(SubStr(1:Iend-1));

Str = webread(URL_up);
SearchStr = 'barns/atom) = [<i>&mu;</i>/<i>&rho;</i>](cm<sup>2</sup> g<sup>-1</sup>) &#160;&#215;&#160; ';
%D=regexp(Str,'barns/atom) = [<i>&mu;</i>/<i>&rho;</i>](cm<sup>2</sup> g<sup>-1</sup>) &#160;&#215;&#160; (?<val>\w+) ','names');
I = strfind(Str,SearchStr);
SubStr = Str(I+numel(SearchStr):I+numel(SearchStr)+20);
Iend = strfind(SubStr,'<dd>');
ConvUp = str2double(SubStr(1:Iend-1));



C_low = www.parse_html_table(URL_low);
C_low{2}=regexprep(C_low{2},'&nbsp;&nbsp;','');

Low = [str2double(C_low{1}),str2double(C_low{2})];

C_up = www.parse_html_table(URL_up);
C_up{2}=regexprep(C_up{2},'&nbsp;&nbsp;','');

Up  = [str2double(C_up{1}),str2double(C_up{2})];

Kappa      = [Low; Up];        % (energy [keV], opacity [cm^2 g^-1])
Data.E     = Kappa(:,1);       % [keV]
Data.Kappa = Kappa(:,2);       % cm^2 gr
Data.Sigma = [ConvLow.*Low(:,2); ConvUp.*Up(:,2)].*1e-24;  % cm^-2
