function [Color,Err]=star_sptype_color(SpType,LumClass,Fam1,Filt1,Sys1,Fam2,Filt2,Sys2)
%--------------------------------------------------------------------------
% star_sptype_color function                                     AstroSpec
% Description: Given a star spectral type and luminosity class, get the star
%              color between any two filters.
% Input  : - Spectral type (e.g., 'A4').
%          - Luminosity class (e.g., 'V').
%          - The familiy name for the first filter (e.g., 'SDSS').
%          - The filter name for the first filter (e.g., 'g').
%          - Magnitude system for the first filter (e.g., 'AB').
%          - The familiy name for the second filter (e.g., 'SDSS').
%          - The filter name for the second filter (e.g., 'r').
%          - Magnitude system for the second filter (e.g., 'AB').
% Output : - Color.
%          - Rough interpolation error.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Oct 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [C,E]=star_sptype_color('A3','IV','SDSS','g','AB','SDSS','r','AB');
%          [C,E]=star_sptype_color('A3','IV','Johnson','B','Vega','Johnson','V','Vega');
% Reliable: 2
%--------------------------------------------------------------------------


Letter    = upper(SpType(1));
SubLetter = str2double(SpType(2:end));
LumClass  = upper(LumClass);

Dic.Class = {'O','B','A','F','G','K','M'};
Dic.Num   = [0   10  20  30  40  50  60];

% convert spectral type to running index
SpTypeInddex = Dic.Num(find(isempty_cell(strfind(Dic.Class,Letter))==0))+SubLetter;


% Dictionary of existing spectra:
I = 1;
Exist.V.SpType{I}='O5'; Exist.V.Num(I)=  5;
I = I + 1;
Exist.V.SpType{I}='O9'; Exist.V.Num(I)=  9;
I = I + 1;
Exist.V.SpType{I}='B0'; Exist.V.Num(I)= 10;
I = I + 1;
Exist.V.SpType{I}='B1'; Exist.V.Num(I)= 11;
I = I + 1;
Exist.V.SpType{I}='B3'; Exist.V.Num(I)= 13;
I = I + 1;
Exist.V.SpType{I}='B8'; Exist.V.Num(I)= 18;
I = I + 1;
Exist.V.SpType{I}='B9'; Exist.V.Num(I)= 19;
I = I + 1;
Exist.V.SpType{I}='A0'; Exist.V.Num(I)= 20;
I = I + 1;
Exist.V.SpType{I}='A2'; Exist.V.Num(I)= 22;
I = I + 1;
Exist.V.SpType{I}='A3'; Exist.V.Num(I)= 23;
I = I + 1;
Exist.V.SpType{I}='A5'; Exist.V.Num(I)= 25;
I = I + 1;
Exist.V.SpType{I}='A7'; Exist.V.Num(I)= 27;
I = I + 1;
Exist.V.SpType{I}='F0'; Exist.V.Num(I)= 30;
I = I + 1;
Exist.V.SpType{I}='F2'; Exist.V.Num(I)= 32;
I = I + 1;
Exist.V.SpType{I}='F5'; Exist.V.Num(I)= 35;
I = I + 1;
Exist.V.SpType{I}='F6'; Exist.V.Num(I)= 36;
I = I + 1;
Exist.V.SpType{I}='F8'; Exist.V.Num(I)= 38;
I = I + 1;
Exist.V.SpType{I}='G0'; Exist.V.Num(I)= 40;
I = I + 1;
Exist.V.SpType{I}='G2'; Exist.V.Num(I)= 42;
I = I + 1;
Exist.V.SpType{I}='G5'; Exist.V.Num(I)= 45;
I = I + 1;
Exist.V.SpType{I}='G8'; Exist.V.Num(I)= 48;
I = I + 1;
Exist.V.SpType{I}='K0'; Exist.V.Num(I)= 50;
I = I + 1;
Exist.V.SpType{I}='K2'; Exist.V.Num(I)= 52;
I = I + 1;
Exist.V.SpType{I}='K3'; Exist.V.Num(I)= 53;
I = I + 1;
Exist.V.SpType{I}='K4'; Exist.V.Num(I)= 54;
I = I + 1;
Exist.V.SpType{I}='K5'; Exist.V.Num(I)= 55;
I = I + 1;
Exist.V.SpType{I}='K7'; Exist.V.Num(I)= 57;
I = I + 1;
Exist.V.SpType{I}='M0'; Exist.V.Num(I)= 60;
I = I + 1;
Exist.V.SpType{I}='M2'; Exist.V.Num(I)= 62;
I = I + 1;
Exist.V.SpType{I}='M3'; Exist.V.Num(I)= 63;
I = I + 1;
Exist.V.SpType{I}='M4'; Exist.V.Num(I)= 64;
I = I + 1;
Exist.V.SpType{I}='M5'; Exist.V.Num(I)= 65;
I = I + 1;
Exist.V.SpType{I}='M6'; Exist.V.Num(I)= 66;

I = 1;
Exist.IV.SpType{I}='B2'; Exist.IV.Num(I)= 12;
I = I + 1;
Exist.IV.SpType{I}='B6'; Exist.IV.Num(I)= 16;
I = I + 1;
Exist.IV.SpType{I}='A0'; Exist.IV.Num(I)= 20;
I = I + 1;
Exist.IV.SpType{I}='F5'; Exist.IV.Num(I)= 35;
I = I + 1;
Exist.IV.SpType{I}='F8'; Exist.IV.Num(I)= 38;
I = I + 1;
Exist.IV.SpType{I}='G0'; Exist.IV.Num(I)= 40;
I = I + 1;
Exist.IV.SpType{I}='G2'; Exist.IV.Num(I)= 42;
I = I + 1;
Exist.IV.SpType{I}='G5'; Exist.IV.Num(I)= 45;
I = I + 1;
Exist.IV.SpType{I}='G8'; Exist.IV.Num(I)= 48;
I = I + 1;
Exist.IV.SpType{I}='K0'; Exist.IV.Num(I)= 50;
I = I + 1;
Exist.IV.SpType{I}='K1'; Exist.IV.Num(I)= 51;
I = I + 1;
Exist.IV.SpType{I}='K3'; Exist.IV.Num(I)= 53;



I = 1;
Exist.III.SpType{I}='O8'; Exist.III.Num(I)=  8;
I = I + 1;
Exist.III.SpType{I}='B3'; Exist.III.Num(I)= 13;
I = I + 1;
Exist.III.SpType{I}='B5'; Exist.III.Num(I)= 15;
I = I + 1;
Exist.III.SpType{I}='B9'; Exist.III.Num(I)= 19;
I = I + 1;
Exist.III.SpType{I}='A0'; Exist.III.Num(I)= 20;
I = I + 1;
Exist.III.SpType{I}='A3'; Exist.III.Num(I)= 23;
I = I + 1;
Exist.III.SpType{I}='A5'; Exist.III.Num(I)= 25;
I = I + 1;
Exist.III.SpType{I}='A7'; Exist.III.Num(I)= 27;
I = I + 1;
Exist.III.SpType{I}='F0'; Exist.III.Num(I)= 30;
I = I + 1;
Exist.III.SpType{I}='F2'; Exist.III.Num(I)= 32;
I = I + 1;
Exist.III.SpType{I}='F5'; Exist.III.Num(I)= 35;
I = I + 1;
Exist.III.SpType{I}='G0'; Exist.III.Num(I)= 40;
I = I + 1;
Exist.III.SpType{I}='G5'; Exist.III.Num(I)= 45;
I = I + 1;
Exist.III.SpType{I}='G8'; Exist.III.Num(I)= 48;
I = I + 1;
Exist.III.SpType{I}='K0'; Exist.III.Num(I)= 50;
I = I + 1;
Exist.III.SpType{I}='K1'; Exist.III.Num(I)= 51;
I = I + 1;
Exist.III.SpType{I}='K2'; Exist.III.Num(I)= 52;
I = I + 1;
Exist.III.SpType{I}='K3'; Exist.III.Num(I)= 53;
I = I + 1;
Exist.III.SpType{I}='K4'; Exist.III.Num(I)= 54;
I = I + 1;
Exist.III.SpType{I}='K5'; Exist.III.Num(I)= 55;
I = I + 1;
Exist.III.SpType{I}='M0'; Exist.III.Num(I)= 60;
I = I + 1;
Exist.III.SpType{I}='M1'; Exist.III.Num(I)= 61;
I = I + 1;
Exist.III.SpType{I}='M2'; Exist.III.Num(I)= 62;
I = I + 1;
Exist.III.SpType{I}='M3'; Exist.III.Num(I)= 63;
I = I + 1;
Exist.III.SpType{I}='M4'; Exist.III.Num(I)= 64;
I = I + 1;
Exist.III.SpType{I}='M5'; Exist.III.Num(I)= 65;
I = I + 1;
Exist.III.SpType{I}='M6'; Exist.III.Num(I)= 66;
I = I + 1;
Exist.III.SpType{I}='M7'; Exist.III.Num(I)= 67;
I = I + 1;
Exist.III.SpType{I}='M8'; Exist.III.Num(I)= 68;
I = I + 1;
Exist.III.SpType{I}='M9'; Exist.III.Num(I)= 69;





I = 1;
Exist.II.SpType{I}='B2'; Exist.II.Num(I)= 12;
I = I + 1;
Exist.II.SpType{I}='B5'; Exist.II.Num(I)= 15;
I = I + 1;
Exist.II.SpType{I}='F0'; Exist.II.Num(I)= 30;
I = I + 1;
Exist.II.SpType{I}='F2'; Exist.II.Num(I)= 32;
I = I + 1;
Exist.II.SpType{I}='G5'; Exist.II.Num(I)= 45;
I = I + 1;
Exist.II.SpType{I}='M3'; Exist.II.Num(I)= 63;



I = 1;
Exist.I.SpType{I}='B0'; Exist.I.Num(I)= 10;
I = I + 1;
Exist.I.SpType{I}='B1'; Exist.I.Num(I)= 11;
I = I + 1;
Exist.I.SpType{I}='B3'; Exist.I.Num(I)= 13;
I = I + 1;
Exist.I.SpType{I}='B5'; Exist.I.Num(I)= 15;
I = I + 1;
Exist.I.SpType{I}='B8'; Exist.I.Num(I)= 18;
I = I + 1;
Exist.I.SpType{I}='A0'; Exist.I.Num(I)= 20;
I = I + 1;
Exist.I.SpType{I}='A2'; Exist.I.Num(I)= 22;
I = I + 1;
Exist.I.SpType{I}='F0'; Exist.I.Num(I)= 30;
I = I + 1;
Exist.I.SpType{I}='F5'; Exist.I.Num(I)= 35;
I = I + 1;
Exist.I.SpType{I}='F8'; Exist.I.Num(I)= 38;
I = I + 1;
Exist.I.SpType{I}='G0'; Exist.I.Num(I)= 40;
I = I + 1;
Exist.I.SpType{I}='G2'; Exist.I.Num(I)= 42;
I = I + 1;
Exist.I.SpType{I}='G5'; Exist.I.Num(I)= 45;
I = I + 1;
Exist.I.SpType{I}='G8'; Exist.I.Num(I)= 48;
I = I + 1;
Exist.I.SpType{I}='K2'; Exist.I.Num(I)= 52;
I = I + 1;
Exist.I.SpType{I}='K3'; Exist.I.Num(I)= 53;
I = I + 1;
Exist.I.SpType{I}='K4'; Exist.I.Num(I)= 54;
I = I + 1;
Exist.I.SpType{I}='M2'; Exist.I.Num(I)= 62;



I1=find(Exist.(LumClass).Num>=SpTypeInddex,1);
I2=find(Exist.(LumClass).Num<=SpTypeInddex,1,'last');

if (isempty(I1) || isempty(I2)),
   Color = NaN;
   Err   = NaN;
else
   if (I1==I2),
      % no interpolation needed
      SpName = sprintf('uk%s%s.mat',lower(Exist.(LumClass).SpType{I1}),lower(LumClass));
      Spec   = load2(SpName);
      Color  = synphot(Spec,Fam1,Filt1,Sys1) - synphot(Spec,Fam2,Filt2,Sys2);
      Err    = 0;
   else
      SpName1 = sprintf('uk%s%s.mat',lower(Exist.(LumClass).SpType{I1}),lower(LumClass));
      SpName2 = sprintf('uk%s%s.mat',lower(Exist.(LumClass).SpType{I2}),lower(LumClass));
      % weights
      D1 = abs(SpTypeInddex - Exist.(LumClass).Num(I1));
      D2 = abs(SpTypeInddex - Exist.(LumClass).Num(I2));
      W1 = (D1+D2)./D1;
      W2 = (D1+D2)./D2;
   
      Spec1   = load2(SpName1);
      Spec2   = load2(SpName2);
      % linearly interpolate color
      Color1  = synphot(Spec1,Fam1,Filt1,Sys1) - synphot(Spec1,Fam2,Filt2,Sys2);
      Color2  = synphot(Spec2,Fam1,Filt1,Sys1) - synphot(Spec2,Fam2,Filt2,Sys2);
      Color   = (Color1.*W1 + Color2.*W2)./(W1+W2);
      Err     = abs(Color1-Color2).*0.5; 
   
   end
end
