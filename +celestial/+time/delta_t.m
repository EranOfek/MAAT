function DeltaT=delta_t(JD)
% Return \Delta{T}
% Package: celestial.time
% Description: Return \Delta{T} at a vector of Julian days.
%              DeltaT is defined as ET-UT prior to 1984, and 
%              TT-UT1 after 1984 (= 32.184+(TAI-UTC)-(UT1-UTC)).
% Input  : - Vector of JD.
% Output : - DeltaT [seconds].
%            Return NaN if DeltaT is not available (out of range).
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: DeltaT=celestial.time.delta_t([2451545;celestial.time.julday])
% Reliable: 2
%--------------------------------------------------------------------------


Ranges = [-390 948 1620 1800 1900 1961].';
YearJD = convert.time(Ranges,'1','JD');  % 1st Jan of the year

%RangeInd = assoc_range(JD,YearJD);
DeltaT = zeros(size(JD)).*NaN;

%--- 390BC to 948AD ---
Flag = JD>YearJD(1) & JD<=YearJD(2);
Theta = (JD(Flag)-convert.date2jd([1 1 1800]))./36525;
DeltaT(Flag) = 1360 + 320.*Theta + 44.3.*Theta.^2;

%--- 948 to 1620 ---
Flag = JD>YearJD(2) & JD<=YearJD(3);
Theta = (JD(Flag)-convert.date2jd([1 1 1800]))./36525;
DeltaT(Flag) = 25.5.*Theta.^2;

%--- 1620 to 1800 ---
% max error 0.76s
% polynomial fit to Table 10.A in Meeus (1991)
Flag = JD>YearJD(3) & JD<=YearJD(4);
Theta = (JD(Flag)-convert.date2jd([1 1 1700]))./36525;
PolyPar1620_1800 = [3592.368  -4565.475  -6894.487   9653.180   4766.756  -7916.634  -1138.191   2991.339    -63.892   -551.212    124.603      7.595      7.115];
DeltaT(Flag) = polyval(PolyPar1620_1800,Theta);

%--- 1800 to 1899 ---
% max error 0.9s
% Meeus (1991)
Flag = JD>YearJD(4) & JD<=YearJD(5);
Theta = (JD(Flag)-convert.date2jd([1 1 1900]))./36525;
PolyPar1800_1899 = [123563.95 727058.63 1818961.41 2513807.78 2087298.89 1061660.75 324011.78 56282.84 5218.61 228.95 -2.50];
DeltaT(Flag) = polyval(PolyPar1800_1899,Theta);


%--- 1900 to 1997 ---
% max error 0.9s
% Meeus (1991)
Flag = JD>YearJD(5) & JD<=YearJD(6);
Theta = (JD(Flag)-convert.date2jd([1 1 1900]))./36525;
PolyPar1900_1997 = [58353.42 -232424.66 372919.88 -303191.19 124906.15 -18756.33 -2637.80 815.20 87.24 -2.44];
DeltaT(Flag) = polyval(PolyPar1900_1997,Theta);

%--- after 1961 ---
% exact
Flag = JD>YearJD(6);
TAImUTC = celestial.time.tai_utc(JD(Flag));
UT1mUTC = celestial.time.ut1_utc(JD(Flag));
DeltaT(Flag) = 32.184 + TAImUTC - UT1mUTC;
