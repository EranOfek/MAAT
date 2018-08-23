function Map=catalog_mapping(CatName,MapName)
% Mapping of VizieR catalogs columns
% Package: VO.VizieR
% Description: Get the VizieR catalogs mapping indicating the format
%              (column names, units and location) for each catalog.
% Input  : - Catalog program name. If not provided than return all
%            catalogs.
%          - Map name. If more than one mapping exist per catalog than this
%             can be used to identify the specific map. Default is empty.
% Output : - A structure array containing the column mapping information
%            for the all catalogs.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Map=VO.VizieR.catalog_mapping
% Reliable: 2



if (nargin<2)
    MapName = [];
end

I = 0;
%--- UCAC4 ---
I = I + 1;
Map(I).CatName = 'finducac4';
Map(I).MapName = 'regular';
Map(I).Cols    = {'Name','RA',  'Dec',   'RAErr','DecErr','Err',  'EpochRA','EpochDec','Fmag', 'Amag', 'MagErr','muRA',  'muDec',  'muRAErr', 'muDecErr'};
Map(I).Pos     = {[1 10],[12 22],[23 33],[35 37],[39 41], [43 46],[48 54],  [56 62],   [64 69],[71 76],[78 81], [102 109],[111 118],[120 123],[125 128]'};
Map(I).Format  = {'%s',  '%f',   '%f',   '%f',   '%f',    '%f',   '%f',     '%f',      '%f',   '%f',   '%f',   '%f',     '%f',     '%f',      '%f'};
Map(I).Units   = {''    ,'deg',  'deg',  'mas',  'mas',   'mas',  'jyear',  'jyear',   'mag',  'mag',  'mag',   'mas/yr','mas/yr', 'mas/yr',  'mas/yr'};
   
%#UCAC4    |    RA  (ICRS) Dec     +/- +/-  mas  EpRA    EpDE  | f.mag  a.mag  +/- |of db| Na  Nu  Cu| pmRA(mas/yr)pmDE  +/-  +/-|MPOS1      UCAC2      Tycho-2    |     2Mkey   Jmag  +/-:cq   Hmag  +/-:cq   Kmag  +/-:cq|  Bmag:1s.   Vmag:1s.   gmag:1s.   rmag:1s.   imag:1s.|gc.HAbhZBLNS|LED 2MX|;     r(")
%          1        2          3         4        5         6         7         8          9        10        11        12
%123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
%374-010653|089.8424421-15.2051623  63  60  134 1982.90 1980.76|15.963 15.875  --- | 0  0|  1   1   3|    -1.7     -8.6  4.0  4.0|084056380                        |1326243385 13.972 0.03:05 13.336 0.03:05 13.173 0.04:05|  --- :      --- :      --- :      --- :      --- :   |04.000000011|  0   0|;    547.67
    
%--- 2MASS ---
I = I + 1;
Map(I).CatName = 'find2mass';
Map(I).MapName = 'regular';
Map(I).Cols    = {'RA',  'Dec',  'MagJ', 'MagJErr', 'MagH', 'MagHErr','MagK',  'MagKErr'};
Map(I).Pos     = {[1 10],[12 21],[55 60],[62 66],   [85 90],[92 96],  [115 120],[122 126]};
Map(I).Format  = {'%f',  '%f',   '%f',   '%f',      '%f',   '%f',     '%f',     '%f'};
Map(I).Units   = {'deg', 'deg',  'mag',  'mag',     'mag',   'mag',   'mag',   'mag'};
   
%#.RAdeg      DEdeg   |errM errm ePA|   2MASS_name    | Jmag   +/-  m.e.        S/N | Hmag   +/-  m.e.        S/N | Kmag   +/-  m.e.        S/N |Qfl Rfl Bfl Cfl nJnHnK| prox  PA| proxCntr  X A      Cntr  H|  ObsDate Scan|  Glon    Glat | Xscan MeasureJD    | Jchi   Hchi   Kchi| Japm  +/-  | Hapm  +/-  | Kapm  +/-  |edgeN  EW ..|d u|O rO"  PA| Bmag  Rmag N| Extnd# Scan#  CoAdd#  ##
%         1         2         3         4         5         6         7
%123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
%179.918265 -00.139934|0.10 0.09  85|11594038-0008237 |15.975 0.066 0.067       16.9|15.426 0.095 0.096       13.8|15.023 0.141 0.141        8.0|AAB 222 111 000 361606| 61.2 164| 222867512 0 0  222867542 s|1999-01-14  93|2

%--- USNO-B1 ---
I = I + 1;
Map(I).CatName = 'findusnob1';
Map(I).MapName = 'regular';
Map(I).Cols    = {'RA',   'Dec',  'RAErr', 'DecErr', 'Epoch', 'muRA', 'muDec', 'MagB1', 'MagR1',  'MagB2',  'MagR2',  'MagI'};
Map(I).Pos     = {[27 36],[37 46],[48 50], [52 54],  [56 61], [63 68],[70 75], [97 102],[129 133],[160 164],[191 195],[222 226]};
Map(I).Format  = {'%f',  '%f',   '%f',   '%f',      '%f',   '%f',     '%f',     '%f',   '%f',      '%f',    '%f',     '%f'};
Map(I).Units   = {'deg', 'deg',  'mas',  'mas',     'jyear',   'mas/yr','mas/yr','mag',  'mag',   'mag',    'mag',    'mag'};
   
%#USNO-B1.0   Tycho-2        RA  (J2000) Dec    sRA sDE  Epoch   pmRA   pmDE P spA spD Fit N MsY| Bmag1 C Surv. cl <-xi-><-eta>| Rmag1 C Surv. cl <-xi-><-eta>| Bmag2 C Surv. cl <-xi-><-eta>| Rmag2 C Surv. cl <-xi-><-eta>|  Imag C Surv. cl <-xi-><-eta>| ;     r(")
%          1        2          3         4        5         6         7          8          9        10        11        12       13        14         15
%123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
%0899-0206737              179.840803-00.005959 999 590 1975.8     +0     +0 0   0   0 9 8 2 ...|  ---  - --    --             | 18.75 3 1-615  8 +02.22+00.59| 21.90 2 2-859  1 -02.23-00.60|  ---  - --    --             |  ---  - --    --             | ;    573.51
    
%--- WISE ---
I = I + 1;
Map(I).CatName = 'findwise';
Map(I).MapName = 'regular';
Map(I).Cols    = {'RA',   'Dec',  'W1mag', 'W2mag', 'W3mag',   'W4mag',  'e_W1mag', 'e_W2mag', 'e_W3mag',  'e_W4mag'};
Map(I).Pos     = {[24 34],[41 51],[81 86], [89 94],  [97 102], [104 109],[137 141], [145 149],[153 157],[161 165]};
Map(I).Format  = {'%f',  '%f',   '%f',   '%f',      '%f',   '%f',     '%f',     '%f',   '%f',      '%f'};
Map(I).Units   = {'deg', 'deg',  'mag',  'mag',     'mag',   'mag',   'mag',     'mag',     'mag',      'mag'};
   
%#WISE_ALLSKY    RAdeg   DEdeg   errHalfMaj      errHalfMin      errPosAng       W1mag   W2mag   W3mag   W4mag   Jmag    Hmag    Kmag    e_W1mag e_W2mag e_W3mag e_W4mag e_Jmag  e_Hmag  e_Kmag  ID      snr1    chi2W1  snr2    chi2W2  snr3    chi2W3  snr4    chi2W4  nb      na      sat1    sat2    sat3    sat4    ccf     ex      var     nW1     mW1     nW2     mW2     nW3     mW3     nW4     mW4     2Mkey   d2M
%          1        2          3         4        5         6         7          8          9        10        11        12       13        14         15
%123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
%J060004.30-150101.9     90.017928       -15.017205      0.622   0.566   100     16.568  16.833  12.090  9.253                           0.105   0.398   0.247   0.513                           907015201241025397      10.3    9.962e-01       2.7     8.545e-01       4.4     1.041e+00       2.1     9.008e-01       1       0       0.000   0.000   0.000   0.000   0000    0       2nnn    0       7       15      0       15      0       14      0       0
%J060000.89-150120.6     90.003726       -15.022402      0.831   0.755   100     16.883  16.888  12.743  9.286                           0.129   0.466                                           907015201241025540      8.4     9.264e-01       2.3     7.365e-01       0.2     8.901e-01       -0.9    1.065e+00       1       0       0.000   0.000   0.000   0.000   0000    0       nnnn    0       3       14      1       14      1       13      0       0
%J034911.52-571833.6     57.2980194      -57.3093561  
%J034911.52-571833.6     57.2980194      -57.3093561     0.2
%J034911.52-571833.6     57.2980194      -57.3093561     0.2374  

% search catalog
if (nargin>0)
    if (~isempty(MapName))
        Icat = find(strcmp({Map.CatName},CatName) & strcmp({Map.MapName},MapName));
    else
        Icat = find(strcmp({Map.CatName},CatName));
    end
    if (isempty(Icat))
        error('Catalog %s not found - you may need to add the catalog mapping to function',CatName);
    end
    Map  = Map(Icat);
end
