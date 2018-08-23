function DB=meteors_db
% Return a meteor shower database (incomplete)
% Package: celestial.meteors
% Description: Return a meteor showers database (not complete).
% Input  : null
% Output : - A structure array containing the meteor showers DB.
%            The following fields are available.
%            .Name
%            .Code
%            .CommentZHR
%            .Comet
%            .Data - a vector containing
%                    begin date day
%                    begin date month
%                    end date day
%                    end date day
%                    peak date day
%                    peak date day
%                    solar long. [deg]
%                    J2000 RA [deg]
%                    J2000 Dec [deg]
%                    Velocity [km/s]
%                    r
%                    ZHR
%                    radiant daily motion RA
%                    radiant daily motion Dec
% Tested : Matlab 7.13
%     By : Eran Ofek                       Aug 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: DB=celestial.meteors.meteors_db;
% Reliable: 2
%--------------------------------------------------------------------------

I=0;

I=I+1;
DB(I).Name='Quadrantids          ';
DB(I).Code='QUA';
DB(I).CommentZHR='vary 60-200';
DB(I).Comet='';
DB(I).Data=[ 01 01 01 05  01 03 283.16 230. +49. 41 2.1 120   3.2 -0.2];

I=I+1;
DB(I).Name='delta Cancrids       ';
DB(I).Code='DCA';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 01 01 01 24  01 17 297.   130. +20. 28 3.0 4    NaN NaN];

I=I+1;
DB(I).Name='alpha Centaurids     ';
DB(I).Code='ACE';
DB(I).CommentZHR='may reach 25+';
DB(I).Comet='';
DB(I).Data=[ 01 28 02 21  02 07 319.2  210. -59. 56 2.0 6   NaN  NaN];

I=I+1;
DB(I).Name='delta Leonids        ';
DB(I).Code='DLE';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 02 15 03 10  02 24 336.   168. +16. 23 3.0 2   NaN NaN];

I=I+1;
DB(I).Name='gamma Normids        ';
DB(I).Code='GNO';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 02 25 03 22  03 13 353.   249. -51. 56 2.4 8  NaN NaN];

I=I+1;
DB(I).Name='Virginids            ';
DB(I).Code='VIR';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 01 25 04 15  03 24   4.   195. -04. 30 3.0 5  +2.0 -0.3];

I=I+1;
DB(I).Name='Lyrids               ';
DB(I).Code='LYR';
DB(I).CommentZHR='up to 90';
DB(I).Comet='';
DB(I).Data=[ 04 16 04 25  04 22 032.32 271. +34. 49 2.1 18  +4.4 0.0];

I=I+1;
DB(I).Name='pi-Puppids           ';
DB(I).Code='PPU';
DB(I).CommentZHR='13-40';
DB(I).Comet='26P/Grigg-Skjellerup';
DB(I).Data=[ 04 15 04 28  04 23 033.5  110. -45. 18 2.0 2   NaN NaN];

I=I+1;
DB(I).Name='eta Aquarids         ';
DB(I).Code='ETA';
DB(I).CommentZHR='40-85';
DB(I).Comet='1P/Halley';
DB(I).Data=[ 04 19 05 28  05 05 045.5  338. -01. 66 2.4 60    3.6 +0.4];

I=I+1;
DB(I).Name='Sagittarids          ';
DB(I).Code='SAG';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 04 15 07 15  05 19  59.   247. -22. 30 2.5 5   NaN NaN];

I=I+1;
DB(I).Name='June Bootids         ';
DB(I).Code='JBO';
DB(I).CommentZHR='0-100+';
DB(I).Comet='7P/Pons-Winnecke';
DB(I).Data=[ 06 26 07 02  06 27 095.7  224. +48. 18 2.2  0  1.6 -0.2];

I=I+1;
DB(I).Name='Pegasids             ';
DB(I).Code='JPE';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 07 07 07 13  07 09 107.5  340. +15. 70 3.0 3   +3.2 +0.2];

I=I+1;
DB(I).Name='July Phoenicids      ';
DB(I).Code='PHE';
DB(I).CommentZHR='3-10';
DB(I).Comet='';
DB(I).Data=[ 07 10 07 16  07 13 111.    32. -48. 47 3.0 3   NaN NaN];

I=I+1;
DB(I).Name='Piscis Austrinids    ';
DB(I).Code='PAU';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 07 15 08 10  07 28 125.   341. -30. 35 3.2 5   3.4 0.4];

I=I+1;
DB(I).Name='South delta Aquarids ';
DB(I).Code='SDA';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 07 12 08 19  07 28 125.   339. -16. 41 3.2 20   3.0 +0.2];

I=I+1;
DB(I).Name='alpha Capricornids   ';
DB(I).Code='CAP';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 07 03 08 15  07 30 127.   307. -10. 23 2.5 4    3.6 +0.3];

I=I+1;
DB(I).Name='South iota Aquarids  ';
DB(I).Code='SIA';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 07 25 08 15  08 04 132.   334. -15. 34 2.9 2   4.3 +0.2];

I=I+1;
DB(I).Name='North delta Aquarids ';
DB(I).Code='NDA';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 07 15 08 25  08 08 136.   335. -05. 42 3.4 4    3.0 +0.2];

I=I+1;
DB(I).Name='Perseids             ';
DB(I).Code='PER';
DB(I).CommentZHR='';
DB(I).Comet='109P/Swift-Tuttle';
DB(I).Data=[ 07 17 08 24  08 12 140.    46. +58. 59 2.6 100    5.6 +0.2];

I=I+1;
DB(I).Name='kappa Cygnids        ';
DB(I).Code='KCG';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 08 03 08 25  08 17 145.   286. +59. 25 3.0 3      0.8 +0.1];

I=I+1;
DB(I).Name='North iota Aquarids  ';
DB(I).Code='NIA';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 08 11 08 31  08 19 147.   327. -06. 31 3.2 3    4.1 +0.1];

I=I+1;
DB(I).Name='alpha Aurigids       ';
DB(I).Code='AUR';
DB(I).CommentZHR='up to 40';
DB(I).Comet='';
DB(I).Data=[ 08 25 09 08  09 01 158.6   84. +42. 66 2.6 7   4.4 +0.0];

I=I+1;
DB(I).Name='delta aurigids       ';
DB(I).Code='DAU';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 09 05 10 10  09 09 166.7   60. +47. 64 2.9 5   4.0 +0.1];

I=I+1;
DB(I).Name='Piscids              ';
DB(I).Code='SPI';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 09 01 09 30  09 19 177.     5. -01. 26 3.0 3   3.6 +0.2];

I=I+1;
DB(I).Name='Draconids            ';
DB(I).Code='GIA';
DB(I).CommentZHR='20-500+';
DB(I).Comet='21P/Giacobini-Zinner';
DB(I).Data=[ 10 06 10 10  10 08 195.4  262. +54. 20 2.6 20   NaN  NaN];

I=I+1;
DB(I).Name='epsilon Geminids     ';
DB(I).Code='EGE';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 10 14 10 27  10 18 205.   102. +27. 70 3.0 2   NaN  NaN];

I=I+1;
DB(I).Name='Orionids             ';
DB(I).Code='ORI';
DB(I).CommentZHR='14-31';
DB(I).Comet='1P/Halley';
DB(I).Data=[ 10 02 11 07  10 21 208.    95. +16. 66 2.5 23  2.6 +0.1];

I=I+1;
DB(I).Name='Southern Taurids     ';
DB(I).Code='STA';
DB(I).CommentZHR='';
DB(I).Comet='2P/Encke';
DB(I).Data=[ 10 01 11 25  11 05 223.    52. +13. 27 2.3 5   3.2 +0.1];

I=I+1;
DB(I).Name='Northern Taurids     ';
DB(I).Code='NTA';
DB(I).CommentZHR='';
DB(I).Comet='2P/Encke';
DB(I).Data=[ 10 01 11 25  11 12 230.    58. +22. 29 2.3 5   3.0 +0.1];

I=I+1;
DB(I).Name='Leonids              ';
DB(I).Code='LEO';
DB(I).CommentZHR='10-50000';
DB(I).Comet='55P/Tempel-Tuttle';
DB(I).Data=[ 11 14 11 21  11 17 235.27 153. +22. 71 2.5 20   2.8 -0.4];
%Data(I,:)=[ 11 14 11 21  11 17 236.709 153. +22. 71 2.5 20   2.8 -0.4];  % 2006

I=I+1;
DB(I).Name='alpha Monocerotids   ';
DB(I).Code='AMO';
DB(I).CommentZHR='5-400+';
DB(I).Comet='';
DB(I).Data=[ 11 15 11 25  11 21 239.32 117. +01. 65 2.4 5   3.6 -0.2];

I=I+1;
DB(I).Name='Chi Orionids         ';
DB(I).Code='XOR';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 11 26 12 15  12 02 250.    82. +23. 28 3.0 3   NaN NaN];

I=I+1;
DB(I).Name='Phoenicids           ';
DB(I).Code='PHO';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 11 28 12 09  12 06 254.25  18. -53. 18 2.8 10   3.2 -0.2];

I=I+1;
DB(I).Name='Puppid/Velids        ';
DB(I).Code='PUP';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 12 01 12 15  12 07 255.   123. -45. 40 2.9 10   0.8 0.0];

I=I+1;
DB(I).Name='Monocerotids         ';
DB(I).Code='MON';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 11 27 12 17  12 09 257.   100. +08. 42 3.0 3  3.2 0.0];

I=I+1;
DB(I).Name='Hydrids              ';
DB(I).Code='HYD';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 12 03 12 15  12 12 260.   127. +02. 58 3.0 2   3.2 -0.2];

I=I+1;
DB(I).Name='Geminids             ';
DB(I).Code='GEM';
DB(I).CommentZHR='';
DB(I).Comet='3200 Phaethon';
DB(I).Data=[ 12 07 12 17  12 14 262.2  112. +33. 35 2.6 120   3.9 -0.1]; 

I=I+1;
DB(I).Name='Coma Berenicids      ';
DB(I).Code='COM';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 12 12 01 23  12 19 268.   175. +25. 65 3.0 5   3.6 -0.3];

I=I+1;
DB(I).Name='Ursids               ';
DB(I).Code='URS';
DB(I).CommentZHR='';
DB(I).Comet='';
DB(I).Data=[ 12 17 12 26  12 22 270.7  217. +76. 33 3.0 10    0.0 -0.4];

