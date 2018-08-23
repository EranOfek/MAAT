function [LC,Col,Head]=read_lc_file(FileName)
% Read Kepler light curve FITS file.
% Package: VO.Kepler
% Description: Read kepler light curve fits file.
% Input  : - String containing file name.
% Output : - Matrix containing light curve information.
%          - Structure containing columns list:
%            .FieldsCell - cell array of field names.
%            .FieldsStruct - structure containing column index of fields.
%          - Structure containing header information.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Jan 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2 (Modifications were not tested)
%------------------------------------------------------------------------------

Data = fitsread(FileName,'BinTable');
LC = cell2mat(Data);
Col.FieldsStruct.Time          = 1;   % Barycentric time  % [days] MJD = JD - 2400000.5
Col.FieldsStruct.BaryCorr      = 2;   % barycentric corrections [s]
Col.FieldsStruct.CadenceNumber = 3;   % cadence point number
Col.FieldsStruct.Row           = 4;   % row position [pix]
Col.FieldsStruct.SigmaRow      = 5;   % row position error [sd pix]
Col.FieldsStruct.Col           = 6;   % column position [pix]
Col.FieldsStruct.SigmaCol      = 7;   % column position error [sd pix]
Col.FieldsStruct.RawFlux       = 8;   % Aperture photometry raw flux [e-/cadence]
Col.FieldsStruct.SigmaRawFlux  = 9;   % Error in raw flux [e-/cadence]
Col.FieldsStruct.CalFlux       = 10;  % Corrected photometry [e-/cadence]
Col.FieldsStruct.SigmaCalFlux  = 11;  % Error in Corrected photometry [e-/cadence]
Col.FieldsStruct.InsMag        = 12;
Col.FieldsStruct.SigmaInsMag   = 13;
Col.FieldsStruct.D_RawFlux     = 14;
Col.FieldsStruct.D_SigmaRawFlux= 15;
Col.FieldsStruct.D_CalFlux     = 16;
Col.FieldsStruct.D_SigmaCalFlux= 17;
Col.FieldsStruct.D_Insmag      = 18;
Col.FieldsStruct.D_SigmaInsMag = 19;
Col.FieldsCell = fieldnames(Col.FieldsStruct);

%-------------------
%--- Read Header ---
%-------------------
% Time in MJD
Keys = {'LC_START','LC_END','KEPLERID','RA','DEC','GMAG','RMAG','IMAG','ZMAG','D51MAG','JMAG','HMAG','KMAG','KEPMAG','GALAXY','BLEND','TEFF','LOGG','FEH','AV','EBMINUSV','RADIUS'};

Vals = FITS.get_keys(FileName,Keys);
Head = cell2struct(Vals,Keys,1);
