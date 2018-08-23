function Dir=files_by_date(Files,MinDate,MaxDate)
% Select files by date
% Package: Util.files
% Input  : - Files string
%          - Minimum date in any format acceptable by
%            celestial.time.julday.
%          - Maximum date.
% Output : - Selected files
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Dir=Util.files.files_by_date('*.hdf5',[1 1 2000],[2 2 2017]);
% Reliable: 

Dir = dir(Files);
DateV = datevec([Dir.datenum]);
JD    = celestial.time.julday(DateV(:,[3 2 1 4 5 6]));
MinJD = celestial.time.julday(MinDate);
MaxJD = celestial.time.julday(MaxDate);

Flag  = JD>MinJD & JD<MaxJD;
Dir   = Dir(Flag);