function wget_all
% wget all Chandra observations in cats.X.ChandraObs
% Package: VO
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.Chandra.wget_all
% Reliable: 
%--------------------------------------------------------------------------


%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Cat = cats.X.ChandraObs;

Nid = numel(Cat.Cat.ObsID);
for Iid=1:1:Nid
    ObsID = Cat.Cat.ObsID(Iid);
    
    fprintf('---------------------------------\n');
    fprintf('wget ObsID=%d  (dir %d out of %d)\n',ObsID,Iid,Nid);
    fprintf('---------------------------------\n');
    
    VO.Chandra.wget_obsid(ObsID)
    
    pause(1);
    
end
