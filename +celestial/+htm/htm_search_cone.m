function ID=htm_search_cone(HTM,Long,Lat,Radius,Ind)
% Search for all HTM leafs interscting a small circle (cone search)
% Package: celestial.htm
% Description: Search for all HTM leafs interscting a small circle
%              (i.e., cone search).
% Input  : - HTM structure. See celestial.htm.htm_build
%          - Longitude [radians] to search.
%          - Latitude [radians] to search.
%          - Search radius [radians].
% Example:  [HTM,LevList]=celestial.htm.htm_build(4);
%           ID=celestial.htm.htm_search_cone(HTM,1,1,0.0001)
% Reliable : 2

if (nargin<5)
    Ind = [];
end


if isempty(Ind)
    % first iteration
    Sons = (1:1:8);
else
    Sons = Ind;
end

ID = [];
Nsons = numel(Sons);
PolesLong = zeros(3,Nsons);
PolesLat  = zeros(3,Nsons);

for Isons=1:1:Nsons
    %CSon = Sons(Isons);
    PolesLong(:,Isons) = HTM(Sons(Isons)).PolesCoo(:,1);
    PolesLat(:,Isons)  = HTM(Sons(Isons)).PolesCoo(:,2);
end
Flag = celestial.htm.cone_in_polysphere(PolesLong,PolesLat,Long,Lat,Radius);


for Isons=1:1:Nsons
    if (Flag(Isons))
        % cone overlap HTM
        CSon = Sons(Isons);
        if isempty(HTM(CSon).son)
            % arrived at last leaf
            % return ID
            ID = [ID, CSon];
        else
            Ind = HTM(CSon).son;
            ID  = [ID, celestial.htm.htm_search_cone(HTM,Long,Lat,Radius,Ind)];
            %ID = cat(2,ID,celestial.htm.htm_search_cone(HTM,Long,Lat,Radius,Ind));
        end
    end
end      