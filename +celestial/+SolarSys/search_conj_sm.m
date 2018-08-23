function [Conj,IndMinSep]=search_conj_sm(JD,List1,Col1,AngFlag1,List2,Col2,AngFlag2,MinSep)
% Celestial conjunctions on the between moving and stationary objects
% Package: celestial.SolarSys
% Description: Search for conjunctions on the celestial sphere between a
%              list of stationary points and a moving object given the
%              coordinates of the moving object as a function of time.
% Input  : - Column vector of JD, refering to the moving object.
%          - First list of coordinates [radians] and other properties
%            for moving source.
%          - Columns indices of RA and Dec in first source list.
%            If empty matrix (i.e., []), then set it to [1 2].
%          - Vector of flaga indicating the type of column in the first
%            source list: 0 - non angle, or 1 - angle property.
%            Angle properties (0..2pi) are interpolated using
%            interp_diff_ang.m. If empty matrix (i.e., []), then
%            assumes onlt the RA column is an angle property.
%          - List of coordinates [radians] and other properties for
%            stationary points.
%          - Columns indices of RA and Dec in second source list.
%            If empty matrix (i.e., []), then set it to [1 2].
%          - Vector of flaga indicating the type of column in the first
%            source list: 0 - non angle, or 1 - angle property.
%            Angle properties (0..2pi) are interpolated using
%            interp_diff_ang.m. If empty matrix (i.e., []), then
%            assumes onlt the RA column is an angle property.
%          - Minimum seperation to look for [radians].
% Output : - Structure of all the conjunction between the bodies in
%            the two lists. The structure contains the following fields:
%            .MinJD    - Vector of JD of conjunctions.
%            .MinDist  - Vector of minimum separation [radians] in
%                        each conjunction.
%            .MinPA    - Vector of PA [radians] in each conjunction
%                        (at .MinJD)
%            .List1    - Matrix of List1 interpolated to the time
%                        of each conjunction (line per conjunction).
%            .List2    - Matrix of List2 interpolated to the time
%                        of each conjunction (line per conjunction).
%          - Indices of conjunctions that their .MinDist is smaller
%            than the minimum required separation. Note that the
%            conjunction structure may include some events with
%            minimum distance larger than the required minimum
%            separation.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Web example: http://astroclub.tau.ac.il/ephem/PlanetsGC/
% Example: %celestial.SolarSys.search_conj(JD,List1,[1 2],[],List2,[1 2],[],1./RAD);
% Reliable: 2
%------------------------------------------------------------------------------
RAD   = 180./pi;

NColList1 = size(List1,2);
NColList2 = size(List2,2);
Nstat     = size(List2,1);

if (isempty(Col1)==1),
   Col1 = [1 2];
end
if (isempty(Col2)==1),
   Col2 = [1 2];
end
if (isempty(AngFlag1)==1),
  AngFlag1 = zeros(NColList1,1);
  AngFlag1(Col1(1)) = 1;
end
if (isempty(AngFlag2)==1),
  AngFlag2 = zeros(NColList2,1);
  AngFlag2(Col2(1)) = 1;
end


Rate1 = celestial.coo.sphere_dist(List1(1:end-1,Col1(1)), List1(1:end-1,Col1(2)), List1(2:end,Col1(1)), List1(2:end,Col1(2))); 
Rate2 = 0;



MaxRate = max(Rate1 + Rate2);

% set MinSep to account for the fact the object are moving
% between tabulated points.
SearchMinSep = MinSep + MaxRate;

Conj.MinJD   = zeros(0,1);
Conj.MinDist = zeros(0,1);
Conj.MinPA   = zeros(0,1);
Conj.List1   = zeros(0,NColList1);
Conj.List2   = zeros(0,NColList2);

ConjInd = 0;
for Istat=1:1:Nstat,
    % find distances between current stationary object and all points
    % in moving object
   [Dist,PA] = celestial.coo.sphere_dist(List1(:,Col1(1)),List1(:,Col1(2)),List2(Istat,Col2(1)),List2(Istat,Col2(2)));

   Flag = double(Dist<=SearchMinSep);    % nearby flag
   % convolve Flag with [1 .. 1] to ensure that there are enough points
   % for interpolation:
   FlagConvLargeWin = conv(Flag,ones(20,1),'same');
   FlagConvSmallWin = conv(Flag,ones(3,1),'same');

   I    = find(FlagConvLargeWin>0.5);    % select all Flaged
   I    = I(find(I<=length(Dist)));
   if (length(I)<10),
      % skip - conjunction near end/begining of time in moving object list 
   else
      % search for local extramum in flaged list
      Extram = find_local_extramum(JD(I),Dist(I));

      if(isempty(Extram)==1),
         % skip
      else
         % select only minima:
         Im = find(Extram(:,3)>0);
         Minima = Extram(Im,:);
      
         % go over Minima
         Nm = size(Minima,1);
      
         for Im=1:1:Nm,
            ConjInd = ConjInd + 1;
      
            MinJD    = Minima(Im,1);
            MinDist  = Minima(Im,2);
            % check if its a real minima or artefact due to removal
            % of gaps between near approaches
            [Min,MinInd] = min(abs(MinJD-JD));
            if (FlagConvSmallWin(MinInd)>0.5),
               % real minima
               % [MinJD-2450000,MinDist.*RAD]
         
               % Interpolate data in both tables for MinJD
               InterpList1 = zeros(1,NColList1);
               for Ic1=1:1:NColList1,
                  if (AngFlag1(Ic1)>0.5),
                     % angular proprty
                     IntVal = interp_diff_ang(JD,List1(:,Ic1),MinJD);
                  else
                     % NON-angular proprty
                     IntVal = interp_diff(JD,List1(:,Ic1),MinJD);
                  end
                  InterpList1(Ic1) = IntVal;
               end

               % calculate PA at MinJD
               [D,PA] = celestial.coo.sphere_dist(InterpList1(Col1(1)),InterpList1(Col1(2)), List2(Istat,Col2(1)),List2(Istat,Col2(2)));
         
               Conj.MinJD(ConjInd,1)   = MinJD;
               Conj.MinDist(ConjInd,1) = MinDist;
               Conj.MinPA(ConjInd,1)   = PA;
               Conj.List1(ConjInd,:)   = InterpList1;
               Conj.List2(ConjInd,:)   = List2(Istat,:);
         
            else
               % artfact - skip
            end
         end
      end
   end
end   

% sort by time
[Conj.MinJD,SortI] = sort(Conj.MinJD);
Conj.MinDist       = Conj.MinDist(SortI);
Conj.MinPA         = Conj.MinPA(SortI);
Conj.List1         = Conj.List1(SortI,:);
Conj.List2         = Conj.List2(SortI,:);

% indices of conjunctios with MinSep or less...
IndMinSep = find(Conj.MinDist<=MinSep);
