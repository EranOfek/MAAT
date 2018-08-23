function [Conj,IndMinSep]=search_conj(JD,List1,Col1,AngFlag1,List2,Col2,AngFlag2,MinSep)
% Celestial conjunctions between moving objects
% Package: celestial.SolarSys
% Description: Search for conjunctions on the celestial sphere between
%              two moving objects given thier coordinates as a function
%              of time.
% Input  : - Column vector of JD.
%          - First list of coordinates [radians] and other properties
%            for moving source.
%          - Columns indices of RA and Dec in first source list.
%            If empty matrix (i.e., []), then set it to [1 2].
%          - Vector of flags indicating the type of column in the first
%            source list: 0 - non angle, or 1 - angle property.
%            Angle properties (0..2pi) are interpolated using
%            interp_diff_ang.m. If empty matrix (i.e., []), then
%            assumes onlt the RA column is an angle property.
%          - List of coordinates [radians] and other properties for
%            second moving source.
%          - Columns indices of RA and Dec in second source list.
%            If empty matrix (i.e., []), then set it to [1 2].
%          - Vector of flags indicating the type of column in the second
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
% Web example: http://astroclub.tau.ac.il/ephem/PlanetsConj/
%              http://astroclub.tau.ac.il/ephem/AsteroidsPlanetsConj/
%              http://astroclub.tau.ac.il/ephem/LunarOcc/PlanetsConj/
%              http://astroclub.tau.ac.il/ephem/LunarOcc/Planets/
% Example: %search_conj(JD,List1,[1 2],[],List2,[1 2],[],1./RAD);
% Reliable: 2
%------------------------------------------------------------------------------
RAD   = 180./pi;

NColList1 = size(List1,2);
NColList2 = size(List2,2);

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
Rate2 = celestial.coo.sphere_dist(List2(1:end-1,Col2(1)), List2(1:end-1,Col2(2)), List2(2:end,Col2(1)), List2(2:end,Col2(2))); 

MaxRate = max(Rate1 + Rate2);

% set MinSep to account for the fact the object are moving
% between tabulated points.
SearchMinSep = MinSep + MaxRate;

[Dist,PA] = celestial.coo.sphere_dist(List1(:,Col1(1)),List1(:,Col1(2)),List2(:,Col2(1)),List2(:,Col2(2)));

Flag = double(Dist<=SearchMinSep);    % nearby flag
% convolve Flag with [1 .. 1] to ensure that there are enough points
% for interpolation:
FlagConvLargeWin = conv(Flag,ones(20,1),'same');
FlagConvSmallWin = conv(Flag,ones(3,1),'same');

I    = find(FlagConvLargeWin>0.5);    % select all Flaged
if (isempty(I)==1),
   % no conjunctions
   Conj = [];
   IndMinSep = [];
else
   
   % search for local extramum in flaged list
   size(I)
   size(Dist)
   size(JD)
   max(I)
   min(I)

   Extram = find_local_extramum(JD(I),Dist(I));

   if (isempty(Extram)==1),
      % skip
      Conj = [];
      IndMinSep = [];
   else
      Im = find(Extram(:,3)>0);
      Minima = Extram(Im,:);
   
      % go over Minima
      Nm = size(Minima,1);

      Conj.MinJD   = zeros(Nm,1);
      Conj.MinDist = zeros(Nm,1);
      Conj.MinPA   = zeros(Nm,1);
      Conj.List1   = zeros(Nm,NColList1);
      Conj.List2   = zeros(Nm,NColList2);
      for Im=1:1:Nm,
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
            % and List2...
            InterpList2 = zeros(1,NColList2);
            for Ic2=1:1:NColList2,
               if (AngFlag2(Ic2)>0.5),
                  % angular proprty
                  IntVal = interp_diff_ang(JD,List2(:,Ic2),MinJD);
               else
                  % NON-angular proprty
                  IntVal = interp_diff(JD,List2(:,Ic2),MinJD);
               end
               InterpList2(Ic2) = IntVal;
            end
      
            % calculate PA at MinJD
            [D,PA] = celestial.coo.sphere_dist(InterpList1(Col1(1)),InterpList1(Col1(2)), InterpList2(Col2(1)),InterpList2(Col2(2)));
      
            Conj.MinJD(Im)   = MinJD;
            Conj.MinDist(Im) = MinDist;
            Conj.MinPA(Im)   = PA;
            Conj.List1(Im,:) = InterpList1;
            Conj.List2(Im,:) = InterpList2;
      
         else
            % artfact - skip
         end
      end
      % indices of conjunctios with MinSep or less...
      IndMinSep = find(Conj.MinDist<=MinSep);
   end   
end
