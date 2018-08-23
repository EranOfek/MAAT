function [HTM,LevList]=htm_build_son(HTM,LevList,Level,Ind)
% An auxilary function for htm_build
% Package: celestial.htm
% Description: An auxilary function for htm_build.m for building
%              Hierarchical Triangular Mesh (HTM) structure.
% Input  : - The HTM structure array (see htm_build.m).
%          - The LevList structure array (see htm_build.m).
%          - Number of levels required.
%          - The last populated index
% Output : - The HTM structure array (see htm_build.m).
%          - The LevList structure array (see htm_build.m).
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jul 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [HTM,LevList]=celestial.htm.htm_build(4);
%          [HTM,LevList]=htm_build_son(HTM,LevList,2)
% Reliable: 2
%------------------------------------------------------------------------------

%Ind        = numel(HTM);
LevelDepth = numel(LevList);

if (Level>LevelDepth)
   % didn't reach yet to the lowest level - continuue recusrion

   %LevList(LevelDepth+1).ptr   = []; 
   LevList(LevelDepth+1).ptr   = zeros(1,8.*4.^(LevelDepth));
   Nleaf = numel(LevList(LevelDepth).ptr);
   K = 0;
   for Ileaf=1:1:Nleaf
       % for each leaf

       FatherPtr   = LevList(LevelDepth).ptr(Ileaf);
       %FatherPtr
       FatherLevel = HTM(FatherPtr).level;
       % unit vector of verteces
       Vert1 = HTM(FatherPtr).cosd(1,:);
       Vert2 = HTM(FatherPtr).cosd(2,:);
       Vert3 = HTM(FatherPtr).cosd(3,:);

       % unit vector of mid points in between verteces
       %Cen12 = gc_mid_section(Vert1,Vert2);
       %Cen23 = gc_mid_section(Vert2,Vert3);
       %Cen31 = gc_mid_section(Vert3,Vert1);

       Cen = celestial.htm.gc_mid_section([Vert1;Vert2;Vert3],[Vert2;Vert3;Vert1]);
       %Cen12 = Cen(1,:);
       %Cen23 = Cen(2,:);
       %Cen31 = Cen(3,:);
       
       % sub traiangle corners
       % Vert1 Cen12 Cen31
       % Vert2 Cen23 Cen12
       % Vert3 Cen31 Cen23
       % Cen12 Cen23 Cen31

       for Il=1:1:4
          Ind = Ind + 1;
          HTM(Ind).level    = FatherLevel + 1;

          switch Il
           case 1
              %HTM(Ind).cosd     = [Vert1; Cen12; Cen31];  
              HTM(Ind).cosd     = [Vert1; Cen(1,:); Cen(3,:)];  
           case 2
              %HTM(Ind).cosd     = [Vert2; Cen23; Cen12];  
              HTM(Ind).cosd     = [Vert2; Cen(2,:); Cen(1,:)];  
           case 3
              %HTM(Ind).cosd     = [Vert3; Cen31; Cen23];  
              HTM(Ind).cosd     = [Vert3; Cen(3,:); Cen(2,:)];  
           case 4
              %HTM(Ind).cosd     = [Cen12; Cen23; Cen31];
              HTM(Ind).cosd     = [Cen(1,:); Cen(2,:); Cen(3,:)];
              
           otherwise
              error('Illegal Il option');
          end

          [Long, Lat]       = celestial.coo.cosined2coo(HTM(Ind).cosd(:,1),HTM(Ind).cosd(:,2),HTM(Ind).cosd(:,3));
          HTM(Ind).coo      = [Long, Lat];
          HTM(Ind).id       = [HTM(FatherPtr).id, Il-1];
          HTM(Ind).father   = FatherPtr;
          HTM(Ind).son      = [];
          HTM(Ind).cat      = [];
          HTM(FatherPtr).son = [HTM(FatherPtr).son, Ind];
          %HTM(Ind).center_cosd = mean(HTM(Ind).cosd);
          %[LonC,LatC]       = celestial.coo.cosined2coo(HTM(Ind).center_cosd(1),HTM(Ind).center_cosd(2),HTM(Ind).center_cosd(3));
          %HTM(Ind).center_coo = [LonC, LatC];
          [PoleLong,PoleLat]  =                 celestial.htm.polysphere_poles(HTM(Ind).coo(:,1), HTM(Ind).coo(:,2));
          HTM(Ind).PolesCoo   = [PoleLong(:), PoleLat(:)];
   
          K = K + 1;
          %[Ind, K,numel(LevList(LevelDepth+1).ptr)]
          %LevList(LevelDepth+1).ptr   = [LevList(LevelDepth+1).ptr, Ind];
          LevList(LevelDepth+1).ptr(K)   = Ind; 
      end
      LevList(LevelDepth+1).level = LevList(LevelDepth).level + 1;
      LevList(LevelDepth+1).side  = LevList(LevelDepth).side.*0.5;
   end

   if (Level<(LevelDepth+1))
      % stop
   else
      % build HTM recursively
      [HTM,LevList] = celestial.htm.htm_build_son(HTM,LevList,Level,Ind);
   end

else
   % stop - desired level reached
end







 
