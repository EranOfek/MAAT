function prep_unWISE_htm(varargin)
% SHORT DESCRIPTION HERE
% Package: VO
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.prep.prep_unWISE_htm
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.Verbose             = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


[~,List] = Util.files.create_list('*.cat.fits');
N        = numel(List);
Tmp = regexp(List,'(p)|(m)|(\.)','split');
RA  = nan(N,1);
Dec = nan(N,1);
Sign = nan(N,1);
for I=1:1:N
    switch List{I}(5)
        case 'm'
            Sign(I) = -1;
        case 'p'
            Sign(I) = 1;
    end
    
    % because dec near zero is either 0 or 1.5 deg - can multiply by sign
    
    RA(I)  = str2double(Tmp{I}{1});
    Dec(I) = Sign(I).*str2double(Tmp{I}{2});
    
end

RA = RA./10;
Dec = Dec./10;

DecBuffer = 1;
DecStep = 5;
DecRange = (-90:DecStep:90)';
Ndec     = numel(DecRange);

for Idec=1:1:Ndec-1
    D1 = DecRange(Idec);
    D2 = DecRange(Idec+1);
    
     if (InPar.Verbose)
         fprintf('Idec=%d out of %d\n',Idec,Ndec-1);
     end
    
    Ilist = find(Dec>(D1-DecBuffer) & Dec<(D2+DecBuffer));
    
    SubCol = {'RA','Dec','MagAB_w1','MagAB_w2','ErrMag_w1','ErrMag_w2','SkyMag_w1','SkyMag_w2','N_w1','N_w2','Flags_w1','Flags_w2','FlagsAdd_w1','FlagsAdd_w2'};
    
    % read all files in FlagList
    Nf = numel(Ilist);
    for If=1:1:Nf
        AC(If) = FITS.read_table(List{Ilist(If)});
        
        MagABw1 = 22.5 - 2.5.*log10(AC(If).Cat(:,19)) + 2.699;
        MagABw2 = 22.5 - 2.5.*log10(AC(If).Cat(:,20)) + 3.339;
        ErrMagABw1 = 1.086.*AC(If).Cat(:,21)./AC(If).Cat(:,19);
        ErrMagABw2 = 1.086.*AC(If).Cat(:,22)./AC(If).Cat(:,20);
        SkyABw1    = 22.5 - 2.5.*log10(AC(If).Cat(:,29)) + 2.699;
        SkyABw2    = 22.5 - 2.5.*log10(AC(If).Cat(:,30)) + 3.339;

        SubCat = [AC(If).Cat(:,69)./RAD, AC(If).Cat(:,70)./RAD, MagABw1, MagABw2, ErrMagABw1, ErrMagABw2,  SkyABw1, SkyABw2, AC(If).Cat(:,[61 62 65 66 67 68])];

        AC(If).Cat = SubCat;
        AC(If).ColCell = SubCol;
        AC(If)  = colcell2col(AC(If));
    end
    MAC = merge(AC);
    clear AC
    
    % save only necessey columns
    % replace 0 RA/Dec with nan
    
    
    
%     fluxlbs_1_   19
%     fluxlbs_2_   20
%     dfluxlbs_1_  21
%     dfluxlbs_2_  22
%     sky_1_       29 
%     sky_2_       30
%     ra 69
%     dec 70
%     nm_1_      61
%     nm_2_      62
%     flags_unwise_1_  65
%                      66
%     flags_info_1_    67
%                      68

%     MagABw1 = 22.5 - 2.5.*log10(MAC.Cat(:,19)) + 2.699;
%     MagABw2 = 22.5 - 2.5.*log10(MAC.Cat(:,20)) + 3.339;
%     ErrMagABw1 = 1.086.*MAC.Cat(:,21)./MAC.Cat(:,19);
%     ErrMagABw2 = 1.086.*MAC.Cat(:,22)./MAC.Cat(:,20);
%     SkyABw1    = 22.5 - 2.5.*log10(MAC.Cat(:,29)) + 2.699;
%     SkyABw2 = 22.5 - 2.5.*log10(MAC.Cat(:,30)) + 3.339;
%     
%     SubCat = [MAC.Cat(:,69)./RAD, MAC.Cat(:,70)./RAD, MagABw1, MagABw2, ErrMagABw1, ErrMagABw2,  SkyABw1, SkyABw2, MAC.Cat(:,[61 62 65 66 67 68])];

    
%     Cat = AstCat;
%     Cat.Cat = SubCat;
%     Cat.ColCell = SubCol;
%     Cat     = colcell2col(Cat);
    MAC     = sortrows(MAC,2);
    
    %call VO.prep.build_htm_catalog
    Nsrc = VO.prep.build_htm_catalog(MAC,'CatName','unWISE','HTM_Level',9,'DecRange',[D1 D2]./RAD);
    clear MAC;
    
   
        
end