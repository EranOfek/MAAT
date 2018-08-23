function SimT=transform(Sim,Tran,varargin)
% Apply spatial transformations to SIM images
% Class  : class/@SIM
% Description: Apply a set of spatial transformations to a set of SIM
%              images. The transformations may be of various types.
%              To find the best fit transformation and apply it use
%              SIM/align.m.
% Input  : - A SIM object.
%          - Transformation for each SIM image.
%            The following transformation are supported:
%            Option (1):
%            A matrix in which each line represent a transformation of the
%            corresponding image.
%            The matrix may have 2, 3 or 6 columns.
%            2 column matrix represent [X, Y] shifts.
%            3 column matrix represent [Xshift, Yshift, rotation(deg)]
%            6 col matrix represent affine transformation of the form:
%            [SX,SY, a,b,c,d] where the affine 2-D matrix is:
%            [a b 0;
%             c d 0;
%             SX, SY 1].
%            Option (2):
%            An affine2d transformation object that represent the affine
%            transformation per image.
%            Option (3):
%            An images.geotrans.PolynomialTransformation2D transformation
%            object that represent the 2D polynomial transformation per
%            image.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - Cell array of SIM field names on which to
%                        execute the transformation. Default is {'Im'}.
%            'TransFun'  - Transformation function name.
%                          'imwarp' - use @imwarp.
%            'InterpMethod' - Interpolation method.
%                        Default is 'cubic'.
%            'InterpMethodMask' - Interpolation method for the MASK image
%                        Default is 'nearest'.
%            'FillValues' - Fill values. Default is NaN.
%            'SmoothEdges' - imwarp SmoothEdges option. Default is false.
%            'Verbose'    - Verbose. Default is true.
% Output : - A SIM object with the transformed images.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=FITS.read2sim('test_PTF_*f02_p004162_c08.fits');
%          S=mextractor(S);
%          [AstM]=match(S);
%          [MatchedCat]=astcat2matched_array(AstM,{'XWIN_IMAGE','YWIN_IMAGE','MAG_PSF'});
%          TranC = TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunTiltXp,[],@FunTiltXn,[]}, {@FunOne, [],@FunX,[],@FunY,[],@FunTiltYp,[],@FunTiltYn,[] })
%          Res=fit_transform(TranC,MatchedCat)
%          SimT=transform(S(2),Res(2).TranC);
%          %ds9(S(1),1); ds9(SimT,2);
%          % Shift an image
%          SimT=transform(S(1),[10 10]);
% Reliable: 2
%--------------------------------------------------------------------------

MaskField = MASK.MaskField;


DefV.ExecField            = {SIM.ImageField};
DefV.TransFun             = 'imwarp';
DefV.InterpMethod         = 'cubic';
DefV.InterpMethodMask     = 'nearest';
DefV.FillValues           = NaN;
DefV.SmoothEdges          = false;
DefV.Verbose              = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);




% Deal with transformation types
if (isnumeric(Tran))
    % assume Tran is a matrix of shift, rotation or affine 2d
    % transformation
    
    % assume each line in Tran contains a transformation:
    [Ntran,Ncol] = size(Tran);
    for Itran=1:1:Ntran
        % 
        % 3-by-3 double-precision, floating point matrix that defines the 2-D forward affine transformation
        % [ShiftX, ShiftY, a b c d]
        % The matrix T uses the convention:
        % [x y 1] = [u v 1] * T
        % where T has the form:
        %[a b 0;
        % c d 0;
        % e f 1];
        % 10 deg counter clock wise:
        % tform = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1])
        if (Ncol==2)
            % [ShiftX, ShiftY]
            Tform(Itran) = affine2d([1 0 0;0 1 0; Tran(Itran,1:2) 1]);
        elseif (Ncol==3)
            % [ShiftX, ShiftY, Theta]
            Theta = Tran(Itran,3);
            Tform(Itran) = affine2d([cosd(Theta) -sind(Theta) 0; sind(Theta) cosd(Theta) 0; Tran(Itran,1:2) 1]);
        elseif (Ncol==6)
            % [ShiftX, ShiftY, a b c d]
            Tform(Itran) = affine2d([Tran(Itran,3:4) 0; Tran(Itran,5:6) 0; Tran(Itran,1:2) 1]);
        else
            error('Unknown Affine transformation');
        end
            
    end
elseif (isa(Tran,'TranClass'))
    % Tran is a TranClass transformation object
    % convert to images.geotrans.PolynomialTransformation2D
    Tform = tranclass2PolynomialTransformation2D(Tran,true);
    
elseif (isa(Tran,'affine2d'))
    % Tran is affine2d transformation
    Tform = Tran;
    
elseif (isa(Tran,'images.geotrans.PolynomialTransformation2D'))
    % Tran is images.geotrans.PolynomialTransformation2D transformation
    Tform = Tran;
   
else
    error('Unsupported transformation class');
end
Ntran = numel(Tform);


if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

SimT = Sim;
Nsim = numel(Sim);
for Isim=1:1:Nsim
    % for each SIM element
    for If=1:1:Nf
        % for each field
        Itran = min(Ntran,If);
        switch lower(InPar.TransFun)
            case 'imwarp'
                ImageSize = size(Sim(Isim).(InPar.ExecField{If}));
                RI    = imref2d(ImageSize);
                
                if (strcmp(InPar.ExecField{If}, MaskField))
                    % treat the interpolation of mask differently
                    InterpMethod = InPar.InterpMethodMask;
                else
                    InterpMethod = InPar.InterpMethod;
                end
                
                % transform the image
                SimT(Isim).(InPar.ExecField{If}) = imwarp(Sim(Isim).(InPar.ExecField{If}), Tform(Itran),...
                                                          InterpMethod,...
                                                          'FillValues',InPar.FillValues,...
                                                          'SmoothEdges',InPar.SmoothEdges,...
                                                          'OutputView',RI);
            otherwise
                error('Unknown TransFun option');
        end
    end
    % correct WCS
    if (InPar.Verbose)
        fprintf('SIM/transform: WCS was not corrected\n');
    end
end

        
        
        
        