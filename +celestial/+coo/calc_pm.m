function [PM_RA,PM_Dec,PM_Tot,PM_PA]=calc_pm(PosList,ErrMat,TimeUnit)
%--------------------------------------------------------------------------
% calc_pm function     
% Description: Calculate the proper motion of a star from a set of
%              measurments.
%              OBSOLETE: 
% Input  : - List of positions [Time, RA, Dec] per line
%            (coordinates in radians).
%          - Error matrix: [RA_Err, Dec_Err] per line (in radians).
%          - Time units:
%            'j'   - Julian years  (default).
%            'jd'  - Julian days
% Output : - PM_RA (RA*cos(Dec)) [Value, left_err, right_err] ["/year]
%          - PM_Dec              [Value, left_err, right_err] ["/year]
%          - PM_Tot              [Value, left_err, right_err] ["/year]
%          - PM_PA               [Value, left_err, right_err] [deg]
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      January 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------
RAD = 180./pi;
OneSigma = 0.6827;

TimeCol = 1;
RACol   = 2;
DecCol  = 3;

N = length(PosList(:,1));

if (nargin==2)
   TimeUnit = 'j';
end

% convert JD to Julian years:
switch TimeUnit
 case 'j'
    Time = PosList(:,TimeCol);
 case 'jd'
    Time = 2000.0 + (PosList(:,TimeCol)-2451545.5)./365.25;
 otherwise
    error('Unknown TimeUnit');
end

if (N==1)
   errorxy('Impossible to calculate PM from one observation');
else
   RA      = PosList(:,2);
   Dec     = PosList(:,3);
   RA_Err  = ErrMat(:,1);
   Dec_Err = ErrMat(:,2);
   RA      = RA.*cos(Dec);
   RA_Err  = RA_Err.*cos(Dec);
   % Proper Motion Monte-Carlo
   %rand('state',sum(clock*100));

   Nsim = 10000;
   DiffX = zeros(Nsim.*(N-1),1);
   DiffY = zeros(Nsim.*(N-1),1);
   Rate  = zeros(Nsim.*(N-1),1);
   PA    = zeros(Nsim.*(N-1),1);

   for I=1:(N-1):(Nsim-1).*(N-1)

      NewRA  = RA + RA_Err.*randn(N,1);
      NewDec = Dec + Dec_Err.*randn(N,1);

      for J=1:1:N-1
         [Dist,PAs] = sphere_dist(NewRA(J),NewDec(J),NewRA(J+1),NewDec(J+1));
         Rate(I+J-1) = Dist;
         PA(I+J-1)   = PAs;
      end

      DiffX(I:I+N-2) = diff(NewRA);
      DiffY(I:I+N-2) = diff(NewDec);
   end

   %Rate = sqrt(DiffX.^2 + DiffY.^2);

   DiffTime = repmat(diff(Time),Nsim,1);

   DiffX = DiffX./DiffTime;
   DiffY = DiffY./DiffTime;
   Rate  = Rate./DiffTime;

   %hist(DiffX.*RAD.*3600,100)

   % errors
   X_CI   = err_cl(DiffX,OneSigma);
   Y_CI   = err_cl(DiffY,OneSigma);
   R_CI   = err_cl(Rate,OneSigma);
   PA_CI  = err_cl(PA,OneSigma);
   PM_RA  = median(DiffX);
   PM_Dec = median(DiffY);
   PM_Tot = median(Rate);
   PM_PA  = median(PA);
   PM_RA  = [PM_RA, PM_RA-X_CI(1), X_CI(2)-PM_RA];
   PM_Dec = [PM_Dec, PM_Dec-Y_CI(1), Y_CI(2)-PM_Dec];
   PM_Tot = [PM_Tot, PM_Tot-R_CI(1), R_CI(2)-PM_Tot];
   PM_PA  = [PM_PA, PM_PA-PA_CI(1), PA_CI(2)-PM_PA];

   PM_RA  = PM_RA.*RAD.*3600;
   PM_Dec = PM_Dec.*RAD.*3600;
   PM_Tot = PM_Tot.*RAD.*3600;
   PM_PA  = PM_PA.*RAD;

end
