function [RT]=fdrt(Image)
% Fast discrete radio transform
% Package: radio.fdmt
% Description: Calculate the fast discrete radio transform
% Input  : - 2 dimensional image, both dimensions should be powers of 2,
%            not necessarily the same power
% Output : - Radon transform of the image
%            each coordinate is a sum of an integral of the input over a
%            line. only lines between -45<theta<45 degrees that intersect 
%            the vertical line x=0 are represented. 
%            For the full Radon transform, apply also for the transpose of the data  
%            This function will initiate data structure that has: 
%            X position - the position of the radon transform, starts with a regular x
%                position and after i iterations, only the k*2^i +1 for k=1:(x_len/2^i) are legitimate
%            Y position - the Y position sub-index of the radon transform
%            deltaY position - The third coordinate will represent the 
%                theta by the formula sin(theta) = (deltaY/2**i) in the i'th iteration
%                deltaY starts with 0 as the only legitimate value (change, this is good
% only for size(input,1) == size(input,2) == 2^n) and changes each
% iteration to contain all the values between -deltaX:deltaX, where deltaX
% = 2^i (if i is the iteration)
% Tested : Matlab R2014a
%     By : Barak Zackay                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

[N,M] = size(Image);
[n,m] = deal(log2(N),log2(M));

State1 = FDRT_initialization(Image);
State2 = FDRT_initialization(Image);

% the first dimension is taken care by the initialization ...
%for i_n = 1:n,
%    State = FDRT_iteration(State,1);
%end

for i_m = 1:m,
    %tic;
    State1 = FDRT_iteration_Fast(State1,i_m);
    %toc
end


%tic;
%for i_m = 1:m,
%    State2 = FDRT_iteration_Slow(State2,i_m);
%end
%toc
[X,Y,dY] = size(State1);

RT = reshape(State1, Y,dY);
%RT2 = reshape(State2,Y,dY);
end

function [Output] = FDRT_iteration_Slow(Input,iteration_num)
% internal function that does one iteration of radon transform on the Input
% with dimension "dim" shrinking and the result dimensions getting larger
%iteration_num
input_dims = size(Input);
output_dims = input_dims;
% x position is being shrinked
% deltaX is the assumed true distance after the iteration is ended
deltaX = 2^(iteration_num);

output_dims(1) = output_dims(1)/2;

% You might want to extend the Y axis of the input by a factor of 3, and by
% that vectorize the dY axis operations (along with the already vectorized X axis)

% the deltaY is limited to the range -deltaX:deltaX including (which maybe
% is a mistake, maybe I should not do including...)
output_dims(3) = 2*deltaX + 1; 
Output = zeros(output_dims);
RadonShiftOutput = deltaX+1;
RadonShiftInput = floor(deltaX/2)+1;
Y = output_dims(2);
for i_Y = 1:Y,
    for i_dY = -deltaX:deltaX,
        DY = round(i_dY/2) + RadonShiftInput;
        if DY == 0,
            % this case happens only at the lower limit i_dY = -deltaX, and
            % then, 1 is the correct value...
            DY = 1;
        end
        if DY >= size(Input,3),
            DY = size(Input,3);
        end
        if (i_Y + round(i_dY/2) < Y) && (i_Y + round(i_dY/2))>0,
            Output(:,i_Y,i_dY + RadonShiftOutput) = Input(1:2:input_dims(1),i_Y, DY) + Input(2:2:input_dims(1),i_Y + round(i_dY/2), DY);
        else 
            Output(:,i_Y,i_dY + RadonShiftOutput) = Input(1:2:input_dims(1),i_Y, DY);
        end
        % MUST attend the ELSE IF! 
    end
end

end

function [Output] = FDRT_iteration_Fast(Input,iteration_num)
% internal function that does one iteration of radon transform on the Input
% with dimension "dim" shrinking and the result dimensions getting larger
%iteration_num
input_dims = size(Input);
output_dims = input_dims;
% x position is being shrinked
% deltaX is the assumed true distance after the iteration is ended
if numel(input_dims) == 2, 
    input_dims(3)=1;
end
deltaX = 2^(iteration_num);
ExtendedInput = zeros(input_dims(1),input_dims(2) + 2*deltaX+1,input_dims(3));
ExtendedInput(:,[1+deltaX:deltaX + input_dims(2)],:) = Input;
InputYShift = 1+deltaX;
output_dims(1) = output_dims(1)/2;

output_dims(3) = 2*deltaX + 1; 
Output = zeros(output_dims);
RadonShiftOutput = deltaX+1;
RadonShiftInput = floor(deltaX/2)+1;
Y = output_dims(2);

% this array is the array of Delta Y's, regular almost all the way, with
% the exception of the boundaries taken care by the ifs
DY = zeros(1,2*deltaX+1);
for i_dY = -deltaX:deltaX,
        DY(i_dY + deltaX+1) = round(i_dY/2) + RadonShiftInput;
        if DY(i_dY + deltaX+1) == 0,
            % this case happens only at the lower limit i_dY = -deltaX, and
            % then, 1 is the correct value...
            DY(i_dY + deltaX+1) = 1;
        end
        if DY(i_dY + deltaX+1) >= size(Input,3),
            DY(i_dY + deltaX+1) = size(Input,3);
        end
end
i_dY = -deltaX:deltaX;

% all logic of the function is here...
for i_Y = 1:Y,
    v1 = InputYShift+i_Y + round(i_dY/2);
    IndsList = v1 + (DY - 1).*size(ExtendedInput,2);
    %IndsList = sub2ind([size(ExtendedInput,2),size(ExtendedInput,3)],InputYShift+i_Y + round(i_dY/2), DY);
    A = ExtendedInput(1:2:input_dims(1),InputYShift+i_Y, DY);
    B = ExtendedInput(2:2:input_dims(1),IndsList);
    Output(:,i_Y,i_dY + RadonShiftOutput) = A + reshape(B,size(A)); 
end
end


function [Output] = FDRT_initialization(Image)
% Builds the third dimension, copies the data
[X,Y] = size(Image);
Output = zeros(X,Y,1);
Output(:,:,1) = Image;
end
%Gaps = cellfun(@isempty, varargin);
%DefArgs = {InPar1Def InPar2Def InPar3Def};    % default input arguments
%Suboptargs = DefArgs(1 : NumVarargs);
%varargin(Gaps) = Suboptargs(Gaps);
%DefArgs(1 : NumVarargs) = varargin;
%[Par1 Par2 Par3] = DefArgs{:}


%DefV. = 
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
