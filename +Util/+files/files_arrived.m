function files_arrived(Files)
% Check if all files in a list arrived to disk (i.e., size not increasing).
% Package: Util.files
% Description: Check if all files in a list arrived to disk. This is done
%              by checking that the file size does not increase with time.
% Input  : - Cell array of files to check. 
% Output : null. Return when file sizes converged.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Util.files.files_arrived({'files_arrived'});
% Reliable: 2
%--------------------------------------------------------------------------

Arrived = false;
Dir1 = Util.files.dir_cell(Files);
while ~Arrived    
    % check if all files arrived
    pause(0.2);
    Dir2 = Util.files.dir_cell(Files);

    try
        if (all(([Dir1.bytes] - [Dir2.bytes])==0))
            Arrived = true;
        end
    catch
        fprintf('May be a problem Dir1 vs. Dir2');
        pause(1);
    end

    Dir1 = Dir2;
end
pause(0.2);

    