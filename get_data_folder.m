function dirname = get_data_folder
% Usage: dirname = get_data_folder
% Find the folder with all the data on any system
% To get this to work you must define the environmental variable:
%     DATA = path/to/data/folder

    dirname = getenv('DATA');
    
    if isempty(dirname) % if no environmental variable is found, use Eran's default hard-coded location
        dirname = '/home/eran/matlab/data/';
    end
    
    if ~exist(dirname, 'dir')
        warning('Cannot find the directory "%s". Please define the environmental variable "DATA" in your system!', dirname);
    end

end