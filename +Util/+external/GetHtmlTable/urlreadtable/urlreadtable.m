function [mytable] = urlreadtable(url)
% "Copyright 2016 The MathWorks, Inc."

mypage = urlread(url);

startind = regexp(mypage, '<table');
endind =   regexp(mypage, '</table>');
endind = endind+7;
if length(startind)~=length(endind)
    disp('could not find start and end of tables');
    
end

for z = 1:length(startind)
    mytable{z} = table;
    web = [];    
    htmltable = mypage(startind(z):endind(z));     
    %Establish a row index
    rowind = 0;  
    % Build cell aray of table data
    try
        rows = regexpi(htmltable, '<tr.*?>(.*?)</tr>', 'tokens');
        for i = 1:numel(rows)
            colind = 0;
            if (isempty(regexprep(rows{i}{1}, '<.*?>', '')))
                continue
            else
                rowind = rowind + 1;
            end            
            headers = regexpi(rows{i}{1}, '<th.*?>(.*?)</th>', 'tokens');
            if ~isempty(headers)
                for j = 1:numel(headers)
                    colind = colind + 1;
                    data = regexprep(headers{j}{1}, '<.*?>', '');
                    if (~strcmpi(data,'&nbsp;'))
                        web{rowind,colind} = strtrim(data);
                    end
                end
                continue
            end
            cols = regexpi(rows{i}{1}, '<td.*?>(.*?)</td>', 'tokens');          
            for j = 1:numel(cols)
                if(rowind==1)
                    if(isempty(cols{j}{1}))
                        continue
                    else
                        colind = colind + 1;
                    end
                else
                    colind = j;
                end
                data = regexprep(cols{j}{1}, '&nbsp;', ' ');
                data = regexprep(data, '<.*?>', '');               
                if (~isempty(data))
                    web{rowind,colind} = strtrim(data) ;
                end
            end
        end
    catch
        disp('Could not parse this table');
    end
    mytable{z} = cell2table(web);
    aux = (regexprep(mytable{z}{1,:},'[ \[\]\-:=,\(\)\.;&#(0-9)\\^\{\}+''//\n]',''));
    for y = 1:length(aux)
        if length(char(aux{y}))>27
            aux{y} = aux{y}(1:27);
        end       
        mytable{z}.Properties.VariableNames{y}  = [char(mytable{z}.Properties.VariableNames{y}) char(aux{y})];
    end
    mytable{z}(1,:)=[];  
end
