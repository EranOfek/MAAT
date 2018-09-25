function [Im, metadata_struct] = getImageRGBG(filename)
    
    I = dslr.open_raw(filename);
    
    Ir = squeeze(I(1,:,:));
    Ir = Ir(1:2:end,1:2:end);
    
    Ig1 = squeeze(I(2,:,:));
    Ig1 = Ig1(2:2:end,1:2:end);
    
    Ig2 = squeeze(I(4,:,:));
    Ig2 = Ig2(1:2:end,2:2:end);
    
        
    Ib = squeeze(I(3,:,:));
    Ib = Ib(2:2:end,2:2:end);
        
    Im = cat(3,Ir,Ig1,Ib,Ig2);       

    Im = flip(Im,1);            %flip image
    
    Im = rot90(Im,3);           %rotate image
        
end