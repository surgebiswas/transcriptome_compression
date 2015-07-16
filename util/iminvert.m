function iminvert( imfile )

imwrite(imcomplement(imread(imfile)), [imfile,'_inverted.png'], 'PNG');


end

