function mask = convhull2mask(imgi)
    imga = squeeze(bdilation(imgi,1,3,0));
    [r,c,v] = ind2sub(size(imga),find(imga ==1 ));
    img=imga(min(r)-1:max(r)-1,min(c)-1:max(c)-1,min(v)-1:max(v)-1);
    imglab = label(~logical(img));
    CH = bwconvhull(logical(img(:,:,ceil(size(img,3)/2))));
    for i=1:max(imglab)
        sz(i) = sum(CH & (imglab(:,:,ceil(size(img,3)/2)) == i));
    end
    [val in] = max(sz);  
    ff = newim(imgi,'bin');
    ff(min(r)-1:max(r)-1,min(c)-1:max(c)-1,min(v)-1:max(v)-1) = (imglab == in);
    mask = berosion((ff > 0 | imga),1,3,0);
end