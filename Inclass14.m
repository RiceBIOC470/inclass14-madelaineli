%Inclass 14

%Work with the image stemcells_dapi.tif in this folder
% (1) Make a binary mask by thresholding as best you can
img = ('stemcells_dapi.tif');
img = imread(img);
figure(1)
title('original')
imshow(imadjust(img))

mean = sum(img(:))/length(img(:));
mask = img>(mean+30);
figure(2)
title('binary mask');
imshow(mask)

% (2) Try to separate touching objects using watershed. Use two different
% ways to define the basins. (A) With erosion of the mask (B) with a
% distance transform. Which works better in this case?
%%
%figure out which objects potentially contain more than 1 nucleus
CC = bwconncomp(mask);
stats = regionprops(CC,'Area');
area = [stats.Area];
mean = sum(area)/length(area);
fusedCandidates = area > (mean + std(area));
sublist = CC.PixelIdxList(fusedCandidates);
sublist = cat(1,sublist{:});
fusedMask = false(size(mask));
fusedMask(sublist) = 1;

%%
%distance transform
D = bwdist(~fusedMask);
D = -D;
D(~fusedMask) = -Inf;
L = watershed(D);
rgb = label2rgb(L,'jet',[.5,.5,.5]);
figure(3)
imshow(rgb);
title('watershed transform of D');
newNuclearMask = L>1|(mask - fusedMask);
figure(4)
imshow(newNuclearMask)
title('final result with distance transform')

%%
%erode
s = round(1.2*sqrt(mean)/pi);
nucmin = imerode(fusedMask,strel('disk',s));
outside = ~imdilate(fusedMask,strel('disk',1));
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin,nucmin|outside);
pcolor(basin);shading flat;
L_2 = watershed(basin);
figure(5)
imshow(L_2,[]); colormap('jet');caxis([0 20]);
newNuclearMask_2 = L_2>1|(mask - fusedMask);
figure(6)
imshow(newNuclearMask_2)
title('final result with erosion')