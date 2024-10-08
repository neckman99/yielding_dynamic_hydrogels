close all; clear;
file_name = 'Control-2.png';
imdata = imread(file_name);
imshow(imdata);
r = drawrectangle();
verts = uint32(r.Vertices);
imcrop = imdata((verts(1,2):verts(2,2)),(verts(1,1):verts(4,1)),:);
verts = double(verts);
area = abs(verts(1,2)-verts(2,2))*abs(verts(1,1)-verts(4,1));
imcrop = imresize(imcrop,4);
%figure;
%imshow(imcrop);
[centers, radii, metric] = imfindcircles(imcrop,[6,30],'Sensitivity',0.88);
figure; imshow(imcrop); hold on;
viscircles(centers,radii,'EdgeColor','w');
disp(['Cells per 1000 px: ',num2str(numel(centers)/(area/1000))]);