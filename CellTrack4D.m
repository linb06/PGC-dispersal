%%
%Use for tracking individual cells from a 4D image sequence
%This code assumes the cell is well separate in XY and Z to permit use of a
%single ROI for segmentation.

%Written by Benjamin Lin 2020 in Matlab 2016a

%Required m files from file exchange
    %1. sort_Nat
    %2. findFiles
    %3. ReadImageJROI
    %4. func_threshold
    

%Set dir1 to directory where the ROIs are located
%Set dir2 to the directory where the tiff stacks are located. Tiff stacks
%are assumed to be in a temporal sequence in numerical order.
ROI = sort_nat(findFiles('.roi',dir1));
MyoIIstack = sort_nat(findFiles('.tif',dir2));
Roinum = size(ROI,2);

%get image size for processing the Rois
IM = imread(MyoIIstack{1},1);
IM_s1 = size(IM,1);
IM_s2 = size(IM,2);

%find the size of the stack
Stackinfo = imfinfo(MyoIIstack{1});
Finstack = size(Stackinfo,1);

%set up green filter for viewing segmentation
green = cat(3, zeros(IM_s1, IM_s2),ones(IM_s1,IM_s2),  zeros(IM_s1,IM_s2));

%PreAllocate Centroid matrix for Cell center and for Patch center
CellCentMat = zeros(Roinum,3);

%Structuring element for getting cell perimeter
S = strel('disk',8);

%Start = initial z position -1
Zstart = 0;

%Define ZT as the top Z slice position over time
for k = 1:Roinum
    [sROI] = ReadImageJROI(ROI{k});
    Firstslice = sROI.vnPosition(2);
    
    %The assumption here is that the ROI was defined in the top slice where
    %the cell was present. Another assumption is that the cell of interest
    %is 5 slices in depth
    Lastslice = Firstslice + 4;
    if Lastslice > Finstack
        Lastslice = Finstack;
    end
    
    Imholder = zeros(IM_s1, IM_s2, 30);
    counter = 1;
    figure;
    ha = tight_subplot(1,5,0,0,0);
    
    %Segment cell in each Z slice
    for p = Firstslice:Lastslice
        
        TotCoord = size(sROI.mnCoordinates,1);
        bw = poly2mask(sROI.mnCoordinates(1:TotCoord), sROI.mnCoordinates(TotCoord+1:TotCoord*2), IM_s1, IM_s2);
        IM2 = double(imread(MyoIIstack{k+Zstart},1,p));
        IM3 = double(bw).*IM2;
        
        %Refine segmentation
        Thresh = func_threshold(nonzeros(IM3));
        SegCell = IM3 > Thresh;
        SegCell2 = bwareaopen(SegCell,50);
        
        SegCell3 = imfill(SegCell2, 'holes');
        SegCell3 = imdilate(SegCell3,S);
        SegCell3 = imerode(SegCell3,S);
        
        %Evaluate the segmentation
        axes(ha(counter));
        imshow(IM2,[], 'border', 'tight', 'InitialMagnification', 200);
        hold on
        h = imshow(green);
        hold off
        set(h, 'AlphaData',SegCell3.*.3);
        ax = gca;
        

        Imholder(:,:,p) = SegCell3;
        counter = counter + 1;
        
    end
    %Export the segmentation for quality control
    name = [dir1, '\', '4dseg_', num2str(k), '_', num2str(p), '.tif'];
    saveas(gcf,name);
    close all;
    
    %Calculate centroid based on 3D segmented cell
    RowS = 0;
    ColS = 0;
    Zpos = 0;
    TotPix = 0;
    for m = Firstslice:Lastslice
        Coords = find(Imholder(:,:,m));
        [Row, Col] = ind2sub([IM_s1, IM_s2], Coords);
        RowS = RowS + sum(Row);
        ColS = ColS + sum(Col);
        Zpos = Zpos + size(Coords,1)*m;
        TotPix = TotPix + size(Coords,1);
    end
    
    %Calculate centroid and save it
    Cent = [RowS/TotPix ; ColS/TotPix; Zpos/TotPix];
    CellCentMat(k,:) = Cent;
    

   
end

%Pixel to micron conversion depending on the system
xpix = .389;
ypix = .389;
Zpix  = 2;

%Calculate cell speed
DistMat = zeros(Roinum-1,1);
for r = 2:Roinum;
    Xdif = abs((CellCentMat(r,1) - CellCentMat(r-1,1))*xpix);
    Ydif = abs((CellCentMat(r,2) - CellCentMat(r-1,2))*ypix);
    Zdif =  abs((CellCentMat(r,3) - CellCentMat(r-1,3))*zpix);
    DistMat(r-1,1) = Xdif + Ydif + Zdif;
end

%Calculate cell persistence
Persistence = zeros(3,1);
Persistence(1) = sum(DistMat);
 Xdif = abs((CellCentMat(Roinum,1) - CellCentMat(1,1))*xpix);
    Ydif = abs((CellCentMat(Roinum,2) - CellCentMat(1,2))*ypix);
    Zdif =  abs((CellCentMat(Roinum,3) - CellCentMat(1,3))*zpix);
    Persistence(2) = Xdif + Ydif + Zdif;
    Persistence(3) = Persistence(2)/Persistence(1);

%Export relevant parameters
name = [dir1 '\' '4Dtrack.xls'];
xlswrite(name,CellCentMat, 'Cell_Centroid');
xlswrite(name,DistMat, 'Distance');
xlswrite(name,Persistence, 'Persistence');
