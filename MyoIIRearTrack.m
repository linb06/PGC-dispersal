%%
%Code tracking myosin II intensity at the rear of a given cell over time.
%The rear is defined here as the first 20% of rows defining the segmented
%cell. Cells were rotated to be vertical. Cells maintained a consistent
%orientation so a single angle was used. 
%The user supplies ROIs defined in ImageJ that segment a cell of interest
%in a given Z slice.

%Written by Benjamin Lin 2020 in Matlab 2016a.


%Required m files from file exchange
    %1. sort_Nat
    %2. findFiles
    %3. ReadImageJROI

%dir1 is set to directory where Rois are stored. 
%dir2 is set to directory where each image stack is located. Each stack is
%assumed to be one time point in numerical order
ROI = sort_nat(findFiles('.roi',dir1));
Fstack = sort_nat(findFiles('.tif',dir2));
Roinum = size(ROI,2);

%get image size for processing the Rois
IM = imread(Fstack{1},1);
IM_s1 = size(IM,1);
IM_s2 = size(IM,2);

%find the number of Z slices in the stack
Stackinfo = imfinfo(Fstack{1});
Finstack = size(Stackinfo,1);

%set up green and red filter for viewing segmentation of both the cell and
%posterior region
green = cat(3, zeros(IM_s1, IM_s2),ones(IM_s1,IM_s2),  zeros(IM_s1,IM_s2));
red = cat(3, ones(IM_s1, IM_s2),zeros(IM_s1,IM_s2),  zeros(IM_s1,IM_s2));




%Initialize variable to contain rear mean intensities
RearMean = zeros(40,1);

for k = 1:Roinum
    [sROI] = ReadImageJROI(ROI{k});
    Slice = sROI.vnPosition(2);
    Tstack = sROI.vnPosition(3);
    
    TotCoord = size(sROI.mnCoordinates,1);
    bw = poly2mask(sROI.mnCoordinates(1:TotCoord), sROI.mnCoordinates(TotCoord+1:TotCoord*2), IM_s1, IM_s2);
    IM2 = double(imread(Fstack{Tstack},1,Slice));
    
    
    %RotVal needs to be predefined as an array that specifies the angle to
    %rotate the given cell to a vertical position for each time point
    bwRot = double(imrotate(bw,RotVal(k), 'bicubic', 'crop'));
    IM2Rot = imrotate(IM2,RotVal(k), 'bicubic', 'crop');
    greenRot = imrotate(green,RotVal(k), 'bicubic', 'crop');
    redRot = imrotate(red, RotVal(k), 'bicubic', 'crop');
    CatS1 = size(bwRot,1);
    CatS2 = size(bwRot,2);
    
    CellInd = find(bwRot);
    CellInd2 = zeros(size(CellInd,1), 2);
    [CellInd2(:,1),CellInd2(:,2)] = ind2sub([CatS1,CatS2],CellInd);
    
    %Calculate the bottom 20% of rows to use for rear segmentation
    MaxHeight = max(CellInd2(:,1));
    MinHeight = min(CellInd2(:,1));
    RearThresh = MaxHeight - round(.2*(MaxHeight-MinHeight));
    RearInd = CellInd2(:,1) > RearThresh;
    RearInd2 = [RearInd RearInd];
    RearInd3 = RearInd2 .* CellInd2;
    RearInd4 = [nonzeros(RearInd3(:,1)) nonzeros(RearInd3(:,2))];
    ValSize = size(RearInd4,1);
    
    RearSeg = zeros(size(bwRot,1), size(bwRot,2));
    
    for p = 1:ValSize
        RearSeg(RearInd4(p,1), RearInd4(p,2)) = 1;
    end

    %Show segmentation for quality control 
    figure;
    imshow(IM2Rot,[], 'border', 'tight', 'InitialMagnification', 200);
    hold on
    h = imshow(greenRot);
    hold off
    set(h, 'AlphaData',bwRot*.3);
    hold on
    h = imshow(redRot);
    hold off
    set (h, 'AlphaData', RearSeg*.3);
    
    name = [dir1, '\', 'SegmentedCell_', num2str(k), '.tif'];
    saveas(gcf,name);
    close all;
    
    %Calculate mean intensity of posterior
    MyoIIVal = double(RearSeg).*IM2Rot;
    RearMean(k) = mean2(nonzeros(MyoIIVal));
    

 
end

%Calculate normalized mean to initial time point
RearMeanNorm = RearMean./(RearMean(1));
RearMeanFin = [RearMean RearMeanNorm];


close all;
            
            
            
            
            
            
            
            