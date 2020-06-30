%
%Code for tracking membrane to cytoplasm ratio in a given cell over time.
%This was used to track the Ecadherin membrane to cytoplasm ratio over
%time. Input into this code are Rois from ImageJ which segment the cell of
%interest in a given Z slice.

%Written by Benjamin Lin 2020 in Matlab 2016a.

%Required m files from file exchange
    %1. sort_Nat
    %2. findFiles
    %3. ReadImageJROI

%set dir1 to directory where the Rois from ImageJ are located
%set dir2 to directory where the tiff file stacks are located. Assuming
%they are in numerical order.
%%
ROI = sort_nat(findFiles('.roi',dir1));
Ecad = sort_nat(findFiles('.tif',dir2));
Roinum = size(ROI,2);

%get image size for processing the Rois
IM = imread(Ecad{1},1);
IM_s1 = size(IM,1);
IM_s2 = size(IM,2);

%set up green filter for viewing segmentation
green = cat(3, zeros(IM_s1, IM_s2),ones(IM_s1,IM_s2),  zeros(IM_s1,IM_s2));

%The define the number of timepoints to analyze
T = 16;
%PreAllocate Centroid matrix for Cell center and for Patch center
CellCentMat = zeros(T,2);
EcadRatio = zeros(T,1);
MemMeanMat = zeros(T,1);
CytoMeanMat = zeros(T,1);


%Structuring element for getting cell perimeter
S = strel('disk',5);
S2 = strel('disk',1);

for k = 1:Roinum
    [sROI] = ReadImageJROI(ROI{k});
    TotCoord = size(sROI.mnCoordinates,1);
    bw = poly2mask(sROI.mnCoordinates(1:TotCoord), sROI.mnCoordinates(TotCoord+1:TotCoord*2), IM_s1, IM_s2);
    IM2 = double(imread(Ecad{k},sROI.vnPosition(2)));
    
    %Smooth the segmentation
    Mod1 = imerode(bw,S2);
    SegCell3 = double(Mod1);

    %Show how cell the segmentation worked
    fig = figure; imshow(IM2,[], 'border', 'tight', 'InitialMagnification', 200);
    hold on
    h = imshow(green);
    hold off
    set(h, 'AlphaData',SegCell3.*.3);
    ax = gca;
    
    
    %Export the segmentation for quality control
    name = [dir1, '\', 'SegmentedCell_', num2str(k), '.tif'];
    saveas(gcf,name);
    close all;
    
    %Extract Centroid of cell and save it
    CellCent = regionprops(SegCell3, 'centroid');
    CellCentMat(k,1) = CellCent.Centroid(1);
    CellCentMat(k,2) = CellCent.Centroid(2);
    
    %Segment cell periphery and cytoplasm and take the mean values
    Mem1 = imerode(SegCell3,S);
    CytoMean = mean2(nonzeros(double(Mem1).* IM2));
    Mem2 = SegCell3 - double(Mem1);
    MemMean = mean2(nonzeros(double(Mem2).* IM2));
    
    %Show how cell the segmentation worked
    fig = figure; imshow(IM2,[], 'border', 'tight', 'InitialMagnification', 200);
    hold on
    h = imshow(green);
    hold off
    set(h, 'AlphaData',Mem2.*.3);
    
    %Export the segmentation for quality control
    name = [dir1, '\', 'Mem_', num2str(k), '.tif'];
    saveas(gcf,name);
    close all;
    
    EcadRatio(k) = MemMean/CytoMean;
    MemMeanMat(k) = MemMean;
    CytoMeanMat(k) = CytoMean;
       
end

%Export relevant stats to excel file.
name = [dir1 '\' 'EcadStats_Z.xls'];
xlswrite(name,CellCentMat, 'CellCentroid');
xlswrite(name,EcadRatio,'EcadMemCyto');
xlswrite(name,MemMeanMat, 'Mem_Mean');
xlswrite(name, CytoMeanMat, 'Cyto_Mean');

