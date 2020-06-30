%
%Code for analyzing myosin II hot spots in single cells over time. Calculates presence of spots and orientation relative to a fixed point. 
%In this case, the center of endoderm was used as the reference point. The
%centroid of the cell of interest is also determined here for cell tracking
%for quantifications of cell speed and persistence

%Written by Benjamin Lin 2020 in Matlab 2016a.
%Input into this code is a set of ROIs which roughly outline a cell of
%interest. The segmentation is refined by thresholding. 

%Required m files from file exchange
    %1. sort_Nat
    %2. findFiles
    %3. ReadImageJROI
    %4. func_threshold

% Set dir1 to the directory where the Rois from ImageJ are located
%Set dir2 to the directory where the tif file where Cell images are. In
%this analysis I utilzied max intensity projections from cells expressing
%Myo-II GFP.
%%
ROI = sort_nat(findFiles('.roi',dir1));
MaxMyoII = findFiles('.tif',dir2);
Roinum = size(ROI,2);

%get image size for processing the Rois
IM = imread(MaxMyoII{1},1);
IM_s1 = size(IM,1);
IM_s2 = size(IM,2);

%set up green filter for viewing segmentation
green = cat(3, zeros(IM_s1, IM_s2),ones(IM_s1,IM_s2),  zeros(IM_s1,IM_s2));


%PreAllocate Centroid matrix for Cell center and for Patch center
CellCentMat = zeros(Roinum,2);
MyoIIPatchArea = zeros(Roinum,1);
MyoIICentMat  = zeros(Roinum,2);
MyoIIAngle = ones(Roinum,2).*10;
NumPatches = zeros(Roinum,1);

%Initialize structuring element for getting cell perimeter
S = strel('disk',6);


for k = 1:Roinum
        [sROI] = ReadImageJROI(ROI{k});
        TotCoord = size(sROI.mnCoordinates,1);
        bw = poly2mask(sROI.mnCoordinates(1:TotCoord), sROI.mnCoordinates(TotCoord+1:TotCoord*2), IM_s1, IM_s2);
        IM2 = double(imread(MaxMyoII{1},sROI.nPosition));
        IM3 = double(bw).*IM2;
    
    %Refine segmeantation
        Thresh = func_threshold(nonzeros(IM3));
        SegCell = IM3 > Thresh;
        SegCell2 = bwareaopen(SegCell,50);
        SegCell3 = imfill(SegCell2, 'holes');
        SegCell3 = imdilate(SegCell3,S);
        SegCell3 = imerode(SegCell3,S);

    
    %Show how well the segmentation worked
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
   
    %Extract Centroid of cell and save it. This is used as a reference
    %point to calculate orientation of the myosin II hot spot.
         CellCent = regionprops(SegCell3, 'centroid');
         CellCentMat(k,1) = CellCent.Centroid(1);
         CellCentMat(k,2) = CellCent.Centroid(2);
         
    %Segment cell periphery and threshold to find peripheral hotspots 
         Mem1 = imerode(SegCell3,S);
         Mem2 = SegCell3 - double(Mem1);
         CellMem1 = Mem2 .* IM3;
         Patch1 = CellMem1 > (mean2(nonzeros(CellMem1)) + 1.5*std(nonzeros(CellMem1)));
         
         %Remove patches with area less than 10 pixels to remove noise
         Patch2 = bwareaopen(Patch1,10);
         
         
    %Show how cell the segmentation worked
         fig = figure; imshow(IM2,[], 'border', 'tight', 'InitialMagnification', 200);
         hold on
         h = imshow(green);
         hold off
         set(h, 'AlphaData',Patch2.*.3);
         
    %Export the segmentation for quality control
         name = [dir1, '\', 'MyoIIPatch_', num2str(k), '.tif'];
         saveas(gcf,name);
         close all;
         
    %If a hotspot was detected, calculate the angle relative to the cell centroid and a fixed point.
    %Cent is defined as the X,Y,Z coordinates of that fixed point.
         if sum(sum(Patch2)) > 0
       
             MyoPatch = regionprops(Patch2);
             NumPatches(k,1) = size(MyoPatch,1);
     
            %If multiple hot spots are detected, the largest hotspot was
            %used for analysis
             [M,I] = max(extractfield(MyoPatch, 'Area'));
             MyoIIPatchArea(k) = MyoPatch(I).Area;
             MyoIICentMat(k,1) = MyoPatch(I).Centroid(1);
             MyoIICentMat(k,2) = MyoPatch(I).Centroid(2);
           
            %Calculate angle the myosin II hotspot relative to cell centroid and a fixed point (cent) 
             a = [CellCent.Centroid(1) - MyoPatch(I).Centroid(1), CellCent.Centroid(2) - MyoPatch(I).Centroid(2)];
             b = [CellCent.Centroid(1) - Cent(1), CellCent.Centroid(2) - Cent(2)];
             Angles = acos(min(1,max(-1, a(:).' * b(:) / norm(a) / norm(b) )));
             CosAngles = cos(Angles);
             
             MyoIIAngle(k,1) = Angles;
             MyoIIAngle(k,2) = CosAngles;
         end
         


   
end

%Export relative information to an excel file
name = [dir1 '\' 'Myo2PatchStats.xls'];
xlswrite(name,CellCentMat, 'Cell_Centroid');
xlswrite(name,MyoIIPatchArea,'MyoII_Patch_Area');
xlswrite(name,MyoIICentMat, 'MyoII_Centroid');
xlswrite(name, MyoIIAngle, 'MyoII_Angles');
xlswrite(name, NumPatches, 'NumPatches');
