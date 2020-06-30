%%
%Tracking Factin intensity at posterior of PGCs 
%Written by Benjamin Lin 2020 in Matlab 2016a
%Analyzes intensity of F-actin based on lifeact signal
%Assumes data is in the form of ROIs from imageJ and the series of Z stacks
%from a timelapse experiment


%Required m files from file exchange
    %1. sort_Nat
    %2. findFiles
    %3. ReadImageJROI
    
%Set dir1 to directory where Rois are held
%Set dir2 to directory where Tiff are held
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

%set up green filter for viewing segmentation
%Variables were set up for analyzing 10 timepoints.
green = cat(3, zeros(IM_s1, IM_s2),ones(IM_s1,IM_s2),  zeros(IM_s1,IM_s2));
Mint = zeros(10,1);
AreaROI = zeros(10,1);

for k = 1:Roinum
    [sROI] = ReadImageJROI(ROI{k});
    Slice = sROI.vnPosition(2);
    Tstack = sROI.vnPosition(3);

        TotCoord = size(sROI.mnCoordinates,1);
        bw = poly2mask(sROI.mnCoordinates(1:TotCoord), sROI.mnCoordinates(TotCoord+1:TotCoord*2), IM_s1, IM_s2);
        %IM2 = double(imread(MyoIIstack{k},1,p));
        IM2 = double(imread(Fstack{Tstack},1,Slice));
        
        %Display segmentation in transparent green mask with cell
        figure; 
        imshow(IM2,[], 'border', 'tight', 'InitialMagnification', 200);
                hold on
                h = imshow(green);
                hold off
                set(h, 'AlphaData',bw*.3);

        %Export visualization of segmentation for checking
        name = [dir1, '\', 'SegmentedCell_', num2str(k), '.tif']; 
        saveas(gcf,name);
        close all;
       
        %Calculate mean intensity
        IM3 = double(bw).*IM2;
        Mint(k) = mean2(nonzeros(IM3));
 
        %Rescale based on initial value and calculate total area of the
        %segementation based on the ROI.
        Mint2 = Mint./(Mint(1));
        AreaROI(k) = sum(sum(bw));
end

close all;
            
            
            
            
            
            
            
            