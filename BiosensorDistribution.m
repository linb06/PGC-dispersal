%%
%Code to plot intensity vs. distance from centroid of a segmented 3D cell
%group
%by Benjamin Lin 2020 in Matlab 2016a

%Input into this code are ROIs defined in imageJ with cell segmentation 
%and a background ROI for each slice for normalization. The background ROI used in this analysis was a 10x10 pixel box. 
%A single Z stack comprising the cell group to be analyzed is expected to be in the same folder

%Required m files from file exchange
    %1. sort_Nat
    %2. findFiles
    %3. ReadImageJROI

%define dir1 to the directory where the ROIs are stored
ROI_list = sort_nat(findFiles('.roi',dir1));
Istack = findFiles('.tif',dir1);

%define the pixel to micron conversion. This will vary based on the system.
xpix = .3;
ypix = .3;

%Get size of images in the stack and number of ROIs
IM = imread(Istack{1});
IM_s1 = size(IM,1);
IM_s2 = size(IM,2);
Roinum = size(ROI_list,2);

%Pre-initialize variables
Distdat = zeros(Roinum/2,40000);
Intdat = zeros(Roinum/2,40000);
Backdat = zeros(Roinum/2,1);
MaxVal = zeros(Roinum/2, 1);


%calculate the 3d centroid of the segmentation
%this assumes the stack is not over 20 slices 
IMmat = zeros(IM_s1,IM_s2,20);
for k = 1:Roinum
[sROI] = ReadImageJROI(ROI_list{k});
    % for even numbers- these ROIs with the cell segmentations
    if mod(k,2) == 1
    TotCoord = size(sROI.mnCoordinates,1);
    bw = poly2mask(sROI.mnCoordinates(1:TotCoord), sROI.mnCoordinates(TotCoord+1:TotCoord*2), IM_s1, IM_s2);
    IMmat(:,:,sROI.nPosition) = bw;  %Correct for Z drift
    else
    continue
    end
end

RowS = 0;
ColS = 0;
Zpos = 0;
TotPix = 0;
for p = 1:Roinum/2
    Coords = find(IMmat(:,:,p));
    [Row, Col] = ind2sub([IM_s1, IM_s2], Coords);
    RowS = RowS + sum(Row);
    ColS = ColS + sum(Col);
    Zpos = Zpos + size(Coords,1)*p;
    TotPix = TotPix + size(Coords,1);
end

Cent = [RowS/TotPix ; ColS/TotPix; Zpos/TotPix];


%  Roiholder = zeros(Roinum,2);

counter = 1;
counter2 = 1;
for k = 1:Roinum
[sROI] = ReadImageJROI(ROI_list{k});
    
%      Roiholder(k,1) = sROI.nPosition;
%      Roiholder(k,2) = sROI.nPosition;
    if mod(k,2) == 1
        TotCoord = size(sROI.mnCoordinates,1);
        bw = poly2mask(sROI.mnCoordinates(1:TotCoord), sROI.mnCoordinates(TotCoord+1:TotCoord*2), IM_s1, IM_s2);
        IM2 = double(imread(Istack{1},sROI.nPosition));  %Correct for Z drift
        Pix = find(bw);
        [Row, Col] = ind2sub([IM_s1, IM_s2], Pix);
   
    %Calculate the distance of each pixel from the centroid
        % euclideanDistance = sqrt((x2-x1)^2+(y2-y1)^2);
        %xpix and ypix are the pixels to micron conversion
        for p = 1:size(Pix,1)
            Distdat(counter,p) = sqrt(xpix*(Row(p)-Cent(1))^2 + ypix*(Col(p) - Cent(2))^2 + 2*(sROI.nPosition - Cent(3))^2);
            Intdat(counter,p) = IM2(Row(p), Col(p));
            
%             if Intdat(counter,p) > 0 && IntdatG(counter,p)== 0
%                IntdatG(counter,p) = Intdat(counter,p);
%             elseif IntdatG(counter,p) > 0 && Intdat(counter,p)== 0
%                Intdat(counter,p) = IntdatG(counter,p);
%             end     
        end
        
        MaxVal(counter,1) = size(Pix,1);   
        counter = counter +1;
    else
        %Find edges of the background ROI and convert back 10x10 pixel box
        %was used.
        RectCoord = sROI.vnRectBounds;
        bw2 = zeros(IM_s1, IM_s2);   
            for p = 1:11
                bw2(RectCoord(1):RectCoord(3), RectCoord(2) + p-1)= 1;
            end
        IM2 = double(imread(Istack{1},sROI.nPosition)); 
        IM3 = double(bw2).*IM2;
        Backdat(counter2,1) = mean2(nonzeros(IM3));
        counter2 = counter2 + 1;

    end

end

%Truncate zeros off initial data matrix
Distdat2 = Distdat(:,1:max(MaxVal));
Distdat3 = Distdat2;

Intdat2 = Intdat(:,1:max(MaxVal));
Intdat3 = Intdat2;



%Sort the nonzero values to be from min to max distance from the centroid
for kk = 1:Roinum/2    
    [Vals1, idx1] = sort(nonzeros(Distdat2(kk,:)), 'ascend');
    
    for pp = 1:size(Vals1,1)
        Distdat3(kk,pp) = Distdat2(kk,idx1(pp)) ;
        Intdat3 (kk,pp) = Intdat2 (kk, idx1(pp));
    end
    
end

%Normalize each row by the ratio to background at each Z slice
Intdat4 = Intdat3;

for i = 1:Roinum/2
    Intdat4(i,:) = Intdat3(i,:)/Backdat(i);
end


%Respape the distance and intensity matrix to an array and remove 0s. 
DistArray = reshape(Distdat3.',1,[]);
IntArray = reshape(Intdat4.', 1,[]);
DistArray(DistArray == 0) = [];
IntArray(IntArray == 0) = [];


%Resort the distance array from min to max from centroid
DistArray2 = DistArray;
IntArray2 = IntArray;
[Vals1, idx1] = sort(DistArray, 'ascend');

     for pp = 1:size(Vals1,2)
        
        DistArray2(pp) = DistArray(idx1(pp)) ;
        IntArray2(pp) = IntArray(idx1(pp));

     end

%Rescale the distance to 1
DistArray3 = DistArray2/(max(DistArray2));


%Bin the intensity data based on distance from centroid into 50 bins
counter = 0;
bins = zeros(50,1);

for i = 1:50
    binVals = find(DistArray3>counter & DistArray3<counter+.02);
    bins(i) = mean(IntArray2(binVals));
    counter = counter + .02;
end

%Normalized bin values
nbins = bins/max(bins);

%Visualize the distribution
plot(bins);




