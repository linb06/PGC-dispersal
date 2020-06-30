%Code to calculate the centroid of a 3D cell cluster based on ROIs defining
%the stack in each Z slice from ImageJ.
%Written by Benjamin Lin 2020 in Matlab 2016a.

%Required m files from file exchange
    %1. sort_Nat
    %2. findFiles
    %3. ReadImageJROI
    
%%set dir1 to directory where the ROIs are present.
%The code assumes the Z stack tiff is in the same folder
ROI = sort_nat(findFiles('.roi',dir1));
Stack = findFiles('.tif',dir1);

%Find size of image for reconstructing the ROI
IM = imread(Stack{1});
IM_s1 = size(IM,1);
IM_s2 = size(IM,2);
Roinum = size(ROI,2);


%Intialize matrix that assumes there are no more than 30 slices to this
%stack
IMmat = zeros(IM_s1,IM_s2,30);


for k = 1:Roinum
    [sROI] = ReadImageJROI(ROI{k});
    TotCoord = size(sROI.mnCoordinates,1);
    bw = poly2mask(sROI.mnCoordinates(1:TotCoord), sROI.mnCoordinates(TotCoord+1:TotCoord*2), IM_s1, IM_s2);
    IMmat(:,:,sROI.nPosition) = bw;
end

RowS = 0;
ColS = 0;
Zpos = 0;
TotPix = 0;
for p = 1:30
    Coords = find(IMmat(:,:,p));
    [Row, Col] = ind2sub([IM_s1, IM_s2], Coords);
    RowS = RowS + sum(Row);
    ColS = ColS + sum(Col);
    Zpos = Zpos + size(Coords,1)*p;
    TotPix = TotPix + size(Coords,1);
end

Cent = [ColS/TotPix;RowS/TotPix; Zpos/TotPix];

%%