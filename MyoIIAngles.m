%Find angle between two lines with a common origin
%The code was used to analyze the angle of myosin II crescents relative to
%the center of a cell cluster.
%Written by Benjamin Lin 2020 in Matlab 2016a.

%Angle is calculated using the formula below
%%theta = acos(min(1,max(-1, a(:).' * b(:) / norm(a) / norm(b) )));
%Where a and b are the vectors obtained by subtracting the positions.

%Assuming there is a varible defined as Aset (angle set) where consecutive coordinates
%of a given Myosin II crescent and cell nucleus were defined in ImageJ
%using the mark points function. 

%The CalCent.m was used to calculate the centroid reference point based on
%segmentation of the whole cell stack from ROIs defined in imageJ. The
%variable Cent used here is the X,Y,Z coordinates of that position.

tot = size(Aset,1);
Angles = zeros(size(Aset,1)/2,1);
CosAngles = zeros(size(Aset,1)/2,1);
counter = 0;
for i = 1:size(Aset,1)/2
a = [Aset(i+counter,1) - Aset(i+counter+1,1), Aset(i+counter,2) - Aset(i+counter+1,2)];
b = [Aset(i+counter,1)-Cent(1), Aset(i+counter,2)-Cent(2)];
Angles(i) = acos(min(1,max(-1, a(:).' * b(:) / norm(a) / norm(b) )));
CosAngles(i) = cos(Angles(i));
counter = counter +1;
end





