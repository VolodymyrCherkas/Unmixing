
% This scripts uses D array and OME structure to make red-green picture 
% and visualise it with modified version of matVis

% You must have D array and OME structure in your workspace and
% matVis_mod.m file in your current folder

%%  Specify your parameters before executing!

avNst = 3; % number of frames to average first
spN = 3; % space between an averaged frames
avNfin = 3; % number of frames to average last 

%%
RG = zeros(480,320);

for j=1:OME.nC
    for i = 1:OME.nT-avNfin-spN-avNst+1
        RG(:,:,j,i) = mean(D(:,:,j,i+avNst+spN-1:i+avNst+spN+avNfin-1),4) - mean(D(:,:,j,i:i+avNst-1),4) ;
    end
end

diff_par = [avNst spN avNfin]; %save parameter values to array
n_fr_RG = i; % number of frames in diff movie

clear i j spN avNst avNfin;

%%
% use averaging filter to get less noise

h=ones(3,3)/9;
RG = imfilter(RG,h);

clear h;

%%
% visualise with matVis

matVis_mod(RG);
