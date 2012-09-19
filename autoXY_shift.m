function tform = autoXY_shift(C0,C1,fun)  
%% Function to calculate precise XY transformation fully automatically
% tform = autoXY_shift(C0,C1,fun)  finds a transformation function C0 to C1
% fun = 'nonreflective similarity', 'similarity', 'affine', 'projective',
% 'polynomial', 'piecewise linear', 'lwm'
% 
% should work correctly with uint8, uint16, single, double, speed independent, so good reason to use double
% would appreciate having background removed from input images for faster
% 
% processing and easier thresholding% is not initially designed to work with shift of more then 5 pixels (can hardly be
% corrected, but speed will dramatically degrade), for more then 5 pixels
% shift can be be used twice: first with binned images (e.g. 3 by 3), then
% calculate transformation function 1,  shift data, and apply second
% loop of transformation in normal way (also possible to use manual "preshift") 
%
% Minor bug fixes
% Version 2.1 (04.04.2011, by Volodymyr Cherkas)


%%
% create a grid of pixels to process correlation optimization (checked,
% works correctly even with M by N nonsquare images

nw(1)=floor((size(C0,1)-20)/5); % n in web (x y)
nw(2)=floor((size(C0,2)-20)/5);

ii(1,:)=0:(nw(1)*nw(2)-1); %max number of pixels in web
grid(ii+1,2)=floor(ii/nw(1))*5+10;
grid(ii+1,1)=mod(ii,nw(1))*5+10;

% Im=zeros(size(C0,1),size(C0,2));
% Im(grid(ii+1,1),grid(ii+1,2))=1; figure; imshow(Im);
%%
% optimize web with local maximum of track2 (+/-2 pixels to cover complete area)
Threshold=max(max(C1))/150; % shoud sometimes think about how to set depending on image received
web=zeros(length(ii),2);
for i=1:length(ii);  
    [co(:,1) co(:,2)]=max(C1((grid(i,1)-2):(grid(i,1)+2), (grid(i,2)-2):(grid(i,2)+2)),[],2); % y first, x second
    [ci(1) ci(2)]=max(co(:,1));
    web(i,1)=double(ci(1)>Threshold).*(grid(i,1)+ci(2)-3);
    web(i,2)=double(ci(1)>Threshold).*(grid(i,2)+co(ci(2),2)-3);
end
web((web(:,2)==0),:)=[]; % Thanks a lot to Andre!!!
web=flipdim(web,2);
webcorr=cpcorr(web,web,C0,C1); % maxumum correlation function
webcorr(:,3)=webcorr(:,1)-web(:,1);
webcorr(:,4)=web(:,1);
webcorr(:,5)=web(:,2);
webcorr((webcorr(:,3)==0),:)=[];
tform= cp2tform(webcorr(:,1:2), webcorr(:,4:5), fun); % transformation function creation
% C0c = imtransform(C0, tform, 'bicubic','XData',[1 size(C0,2)],'YData',[1 size(C0,1)]);