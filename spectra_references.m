% this is to read pure Donor and Acceptor files exported directly from
% TillVision without any metadata (all necessary information is to be
% provided manually)
%%
[filename, path]=uigetfile({'*.tif';'*.*';'*.stk'},'Please select 435/505/575nm channel files to process','MultiSelect', 'on');
fn=strcat(path, filename);
if iscell(fn)
    fnc=fn;
else fnc{1}=fn;
end
% [path, name, ext] = fileparts(filename); %split filename to path, name and extension
nf=length(fnc); %count files
for i=1:nf
    info{i} = imfinfo(fnc{1});
    num_images{i} = numel(info{i}); %count the number of images in each file
end
% one should think about reading xml metadata before continue (to have
% precise time for each frame (this can also be used to separate channels)
% and exposure time.

% 5/5/5ms exposures from TillVision for HP-YFP

FT=double(imread(fnc{1},1)); % read the first frame of the first file to create split function (for DV)
x=size(FT,2); % should take x and y from info{:}(:)
y=size(FT,1);
C0=FT(1:y, 1:x/2);
C1=FT(1:y, (x/2+1):x);
tform=autoXY_shift(C0, C1,'nonreflective similarity'); % registration structure
% C0c = imtransform(C0, tform, 'bicubic','XData',[1 size(C0,2)],'YData',[1 size(C0,1)]); % test for verbose mode
% F1=double(imread(Acceptor_fn,1));
% C1=FT(1:480, 321:640);
clear FT C0 C1;
% now should have list of files with image size, channel number, time
% intervals, exposure, wavelength, etc. It would be great to have all of this...
par.verbose=1;
par.nchannels=3; %excitation
par.cam_offset=50; % for PCO Sensicam VGA
%% read data
% desired data type: F(y,x,channel(x1:em0-em1,x2:em0-em1,x3:em0-em1),frame)
clear F f ex frame C0 C0c;
% F{1,1,1}=zeros(y,x); %filenumber, channel, frame
f=1;
for fi=1:nf
    for i=1:num_images{fi}
%         FT=double(imread(fn{fi},i,'Info',info{fi})); % Comment to TILL - bad TIFF offset field!!!
        FT=double(imread(fnc{fi},i));
        C0=FT(1:y, 1:x/2);
        C0c = imtransform(C0, tform, 'bicubic','XData',[1 size(C0,2)],'YData',[1 size(C0,1)]);
%         F{fi,1,i}=C0c;
%         F{fi,2,i}=FT(1:y, (x/2+1):x); %now have F{filenumber, emission, index}(y,x) index is excitation-time
        
%         ex=mod((f+par.nchannels-1), par.nchannels)+1;
%         frame=floor((f+par.nchannels-1)/par.nchannels);
        F(:,:,1, fi, i)=C0c;           % (y,x,ch_em,ch_ex, time) fi (file index) is an excitation index for separate files, frame is an index i
        F(:,:,2, fi, i)=FT(1:y, (x/2+1):x);
%         f=f+1;
    end
end
% F(y,x,L,l,t);
% I would convert data now to F(y,x,em,ex,frame)
% Sure I need time stamps for each frame.

% Ftemp=F{1,1,1};
% Ftemp(:,:,2)=F{1,2,1};
% matVis(Ftemp);
sz=size(F);
F=permute(reshape(F,[sz(1) sz(2) sz(3)*sz(4) 1 sz(5)]),[1 2 3 5 4]);


%%
% AutoBackground as maximum of histogram with offset +/-20%
sz=size(F);
nch=sz(3); % duplicate size variables for convenience; nch - number of channels
BG_lim=[0.8 1.2];
BG=zeros(nch,1);
Offset=par.cam_offset;
nbins=1000;
for j=1:nch; %number of L channels
    [Dat_H, Dat_C]=diphist(F(:,:,j,:), [Offset*BG_lim(1) Offset*BG_lim(2)],nbins); %control of maximum is required, fitting, errors, etc
    [M(1), M(2)]=max(Dat_H);
    BG(j)=Dat_C(M(2));
end
clear BG_* M Dat_*
% I should find clipping points before BG removal!!!
F(F>4094)=0;
%% Remove background
parfor j=1:nch;
    F(:,:,j,:)=F(:,:,j,:)-BG(j); %background removal
end

%% set calculation region limits
Fm=mean(F,4);
% mask should be equal weight for both channels
Max=permute(max(max(Fm,[],1),[],2),[3 1 2]);
Fmask=zeros(sz(1), sz(2));
for i=1:nch
    Fmask=Fmask+Fm(:,:,i)/Max(i);
end


