%% Main processing script for FRET and "concentrations" calculations
%% set parameters to par.
par.calculateref=0;     % calculate reference spectras (need additional input to select references)
% "0" will load precalculated spectra, "1" will ask for references files
par.checkref=1;         % visualise spectras and check (plot references after calculations)
% par.processall=1;       % procass all files in subfolders of current
% folder (1 level) useful for in-line tree parsing
par.verbose=0;          % display intermediate steps of computation
par.disablewarnings=1;  % disable warnings
par.autoXY_shift=1;     % perform automatic image alignment
par.savemat=1;          % save sequences of images to mat: EfD, EfA, xD, Rt
par.saveFig=1;          % save firures to Images/ folder
par.calHist=1;          % displays histograms
par.createhtml=1;       % create html with summary of results
par.publishonline=0;    % put html to public folder to be online available
par.range=[0.5 0.95];   % Range of intensities in max-z, sum-L for determination of useful area for spectra calculation
par.DV=1;               % "FRET mode" in L.A. (frames split into sequential)
par.emulateDV=0;        % DV on, but L.A. in "EPI-mode" (full frame)
par.verbose=1;
par.nchannels=3;        % excitation
par.cam_offset=50;      % for PCO Sensicam VGA

% obsolete parameters
% par.Qd=0.40;            % Quantum efficiency of Donor
% par.Qa=0.61;            % Quantum efficiency of acceptor


%% References analysis: Donor, Acceptor, Tandem1, Tandem2 (alternative Tandem+Tandem bleaching)
% DONOR
[filename, path]=uigetfile({'*.tif';'*.*';'*.stk'},'Select all DONOR files you want to process as one dataset ','MultiSelect', 'on');
fn=strcat(path, filename);

if iscell(fn)
    fnd=fn;             % fnd == file names for Donor
else fnd{1}=fn;
end
% [path, name, ext] = fileparts(filename); %split filename to path, name and extension
nf=length(fnd);         % count files
for i=1:nf
    info_d{i} = imfinfo(fnd{1});
    num_images_d{i} = numel(info_d{i}); %count the number of images in each file - does not work properly!!!
end
if par.verbose
    tStart = tic;
    t = output([nf ' DONOR files received'], whos, toc(tStart), 0);
end
clear filename fn i nf num path

%% ACCEPTOR
[filename, path]=uigetfile({'*.tif';'*.*';'*.stk'},'Select all ACCEPTOR files you want to process as one dataset ','MultiSelect', 'on');
fn=strcat(path, filename);

if iscell(fn)
    fna=fn;             % fna == file names for Acceptor
else fna{1}=fn;
end
% [path, name, ext] = fileparts(filename); %split filename to path, name and extension
nf=length(fna);         % count files
for i=1:nf
    info_a{i} = imfinfo(fna{1});
    num_images_a{i} = numel(info_a{i}); %count the number of images in each file - does not work properly!!!
end
if par.verbose
    t = output([nf ' ACCEPTOR files received'], whos, toc(tStart), t1);
end
clear filename fn i nf num path

%% TANDEM DONOR-ACCEPTOR v1
[filename, path]=uigetfile({'*.tif';'*.*';'*.stk'},'Select all TANDEM v1 files you want to process as one dataset ','MultiSelect', 'on');
fn=strcat(path, filename);

if iscell(fn)
    fnt1=fn;             % fnt1 == file names for Tandem v1
else fnt1{1}=fn;
end
% [path, name, ext] = fileparts(filename); %split filename to path, name and extension
nf=length(fnt1);         % count files
for i=1:nf
    info_t1{i} = imfinfo(fnt1{1});
    num_images_t1{i} = numel(info_t1{i}); %count the number of images in each file - does not work properly!!!
end
if par.verbose
    t = output([nf ' TANDEM v1 files received'], whos, toc(tStart), t1);
end
clear filename fn i nf num path

%% 

% one should think about reading xml metadata before continue (to have
% precise time for each frame (this should also be used to separate channels)
% and exposure time.

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
        
        ex=mod((f+par.nchannels-1), par.nchannels)+1;
        frame=floor((f+par.nchannels-1)/par.nchannels);
        F(:,:,1, ex, frame)=C0c;           % (y,x,ch_em,ch_ex, time)
        F(:,:,2, ex, frame)=FT(1:y, (x/2+1):x);
        f=f+1;
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
