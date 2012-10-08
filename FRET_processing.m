% References calculation
% GUI to set reference name and parameters that are not set in metadata
% Set filenames
% 1. Donor
% 2. Acceptor
% 3. Tandem
% load data, load metadata
% +/- check for saturated pixels
% - compensate shading correction
% subtract background
% allign XY
% set region for averaging
% calculate spectra
% normalize by exposure
% store to a database
% for tandem make before/after photobleaching spectras separately
% % at this moment should have 4 spectras
% calculate FRET positive-negtive spectra
% calculate G from it
%%
par.verbose=0;          % display intermediate steps of computation
par.disablewarnings=0;  % disable warnings - not implemented
par.emulateDV=0;        % DV on, but L.A. in "EPI-mode" (full frame)
par.shading=0;          % Shading correction processing - not implemented yet
par.saturated=0;        % Mark saturated pixels - not implemented
par.bgcorr=1;           % search for and subtract background
par.cam_offset=50;      % for PCO Sensicam VGA
par.BG_range=[0.8 1.2]; % BG to be searched within par.cam_offset*par.BG_range
par.XY_shift.enable=1;  % perform automatic image alignment
par.XY_shift.channel=0; % select channel to shift
par.XY_shift.type='nonreflective similarity';   % transformation type
par.XY_shift.interp='bicubic';                  % interpolation type
par.XY_shift.exch=1;    % excitation channel to be used for transformation function creation
% par.spectra.time1=0;    % frames to take average over time (spectra 1), 0 means take all frames
% par.spectra.time2=0;    % frames to take average over time (spectra 2), 0 means take all frames
par.spectra.range=[0.5 0.95];   % Range of intensities in max-z, sum-L for determination of useful area for spectra calculation
par.savemat=1;          % save sequences of images to mat: EfD, EfA, xD, Rt
par.saveFig=1;          % save firures to Images/ folder
par.calHist=1;          % displays histograms
par.createhtml=1;       % create html with summary of results
par.publishonline=0;    % put html to public folder to be online available


%% Read files into cells with data and OME
[filename, path]=uigetfile({'*.tif';'*.*';'*.stk'},'Select all files you want to process','MultiSelect', 'on');
fn=strcat(path, filename);

if iscell(fn)
    fnd=fn;             % fnd == file names for Donor
else fnd{1}=fn;
end
% [path, name, ext] = fileparts(filename); %split filename to path, name and extension
nf=length(fnd);         % count files
for i=1:nf
    [D{i}, OME{i}] = OME_read(fnd{i});    % donor Data for OME
end

if par.verbose
    tStart = tic;
    t = output([nf ' files received successfully'], whos, toc(tStart), 0);
end
clear filename fn i num path
%% this is a temporary part to split data into DV when FRET was not selected in software
if par.emulateDV
    for i=1:nf 
        sz=size(D{i});
        D1{i}(1:sz(1),1:sz(2)/2,1:2:(sz(3)*2-1),1:sz(4))=D{i}(1:sz(1),1:sz(2)/2,1:sz(3),1:sz(4));
        D1{i}(1:sz(1),1:sz(2)/2,2:2:(sz(3)*2),1:sz(4))=D{i}(1:sz(1),sz(2)/2+1:sz(2),1:sz(3),1:sz(4));
    end
    clear D
    D=D1;
    clear D1
end
%% shading correction
if par.shading
    ;
end
%% calculate background with maximum of the histogram
if par.bgcorr
    nbins=100;
    for i=1:nf
        sz=size(D{i});
        BG{i}=zeros(sz(3),1);
        for j=1:sz(3)
            [Dat_H, Dat_C]=diphist(D{i}(:,:,j,:), [par.cam_offset*par.BG_range(1) par.cam_offset*par.BG_range(2)],nbins); %control of maximum is required, fitting, errors, etc
            [M(1), M(2)]=max(Dat_H);
            BG{i}(j)=Dat_C(M(2));
            D{i}(:,:,j,:)=D{i}(:,:,j,:)-BG{i}(j);
        end
    end
end
%% check for saturation - need to think what to do with saturated pixels
if par.saturated
    ; % currently discarded 5% of top intensity pixels to avoid this
end
%% Automatic compensation of channel missalignment
if par.XY_shift.enable
    for i=1:nf
        Dt=sum(D{i},4);
        sz=size(D{i});
        C0=Dt(:,:,(2*par.XY_shift.exch)-1);
        C1=Dt(:,:,2*par.XY_shift.exch);
        if ~par.XY_shift.channel
            tform=autoXY_shift(C0, C1, par.XY_shift.type);
            for j=1:sz(4)
                for k=1:2:sz(3)
                    D{i}(:,:,k,j) = imtransform(D{i}(:,:,k,j), tform, 'bicubic','XData',[1 sz(2)],'YData',[1 sz(1)]);
                end
            end
        else 
            tform=autoXY_shift(C1, C0, par.XY_shift.type);
            for j=1:sz(4)
                for k=2:2:sz(3)
                    D{i}(:,:,k,j) = imtransform(D{i}(:,:,k,j), tform, 'bicubic','XData',[1 sz(2)],'YData',[1 sz(1)]);
                end
            end
        end
    end
end
%% up to this everything's OK
% now many problems started with indexes (different protocols)
% have to unify processing to get R{i} independently on 4 or 6 channels
%% load references
load('Refs.mat');
Sa=[Sa_0' Sa_1']';
Sd=[Sd_0' Sd_1']';
Saf=[Saf_0' Saf_1']';
% size(D)=XYCT
Sa(6)=[];Sa(3)=[];
Sd(6)=[];Sd(3)=[];
Sf.spectrum(6)=[];Sf.spectrum(3)=[];

for i=1:nf
    ic=[1 3 2 4];
    sz=size(D{i});
    for j=1:sz(4)
        R{i}(:,:,:,j)=nonnegative_unmix([Sd Sa Sf.spectrum'], D{i}(:,:,ic,j));
    end
end
        
