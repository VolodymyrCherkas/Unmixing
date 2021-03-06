%% new references calculator (started 2012.11.21)
% the main difference from the previous version is to use autofluorescence
% as a separate component for unmixing. This should cause higher precision
% and lower errors in low-signal recordings.

% GUI to set reference name and parameters that are not set in metadata
% Set filenames
% 1. Donor
% 2. Acceptor
% 3. Tandem
% load data, load metadata
% read channel configuration;
% reassign channels orger;
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
par.emulateDV=1;        % DV on, but L.A. in "EPI-mode" (full frame)
par.wlsort='lambda';    % sort excitation wawelenghts with 'lambda', 'time'
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
[filename, path]=uigetfile({'*.tif';'*.*';'*.stk'},'Select all DONOR files you want to process as one dataset ','MultiSelect', 'on');
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
    t = output([nf ' DONOR files received'], whos, toc(tStart), 0);
end
clear filename fn i num path
%% this is a temporary part to split data into DV when FRET was not selected in software
for i=1:nf
    if par.emulateDV&&(~OME{i}.fret)
        sz=size(D{i});
        D1{i}(1:sz(1),1:sz(2)/2,1:2:(sz(3)*2-1),1:sz(4))=D{i}(1:sz(1),1:sz(2)/2,1:sz(3),1:sz(4));
        D1{i}(1:sz(1),1:sz(2)/2,2:2:(sz(3)*2),1:sz(4))=D{i}(1:sz(1),sz(2)/2+1:sz(2),1:sz(3),1:sz(4));
        % add changes to OME to fit changes in dataset D
    else
        D1{i}=D{i}; % think about optimization to skip unnecessary steps
    end
end
clear D
D=D1;
clear D1
%% sort wavelenghts
if strcmp(par.wlsort,'lambda')
    for i=1:nf
        [OME{i}.wlconf, reallocate]=sort(OME{i}.wlconf);
        OME{i}.time_ch=OME{i}.time_ch(:,reallocate);
        OME{i}.expconf=OME{i}.expconf(reallocate);
        D{i}(:,:,:,:)=D{i}(:,:,reallocate,:);
    end
end
%% shading correction
if par.shading
    ; % shading correction may increase moise, but make better offset and background correction
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
        Dtt=sum(D{i},4);
        Dt(:,:,1)=sum(Dtt(:,:,2:2:OME{i}.nC),3);   % sum all excitation channels to have universal alignment mechanism
        Dt(:,:,2)=sum(Dtt(:,:,1:2:OME{i}.nC),3);
        clear Dtt;
        sz=size(D{i});
        C0=Dt(:,:,2); % (2*par.XY_shift.exch)-1);
        C1=Dt(:,:,1); % 2*par.XY_shift.exch);
        if ~par.XY_shift.channel
            tform=autoXY_shift(C0, C1, par.XY_shift.type);
            for j=1:sz(4)
                for k=1:2:sz(3)
                    D{i}(:,:,k,j) = imtransform(D{i}(:,:,k,j), tform, par.XY_shift.interp,'XData',[1 sz(2)],'YData',[1 sz(1)]);
                end
            end
        else 
            tform=autoXY_shift(C1, C0, par.XY_shift.type);
            for j=1:sz(4)
                for k=2:2:sz(3)
                    D{i}(:,:,k,j) = imtransform(D{i}(:,:,k,j), tform, par.XY_shift.interp,'XData',[1 sz(2)],'YData',[1 sz(1)]);
                end
            end
        end
    end
end
%% set region for averaging
for i=1:nf
    M=sum(sum(D{i},4),3);   % sum time axis, sum channels (this will make unequal weights of channels)
    Mask{i}=(M>(max(M(:))*par.spectra.range(1))).*(M<(max(M(:))*par.spectra.range(2))); % set Mask using par.spectra.range
    sz=size(D{i});
    for j=1:sz(4)
        for k=1:sz(3)
            S{i}(k,j)=sum(sum(D{i}(:,:,k,j).*Mask{i}))/sum(Mask{i}(:)); % calculate Lambda(time)
        end
    end
    Spectra{i}=mean(S{i},2);    % average over time to have good S/N and low photobleaching
end


