% References calculation
% GUI to set reference name and parameters that are not set in metadata
% Set filenames
% 1. Donor
% 2. Acceptor
% 3. Tandem
% load data, load metadata
% check for saturated pixels
% compensate shading correction
% subtract background
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
par.shading=0;          % Shading correction processing
par.saturated=0;        % Mark saturated pixels
par.bgcorr=1;           % search for and subtract background
par.cam_offset=50;      % for PCO Sensicam VGA
par.BG_range=[0.8 1.2]; % BG to be searched within par.cam_offset*par.BG_range
par.autoXY_shift=1;     % perform automatic image alignment
par.range=[0.5 0.95];   % Range of intensities in max-z, sum-L for determination of useful area for spectra calculation
par.savemat=1;          % save sequences of images to mat: EfD, EfA, xD, Rt
par.saveFig=1;          % save firures to Images/ folder
par.calHist=1;          % displays histograms
par.createhtml=1;       % create html with summary of results
par.publishonline=0;    % put html to public folder to be online available


%% Calculate donor spectra
[filename, path]=uigetfile({'*.tif';'*.*';'*.stk'},'Select all DONOR files you want to process as one dataset ','MultiSelect', 'on');
fn=strcat(path, filename);

if iscell(fn)
    fnd=fn;             % fnd == file names for Donor
else fnd{1}=fn;
end
% [path, name, ext] = fileparts(filename); %split filename to path, name and extension
nf=length(fnd);         % count files
for i=1:nf
    [dD{i}, dOME{i}] = OME_read(fnd{1});    % donor Data for OME
end

if par.verbose
    tStart = tic;
    t = output([nf ' DONOR files received'], whos, toc(tStart), 0);
end
clear filename fn i num path

%% this is a temporary part to split data into DV when FRET was not selected in software
if par.emulateDV
    for i=1:nf 
        sz=size(dD{i});
        dD1{i}(1:sz(1),1:sz(2)/2,1:2:(sz(3)*2-1),sz(4))=dD{i}(1:sz(1),1:sz(2)/2,1:sz(3),sz(4));
        dD1{i}(1:sz(1),1:sz(2)/2,2:2:(sz(3)*2),sz(4))=dD{i}(1:sz(1),sz(2)/2+1:sz(2),1:sz(3),sz(4));
    end
    clear dD
    dD=dD1;
    clear dD1
end
%% 
if par.bgcorr
    nbins=1000;
    for i=1:nf
        sz=size(dD{i});
        BG{i}=zeros(sz(3),1);
        for j=1:sz(3)
            [Dat_H, Dat_C]=diphist(dD{i}(:,:,j,:), [par.cam_offset*par.BG_range(1) par.cam_offset*par.BG_range(2)],nbins); %control of maximum is required, fitting, errors, etc
            [M(1), M(2)]=max(Dat_H);
            BG{i}(j)=Dat_C(M(2));
        end
    end
end
%%



if par.verbose
    figure; plot(1:nch, BG, 1:nch, Offset*BG_lim(1), '.--k', 1:nch, Offset*BG_lim, '.--k') % test how was BG calculated and if there was no limit reached
    % test histogram before I clear temporary variables
end
clear BG_* M Dat_* Offcor

if par.verbose
    t2 = output('autoXY_shift tform calculation', whos, toc(tStart2), t2);
end
