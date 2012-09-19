function [D, OME]=OME_read(filename)
%% image splitter
% OME.filename is <path>\<filename>.xml, unless other specified
[path, name, ext_im] = fileparts(filename); %split filename to path, name and extension
ext='.xml';
OME_filename=strcat(path, '/' , name, ext); % Windows accepts both forward and backward slashes

if ~exist(OME_filename,'file')
    [t_return_status, t_result]=system([cd, '/XMLtools/tiffcomment ',filename,' > ',OME_filename]);
    if t_return_status
        errordlg(strcat('XML_writer:tiffcomment:tiffcoment returned ',num2str(t_return_status)));
        return
    end
end

OME=OME_info(OME_filename); % reads xml to a convenient structure

D=zeros(OME.nY, OME.nX, OME.nC, OME.nT, OME.nZ);    % Z-stack is not implemented yet
for i=1:length(OME.frame)
    D(:,:,OME.exch(i)+1,OME.cycle(i))=imread(filename,i);
end
% i=1:OME.nT;
% for j=1:OME.nC
%     time(i,j)=OME.time((i-1)*OME.nC+j);
% end

%% step-by-step plan
%- 0. calculate and verify references, write "TFP", "YFP", "TFP-YFP1", "TFP-YFP2"
% each of this contains: .wl, .exp, emch, and .spectrum
%+ 1. get input filenames
%+ 2. get optional parameters (optical properties, configuration, etc)
%(GUI) OME_info.m
%- 3. select proper calibration files and load them (GUI-selected spectra)
%+ 4. read data (incl. background calculation and removal)
%+ 5. xy correct
%+ 6. unmix (Donor concentration, Acceptor, and FRET(need to think in donor or acceptor values, currently in Donor))
%+ 7. write tiff's, (or .mat's)
%- 8. plot Efa (Efd) to R (and I) in count plot
%+/- 9. plot intensity-saturation plots for FRET
%- 10. select regions based on count plots, check different features (cluster analysis)
%- 11. plot averaged time traces for regions with similar features
