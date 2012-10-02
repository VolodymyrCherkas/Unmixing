% FUC=FRET Universal Calculator
% tic
%% User Adjustible Parameters
%they start with p_
p_path='/media/Extended/!DATA/Gena/!!Projects/Hippocalcin/Processing_new/PC12/ECFP/08062012/';   %this must end with a /
p_fname='LA_ProtocolCombination_08.06.2012_15_46_16.tif';
p_roi_output_file_name='a_out.3D';
%% init
t_fname=strcat(p_path,p_fname);
if ~exist(t_fname,'file')
    errordlg(strcat('FUC.m: ',t_fname,' file does not exist'));
    return
end
%% extract and parse metadata
%prepare the full .xml filename string to use with extractor and parser
[~,t_name,~]=fileparts(p_fname);
t_xml=strcat(p_path,t_name,'.xml');
if ~exist(t_xml,'file')
    %TODO: platform dependent execution (using the 'computer' function)
    %{
    %Prepare the command line string to be called with the unix()
    %We use the tiffcomment tool from here
    %http://loci.wisc.edu/bio-formats/command-line-tools
    %read the instructions carefully - you need 2 downloads
    %}
    t_unix_tiffcoment_command=strcat('tiffcomment ''',t_fname,''' > ''',t_xml,'''');
    t_return_status=unix(t_unix_tiffcoment_command);
    %     if something goes wrong display the error box and halt
    if t_return_status
        errordlg(strcat('FUC.m:tiffcomment:tiffcoment returned ',num2str(t_return_status)));
        return
    end
end
% if all OK parse the file with the xml2struct function found here
% http://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
t_OME_struct=xml2struct(t_xml);
%% unmarshall

[A,B,C,D,E,F]=f_load(t_fname,'123456');
%% register
%{
find the transform that when applied to the first frame of B will
  give the most correlation with the first frame of A
  needs to be done once per DualView adjustment. (If you are careful
  enough this means once per experiment, this is what I assume here).
Gena Madan

%}

t_transform=autoXY_shift(B(:,:,1),A(:,:,1),'nonreflective similarity');
%     save('t_transform.mat','t_transform');
% now let's apply this transform to every frame in B and D
t_width=size(B,2);
t_height=size(B,1);
%TODO: validate the input (may be not necessary with a good f_load)
parfor k=1:size(B,3)
    B(:,:,k)=imtransform(B(:,:,k),t_transform,'XData',[1 t_width],'YData',[1 t_height]);
end
parfor k=1:size(D,3)
    D(:,:,k)=imtransform(D(:,:,k),t_transform,'XData',[1 t_width],'YData',[1 t_height]);
end
%TODO: show an overlayed image with the registration results
%% Find noise and offset
[t_noise_sigma,t_hwbg]=f_estimate_noise_and_offset(A);

%% TODO: FRET (taketh from Vova)

%% Normalize !!!remove this when FRET is implemented!!!
A=A-t_hwbg;
% A=f_normalize(A);
t_hwbg=0; %<- This is very ugly and should be fixed asap

%% Diff (Pasha's favourite red-green variant is implemented)
p_diff=1; %parameter
p_mid=1; %parameter
AD=filter(ones(1,p_mid),p_mid,f_diff(A,p_diff)); % A differentiated

%% MASK the BG (does not account for out-of-focus fluo)
p_n_sigma_cutoff=3; %how many sigmas to count as noise;
t_threshold=t_hwbg+p_n_sigma_cutoff*t_noise_sigma;
AM(A>t_threshold)=A(A>t_threshold)-t_hwbg;
AM(A<=t_threshold)=0;

%% Make the diff relative to total intensity in the pixel
AT=A(:,:,1:end-p_diff); % AT= A truncated
% ADR(AT>t_threshold)=AD(AT>t_threshold)./AT(AT>t_threshold);  % ADR = A diff relative
% ADR(AT<=t_threshold)=0;
% ADR=reshape(ADR,size(AT)); %indexing with logical expressions makes the output 1D, this line fixes this

%% Filter the relative diff (in x,y)
p_imfilter_averaging_kernel_size=[3 3];
t_filter=fspecial('average',p_imfilter_averaging_kernel_size);
% ADRF=imfilter(ADR,t_filter); %ADR filtered

% %% clean the Absolute DIFF from the Bg pixels
% AD(AT<t_threshold)=0;
% ADF=imfilter(AD,t_filter);

%% Threshold the diff
p_n_sigma_cutoff_roi=3; %how many noise SD's to count as significant
t_threshold_roi=p_n_sigma_cutoff_roi*t_noise_sigma;
ADF_plus=(ADF>t_threshold_roi); %this work ONLY provided the bleaching has been compensated for
ADF_minus=(ADF<-t_threshold_roi);% this too
%% Prepare the second derivative (from scratch)
AD2=f_diff(AD,p_diff);
AD2F=imfilter(AD2,t_filter);

%% Threshold the AD2F;
% Using the sigma estimated above (t_noise_sigma) divided by sqrt(N) of the
% averaging imfilter
t_sqrt_n=sqrt(p_imfilter_averaging_kernel_size(1)^2+p_imfilter_averaging_kernel_size(2)^2);
t_threshold_ad2f_roi=p_n_sigma_cutoff_roi*t_noise_sigma/t_sqrt_n;
AD2F_plus=(AD2F>t_threshold_ad2f_roi);



%% Find 3D ROI's
i_3drois_input=AD2F_plus;
t_rois=bwconncomp(i_3drois_input,6);
% %prefilter we don't need rois of less than 9 voxel
t_roi_area=cell2mat(struct2cell(regionprops(t_rois,'Area')));
% t_roi_index=find(t_roi_area>9);

%we don't need rois that span less then 3 frames
%3 comes from the p_mid and should be changed to this later
t_roi_bbox=cell2mat(struct2cell(regionprops(t_rois,'BoundingBox'))');

delete(p_roi_output_file_name);
t_fid=fopen(p_roi_output_file_name,'w+');
fprintf(t_fid,'x y t roin');
fclose(t_fid);
t_nrois=size(t_rois.PixelIdxList,2);
for k=1:t_nrois;
    if t_roi_area(k)>8 % we don't need rois of less than 27 voxel
        if t_roi_bbox(k,6)>3 %we don't need rois that span less then 3 frames
            if t_roi_area(k)<1000;
            t_ROI=t_rois.PixelIdxList{k};
                    [x,y,z] = ind2sub(size(i_3drois_input),t_ROI);
%             [Fs,Vs,C_slice]=ind2patch(t_ROI,ADF_minus,'v'); %Creating patch data for selection of low voxels
%             t_cmap=colormap(lines(t_nrois));
%             hs=patch('Faces',Fs,'Vertices',Vs,'FaceColor',t_cmap(k,:),'FaceAlpha',0.8);
            %         scatter3(I,J,K,1,repmat(t_cmap(k),size(K)));
            a_out=[x y z k.*ones(size(z))];
             save(p_roi_output_file_name,'a_out','-ascii','-append');
             disp 'found a roi'
            end
        end
    end
end
% daspect([1 1 1]);
% toc
