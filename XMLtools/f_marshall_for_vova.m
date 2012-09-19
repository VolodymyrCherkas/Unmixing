function [] = f_marshall_for_vova( i_folder, i_format)
%f_marshall_for_vova load,register and save the LA output in vova readable
%format
%   As of 30.05.2012 vova readable format means a separate mat file for
%   each channel and accompanying metadata struct.
%
%   naming conventions:
%       underscore separates words
%       most variables (except for indexes) start with a prefix
%       prefixes:
%           i_ - input variable
%           t_ -  a temporary variable
%           p_ - a parameter (constant)
%           o_ - output variable
%% TODO: replace iformat with varargin and assume '3412' if not specified
%% TODO: capability to read a file with matching LA filenameas an formats (useful if I dont need all the files, or have different protocols)
%% pick a file
%TODO: check if i_folder is a valid folder name
%TODO: append a / to the i_folder if necessary
t_dir=dir(i_folder);
t_files={};
for i=1:size(t_dir,1)
    if strfind(t_dir(i).name,'LA_ProtocolCombination');
        if strfind(t_dir(i).name,'.tif');
            t_files{end+1}=t_dir(i).name;
        end
    end
end

for i=1:size(t_files,2)
    %% unmarshall
    % unmarshall is an old word for load
    t_fname=strcat(i_folder,t_files{i});
    % f_load loads the channels in the right order (specified by i_format)
    % A means ch00
    % B means ch01
    % C means ch10
    % D means ch11
    % now only i_format='3412' is implemented 30.05.2012 Gena Madan
    [A,B,C,D]=f_load(t_fname,i_format);
    %% register
    % find the transform that when applied to the first frame of B will
    % give the most correlation with the first frame of A
    % needs to be done once per DualView adjustment. (If you are careful
    % enough this means once per experiment, this is what I assume here).
    % Gena Madan
    
    t_transform=autoXY_shift(B(:,:,1),A(:,:,1),'nonreflective similarity');
    %     save('t_transform.mat','t_transform');
    % now let's apply this transform to every frame in B and D
    t_width=size(B,2);
    t_height=size(B,1);
    parfor k=1:size(B,3)
        B(:,:,k)=imtransform(B(:,:,k),t_transform,'XData',[1 t_width],'YData',[1 t_height]);
    end
    parfor k=1:size(D,3)
        D(:,:,k)=imtransform(D(:,:,k),t_transform,'XData',[1 t_width],'YData',[1 t_height]);
    end
    %TODO: show an overlayed image with the registration results
    %% marshall
    % marshall is an old word for save
    % write video
    p_A_suffix='_A';
    p_B_suffix='_B';
    p_C_suffix='_C';
    p_D_suffix='_D';
    p_out_folder='Vova/'
    mkdir(i_folder,p_out_folder);
    t_write_fname=strcat(i_folder,p_out_folder,t_files{i}(1:end-4))
    save(strcat(t_write_fname,p_A_suffix,'.mat'),'A');
    save(strcat(t_write_fname,p_B_suffix,'.mat'),'B');
    save(strcat(t_write_fname,p_C_suffix,'.mat'),'C');
    save(strcat(t_write_fname,p_D_suffix,'.mat'),'D');
%     %% write metadata (an ad-hoc solution)
%     t_struct.times=(1:size(A,3));
%     t_struct.A.excitation=435;
%     t_struct.B.excitation=435;
%     t_struct.C.excitation=505;
%     t_struct.D.excitation=505;
    
end
