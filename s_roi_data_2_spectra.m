% s_roi_data_2_spectra 
%
% Process group data from 3 rois obtained from matVis for easy
% transformation into the reference spectrum structure
%
% Author Madan Hennadii madang@biph.kiev.ua
% Version: 0.0.0
% Date: 20121004

t_aggr=[roiData_set01_Dim3;roiData_set02_Dim3;roiData_set03_Dim3;roiData_set04_Dim3;roiData_set05_Dim3;roiData_set06_Dim3];
for i=1:3
    t_roi_data{i}=t_aggr(i:3:end,:);
end