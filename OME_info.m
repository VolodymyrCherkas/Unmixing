function OME=OME_info(filename)
% bugfix 13.09.2012

%% OME_info for Till L.A. files
% .frame - starting from 1
% .cycle - starting from 1
% .time  - starting from 0
% .emch  - starting from 0
% .exp   - exact number in s
% .wl    - exact number in nm
% .fret  - 1 is FRET mode with image splitting, 0 - common epifluorescent
% mode
% L.A. current version is 2.2.0.6
% script version 04.07.2012 by Volodymyr Cherkas
% fixed: works with old L.A. 1.6.0.6 files!!!
% 06.07.2012 fixed nEm and nEx were wrong
% 26.11.2012 fixed nC=1 to cause error
% 10.01.2013 added expconf for convenience
%%
XML=xml2struct(filename);
OME.filename=filename;
if isfield(XML.OME, 'CA_colon_CustomAttributes')
    OME.version=XML.OME.CA_colon_CustomAttributes.LAVersion.Attributes.Version;
else OME.version=XML.OME.SA_colon_StructuredAnnotations.SA_colon_XMLAnnotation.SA_colon_Value.LA.OME_colon_LAVersion.Attributes.Version;
end
if ~strcmp(OME.version, '1.6.0.6')
    OME.date=XML.OME.OME_colon_Image.OME_colon_AcquiredDate.Text;  % date and time (for log and reference only)
else OME.date=XML.OME.OME_colon_Image.OME_colon_CreationDate.Text;
end
OME.DimensionOrder=XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.DimensionOrder;  % 'XYCTZ' is default for L.A.
OME.fret=strcmp(XML.OME.OME_colon_Experiment.Attributes.Type,'FRET');
if ~strcmp(OME.version, '1.6.0.6')
    OME.datatype=XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.Type;              % 'uint16' is default for L.A.
else OME.datatype=XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.PixelType;        % version 1.6.0.6
end
OME.nX=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.SizeX);           % size x, pixels
OME.nY=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.SizeY);           % size y, pixels
OME.nC=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.SizeC);           % number of channels
OME.nEx=OME.nC/(OME.fret+1);                                                % number of emission channels
OME.nEm=OME.nC/OME.nEx;                                                     % number of excitation channels
OME.nT=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.SizeT);           % size t, frames
OME.nZ=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.SizeZ);           % size z, z-stack planes
OME.physicalX=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.PhysicalSizeX);    % pixelsize x, um
OME.physicalY=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.PhysicalSizeY);    % pixelsize x, um
OME.physicalZ=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.Attributes.PhysicalSizeZ);    % pixelsize z, um

for i=1:length(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_TiffData)
    OME.frame(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_TiffData{1,i}.Attributes.IFD)+1;
    OME.cycle(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_TiffData{1,i}.Attributes.FirstT)+1;
    OME.exch(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_TiffData{1,i}.Attributes.FirstC);
    if ~strcmp(OME.version, '1.6.0.6')
        OME.time(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_Plane{1,i}.Attributes.DeltaT);
        OME.exp(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_Plane{1,i}.Attributes.ExposureTime);
    else
        OME.time(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_Plane{1,i}.OME_colon_PlaneTiming.Attributes.DeltaT);
        OME.exp(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_Plane{1,i}.OME_colon_PlaneTiming.Attributes.ExposureTime);
    end
end
j=0:OME.nC:length(OME.exch)-1;
OME.wl=zeros(1,length(OME.exch));
if ~strcmp(OME.version, '1.6.0.6')
    if (OME.nC>1)
        for i=1:OME.nC
            OME.wlconf(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_Channel{1,i}.OME_colon_LightSourceSettings.Attributes.Wavelength);
            OME.wl(j+i)=OME.wlconf(i);
        end
    else
        OME.wlconf(i)=str2double(XML.OME.OME_colon_Image.OME_colon_Pixels.OME_colon_Channel.OME_colon_LightSourceSettings.Attributes.Wavelength);
        OME.wl(j+i)=OME.wlconf;
    end
else
    if (OME.nC>1)
        for i=1:OME.nC
            OME.wlconf(i)=str2double(XML.OME.OME_colon_Image.OME_colon_LogicalChannel{1,i}.OME_colon_LightSourceRef.Attributes.Wavelength);
            OME.wl(j+i)=OME.wlconf(i);
        end
    else
        OME.wlconf=str2double(XML.OME.OME_colon_Image.OME_colon_LogicalChannel.OME_colon_LightSourceRef.Attributes.Wavelength);
        OME.wl(j+1)=OME.wlconf;
    end
end
i=1:OME.nT;
for j=1:OME.nC
    OME.time_ch(i,j)=OME.time((i-1)*OME.nC+j);
    OME.expconf(j)=OME.exp(j);
end

