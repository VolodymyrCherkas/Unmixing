function [ varargout ] = f_apply_transform_to_sequence( i_transform, varargin)
%f_apply_transform_to_sequence applies a transform to multiple sequences of
%images 
%
%
%This function accepts a transform as the first argument and 3D
%arrays of images, time in the third dimetion, and applies the transform to
%each image in each array.
%SYNOPSIS:
%function [ varargout ] = f_apply_transform_to_sequence( i_transform, varargin)
%
%Date:20121002
%Version:0.0.0
%Version format: major.minor.revision //major is disruptive,minor
%introduces features,revision fixes bugs

%% init
 if nargin~=nargout+1
     error('f_apply_transform_to_sequence: number of output arguments must be equal to the number of input arguments+1');
 end
%% core:transform
for i=1:(nargin-1)
t_width=size(varargin{i},2);
t_height=size(varargin{i},1);
for k=1:size(varargin{i},3)
    varargout{i}(:,:,k)=imtransform(varargin{i}(:,:,k),i_transform,'XData',[1 t_width],'YData',[1 t_height]);
end
end