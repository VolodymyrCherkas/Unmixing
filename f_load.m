function [ varargout ] = f_load( i_full_fname, i_format)
%f_load Load FRET data from files and output as channels
% 
% Version 0.1 (30.05.2012 by Gena Madan)  
%% MAYBE: Provide facilities to read only set of frames
%% TODO: if cannot open i_full_fname abort here
%% TODO: try to derive i_format if not specified
%% Init
    t_info=imfinfo(i_full_fname);
%%
switch i_format
    % i_format is a string that specifies order in which channels go
    % 1 for donor emission on donor excitation
    % 2 for acceptor emission on donor excitation
    % 3 for donor emission on acceptor excitation
    % 4 for acceptor emmision on acceptor excitation
    % use brackets when LA doesn't separate the emission channels
    case '(12)(34)'
        %This is the old format, without using the LA facility to slice
        %when FRET tab is selected.
        errordlg('f_load:(12)(34) format is not yet implemented' );
        
    case '1234'
        errordlg('f_load:1234 format is not yet implemented' );
    case '3412'
        t_end=size(t_info,1);
        for i=1:floor(t_end/4)
            varargout{1}(:,:,i)=double(imread(i_full_fname,4*(i-1)+3));
            varargout{2}(:,:,i)=double(imread(i_full_fname,4*(i-1)+4));
            varargout{3}(:,:,i)=double(imread(i_full_fname,4*(i-1)+1));
            varargout{4}(:,:,i)=double(imread(i_full_fname,4*(i-1)+2));
        end
        
        % this is the default new format
    case '123456'
        t_end=size(t_info,1);
        for i=1:floor(t_end/6)
            varargout{1}(:,:,i)=double(imread(i_full_fname,6*(i-1)+1));
            varargout{2}(:,:,i)=double(imread(i_full_fname,6*(i-1)+2));
            varargout{3}(:,:,i)=double(imread(i_full_fname,6*(i-1)+3));
            varargout{4}(:,:,i)=double(imread(i_full_fname,6*(i-1)+4));
            varargout{5}(:,:,i)=double(imread(i_full_fname,6*(i-1)+5));
            varargout{6}(:,:,i)=double(imread(i_full_fname,6*(i-1)+6));
        end
    otherwise
        errordlg('f_load:%s: unknown format string',i_format );
        
        
end

