function DotGrayFileWriter(Image, Header, FileName, Folder)
%**************************************************************************
%
% DotGrayFileWriter(Image, Header, FileName, Folder, N)
%
%**************************************************************************
%
% Variables:
%   Image = Image to be written
%   FileName = Name of the file to be converted
%   Header = Header of the dot gray image to be written
%   Folder =  Directory where the file to be converted exsists, Default: see below
% 
%**************************************************************************
%
% This program was written by Cameron Rodriguez
% Last Modified 10/19/2010
%
%**************************************************************************
%% Set Defaults

%Figure out the size of the image
    N = size(Image,1);
    M = size(Image,2);
    %N = str2double(textscan(Header(13:15),'%3c'));
    %M = str2double(textscan(Header(21:23),'%3c'));

       
if ~exist('FileName', 'var')
    FileName = num2str(date);
    if IsOSX
        Folder = '/Users/Cameron/Documents/MATLAB/ImageProcessing/images/';
    else
        Folder = 'C:\Documents and Settings\Cameron\Desktop\Cameron Rodriguez\ImageProcessing\images\';
    end
    DisplayResults = 1;
else
    if ~exist('Folder', 'var')
        if IsOSX
            Folder = '/Users/Cameron/Documents/MATLAB/ImageProcessing/images/';
        else
            Folder = 'C:\Documents and Settings\Cameron\Desktop\Cameron Rodriguez\ImageProcessing\images\';
        end
    end
end

%% Open a file

    SaveFileID = fopen([Folder,FileName,'.gray'], 'w');

%% Write to the file

    fwrite(SaveFileID, Header', 'uchar'); % Write Header

    fwrite(SaveFileID, Image' , 'uchar'); % Write Image

%% Close the file
    fclose(SaveFileID) % That's all she wrote, Ha Ha Ha...
    