function [Image, header] = DotGrayConverter(FileName, Folder)
%**************************************************************************
%
% [Image, header] = DotGrayConverter(FileName, Folder)
%
%**************************************************************************
% Output:
%   Image = Square Image from the converted .gray file
%   header = Header of the .gray file
%
% Variables:
%   FileName = Name of the file to be converted, Default: panic
%   Folder = Directory where the file to be converted exsists, 
%            Default: see below
% 
%**************************************************************************
%
% This program was written by Cameron Rodriguez
% Last Modified 10/19/2010
%
%**************************************************************************

%% Program Settings

SaveAsBMP = 1; % 1 = Save the converted image/figure as a bitmap image, Otherwise don't save

%% Set Defaults

if ~exist('FileName', 'var')
    FileName = 'primate';
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
    DisplayResults = 0;
end

%% Load the .gray file

ImageFileID = fopen([Folder,FileName,'.gray']);
ImageSaveFile = [Folder,FileName,'.bmp'];

% Figure out the image size from the header file
    TotalHeaderLenght = 256; % this seems to be the same for all size images
    header1 = fread(ImageFileID, 25, 'uint8=>char');
    M = str2double(textscan(header1(13:15),'%3c'));
    N = str2double(textscan(header1(21:23),'%3c'));
% Read the rest of the header
    header2 = fread(ImageFileID, (TotalHeaderLenght-25), 'uint8=>char');
% Combine the header reads    
    header = [header1',header2'];

for i=1:N
    Image(i,:) = fread(ImageFileID, M, 'uint8=>float');
end

%% Save image as BMP

if SaveAsBMP == 1
    imwrite(uint8(Image),ImageSaveFile,'bmp')
end

%% Display the Results

if DisplayResults == 1
    disp(header)
    
    imagesc(Image)
    xlabel('N - Pixel #');
    ylabel('M - Pixel #');
    axis image;
    colormap(gray);
    
    Image = []; % suppress the output
end