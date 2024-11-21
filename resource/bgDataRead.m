function [I,bgSet] = bgDataRead(filename,varargin)
% bgDataRead  Read data from .bgData file.
%   [I,bgData] = bgDataRead(filename)
%   [I,bgData] = bgDataRead(filename,"Llimit",Llimit)
%   [I,bgData] = bgDataRead(filename,"Llimit",Llimit,"maxpower",maxpower,"GroupName",GroupName)
%
%   Llimit = [0,1); BeamGage 颜色→z轴刻度下限
%   maxpower = 4000; BeamGage 颜色→z轴刻度峰值
%   GroupName = 1; BeamGage 读取帧序号
%
% Syntax: (这里添加函数的调用格式, `[]`的内容表示可选参数)
%	[I,bgSet] = bgDataRead(filename ...
%							[, 'Llimit', 0 ...
%							 , 'maxpower', 4000...
%							 , 'GroupName', 1]);
%
% Params:
%   - filename  [required]  [char] .bgData文件路径
%   - Llimit    [namevalue]  [numeric; >=0; <1] BeamGage 颜色→z轴刻度下限
%   - maxpower  [namevalue]  [numeric] BeamGage 颜色→z轴刻度峰值
%   - GroupName  [namevalue]  [numeric] BeamGage 读取帧序号
%
% Return:
%   - I 光强矩阵
%   - bgSet .bgSet.numcols .bgSet.numrows .pixelscalexum .pixelscaleyum
%
% Matlab Version: R2023b
%
% Author: oyy
%

maxpower=4000;Llimit=nan;GroupName=1;
if nargin>1
    if nargin>6||~bitget(nargin,1)
        error(sprintf(['bgDataReadError: parameter num error\n' ...
            'bgDataRead(filename)\n' ...
            'bgDataRead(filename,"Llimit",Llimit)\n' ...
            'bgDataRead(filename,"Llimit",Llimit,"maxpower",maxpower,"GroupName",GroupName)\n']))
    end
    for temp = 1:2:nargin-1
        switch varargin{temp}
            case 'maxpower'
                maxpower=varargin{temp+1};
            case 'Llimit'
                Llimit=varargin{temp+1};
            case 'GroupName'
                GroupName=varargin{temp+1};
            otherwise,error('bgDataReadError: unrecognized parameter');
        end
    end
end
h5fileNameBegin = ['/BG_DATA/',num2str(GroupName)];
numcols = h5read(filename,[h5fileNameBegin,'/RAWFRAME/WIDTH']);
numrows = h5read(filename,[h5fileNameBegin,'/RAWFRAME/HEIGHT']);
bgSet.pixelscalexum = double(h5read(filename,[h5fileNameBegin,'/RAWFRAME/PIXELSCALEXUM']));
bgSet.pixelscaleyum = double(h5read(filename,[h5fileNameBegin,'/RAWFRAME/PIXELSCALEYUM']));
power_calibration_multiplier = h5read(filename,[h5fileNameBegin,'/RAWFRAME/ENERGY/POWER_CALIBRATION_MULTIPLIER']);
data = h5read(filename,[h5fileNameBegin,'/DATA']);
encoding = h5read(filename,[h5fileNameBegin,'/RAWFRAME/BITENCODING']);
BitsPerPixel = 32;
if power_calibration_multiplier == 0
    switch lower(encoding)
        case {'l8','l16_8'} % 8 bit
            data = data ./ 2^(BitsPerPixel - 8 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case 'l16_10' % 10 bit
            data = data ./ 2^(BitsPerPixel - 10 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case 'l16_12' % 12 bit
            data = data ./ 2^(BitsPerPixel - 12 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case 'l16_14' % 14 bit
            data = data ./ 2^(BitsPerPixel - 14 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case {'l16_16','l16'} % 16 bit
            data = data ./ 2^(BitsPerPixel - 16 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case {'r8','r16_8'} % 8 bit
            data = data ./ 2^(BitsPerPixel - 8 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case 'r16_10' % 10 bit
            data = data ./ 2^(BitsPerPixel - 10 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case 'r16_12' % 12 bit
            data = data ./ 2^(BitsPerPixel - 12 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case 'r16_14' % 14 bit
            data = data ./ 2^(BitsPerPixel - 14 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case {'r16_16','r16'} % 16 bit
            data = data ./ 2^(BitsPerPixel - 16 - 1);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case 's16_14'  % signed 14 bit
            data = data ./ 2^(BitsPerPixel - 14);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case {'s16_16'} % signed 16 bit
            data = data ./ 2^(BitsPerPixel - 16);
            I = hdf5data_to_matrix(data,numcols,numrows);
        case {'s32'} % signed 32 bit
            I = hdf5data_to_matrix(data,numcols,numrows);
        otherwise
            error('Unknown bitencoding.')
    end
else % Use the Power Calibration Multiplier instead
    I = hdf5data_to_matrix(data,numcols,numrows);
    I = I.*power_calibration_multiplier;
end
I = I/maxpower;
[bgSet.numcols,bgSet.numrows] = deal(double(numcols),double(numrows));
if ~isnan(Llimit),I(I<Llimit)=0;end
end

function [I] = hdf5data_to_matrix(data,width,height)
%  Convert the 1-D array of data into a 2-D matrix, or image.
I = zeros(height,width);
index = 1;
for i = 1:height
    for j = 1:width
        I(i,j) = data(index);
        index = index + 1;
    end
end
end