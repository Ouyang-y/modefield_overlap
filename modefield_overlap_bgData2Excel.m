%% 
% input:
%   - bgDataFolder  .bgData文件夹路径
%   - fiberPath     fiber对应的.bgData路径
%
% Return:
%   在<bgDataFolder>的父文件夹下建立'modefield_overlap.xlsx'
%
%   其1:6列分别为:
%   '文件名'、'WG拟合符合程度(R方)'、'重叠积分'、
%   '重叠积分(dB)'、'拟合重叠积分'、'拟合重叠积分(dB)'
%
%   重叠积分与拟合重叠积分的差距在于:实验or拟合结果进行重叠积分
%
fiberPath = './test/fiber-.bgData';
bgDataFolder = './test/20241107_V';

bgDatas = dir(fullfile(bgDataFolder, '*.bgData'));
num = length(bgDatas);
if ~num,error("There's no .bgData file in selected folder, please check <bgDataFolder>!");end

originalPath = pwd;
T = table('Size',[num,6], ...
    'VariableTypes',{'string','double','double','double','double','double'}, ...
    'VariableNames',{'Name','R^2','Int2','Int2(dB)','Int2_Sim','Int2_Sim(dB)'});
mfilePath = mfilename("fullpath");
addpath([mfilePath(1:end-length(mfilename)),'\resource'])
warning('off','curvefit:prepareFittingData:removingNaNAndInf');
tic
for temp = 1:num
    fprintf('%d/%d:',temp,num)
    WGPath = fullfile(bgDatas(temp).folder,bgDatas(temp).name);
    [Int2,Int2_Sim,gof] = modefield_overlap_bgData(WGPath,fiberPath);
    [~,name,~]=fileparts(WGPath);
    T(temp,:) = {name,gof.rsquare,Int2,-10*log10(Int2),Int2_Sim,-10*log10(Int2_Sim)};
    toc
end
warning('on','curvefit:prepareFittingData:removingNaNAndInf');

cd(bgDataFolder)
cd('..')
writetable(T, 'modefield_overlap.xlsx');
cd(originalPath)
