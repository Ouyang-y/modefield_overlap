# modefield_overlap
对波导模场和光纤模场计算重叠积分
## 文件结构
``` matlab
modefield_overlap
│
│   modefield_overlap_bgData2Excel.m  % 批量处理文件夹至 modefield_overlap.xlsx
│   modefield_overlap_bgDataDebug.m  % 绘图Debug示例
│
├───resource
│       bgDataRead.m  % .bgData读取函数
│       functionSignatures.json  % 自定义代码建议和自动填充
│       gauss2fit.m  % 二维高斯拟合函数
│       modefield_overlap_bgData.m  % 重叠积分计算函数
│       OSI_rainbow.mat  % OSI_rainbow 
└───test  % 测试样例
    │   1.bgData
    │   fiber-.bgData
    │   fiber.bgData
    │   V500KHz-50mW-160μm-13mms.bgData
    │   光纤.bgData
    └───20241107_V  % 批量处理文件夹样例
            11.bgData
            13.bgData
            15.bgData
            17.bgData
            7.bgData
            9.bgData
            fiber-.bgData
            标尺1.bgData
            标尺2.bgData
```
## 更新记录
### V0.3
- 文件结构维护
- 删除2Excel的test文件减小仓库体积
- 遗弃处理图片方式
- 改进函数modefield_overlap_bgData.m
    + 在Debug>figure2中加入R^2输出
    + 更新注释
- 增加modefield_overlap_bgDataDebug.m
- 改进modefield_overlap_bgData2Exce.m内addpath
### V0.2
- 函数化modefield_overlap_bgData.m
- 基于modefield_overlap_bgData.m构建modefield_overlap_bgData2Excel.m
- 更新.m文件注释
### V0.1
- 实现图片重叠积分计算
- 实现.bgData文件重叠积分计算
