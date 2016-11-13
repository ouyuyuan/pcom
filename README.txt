
里面已经包含了一个restr场，初始场，地形场，等，边界场
解压开，把sbcf.data也放到目录里，compil.sh是编译脚本
编译完有个pcom.out，然后直接mpirun 它就可以了

这个由于设置的网格数问题，暂时支持8cpu并行4*1，2*1，3*1，6*1，5*1，10*1，12*1,理论上也可以，我没试过mpirun -n 8 ./pcom.out 就可以了

生成的是sfcdiag.grd文件对应的描述文件是sfcdiag.ctl，可以用grads画画图输出的是月平均的U，V，W，T，S，SSH

topog.data 是地形场
tsobs.data 初始场
restr.data 经过100年spinup的初始场，直接用这个当初始场就不用spinup了


namelist 参数文件

&contrl
imt=362, 纬向格点数
jmt=141, 经向格点数
km=60, 垂直层数
nt=2, 示踪变量数（温度，盐度两个）
dlam=1.0, 纬向分辨率（度）
dphi=1.0, 经向分辨率
phis=-70.0, 起始纬度
runlen=600, 积分时长（月）
restrt=.true., 是否从restrt场开始积分
snbc=0, 盐度自然边界条件是否开启（1=是）
gm90=1, 沿等密度面混合是否开启
implicitvmix=1,  隐式垂直混合是否开启
asselin_b=1, 正压积分过滤是否开启
asselin_c=1, 斜压积分过滤是否开启
asselin_t=1, 温盐积分过滤是否开启
smth=0, 高纬度平滑是否开启
smth_start_nlat=65.0, 高纬度平滑起始纬度（北边界）
smth_start_slat=-65.0,高纬度平滑起始纬度（南边界）
unesco=1, 是否用48项密度计算方案
boussinesq=0, 是否开启boussinesq近似
seadp= 750.00,   2268.41,   3823.64,   5433.89,   7117.36,   8891.88,
10775.23,  12784.64,  14937.31,  17249.66,  19738.11,  22418.13,
25305.19,  28413.60,  31757.69,  35350.49,  39205.04,  43332.94,。。。。 垂直分层（cm）
&mpicontrl
ncpux=8, 纬向cpu数
ncpuy=1, 经向cpu数
 &tsteps
dtts=1.0, 温盐积分步长（小时）
dtuv=0.5, 斜压积分步长（小时）
dtsf=0.5, 正压积分步长（分钟）
&mixing
fam=0, 是否从文件读入am场
fah=0,是否从文件读入ah场
fkm=0,是否从文件读入km场
fkh=0,是否从文件读入kh场
上面设为0则下面有效，下面将这四个参数设为常数
am_c=1.0e7, 侧摩擦系数
ah_c=1.0e7, 侧向混合系数
km_c=40.0, 垂向摩擦系数
kh_c=0.5, 垂向混合系数
cdbot=2.6e-3, 海底拖曳系数（底摩擦）

Zhang Yu
2012

This version has been successfully tested in avatar server
of IAP.
