=========================================================================

Network-Constraint data Integration cancer subtyping Classfication (NCIC)
一种网络约束的数据集成癌症亚型分类方法

=========================================================================

Version 1.0

=========================================================================
代码执行环境：
============

MATLAB version：R2017a
R：R version3.5

=========================================================================
需要加载的包：
=============

survival
flexclust
cluster
CancerSubtypes
Shiny

=========================================================================
系统需求：
========

运行此程序需要一个32GB RAM的计算机，在Windows系统下即可执行。

=========================================================================
安装：
====

1.启动MATLAB
2.通过修改/data目录下的数据文件和名字
3.运行MainNCIC.m
4.得到样本聚类标签

=========================================================================
示例：
====

demo_NCIC_TCGA_KIRC.m 	使用TCGA数据库肾透明细胞癌部分癌症病人的三种组学数据运行NCIC方法
demo_KIRC_get_best_K.R	使用CGA数据库肾透明细胞癌部分癌症病人的三种组学数据运行NCIC方法得到样本标签使用轮廓系数方法得到最佳聚类簇个数K,并且绘制图片
demo_KIRC_KM.R		使用CGA数据库肾透明细胞癌部分癌症病人的三种组学数据运行NCIC方法得到样本标签结合临床数据绘制KM生存曲线，并且基于Cox Log-rank test计算P-value	
demo_KIRC_Accuracy.R	使用CGA数据库肾透明细胞癌部分癌症病人的三种组学数据运行NCIC方法得到样本标签结合金标（如果有）计算NMI、RI和ARI，衡量准确
=======================================================================================

该软件版权所有  ? 2019
西安电子科技大学
可以通过以下方式获得对该软件商业使用的许可：
陕西省西安市雁塔区太白南路2号   邮编：710000
西安电子科技大学
lgao@mail.xidian.edu.cn
17737112575@163.com

=========================================================================
