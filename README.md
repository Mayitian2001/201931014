# 201931014
姓名：马亦恬 学号201931014 日期：2021-11-20

###  实验二: From Fastq data files to Read Count Matrix
---

姓名:李云泽 学号: 201831062 日期: 2020-12-14


### 1. 实验目的

检查测序完成后得到的fastq文件的质量情况，然后使用 STAR 将 reads 比对到参考基因组上，再利用 R 包 DESeq2 进行主成分分析，最后找到差异表达基因，并对差异表达基因进行分析。

### 2. 实验准备

#### 2.1 实验平台

华为云服务器

操作系统：
```
Linux ecs-kc1-large-2-linux-20201119172228 4.15.0-70-generic #79-Ubuntu SMP Tue Nov 12 10:36:10 UTC 2019 aarch64 aarch64 aarch64 GNU/Linux
```

#### 2.2 数据简述

下载自计算生物学课程群

#### 2.3 软件配置

（1）实验中用到的软件及安装方法

（a）STAR：
```
wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
tar -xzf 2.5.3a.tar.gz
cd STAR-2.5.3a
```
（b）fastqc：
```
nohup wget -c http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip &
unzip fastqc_v0.11.5.zip
chmod 755 fastqc
```
（c）R及DESeq2：
```  
apt-get install R

if(!requireNamespace("BiocManager",quietly=TRUE))
install.packages("BiocManager")
```
