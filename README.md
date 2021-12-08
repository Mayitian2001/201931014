
姓名:马亦恬  学号：201931014  日期: 2021-12-8

###  实验二: From Fastq data files to Read Count Matrix


### 1. 实验目的

通过质控即检查测序完成后得到的fastq文件的质量，丰度情况；再使用 STAR 将 reads 比对到参考基因组上得到读长的总长度，平均长度，唯一比对长度，平均比对长度等；接着通过IGV可视化BAM文件观察不同基因表达量的多少；最后利用R包DESeq2 进行差异表达分析。

### 2. 实验准备

#### 2.1 实验平台

Mac OS终端

操作系统：

```
$uname -a
Linux JSvr01 3.10.0-1160.el7.x86_64 #1 SMP Mon Oct 19 16:18:59 UTC 2020 x86_64 x86_64 x86_64 GNU/Linux
```

#### 2.2 数据简述

来源于服务器内/work/cb-data

#### 2.3 软件配置

（1）STAR：

```
cd /work
find . -name STAR
export PATH=$PATH:/work/bin/STAR-2.5.4b/bin/Linux_x86_64
```

（2）fastqc：

```
cd /work
find . -name FastQC
export PATH=$PATH:/work/bin/FastQC
```

（3）R及DESeq2：在linux或者Rstudio完成

```  
apt-get install R
if(!requireNamespace("BiocManager",quietly=TRUE))
install.packages("BiocManager")
```

### 3. 实验内容

实验线路：1，准备工作环境——2，使用fastqc质控——3，使用STAR将读长比对到基因组上——4，使用IGV查看bam文件——5，将命令在脚本中运行——6，主成分分析——7，差异表达分析
`

#### 3.1 Prepare the work directory

创建目录准备工作环境

```
mkdir workdir/myt
cd workdir/myt
ls -al
```

结果如下：

```
(base) S201931014 21:19:44 ~/workdir/myt/test
$ls -al
总用量 500132
drwxrwxr-x. 5 S201931014 S201931014      4096 12月  7 21:13 .
drwxrwxr-x. 7 S201931014 S201931014      4096 12月  7 21:18 ..
-rw-rw-r--. 1 S201931014 S201931014    250600 12月  7 19:21 ERR458493_fastqc.html
-rw-rw-r--. 1 S201931014 S201931014    324047 12月  7 19:21 ERR458493_fastqc.zip
-rw-r--r--. 1 S201931014 S201931014  59532325 12月  7 15:56 ERR458493.fastq.gz
-rw-r--r--. 1 S201931014 S201931014  58566854 12月  7 15:56 ERR458494.fastq.gz
-rw-r--r--. 1 S201931014 S201931014  58114810 12月  7 15:56 ERR458495.fastq.gz
-rw-r--r--. 1 S201931014 S201931014 102201086 12月  7 15:56 ERR458500.fastq.gz
-rw-r--r--. 1 S201931014 S201931014 101222099 12月  7 15:56 ERR458501.fastq.gz
-rw-r--r--. 1 S201931014 S201931014 100585843 12月  7 15:56 ERR458502.fastq.gz
drwxr-xr-x. 8 S201931014 S201931014      4096 10月  4 2018 FastQC
-rw-rw-r--. 1 S201931014 S201931014  10255059 1月  16 2020 fastqc_v0.11.8.zip
drwxrwxr-x. 2 S201931014 S201931014         6 12月  7 21:13 genome
-rw-rw-r--. 1 S201931014 S201931014     13430 12月  7 21:13 Log.out
-rw-------. 1 S201931014 S201931014     16263 12月  7 16:01 nohup.out
-rw-r--r--. 1 S201931014 S201931014  12360704 12月  7 15:56 R64.fa
-rw-r--r--. 1 S201931014 S201931014   8639033 12月  7 15:56 R64.gtf
-rw-r--r--. 1 S201931014 S201931014      1399 12月  7 15:56 runSTAR.sh
-rw-r--r--. 1 S201931014 S201931014        58 12月  7 15:56 samples.txt
drwx------. 2 S201931014 S201931014         6 12月  7 21:13 _STARtmp
```

#### 3.2 Examine the quality of the fastq data files

（1）使用fastqc进行质控

```
fastqc ERR458493.fastq.gz
```

生成了一个html文件ERR458493.fastq.gz,使用FileZilla软件下载并查看如下：

总序列长1093957，GC含量43%

![截屏2021-12-08 上午10.47.55](/Users/mayitian/Desktop/截屏2021-12-08 上午10.47.55.png)

![截屏2021-12-08 上午10.48.00](/Users/mayitian/Desktop/截屏2021-12-08 上午10.48.00.png)

#### 3.3 Run read mapping software

（1）解压缩文件检查数据

```
gunzip -c ERR458493.fastq.gz | head
gunzip -c ERR458493.fastq.gz | wc -l
less R64.fa
less R64.gtf
```

结果如下：

```
(base) S201931014 11:00:07 ~/workdir/myt/test
$gunzip -c ERR458493.fastq.gz | head
@ERR458493.1 DHKW5DQ1:219:D0PT7ACXX:1:1101:1724:2080/1
CGCAAGACAAGGCCCAAACGAGAGATTGAGCCCAATCGGCAGTGTAGTGAA
+
B@@FFFFFHHHGHJJJJJJIJJGIGIIIGI9DGGIIIEIGIIFHHGGHJIB
@ERR458493.2 DHKW5DQ1:219:D0PT7ACXX:1:1101:2179:2231/1
ACTAATCATCAACAAAACAATGCAATTCAAGACCATCGTCGCTGCCTTCGC
+
B@=DDFFFHHHHHJJJJIJJJJJJIJJJJJJJJJJJJJJJJJJJJIJJJJI
@ERR458493.3 DHKW5DQ1:219:D0PT7ACXX:1:1101:2428:2116/1
CTCAAAACGCCTACTTGAAGGCTTCTGGTGCTTTCACCGGTGAAAACTCCG

4375828

>I dna:chromosome chromosome:R64-1-1:I:1:230218:1 REF
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACA
CATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTT
ACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC
CACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCTCCAATTACCCATATC
CAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCACCATAC
TGTTCTTCTACCCACCATATTGAAACGCTAACAAATGATCGTAAATAACACACACGTGCT
TACCCTACCACTTTATACCACCACCACATGCCATACTCACCCTCACTTGTATACTGATTT
TACGTACGCACACGGATGCTACAGTATATACCATCTCAAACTTACCCTACTCTCAGATTC
CACTTCACTCCATGGCCCATCTCTCACTGAATCAGTACCAAATGCACTCACATCATTATG
CACGGCACTTGCCTCAGCGGTCTATACCCTGTGCCATTTACCCATAACGCCCATCATTAT
CCACATTTTGATATCTATATCTCATTCGGCGGTCCCAAATATTGTATAACTGCCCTTAAT
ACATACGTTATACCACTTTTGCACCATATACTTACCACTCCATTTATATACACTTATGTC
AATATTACAGAAAAATCCCCACAAAAATCACCTAAACATAAAAATATTCTACTTTTCAAC
AATAATACATAAACATATTGGCTTGTGGTAGCAACACTATCATGGTATCACTAACGTAAA
AGTTCCTCAATATTGCAATTTGCTTGAACGGATGCTATTTCAGAATATTTCGTACTTACA
CAGGCCATACATTAGAATAATATGTCACATCACTGTCGTAACACTCTTTATTCACCGAGC
AATAATACGGTAGTGGCTCAAACTCATGCGGGTGCTATGATACAATTATATCTTATTTCC
ATTCCCATATGCTAACCGCAATATCCTAAAAGCATAACTGATGCATCTTTAATCTTGTAT
GTGACACTACTCATACGAAGGGACTATATCTAGTCAAGACGATACTGTGATAGGTACGTT
ATTTAATAGGATCTATAACGAAATGTCAAATAATTTTACGGTAATATAACTTATCAGCGG
CGTATACTAAAACGGACGTTACGATATTGTCTCACTTCATCTTACCACCCTCTATCTTAT
TGCTGATAGAACACTAACCCCTCAGCTTTATTTCTAGTTACAGTTACACAAAAAACTATG
CCAACCCAGAAATCTTGATATTTTACGTGTCAAAAAATGAGGGTCTCTAAATGAGAGTTT
GGTACCATGACTTGTAACTCGCACTGCCCTGATCTGCAATCTTGTTCTTAGAAGTGACGC
ATATTCTATACGGCCCGACGCGACGCGCCAAAAAATGAAAAACGAAGCAGCGACTCATTT
TTATTTAAGGACAAAGGTTGCGAAGCCGCACATTTCCAATTTCATTGTTGTTTATTGGAC
ATACACTGTTAGCTTTATTACCGTCCACGTTTTTTCTACAATAGTGTAGAAGTTTCTTTC
TTATGTTCATCGTATTCATAAAATGCTTCACGAACACCGTCATTGATCAAATAGGTCTAT
AATATTAATATACATTTATATAATCTACGGTATTTATATCATCAAAAAAAAGTAGTTTTT
TTATTTTATTTTGTTCGTTAATTTTCAATTTCTATGGAAACCCGTTCGTAAAATTGGCGT

#!genome-build R64-1-1
#!genome-version R64-1-1
#!genome-date 2011-09
#!genome-build-accession GCA_000146045.2
#!genebuild-last-updated 2018-10
IV      sgd     gene    1802    2953    .       +       .       gene_id "YDL248W"; gene_source "sgd"; gene_biotype "protein_coding";
IV      sgd     transcript      1802    2953    .       +       .       gene_id "YDL248W"; transcript_id "YDL248W_mRNA"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_source "sgd"; transcript_biotype "protein_coding";
IV      sgd     exon    1802    2953    .       +       .       gene_id "YDL248W"; transcript_id "YDL248W_mRNA"; exon_number "1"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_source "sgd"; transcript_biotype "protein_coding"; exon_id "YDL248W_mRNA-E1";
IV      sgd     CDS     1802    2950    .       +       0       gene_id "YDL248W"; transcript_id "YDL248W_mRNA"; exon_number "1"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_source "sgd"; transcript_biotype "protein_coding"; protein_id "YDL248W";
IV      sgd     start_codon     1802    1804    .       +       0       gene_id "YDL248W"; transcript_id "YDL248W_mRNA"; exon_number "1"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_source "sgd"; transcript_biotype "protein_coding";
IV      sgd     stop_codon      2951    2953    .       +       0       gene_id "YDL248W"; transcript_id "YDL248W_mRNA"; exon_number "1"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_source "sgd"; transcript_biotype "protein_coding";
IV      sgd     gene    3762    3836    .       +       .       gene_id "YDL247W-A"; gene_source "sgd"; gene_biotype "protein_coding";
IV      sgd     transcript      3762    3836    .       +       .       gene_id "YDL247W-A"; transcript_id "YDL247W-A_mRNA"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_source "sgd"; transcript_biotype "protein_coding";

```

（2）使用STAR进行比对并查看结果

```
mkdir genome
STAR --runMode genomeGenerate --runThreadN 2 --genomeDir genome \ 
--genomeFastaFiles R64.fa --sjdbGTFfile R64.gtf \ 
--sjdbOverhang 50 

STAR --quantMode GeneCounts --genomeDir genome --runThreadN 2 \
--readFilesIn ERR458493.fastq.gz --readFilesCommand zcat \
--outFileNamePrefix wt1_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 \
--outSAMtype BAM SortedByCoordinate

less wt1_Log.final.out
less wt1_ReadsPerGene.out.tab
```

 查看结果如下:

<img src="/Users/mayitian/Desktop/截屏2021-12-08 上午11.20.37.png" alt="截屏2021-12-08 上午11.20.37"  />
输入读长1093957，平均长度为51，唯一比对读长占85.69%，平均比对长度为50.73

<img src="/Users/mayitian/Desktop/截屏2021-12-08 上午11.20.47.png" alt="截屏2021-12-08 上午11.20.47"  />



#### 3.4 Visualize the BAM file with IGV

使用IGV查看bam文件

```
samtools index wt1_Aligned.sortedByCoord.out.bam
```

每一个".bam"生成相对应的".bai"文件，使用filezilla 下载至自己的电脑。

IGV查看结果如下：
![截屏2021-12-08 下午12.16.13](/Users/mayitian/Desktop/截屏2021-12-08 下午12.16.13.png)



#### 3.5 Run the commands in a shell script

在脚本中运行文件，查看运行情况，并将结果合并到gene_count.txt中

```
sh runSTAR.sh >& log &

paste wt1_ReadsPerGene.out.tab wt2_ReadsPerGene.out.tab wt3_ReadsPerGene.out.tab
mu1_ReadsPerGene.out.tab mu2_ReadsPerGene.out.tab mu3_ReadsPerGene.out.tab | \
cut -f1,2,6,10,14,18,22 | \  
tail -n +5 > gene_count.txt
```

<img src="/Users/mayitian/Desktop/截屏2021-12-08 下午12.57.03.png" alt="截屏2021-12-08 下午12.57.03" style="zoom:50%;" />

将得到的gene_count.txt下载,并在excel中打开如下：

<img src="/Users/mayitian/Desktop/截屏2021-12-08 下午1.38.42.png" alt="截屏2021-12-08 下午1.38.42" style="zoom: 33%;" />

#### 3.6 Load the matrix into R and make PCA Plot with DESeq2

在R中运行脚本进行差异表达分析

（1）运行R脚本

```
library(DESeq2)
library(ggplot2)
cts <- as.matrix(read.csv("gene_count.txt", sep="\t", row.names=1, header=FALSE))
coldata <- read.csv("E:/samples.txt", sep="\t", row.names=1)
colnames(cts) <- rownames(coldata)
dds <- DESeqDataSetFromMatrix(countData = cts,colData = coldata,design = ~ Genotype)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Genotype)) +
  geom_point(size=3) +
  xlim(-2.5, 2.5) +
  ylim(-1, 1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)
ggsave("myplot.png")
```

将得到的plot图保存并查看

<img src="/Users/mayitian/Desktop/Rplot01.png" alt="Rplot01"  />
由图可见，wt1、2、3差异较大，而mu1、2、3相似度很高




### 4. 实验总结

#### 4.1 实验结论

（1）通过使用FastQC进行质控操作得到ERR458493.fastq.gz的总序列长1093957，GC含量43% ；

（2）通过使用STAR比对序列操作控制线程数量，读段长度得到输⼊读长1093957，平均长度为51，唯⼀⽐对读长占85.69%，平均⽐对长度为50.73等数据基本信息；

（3）使用IGV将bam文件可视化可知数据比对情况，每一条read的比对详情；

（4）通过将命令全部写入脚本，并运行脚本的方式可以更快地运行完成比对；

（5）通过R包dESeq2将比对得到的基因gene_count.txt进行分析得到ggplot图，可以看出wt1、2、3差异较大，而mu1、2、3相似度很高



#### 4.2 实验收获

​        这次实验是我第一次利用Linux进行相对全面完整的生物信息分析实验，从下载安装软件，获取数据以及进一步质控，比对基因，运用R包分析，IGV可视化等，我对各个软件和实验过程有了更加深入的了解，也学会了基本的操作步骤，对如何比对和分析基因序列有了更多的了解和感悟。实验的过程中遇到了各种问题和困难，老师也不断鼓励我们钻研和探索，同时不要粗心大意。会努力在接下来的学习中进一步了解熟悉生物信息学软件，努力钻研。
