差异基因物种来源网络热图（分组差异检验热图+网络图）

为明确差异基因的物种归属及其在样本组间的表达变异模式，本研究首先提取差异基因所在的重叠群（contig），并采用 Kraken2 联合 GTDB 数据库（或与 NCBI 数据库比对）进行高精度物种注释，以有效规避直接使用 NR/NT 数据库可能引入的比对偏差与物种误判风险；在此基础上，通过整合分组差异检验热图与物种来源网络图，系统评估各基因亚型在不同处理组中的表达丰度差异及其宿主分布特征。该分析旨在将“功能差异”与“物种背景”深度联结，揭示差异基因的物种来源与组间变异之间的耦合关系。这不仅有助于精准解析不同微生物类群在环境或处理响应中的功能贡献，也为后续关键菌株筛选与功能基因溯源提供了坚实的数据支撑，从而进一步深化了对微生物群落功能调控机制的理解。

测试数据及脚本所在连接：
https://github.com/whbeifan/BFRplot/tree/main/heatmap_network

输入数据格式说明：
（1）基因的丰度表（diff_gene2tax.tsv）:第一行为样本的名称，第一列为差异基因的ID，最后一列为基因的名称。
（2）基因对应的物种表（gene_tax.tsv）:第一列为基因的ID，第二列为基因对应的物种名称。
（3）样本的分组文件（group.list）:第一列为样本，第二列为样本对应的组。
运行命令:
Rscript plot_heatmap_network.R diff_gene2tax.tsv --group group.list --dtype gene_tax.tsv --dmethod zscale --prefix diff_gene2tax --width 230 --height 220 --network_width 2.5 --legabel "Zscale(TPM)"

----------------------
usage: plot_heatmap_network.R gene_diff_TPM.tsv --group group.txt --dtype gene_type.tsv

Function:Plot a network heatmap

positional arguments:
  input                 Input gene expression file.	#输入基因表达矩阵或者含量矩阵

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version information.
  -g GROUP, --group GROUP	#输入样本分组进行
                        Input sample group file(group.txt).
  --dtype DTYPE         Input gene source or classification	#输入基因分类信息
                        file(gene_type.tsv).
  --dmethod {zscale,log2,log10,NA}	#设置对数据进行处理的方法
                        Set the processing method for data, default=zscale
  -p PREFIX, --prefix PREFIX	#设置输出结果前缀
                        Set the prefix of the figure, default=out.
  --width WIDTH         Set the figure width, default=200	#设置图片宽度
  --height HEIGHT       Set the figure height,default=100	#设置图片高度
  --font FONT           Set the text fontsize.	#设置文字大小
  -nw NETWORK_WIDTH, --network_width NETWORK_WIDTH	#设置网络图部分的宽度
                        Set network width,default=1.5
  -m {wilcox.test,t.test,anova,kruskal.test}, --method {wilcox.test,t.test,anova,kruskal.test}	#设置进行差异统计的方法（多组比较只有使用anova和kruskal.test）
                        Method for setting up statistical testing of
                        differences, default=t.test
  --padj {none,bonferroni,fdr,holm,sidak}	#设置p值的矫正方法
                        Set P-value correction method, default=none
  -lb LEGABEL, --legabel LEGABEL	#设置热图的图例名称
                        Set legend labels, default=TPM
