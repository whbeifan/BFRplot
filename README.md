# BFRplot
R plotting script(R语言绘图脚本)

后续我们将构建了R的docker镜像，包含常用的R包和部分不常用的R包。需要的请点关注和留言。

## Mantel test（Mantel检验）
运行命令：
```
Rscript plot_mantel.R species_and_gene.tsv --data environmental_factor.tsv --group classified_index.tsv --prefix microbe --width 200 --height 180
Rscript plot_mantel.R species_and_gene.tsv --data environmental_factor.tsv --prefix microbe --width 200 --height 180
```
文件说明：
species_and_gene.tsv为每个样本的物种丰度、基因定量信息和分箱的定量信息等；classified_index.tsv为species_and_gene.tsv对应的分类索引文件，不输入该文件则将species_and_gene.tsv视为一个矩阵表进行分析；environmental_factor.tsv为每个样本的环境影子信息。

绘图结果：

<img width="500" height="500" src="https://github.com/whbeifan/BFRplot/blob/main/mantel_test/microbe.mantel_test.png" />

## Network heatmap of species origins for differential genes（差异基因物种来源网络热图）
运行命令：
```
Rscript plot_heatmap_network.R diff_gene2tax.tsv --group group.list --dtype gene_tax.tsv --dmethod zscale --prefix diff_gene2tax --width 230 --height 220 --network_width 2.5 --legabel "Zscale(TPM)"
```
文件说明：
diff_gene2tax.tsv为基因的丰度表，第一行为样本的名称，第一列为差异基因的ID，最后一列为基因的名称；gene_tax.tsv为基因对应的物种表，第一列为基因的ID，第二列为基因对应的物种名称；group.list为样本的分组文件，第一列为样本，第二列为样本对应的组。

绘图结果：

<img width="500" height="500" src="https://github.com/whbeifan/BFRplot/blob/main/heatmap_network/diff_gene2tax.heatmap_network.png" />
