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
