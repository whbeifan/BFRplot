#!/usr/bin/env Rscript

library("vegan")
library("dplyr")
library("ggcor")
library("ggsci") #Science、Nature等高分期刊配色包
library("ggplot2")
library("argparse") #传参
library("extrafont") #字体

options(bitmapType="cairo") #关闭服务器与界面的互动响应


add_help_args <- function(){

    desc <- "Version: v1.2.1\nAuthor:whbeifan\nEmail: whbeifan@foxmail.com\n"
    hinfo <- "plot_mantel:Run mantel analysis
    URL:https://github.com/whbeifan/BFRplot
    Example:plot_mantel.R species_and_gene.tsv --data environmental_factor.tsv --group classified_index.tsv  --prefix out
        plot_mantel.R species_and_gene.tsv --data environmental_factor.tsv --prefix out
Input species_and_gene.tsv file format:
Sample\tA1\tA2\tA3
ASV_29353\t0\t0\t3529
ASV_15903\t6583\t7238\t3529
ASV_36906\t350\t349\t163
ASV_61291\t0\t0\t3247

Input environmental_factor.tsv file format:
Treat1\tpH\tEh\tEC\tSoil Cr(VI)
A1\t4.76\t180.2\t56.9\t0.68973
A2\t4.46\t181.2\t50.9\t0.68973
A3\t4.36\t189.2\t51.9\t0.68973\n"

 
    parser <- ArgumentParser(description=desc, usage=hinfo)
    parser$add_argument("input", type="character",
        help="Input species abundance file(phylum_abundance.tsv).")
    parser$add_argument("-d", "--data", type="character", required=TRUE,
        help="Input environmental factors table(environmental_factor.tsv).")
    parser$add_argument("-g", "--group", type="character", default="",
        help="Set the group file for input data, default=")
    parser$add_argument("-p", "--prefix", type="character", default="out",
        help="Set input result prefix, default=out")
    parser$add_argument("-w", "--width", type="integer", default=200,
        help="Set the width of the image(mm), default=200")
    parser$add_argument("--height", type="integer", default=160,
        help="Set the height of the image(mm), default=160")
    parser$add_argument("-m", "--method", type="character", default="spearman",
        choices=c("spearman", "pearson"),
        help="Set the method for calculating correlation, default=spearman")
    parser$add_argument("-dt", "--dtype", type="character", default="Microbe",
        help="Set the type of input data, default=Microbe")
    parser$add_argument("--mantel_method", type="character", default="mantel",
	choices=c("mantel", "mantel.randtest", "mantel.rtest", "mantel.partial"),
	help="Set inspection method, default=mantel")
    parser$add_argument("--spec_method", type="character", default="bray",
	choices=c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", 
                  "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", 
                  "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison",
                  "robust.aitchison"),
	help="Set the microbe distance algorithm, default=bray")
    parser$add_argument("--env_method", type="character", default="euclidean",
        choices=c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard",
                  "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial",
                  "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison",
                  "robust.aitchison"),
        help="Set the microbe distance algorithm, default=bray")
    parser$add_argument("--style", type="character", default="upper",
	choices=c("upper", "lower"),
	help="Draw heat map style, upper or lower triangle, default=upper")
    parser$add_argument("--shape", type="character", default="square",
	choices=c("square", "circle", "color", "ellipse", "pie", "star"),
        help="Set the fill shape, default=square")
    parser$add_argument("--high", type="character", default="#DC1623",
	help="Set the high color, default=#DC1623")
    parser$add_argument("--low", type="character", default="#2D6DB1", 
        help="Set the high color, default=#2D6DB1")
    parser$add_argument("--color1", type="character", default="#E64B357F",
        help="Set p value less than 0.01 color, default=#E64B357F")
    parser$add_argument("--color2", type="character", default="#4DBBD57F",
        help="Set p value between 0.01 and 0.05 color, default=#4DBBD57F")
    parser$add_argument("--color3", type="character", default="#00A0877F",
        help="Set p value greater than 0.05, default=#00A0877F")
    parser$add_argument("--typeface", type="character", default="Times New Roman",
        choices=c("Arial", "Times New Roman", "Calibri", "Garamond", "Georgia", "宋体", "仿宋", "黑体"),
        help="Set the typeface, default=Times New Roman")


    args <- parser$parse_args()
    return(args)

}


get_color <- function(color_num, alpha=0.5){

    #color_num为颜色数目
    if(color_num <= 2){ #最多可以获取2种颜色
        colors <- pal_igv("alternating", alpha=alpha)(color_num)
    }else if(color_num <= 10){ #最多可以获取10种颜色
        colors <- pal_npg("nrc", alpha=alpha)(color_num)
    }else if(color_num <= 12){
        colors <- pal_rickandmorty("schwifty", alpha=alpha)(color_num)
    }else if(color_num <= 16){
        colors <- pal_simpsons("springfield", alpha=alpha)(color_num)
    }else if(color_num <= 20){
        colors <- pal_d3("category20", alpha=alpha)(color_num)
    }else if(color_num <= 26){
        colors <- pal_ucscgb("default", alpha=alpha)(color_num)
    }else{ ##最多可以获取51种颜色
         colors <- pal_igv("default")(color_num)
    }

    return(colors)
}


get_group_index <-  function(groups){

    #获取分组的索引
    groupid <- unique(groups$group)
    r <- list()
    for (i in groupid){
        r[[i]] <- which(groups$group==i)
    }

    return(r)
}


plot_mantel <- function(data1, data2, group="", prefix="out", width=200, height=160,
                        method="spearman", dtype="Microbe", mantel_method="mantel", spec_method="bray",
                        env_method="euclidean", style="upper", shape="square", high="#DC1623",
			low="#2D6DB1", color1="#E64B357F", color2="#4DBBD57F", color3="#00A0877F",
                        typeface="Arial"){

    #中文转为英文格式字体
    if (typeface =="宋体"){
        typeface <- "SimSun"
    }else if (typeface =="仿宋"){
        typeface <- "FangSong"
    }else if (typeface =="黑体"){
        typeface <- "SimHei"
    }

    #读取文件格式
    data1 <- read.delim(data1, row.names=1, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, sep="\t")
    data2 <- read.delim(data2, row.names=1, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

    #data2 <- t(data2) #对环境因子进行转置
    data1 <- data1[which(rowSums(data1) > 0),] #过滤所有丰度都为0得行

    sample_id1 <- colnames(data1)
    sample_id2 <- colnames(data2)
    sample_id <- intersect(sample_id1, sample_id2)
    print("分析使用共有的样本信息如下：")
    print(sample_id)

    data1 <- data1[, sample_id] #匹配样本
    data2 <- data2[, sample_id] #匹配样本
    
    if(group!=""){
        group <- read.delim(group, header=FALSE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
        colnames(group) <- c("microbe", "group")
        rownames(group) <- group$microbe
        microbe_id1 <- group$microbe #对数据进行排序，防止样本缺失或者对应关系比对
        microbe_id2 <- rownames(data1)
        microbe_id <- intersect(microbe_id1, microbe_id2)
        group <- group[microbe_id,]
        data1 <- data1[microbe_id,]
    }else{
        microbe_id <- rownames(data1)
        group <- data.frame(sample=microbe_id, group=microbe_id)
        colnames(group) <- c("microbe", "group")
        group$group <- rep(c(dtype),each=length(group$group))
    }
    data2 <- t(data2)

    gindex <- get_group_index(group)
    df_mantel <- mantel_test(t(data1), data2,
        mantel.fun=mantel_method,
        spec.dist.method=spec_method, 
        env.dist.method=env_method, #na.rm=FALSE,
        spec.select=gindex   #输入分组的位置
    )
    
    df_mantel <- df_mantel %>% mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
        labels = c("<0.1", "0.1-0.2", "0.2-0.4", ">=0.4")),#定义Mantel的R值范围标签
        df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
        labels = c("<0.01", "0.01-0.05", ">= 0.05"))
    )
    #df_mantel$spec <- factor(df_mantel$spec, levels=unique(df_mantel$spec))

    write.table(df_mantel,  paste(prefix, ".stat_mantel.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE, na="")
  
    p <- quickcor(data2, method=method, type=style, cor.test=T,
            grid.size=0.5, #网格线条粗细，默认0.25
            grid.colour="#C8C8C8", #网格颜色，默认为"grey50"
            cluster.type="all") +#环境因子之间的相关性热图
        #geom_mark(r=NA, sig.thres=0.05, size=5, colour="black")+#显著性标签
        scale_fill_gradient2(high=high, mid="white", low=low) + #颜色设置
        anno_link(df_mantel, aes(color=df_p, size=df_r), curvature=0.15, label.size=4, label.fontface=1.5)+ #curvature= 0.2,#连接线变为曲线
        scale_size_manual(values=c(0.5, 1, 1.5, 2))+#连线粗细设置
        scale_color_manual(values=c("<0.01"=color1, "0.01-0.05"=color2, ">= 0.05"=color3))+#线条颜色设置
        guides(fill=guide_colorbar(title="Spearman's r", order=1),#图例相关设置
             size=guide_legend(title="Mantel's r", order=2),
             color=guide_legend(title="Mantel's p", order=3, keywidth=unit(0.45, "cm"), override.aes=list(size=4.5)), #设置颜色线的粗细(部分版本使用linewidth)
             linetype="none") + 
        theme(legend.key=element_blank()) #去除图例背景

    if (shape=="square"){
        p <- p + geom_square() +geom_mark(r=NA, sig.thres=0.05, size=5, colour="black") #相关性显示形式
    }else if (shape=="circle"){
        p <- p + geom_circle2()+geom_mark(r=NA, sig.thres=0.05, size=5, colour="black") #显示成圆型
    }else if (shape=="color"){
        p <- p + geom_color() +geom_mark(r=NA, sig.thres=0.05, size=5, colour="black") #全部填充
    }else if (shape=="ellipse"){
        p <- p + geom_ellipse2() +geom_mark(r=NA, sig.thres=0.05, size=5, colour="black")#填充为椭圆形
    }else if (shape=="pie"){
        p <- p + geom_pie2() +geom_mark(r=NA, sig.thres=0.05, size=5, colour="black")#显示成饼图的形式
    }else if (shape =="start"){
        p <- p + geom_star() +geom_mark(r=NA, sig.thres=0.05, size=5, colour="black")#显示成星星的形状
    }
    
    p <- p + theme(text=element_text(family=typeface))

    ggsave(p, filename=paste(prefix, ".mantel_test.png", sep=""), width=width, height=height, units="mm", dpi=350)
    ggsave(p, filename=paste(prefix, ".mantel_test.pdf", sep=""), width=width, height=height, units="mm", device=cairo_pdf, family=typeface)
    
}

args <- add_help_args()
plot_mantel(data1=args$input, data2=args$data, group=args$group, prefix=args$prefix, width=args$width,
            height=args$height, method=args$method, dtype=args$dtype, mantel_method=args$mantel_method,
            spec_method=args$spec_method, env_method=args$env_method, style=args$style, shape=args$shape,
            high=args$high, low=args$low, color1=args$color1, color2=args$color2, color3=args$color3, typeface=args$typeface)

