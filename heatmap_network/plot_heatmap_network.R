#!/usr/bin/env Rscript

library("dplyr")
library("ggsci") #Science、Nature等高分期刊配色包
library("ggplot2")
library("cowplot")  
#library("magrittr") #加载管道
library("patchwork")
library("reshape2")
library("argparse") #传参
library("agricolae") 
library("RColorBrewer")

options(bitmapType="cairo") #关闭服务器与界面的互动响应

add_help_args <- function(){

    version <- "Version: v1.0.0\nAuthor:whbeifan\nEmail:whbeifan@foxmail.com"
    desc <-  "Function:Plot a network heatmap\n"
    usage <- "plot_heatmap_network.R gene_diff_TPM.tsv --group group.txt --dtype gene_type.tsv"

    parser <- ArgumentParser(description=desc, usage=usage)
    parser$add_argument("-v", "--version", action="version", version=version,
        help="Print version information.")
    parser$add_argument("input", type="character",
        help="Input gene expression file.")
    parser$add_argument("-g", "--group", type="character", default="",
        help="Input sample group file(group.txt).")
    parser$add_argument("--dtype", type="character",
        help="Input gene source or classification file(gene_type.tsv).")
    parser$add_argument("--dmethod", type="character", default="zscale",
        choices=c("zscale", "log2", "log10", "NA"),
        help="Set the processing method for data, default=zscale")
    parser$add_argument("-p", "--prefix", type="character", default="out",
        help="Set the prefix of the figure, default=out.")
    parser$add_argument("--width",type="integer", default=200,
        help="Set the figure width, default=200")
    parser$add_argument("--height",type="integer", default=100,
        help="Set the figure height,default=100")
    parser$add_argument("--font", type="integer", default=10,
        help="Set the text fontsize.")
    parser$add_argument("-nw", "--network_width", type="double", default=1.5,
        help="Set network width,default=1.5")
    parser$add_argument("-m", "--method", type="character", default="kruskal.test",
        choices=c("wilcox.test", "t.test", "anova", "kruskal.test"),
        help="Method for setting up statistical testing of differences, default=t.test")
    parser$add_argument("--padj", type="character", default="none",
        choices=c("none", "bonferroni", "fdr", "holm", "sidak"),
        help="Set P-value correction method, default=none")
    parser$add_argument("-lb", "--legabel", type="character", default="TPM",
        help="Set legend labels, default=TPM")

    args <- parser$parse_args()
    return(args)
}    


stat_mean <- function(data, group){

    #计算多组的均值
    gids <- unique(group$Group)
    r <- NULL

    for(i in gids){
        msid <- group$Sample[group$Group == i]  #获取目标组的id
        if(length(msid) <=1){
            mean_col <- data[, msid, drop=FALSE]  # 保持为数据框
        }else{
            temp <- rowMeans(data[, msid, drop=FALSE], na.rm=TRUE)
            mean_col <- matrix(temp, ncol=1, dimnames=list(rownames(data), i))
        }
        if(is.null(r)){
            r <- mean_col
        }else{
            r <- cbind(r, mean_col)
        }
    }

    return(r)
}


hclust_sort <- function(data){
    #对矩阵进行聚类排序
    set.seed(123)
    row_dist <- dist(data, method="euclidean")   #计算距离矩阵
    row_hc <- hclust(row_dist, method="ward.D2") #进行层次聚类
    row_order <- row_hc$order #获取聚类后的顺序
    data <- data[row_order,] #获得排序后的矩阵

    return(data) 
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


diff_stat <- function(data, group, method="anova", padj="none"){

    n <- 0
    for(i in colnames(data)){
        n <- n+1
        temp <- cbind(group, Value=data[, i])
        if(method=="auto"){
            r <- temp %>% group_by(Group) %>% summarise(shapiro_p=shapiro.test(Value)$p.value) #正态性检验,有组内值完全一样会报错
            if(min(r$shapiro_p)>=0.05){
                print(paste("所以组均符合正态分布，进行单因素方差分析", i, sep=":"))
                rtest <- aov(Value ~ Group, data=temp)
                out <- LSD.test(rtest, "Group", p.adj="bonferroni") #p.adj="none" 不进行P值矫正
            }else{
                print(paste("存在部分组不符合正态分布，进行Kruskal-Wallis检验", i, sep=":"))
                print(r)
                out <- kruskal(temp$Value, temp$Group, p.adj=padj)
            } 
        }else if(method=="ANOVA"){
            rtest <- aov(Value ~ Group, data=temp)
            out <- LSD.test(rtest, "Group", p.adj=padj)
        }else{
            out <- kruskal(temp$Value, temp$Group, p.adj=padj)
        }
        rclass <- out$groups[, "groups"] #或者标记的分类
        rclass<- data.frame(Group=rownames(out$groups), Value=rclass)
        colnames(rclass) <- c("Group", i)
        if(n==1){
            dc <- rclass
        }else{
            dc <- merge(dc, rclass)
        }

    }
    return(dc)   
}


plot_heatmap_network <- function(data, group="", dtype="", dmethod="zscale", prefix="out", width=200, height=100, font=10,
                                 network_width=1.5,  method="anova", padj="none", legabel="TPM"){

    data <- read.delim(data, sep="\t", row.names=1, head=TRUE, check.names=FALSE, quote="")
    sample_id <- colnames(data)
    #data <- data[order(data$Gene), ]

    dlabel <- data.frame(Contig=rownames(data), Gene=data$Gene)
    rownames(dlabel) <- dlabel$Contig
    gmaxlen <- max(nchar(as.character(dlabel$Gene)))

    if(group == ""){
        groups <- data.frame(Sample=sample_id, Group=sample_id)
        rownames(groups) <- groups$Sample
        if(length(groups$Sample) >=10){
            show_legend <- FALSE
        }
    }else{
        groups <- read.delim(group, header=FALSE, sep="\t", check.names=F, stringsAsFactors=F)
        colnames(groups) <- c("Sample", "Group")
        rownames(groups) <- groups$Sample
        sample_id <- intersect(sample_id, groups$Sample)
        groups <- groups[sample_id, ]  #对文件进行排序，防止组与样本数据不匹配
        data <- data[, sample_id] #匹配样本
        if(length(unique(groups$Group))==1){
            groups$Group <- groups$Sample
            print("样本只有一个分组，将单个样本视为一个组。")
        }
    }
    print(paste("绘图样本数目", length(sample_id), sep=":"))
    data <- data[rowSums(data, na.rm=TRUE) > 0, ] #移除行和 ≤ 0 的行
    dtype <- read.delim(dtype, header=FALSE, sep="\t", check.names=F, stringsAsFactors=F)
    colnames(dtype) <- c("gene", "type")
    rownames(dtype) <- dtype$gene
    maxlen <- max(nchar(as.character(dtype$type))) #获取最长的分类名称
    print(paste("最长的分类名称长度", maxlen, sep=":"))
    old_data <- data
    data <- stat_mean(data, groups) #计算样本均值
    if(dmethod=="zscale"){   
        data <- t(scale(t(data)))
    }else if(dmethod=="log2"){
        data <- data + 1e-3 #叫一个极小值
        data <- log2(data)
    }else if(dmethod=="log10"){
        data <- data + 1e-6 #叫一个极小值
        data <- log10(data)
    }
    
    data <- as.matrix(hclust_sort(data)) #对基因进行聚类排序
    #匹配顺序
    gene_id <- rownames(data)
    gene_id <- intersect(gene_id, dtype$gene)
    data <- data[gene_id, ]
    dtype <- dtype[gene_id, ]
    print(paste("绘图基因数目", length(gene_id), sep=":"))

    old_data <- old_data[gene_id, ]
    dc <- diff_stat(t(old_data), groups, method, padj)
    rownames(dc) <- dc$Group
    dc$Group <- NULL
    dc <-  melt(t(dc), varnames=c("Gene", "Group"), value.name="Class")
    data_melt <- melt(data, varnames=c("Gene", "Sample"), value.name="Expression")
    data_melt$Class <- dc$Class

    data_melt$Expression <- as.numeric(data_melt$Expression)
    midpoint <- floor(median(data_melt$Expression, na.rm=TRUE)) +1
    print(paste("热图中位数", midpoint, sep=":"))
  
    plot_layout(widths=c(0.25, 1))
    p <- ggplot(data_melt, aes(x=Sample, y=Gene)) +  # 使用ggplot绘图，设置映射
        geom_tile(aes(fill=Expression)) +  # 添加瓷砖图层，用于绘制热图
        geom_text(aes(label=Class), color="black", size=(font-3)/3) +
        geom_text(data=dlabel, aes(x=-Inf, y=Contig, label=Gene),
             hjust=1, vjust=0.5, size=(font-2)/3, color="black", inherit.aes=FALSE) +
        scale_y_discrete(position="right") +  # 设置y轴刻度位置
        scale_fill_gradient2(name=legabel, low="#6495ED", mid="#FFFFFF", high="#FF9999", na.value="grey50") + #midpoint=midpoint
        coord_cartesian(clip="off") +
        #scale_fill_gradientn(colors=colorRampPalette(c("#6495ED", "#FFFFFF", "#FF9999"))(100)) +
        theme(# 设置主题
            panel.background=element_blank(),  # 设置面板背景为空
            plot.background=element_blank(),  # 设置绘图背景为空
            legend.background=element_rect(fill=NA, color=NA),  # 设置图例背景为空
            #plot.margin=margin(20, 200, 20, 20),  # 设置绘图边距
            axis.title=element_blank(),  # 设置轴标题为空
            axis.ticks=element_blank(),  # 设置轴刻度为空
            axis.text.y=element_text(color="black", size=font-2),
            axis.text.x=element_text(angle=-45, color="black", vjust=0, hjust=0.2, size=font-2),  # 设置x轴文字样式
            #legend.position.inside=c(3, 0.8),  # 设置图例位置
            legend.title=element_text(size=font-1, color="black"),
            #legend.title=element_blank(),  # 设置图例标题为空
            plot.margin=margin(5, 5, 5, gmaxlen*4.6)
       )
    plot_build <- ggplot_build(p)
    coords <- plot_build$data[[1]][, c("x", "y")]

    gene_num <- length(dtype$gene)
    dtype$x1 <- rep(1, gene_num)
    dtype$y1 <- seq(1, gene_num)
    dtype$x2 <- rep(5, gene_num)
    group_num <- length(unique(dtype$type))
    gstep <- gene_num/(group_num+1)
    dtype$y2 <- rep(gstep, gene_num)
    dtype$y1 <- dtype$y1+ log2(dtype$y1)*0.1

    n <- 0
    for(i in unique(dtype$type)){
        n <- n +1
        dtype$y2[dtype$type == i] <- dtype$y2[dtype$type == i] * n
    }
    dtype$y2 <- dtype$y2 + gstep/2

    colors1 <- get_color(group_num, alpha=0.5)
    colors2 <- get_color(group_num, alpha=1)
    maxx <- 5+maxlen*0.5

    p2 <- ggplot(dtype)+ geom_segment(aes(x1, y1, xend=x2, yend=y2, color=type), linewidth=0.5)+
        geom_point(aes(x=x1, y=y1), size=2.5, fill="#004eaf", color="#004eaf", stroke=0.5, shape=21)+
        geom_point(aes(x=x2, y=y2, fill=type, color=type), size=3.5, stroke=0.5, shape=21)+
        scale_y_continuous(limits=c(0.5, gene_num+1), expand=c(0,0))+
        scale_color_manual(values=colors1)+
        scale_fill_manual(values=colors2)+
        scale_x_continuous(breaks=c(1, 5), labels=c("Gene", "Host")) +
        theme_void()+ coord_cartesian(xlim=c(1, maxx)) +
        theme(panel.background=element_blank(),  # 设置面板背景为空
            plot.background=element_blank(),  # 设置绘图背景为空
            legend.background=element_blank(),  # 设置图例背景为空
            legend.position="none",
            axis.text.x=element_text(color="black", size=font),
            plot.margin=margin())+
        geom_text(data=subset(dtype, !duplicated(type)), aes(x2+maxlen*0.045, y2, label=type), hjust=0, size=(font-2)*0.45)
    #p2 <- p2 + geom_text(aes(x1-0.5, y1, label=gene))#检查y轴是否对应 
    #p <- p %>% ggdraw() + draw_plot(p2, scale=0.93, x=0.12, y=0.023) 
    p <- p+p2+plot_layout(widths=c(3, network_width), guides="collect") +
        theme(plot.background=element_rect(fill=NA, color=NA),
              panel.background=element_rect(fill=NA, color=NA),
              legend.background=element_rect(fill=NA, color=NA))
   
    print(paste("图片宽度(--width)", width, sep=":"))
    print(paste("图片高度(--height)", height, sep=":")) 
    ggsave(p, filename=paste(prefix, "heatmap_network.png", sep="."), units="mm", width=width, height=height, dpi=300)
    ggsave(p, filename=paste(prefix, "heatmap_network.pdf", sep="."), units="mm", width=width, height=height, device=cairo_pdf)
   
}


args <- add_help_args()
plot_heatmap_network(data=args$input, group=args$group, dtype=args$dtype,
                     dmethod=args$dmethod, prefix=args$prefix, width=args$width,
                     height=args$height, font=args$font, network_width=args$network_width,
                     method=args$method, padj=args$padj, legabel=args$legabel)
