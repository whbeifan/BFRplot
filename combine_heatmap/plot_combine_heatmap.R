#!/usr/bin/env Rscript

library("psych")
library("vegan")
library("ggplot2")
library("reshape2")
library("circlize")
library("argparse") #传参
library("ComplexHeatmap")

options(bitmapType="cairo") #关闭服务器与界面的互动响应


add_help_args <- function(){
  
    version<- "Version:v1.2.2\nAuthor:Whbeifan\nEmail:whbeifan@foxmail.com\n"
    desc <- "plot_combine_heatmap.R: plot a correlation heatmap of two matrices."
    usage <- "plot_combine_heatmap.R genus.tsv --data1 differential_metabolites.tsv -p prefix\n"

    parser <- ArgumentParser(description=desc, usage=usage)
    parser$add_argument("-v", "--version", action="version", version=version,
        help="Print version information.")
    parser$add_argument("data", type="character",
        help="Input the first matrix(genus.tsv)")
    parser$add_argument("-dn", "--data_name", type="character", default="Microbiome",
        help="Input the name of the first matrix, default=Microbiome")
    parser$add_argument("-d1", "--data1", type="character", required=TRUE,
        help="Input the second matrix(differential_metabolites.tsv)")
    parser$add_argument("-dn1", "--data1_name", type="character", default="Metabolome",
        help="Input the name of the second matrix, default=Metabolome")
    parser$add_argument("-mp", "--pmax", type="double", default=0.05,
        help="Filter rows with correlations greater than P value, default=0.05")
    parser$add_argument("-mr", "--rmin", type="double", default=0.5,
        help="Filter rows with correlation less than R value, default=0.6")
    parser$add_argument("-p", "--prefix", type="character", default="out",
        help="Output result file  prefix, defaut=out")
    parser$add_argument("-w", "--width", type="double", default=0,
        help="Set the figure width")
    parser$add_argument("-ht", "--height", type="double", default=0,
        help="Set the figure height")

    args <- parser$parse_args()

    return(args)
}


read_data <- function(data){
  
    #读取文件
    data <- read.table(data, header=T, row.names=1, sep="\t", comment.char="", quote="", check.names=FALSE)
    #check.names=FALSE 保留原始表头
    data <- data[which(rowSums(data) > 0),] #过滤所有丰度都为0得行
    data <- t(data)

    return(data) 
}


filter_pvalue_row <- function(data, pmax=0.05){

    #过滤P值均大于0.05的行
    r <- c()

    for (i in rownames(data)){
        if(min(data[i, ]) <= pmax){ #存在P值小于等于0.05的行
            r <- c(r, i)
        }else{
            print(i)
        }

    }
    return(r)
}


filter_rvalue_row <- function(data, rmin=0.6){

    #过滤R值均小于0.6的行
    r <- c()

    for (i in rownames(data)){
        if(max(abs(data[i, ])) >= rmin){ #存在R值大于等于0.6的行
            r <- c(r, i)
        }else{
            print(i)
        }

    }
    return(r)
}


mark_pvalue <- function(pvalue){
 
    #筛选标注哪些小于0.01 小于0.05并标注 
    if(!is.null(pvalue)){
        site1 <- pvalue<0.01
        pvalue[site1] <- "**"
        site2 <- pvalue > 0.01& pvalue < 0.05
        pvalue[site2] <- "*"
        pvalue[!site2&!site1]<- ""
    } else{
        pvalue <- F
    }
    return(pvalue)
}


filter_pvalue <- function(pvalue, rvalue){

    indexs <- which(apply(pvalue, 1, min, na.rm=TRUE) < 0.05)
    pvalue <- pvalue[indexs,]
    rvalue <- rvalue[indexs,]

    r <- list(pvalue=pvalue, rvalue=rvalue)

    return(r)
}


max_name <- function(samples){

    maxlen <- 0
    for(i in samples){
       if(maxlen <= nchar(i)){
           maxlen <- nchar(i)
       }
    }

    return(maxlen)

}


get_size_picture <- function(sample_num, otu_num, metabolic_num, height=0, width=0){

    if(width<=0){
        width <- metabolic_num*5.5+sample_num*4+18+40
        width <- width*0.06+0.2
    }
    if(height<=0){
        height <- otu_num*3.8+sample_num*3.9+20+40
        height <- height*0.06
    }

    return(c(as.numeric(height), as.numeric(width)))
}


plot_correlation_analysis <- function(data1, data2, data_name="", data1_name="", prefix="out",
                                      pmax=0.05, rmin=0.6, height=0, width=0){
  
    data1 <- read_data(data1) #物种丰度表
    data2 <- read_data(data2)	#代谢丰度表
    sample_id1 <- rownames(data1)
    sample_id2 <- rownames(data2)
    sample_id <- intersect(sample_id1, sample_id2) #获取共有样本名称
    
    data1 <- data1[sample_id, ] #匹配样本
    data2 <- data2[sample_id, ] #匹配样本
    data1 <- data1[,colSums(data1)!=0] #过滤全为0的列
    
    result <- corr.test(data1, data2, method="spearman", adjust="none") #计算相关性 pearson
    rvalue <- result$r
    pvalue <- result$p

    rowpid <- filter_pvalue_row(pvalue, pmax) #获取有P值小于0.05的行
    pvalue <- t(pvalue[rowpid, ])
    colpid <- filter_pvalue_row(pvalue, pmax) #获取有P值小于0.05的列
    rvalue <- rvalue[rowpid, colpid]
    rowrid <- filter_rvalue_row(rvalue, rmin) #获取有R值小于0.6的行
    rvalue <- t(rvalue[rowrid,])
    colrid <- filter_rvalue_row(rvalue, rmin) #R值小于0.6的行
    rvalue <- rvalue[colrid, ]
    pvalue <- pvalue[colrid, rowrid]

    r <- filter_pvalue(pvalue, rvalue)
    rvalue <- r$rvalue
    pvalue <- r$pvalue
  
    result <- melt(rvalue, value.name="cor")
    result$pvalue <- as.vector(pvalue)
    write.table(result, file=paste(prefix, ".correlation_pvalue.xls", sep=""), row.names=FALSE, sep="\t", quote=FALSE, na ="")

    data1 <- scale(data1)
    data2 <- scale(data2)
    rvalue <- t(rvalue)
    pvalue <- mark_pvalue(pvalue)
    pvalue <- t(pvalue)
    data2 <- data2[,match(colnames(rvalue), colnames(data2))]
    data1 <- t(data1[,match(rownames(rvalue), colnames(data1))])

    wh <- get_size_picture(ncol(data1), nrow(data1), ncol(data2), height, width)
    widths <- wh[2]
    heights <- wh[1]
    if(max_name(rownames(data1))>=20){
        widths <- widths+(max_name(rownames(data1))-20)*0.15 + 0.2
    }
    print(paste("widths", widths, sep=":"))
    print(paste("heights", heights, sep=":"))

    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    ha <- rowAnnotation(empty=anno_empty(border=FALSE),  Microbiome=data1, col=list(Microbiome=col_fun), 
        show_legend=FALSE, annotation_label=data_name, foo = anno_text(rownames(data1))
    )

    p1 <- Heatmap(rvalue, name="Correlation", cluster_columns=TRUE, cluster_rows=TRUE, show_row_names=FALSE,
        show_heatmap_legend=TRUE, 
        width=ncol(data2)*unit(5.5,"mm"), height=nrow(data1)*unit(3.8,"mm"),
        column_dend_height=unit(20,"mm"), row_dend_width=unit(18,"mm"), #聚类树的高度和宽度
        cell_fun=function(j, i, x, y, width, height, fill){
            grid.text(sprintf("%s",pvalue[i, j]), x, y, gp=gpar(fontsize=10))
        }, right_annotation=ha)

    p2 <- Heatmap(data2, name=data1_name, column_title=data1_name,
        column_title_side="bottom", cluster_columns=FALSE, show_row_names=TRUE, #show_row_names=T 添加样本名称
        row_names_side="left", cluster_rows=FALSE,
        width=ncol(data2)*unit(5.5,"mm"), height=ncol(data1)*unit(3.9,"mm"))

    lgd.list <- Legend(col_fun=col_fun, title=data_name)
    ht <- p1%v%p2

    pdf(paste(prefix, ".combine_heatmap.pdf", sep=""), width=widths, height=heights)
    ptemp <- dev.cur()
    png(paste(prefix, ".combine_heatmap.png", sep=""), width=widths*70, height=heights*75)
    dev.control("enable")

    draw(ht, annotation_legend_list=lgd.list, column_title="",  #不显示标题column_title="Metabolome" 
        column_title_side="bottom", column_title_gp=gpar(fontface="bold", fontsize=20),
        padding=unit(c(2, 2, 10, 2), "mm"))
    decorate_annotation("Microbiome", { 
        grid.text(data_name, gp=gpar(fontface="bold", fontsize=20) ,
        y=unit(1, "npc") + unit(2,"mm"), just="bottom")}
    )

    #ggsave(ht, file=paste(prefix, ".combine_heatmap.pdf", sep=""), height=heights, width=widths)
    #ggsave(ht, file=paste(prefix, ".combine_heatmap.png", sep=""), height=heights, width=widths, type="cairo-png", bg="white")
    dev.copy(which=ptemp)
    dev.off()
    dev.off()
}


args <- add_help_args()
plot_correlation_analysis(data1=args$data, data2=args$data1, data_name=args$data_name, data1_name=args$data1_name,
                          prefix=args$prefix, pmax=args$pmax, rmin=args$rmin, width=args$width, height=args$height)
