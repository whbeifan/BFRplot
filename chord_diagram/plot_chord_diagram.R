#!/usr/bin/env Rscript

library("ggsci") #Science、Nature等高分期刊配色包
library("statnet")
library("stringr")
library("circlize")
library("argparse")
library("ComplexHeatmap") #可用此包添加图例
library("grid") #可用此包调整画板
library("RColorBrewer")

options(bitmapType="cairo")


add_help_args <- function(){

    version <- "Version: v1.1.0\tAuthor: Beifan\tEmail: whbeifan@foxmail.com"
    desc <- "plot_chord_diagram: plot a chord diagram(绘制和弦图)"
    usage <- "plot_chord_diagram.R species_chord.tsv --prefix out\n
Input species_chord.tsv file format:
    Classify\tGH\tGT\tCBM
    Bacteroides\t1156\t767\t255
    Duncaniella\t909\t888\t233\n"

    parser <- ArgumentParser(description=desc, usage=usage)
    parser$add_argument("-v", "--version", action="version", version=version,
        help="Print version information.")
    parser$add_argument("input", metavar="FILE", type="character",
        help="Input chord file(genus_chord.tsv)")
    parser$add_argument("-rg", "--row_group", metavar="FILE", type="character", default="",
        help="Input row grouping file")
    parser$add_argument("-cl", "--col_group", metavar="FILE", type="character", default="",
        help="Input col grouping file")
    parser$add_argument("-p", "--prefix", metavar="STR", type="character", default="out",
        help="Output result prefix, default=out")
    parser$add_argument("--height", metavar="FLOAT", type="double", default=140,
        help="Set the figure height(mm), default=140")
    parser$add_argument("-w", "--width", metavar="FLOAT", type="double", default=200,
        help="Set the figure width(mm), default=200")
    parser$add_argument("-lm", "--lower_margin", metavar="FLOAT", type="double", default=0.32,
        help="Set the lower margin of the circle(r), default=0.32")
    parser$add_argument("-um", "--upper_margin", metavar="FLOAT", type="double", default=0.1,
        help="Set the upper margin of the circle(r), default=0.1")

    args <- parser$parse_args()
    return(args)
}

COLORS <- c("#B0C4DE", "#8470FF", "#8B5F65", "#FF6A6A", "#8B864E", "#8B4513", "#EEAD0E", "#8B3A62",
            "#FF1493", "#A2CD5A", "#8B008B", "#00008B", "#458B74", "#8B8378", "#8B2323",
            "#00BFFF", "#8A2BE2", "#0000FF", "#FFD700", "#458B00", "#EE3B3B")

get_colors <- function(color_num, alpha=0.5){

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
    }else if(color_num <= 51){ ##最多可以获取51种颜色
        colors <- pal_igv("default")(color_num)
    }else{
        temp <- brewer.pal(n=12, name="Paired")
        pal <- colorRampPalette(temp)
        colors <- pal(color_num)
    }

    return(colors)
}


get_group_color <- function(groups, alpha=0.8){

    group <- unique(unlist(groups$group))
    if(length(group) <= 10){
        colors <- get_colors(length(group), alpha=alpha)
    }else if(length(group) <= 21){
        colors <- COLORS
    }else{
        colors <- get_colors(length(group), alpha=alpha)
    }


    group <- groups$group

    n <- 0
    r <- c()
    temp <- ""

    for(i in group){
        if(i!=temp){
            n <- n + 1
        }
        temp <- i
        r <- c(r, colors[n])

    }

    return(r)
}


polt_chord_diagram <- function(data, row_group, col_group, prefix="out", height=140, width=200, lower_margin=0.32, upper_margin=0.1){

    data <- read.delim(data, sep="\t", stringsAsFactors=FALSE, row.names=1, head=TRUE, check.names=FALSE, quote="")
 
    #对行进行分组
    if(row_group!=""){
        row_group <- read.delim(row_group, sep="\t", stringsAsFactors=FALSE, head=TRUE, check.names=FALSE, quote="")
        colnames(row_group)[1:2] <- c("sample", "group")
        row_group <- row_group[str_order(row_group$group), ]
        rownames(row_group) <- row_group$sample
        sample <- row_group$sample
        sample1 <- colnames(data)
        sample <- intersect(sample, sample1)
        row_group <- row_group[sample, ]
        data <- data[, sample]
        rowg_color <- get_group_color(row_group, alpha=0.5)
        gcdata <- data.frame(group=unique(row_group$group), color=unique(rowg_color))
        row_group <- TRUE
    }

    #对列进行分组
    if(col_group!=""){
       col_group <- read.delim(col_group, sep="\t", stringsAsFactors=FALSE, head=TRUE, check.names=FALSE, quote="")
       colnames(col_group)[1:2] <- c("sample", "group")
       col_group <- col_group[str_order(col_group$group), ]
       rownames(col_group) <- col_group$sample
       sample <- col_group$sample
       sample1 <- rownames(data)
       sample <- intersect(sample, sample1)
       data <- data[sample, ]
       colg_color <- get_group_color(col_group, alpha=0.8)
       gcdata1 <- data.frame(group=unique(unlist(col_group$group)), color=unique(unlist(colg_color)))
       col_group <- TRUE
    }
 
    data <- as.matrix(data)#转换数据格式
    species <- unlist(rownames(data))
    types <- unlist(colnames(data))
    #生成作图数据
    data <- data.frame(from=rep(colnames(data), each=nrow(data)),
                       to=rep(rownames(data), ncol(data)),
                       value=as.vector(data))
    #颜色设定
    color <- NULL
    if(row_group=="" & col_group==""){
        color[types] <- get_colors(length(types), alpha=0.8)
        color[species] <- get_colors(length(species), alpha=0.6)
    }else if(row_group!="" & col_group==""){
        #color[species] <- rep(c("#C0C0C0"), each=length(species))
        color[species] <- get_colors(length(species), alpha=0.6)
        color[types] <- rowg_color
    }else if(col_group!="" & row_group==""){
        color[types] <- get_colors(length(types), alpha=0.8)
        color[species] <- colg_color
    }else{
        color[types] <- rowg_color
        color[species] <- colg_color
    }
 
    #pdf(file=paste(prefix, ".chord_diagram.pdf", sep=""), width=width/25.4, height=height/25.4, onefile=FALSE) #bg="transparent"
    cairo_pdf(file=paste(prefix, ".chord_diagram.pdf", sep=""), width=width/25, height=height/25, onefile=FALSE, bg="transparent")
    fo <- dev.cur()
    png(file=paste(prefix, ".chord_diagram.png", sep=""), width=width, height=height, unit="mm", res=350, bg="transparent")
    dev.control("enable")
    par(mar=c(0.5,0.5,0.5,0.5))

    #species <- species[order(species)]
    #types <- types[order(types)]
    label <- c(species, types)
    circos.par(track.height=0.1, circle.margin=c(0.15, 0.65, lower_margin, upper_margin)) #c(左, 右, 下, 上)
    p <- chordDiagram(data, order=label, grid.col=color, transparency=0.2, link.lwd=0.000001, link.lty=1,
        link.border="black", grid.border="black", small.gap=1, big.gap=10, annotationTrack="grid",
        #directional=1, #表示线条的方向，0代表没有方向，1代表正向，-1代表反向，2代表双向
        diffHeight=mm_h(3), 
        preAllocateTracks=list(track.height=max(strwidth(unlist(dimnames(data))))))

    circos.track(track.index=1, panel.fun=function(x, y){
        circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing="clockwise",
        niceFacing=TRUE, adj=c(0, 0.5)) #cex=1.2
    }, bg.border=NA)  #这个不能缺
    abline(h=0, lty=0, col="#00000080") #设置为水平

    if(row_group!=""){
        legend1 <- Legend(at=gcdata, labels=gcdata$group, labels_gp=gpar(fontsize=12), title="",
            grid_height=unit(0.5, "cm"), grid_width=unit(0.5, "cm"), type="points", pch=NA,
            background=gcdata$color)
        pushViewport(viewport(x=0.85, y=0.70))
        grid.draw(legend1)
        upViewport()
    }
    if(col_group!=""){
        legend2 <- Legend(at=gcdata1, labels=gcdata1$group, labels_gp=gpar(fontsize=12), title="",
            grid_height=unit(0.5, "cm"), grid_width=unit(0.5, "cm"), type="points", pch=NA,
            background=gcdata1$color)
        pushViewport(viewport(x=0.87, y=0.30))
        grid.draw(legend2)
        upViewport()
    }
 
    circos.clear()
    dev.copy(which=fo)
    dev.off()
    dev.off() 
}


args <- add_help_args()
polt_chord_diagram(data=args$input, row_group=args$row_group,
                   col_group=args$col_group, prefix=args$prefix, height=args$height, width=args$width,
                   lower_margin=args$lower_margin, upper_margin=args$upper_margin)
