#!/usr/bin/env Rscript

library("plspm")
library("argparse") #传参

options(bitmapType="cairo") #关闭服务器与界面的互动响应

add_help_args <- function(){

    usage<- "\nVersion: v1.1.0\nAuthor: whbeifan\nEmail: whbeifan@foxmail.com\n"
    desc <- "example: plot_PLS-SEM.R data_sem.tsv --tran --variable_list sem_group.list -p prefix\n
    Input data_sem.tsv file format:
Sample	A1	A2	A3	B1	B2	B3	B4
Lac	30.05770115	24.50620385	24.23156081	32.28039118	25.98926508	27.3777534	31.30899522
AcOH	2.382041205	6.590313708	1.535525539	1.922100998	2.247519665	1.329712621	2.211491444
PA	2.825429147	3.510159526	2.467233239	2.031658301	2.369651521	2.722149721	2.441889107
Time	3	3	3	90	90	90	90
OTC_degrade	0.085730827	0.020423321	0.247927042	0.01773241	0.051466853	0.025669926	0.008779383
DNA_Shannon	4.0234	4.7188	4.3106	6.0324	5.9515	6.1598	3.8772
RNA_Shannon	2.9233	2.9738	2.9738	2.8129	2.8048	2.8069	2.8009
Sum_DNA-ARG	9450.72	12529.81	13653.11	10840.28	7363.55	9986.3	8044.7
Sum_RNA-ARG	1047.19	1016.64	1214.72	417.39	278.07	265.97	2025.24"
    
    parser <- ArgumentParser(description=desc, usage=usage)
    parser$add_argument("-v", "--version", action="store_true", default=FALSE,
        help="Print version information.")
    parser$add_argument("input", type="character",
        help="Input environmental factors and microbial information(environmental_factor.tsv)")
    parser$add_argument("--tran", action="store_true", default=FALSE,
        help="Set to transpose the second matrix, default=FALSE")
    parser$add_argument("-vl", "--variable_list", type="character", required=TRUE,
        help="Input the variable file to be analyzed(variable.list)")
    parser$add_argument("--latent_vars", type="character", default="",
        help="Input corresponding relationship model(latent_vars.list)")
    parser$add_argument("-p", "--prefix", type="character", default="out",
        help="Output result prefix, default=out")
    parser$add_argument("--height",type="integer", default=200,
        help="Set the figure height(mm), default=200")
    parser$add_argument("-w", "--width", type="integer", default=180,
        help="Set the figure width(mm), default=180")

    args <- parser$parse_args()

    return(args)
}


get_group_index <-  function(groups){

    #获取分组的索引
    groupid <- unique(groups$group)

    r <- list() #获取模型的坐标
    for (i in groupid){
        r[[i]] <- groups$microbe[which(groups$group==i)]
    }

    return(r)
}


create_structure_matrix <- function(group, latent_vars){

    group <- unique(group)
    n <- length(group)
    r <- matrix(0, nrow=n, ncol=n, dimnames=list(group, group))
    latent_vars <- read.delim(latent_vars, sep="\t", head=FALSE)

    for(i in 1:nrow(latent_vars)){
        f1 <- latent_vars[i, 1] 
        f2 <- latent_vars[i, 2]
        x <- which(group == f1)[1]
        y <- which(group == f2)[1]
        if(y<x){  #确保为下三角
            r[f1, f2] <- 1
        }else{
            r[f2, f1] <- 1
        }
    }

    return(r)
}


create_lower_tri <-function(names){

    #创建下三角矩阵
    num <- length(names)
    r <- matrix(0, nrow=num, ncol=num)
    r[lower.tri(r, diag=TRUE)] <- 1
    r[row(r) == col(r)] <- 0 #将对角线的设置为0

    colnames(r) <- names
    rownames(r) <- names

    return(r)
}


plot_plot_PLS_SEM <- function(data, tran=FALSE, variable_list="", latent_vars="", prefix="out",
                     height=100, width=600){

      data <- read.delim(data, sep="\t", row.names=1, head=TRUE, check.names=FALSE, quote="")
    if(tran){
        data <- t(data)
    }

    if(variable_list!=""){ #提取需要的变量
        group <- read.delim(variable_list, sep="\t", stringsAsFactors=FALSE, header=FALSE)
        colnames(group) <- c("microbe", "group")
        rownames(group) <- group$microbe
        var_id <- colnames(data)
        var_id <- intersect(group$microbe, var_id) #通过分组文件指定分组顺序和对应关系
        data <- data[, var_id]
        group <- group[var_id,]
    }else{
        var_id <- colnames(data)
        group <- data.frame(microbe=var_id, group=var_id)
    }
    rus_blocks <- get_group_index(group)

    if(latent_vars!=""){ 
        rus_path <- create_structure_matrix(group=unique(group$group), latent_vars=latent_vars) #根据指定保留的路径生成下三角矩阵
    }else{
        rus_path <- create_lower_tri(unique(group$group)) #生成全部路径的下三角矩阵
    }
    print(rus_path)
    modes <- c()
    for(i in names(rus_blocks)){
        if(length(rus_blocks[[i]])==1){
            modes <- c(modes, "B") #单指标用Mode B
        }else{
            modes <- c(modes, "A")
        }
    }

    fit <- plspm(Data=data, #数据加载 
             path_matrix=rus_path,  #内模型（结构模型） 
             blocks=rus_blocks, #外模型（观测模型） 
             modes=modes,  # A=反映型，B=形成型,参数个数与外模型blocks的潜变量数量一致 
             boot.val=TRUE, br=500)  # 启用Bootstrap检验
    print(fit$outer_model) #外模型分析
    print("Loading（载荷）：反映显变量与潜变量相关性，>0.7为优，若>0.5且显著则可保留，<0.5可考虑剔除。") 
    print("--------------------------------------------------------------------------------------------")
    print(fit$unidim) #单维性检验
    print("单维性检验：需满足以下三个指标：")
    print("Cronbach's alpha：>0.7表示内部一致性良好，>0.6可接受。")
    print("Dillon-Goldstein's rho：>0.7表明潜在变量能较好解释观测变量的共同方差。")
    print("第一特征值：需大于1，且第二特征值小于1。")
    print("--------------------------------------------------------------------------------------------")
    print(fit$inner_model) #内模型分析
    print("R²值：反映内生潜变量的解释力度。一般认为：R²>0.75为强解释力，0.5-0.75为中等，<0.25为弱。批注：只有因变量（POLINS）才有R²，自变量（AGRIN, INDEV）不存在R² 。其中Pr(>|t|)为P值")
    print(fit$inner_summary)
    print("--------------------------------------------------------------------------------------------")
    print(fit$effects) #效应值
    print(fit$path_coefs) #直接效应
    print("--------------------------------------------------------------------------------------------")
    print(fit$gof) #模型拟合优度Goodness-of-Fit
    print("GOF>0.7为优秀，0.5-0.7可接受，<0.5需优化模型。")
 
    #print(summary(fit)) 
    png(paste(prefix, "PLS-SEM.png", sep="."), width=width, height=height, units="mm", res=350)
    innerplot(fit, 
          colpos="#CA5023",  # 正向路径红色 
          colneg="#457CC3",  # 负向路径蓝色
          show.values=TRUE,   # 显示系数值
          lcol="black",    #变量名称边框颜色，同步给路径系数相同的颜色 
          box.col="orange",   #变量边框的填充色 
          box.size=0.1,    #变量边框的大小，相对于整个画面的比例 
          box.lwd=0.5,  # 边框粗细
          arr.pos=0.5)   #箭头相对位置，改变该值会同步改变路径系数位置 
    dev.off()

}


args <- add_help_args()
plot_plot_PLS_SEM(data=args$input, tran=args$tran, variable_list=args$variable_list,
    latent_vars=args$latent_vars, prefix=args$prefix, height=args$height, width=args$width)
