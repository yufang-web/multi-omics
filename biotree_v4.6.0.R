# !/usr/bin/env R
# About 关于----"支持自定义对比符号 例如：'-' 或 ' VS '等"--wb
# @ sensichip 2020-07-08
# Encoding in UTF-8

# Dependence 环境设置 ----
#Sys.setenv(R_ZIPCMD = "/Program Files/R/Rtools/bin/zip")
#Sys.setenv(R_ZIPCMD="/Rtools/bin/zip")
#options("repos" = "http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
#options("BioC_mirror" = "https://mirrors.ustc.edu.cn/bioc")
#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
#options(stringsAsFactors = F,timeout = 200)
#Sys.setlocale("LC_ALL","Chinese")
library(httr) #网络连接
library(xml2) #网络连接
library(magrittr)
library(Cairo) #pdf绘图系统
library(dplyr)  # 数据格式整理
library(reshape2)
library(ggplot2)
library(fdrtool)
library(ropls)  # MVDA模型 from BioC
library(plotly)  # 3D PCA
library(htmlwidgets)  # 3D PCA
library(ggrepel)  # 气泡图 from CRAN
library(plotROC)  # ROC曲线 from CRAN
library(openxlsx)
library(pheatmap)
library(treemap)  # 矩形树图
library(rmarkdown) # 生成html报告
library(KEGGgraph) # MA
library(pander) # Rmd
library(rvest) # 编辑html报告

# 自定义S4类 ----
setClass("picDescription",
         slots = list(pcaOutOfCi="numeric",
                      plsDaOutOfCi="numeric",
                      oplsDaOutOfCi="numeric",
                      plsDaClassify="logical",
                      oplsDaClassify="logical",
                      plsDaPermutate="numeric",
                      oplsDaPermutate="numeric"
         ),
         prototype = list(pcaOutOfCi=0,
                          plsDaOutOfCi=0,
                          oplsDaOutOfCi=0,
                          plsDaClassify=c(F,F),
                          oplsDaClassify=c(F,F),# 第一个 T 表示左右能分开，第二个 T 表示 位于竖线两侧
                          plsDaPermutate=c(0,0,0,0,0),
                          oplsDaPermutate=c(0,0,0,0,0)
                          # 第一个 数值 表示R2Y值，
                          # 第二个 数值 表示Q2值
                          # 第三个数值（logical） 表示左侧的点是否都低于右侧的点，
                          # 第四个数值表示残差是否小于零，
                          # 第五个数值表示 斜率是否为正
         )
)


# TODO 仍需修改----
# 好多好多！

# Functions 子函数----
# 数据预处理部分 ----
ArgsVarify <- function(
  type, qc, is, group, sample, add.sample, compare, desc, desc.names, 
  filter, norm, task, ci, pca3D, pls.da, deg.vip, deg.p, deg.named, 
  palette, ellipse, heatmap.ignore, 
  keggOrg, maOrg,kegg.deg.only, kegg.db, kegg.pathview, kegg.ori, kegg.in, 
  ma.treemap, ma.bubble, 
  scaling, transf, input, output, out.temp,keggHitMin,
  statisticalMethod,vsSymbol,anova = anova
) {
  # 确认参数和参数整形，这些参数都是上面赋值过的
  # 
  meta <- GroupInfo(group, sample)
  # 全新的分组信息
  if(!is.na(add.sample)) meta <- AddGroup(meta, add.sample)
  if(!is.na(qc)) meta <- AddQC(meta, qc)
  # 额外分组信息
  list.compare <- Str2Vector(compare, split2 = ':')
  if(is.na(desc)) {
    desc <- switch(type, GC = 6, LC = 8, lipid = 11)
  }
  if(anyNA(desc.names) & length(desc.names) == 1) {
    desc.names <- switch(
      type, 
      GC = c('Peak', 'id', 'Similarity', 'rt', 'Count', 'Mass'), 
      LC = c(
        "id", "MS2 name", "MS2 score", "MS1 name", "MS1 ppm", 
        "type", "mz", "rt"
      ),
      lipid = c(
        'id', 'lipid name', 'LipidIon', 'LipidGroup', 'Class', 'FattyAcid', 'MS2 score',
        'Rt', 'ObsMz', 'CalcMz', 'IonFormula'
      )
    )
  }
  if(desc != 1 & length(desc.names) == 1) {
    desc.names <- Str2Vector(desc.names, split2 = NA)
  }
  if(desc != length(desc.names)) {
    stop("Wrong desc or desc.names @ ArgsVarify")
  }
  task <- Str2Vector(tolower(task), split2 = NA)
  subDir <- sapply(
    1:length(list.compare), function(i) {
      # 为每一个对比组设置二级子文件夹名称
      # 形如S2-S1
      #
      index1 <- names(meta)[list.compare[[i]][1]]
      index2 <- names(meta)[list.compare[[i]][2]]
      dirname <- paste0(index1,vsSymbol, index2)
    }
  )
  deg.exp <- expression(
    df.deg$`P-VALUE` < l.a$deg.p & 
      df.deg$VIP > l.a$deg.vip & 
      if(l.a$deg.named) {
        if(l.a$type == "GC") {
          !grepl("Analyte|unknown", df.deg$Peak)
        } else if(l.a$type == "LC" || l.a$type == "QE") {
          !is.na(df.deg$`MS2 name`)
        } else if(l.a$type == "lipid") {
          !is.na(df.deg$`lipid name`)
        }
      } else {
        T
      }
  )
  scaling = sapply(
    Str2Vector(scaling, split2 = NA), function(x){
      switch(
        x, None = 'none',  Ctr = 'center', 
        UV = 'standard',  Par = 'pareto'
      )
    }
  )
  transf = as.logical(Str2Vector(transf, split2 = NA))
  # 多元变量分析缩放和转换
  if(palette == 'preset1') {
    palette <- c(
      "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99", 
      "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6"
    )
  } else if(palette == 'presetWQ'){
    palette <- c(
      "#ff0000","#0000ff","#800080","#ffa500","#004e00","#FF00FF","#00FFFF"
      #红 蓝 
    )
  } else if(palette == 'preset2'){
    palette <-  c("green","blue","#E41A1C","#984EA3","#FF7F00","#fdb933","#A65628","#E5C494",
                  "#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#FFFF33")
  } else {
    palette <- rep(palette, len = length(sample))
  }
  l.a <- list(
    type = type, qc = qc, is = is, group = meta,  # pol = pol, 
    list.compare = list.compare, desc = desc, desc.names = desc.names, 
    subDir = subDir, filter = filter, norm = norm, task = task, 
    ci = ci, pca3D = pca3D, pls.da = pls.da, 
    deg.vip = deg.vip, deg.p = deg.p, 
    deg.named = deg.named, deg.exp = deg.exp, 
    palette = palette, ellipse = ellipse, 
    scaling = scaling, transf = transf, 
    heatmap.ignore = heatmap.ignore, 
    keggOrg = keggOrg, maOrg = maOrg,kegg.deg.only = kegg.deg.only, 
    kegg.pathview = kegg.pathview, kegg.db = kegg.db, 
    kegg.ori = kegg.ori, kegg.in = kegg.in, 
    ma.treemap = ma.treemap, ma.bubble = ma.bubble, 
    input = input, output = output, out.temp = out.temp,keggHitMin = keggHitMin,
    statisticalMethod = statisticalMethod,vsSymbol = vsSymbol,anova = anova
  )
}
Str2Vector <- function(str, split1 = ' ', split2 = '/') {
  # 用于处理参数列表里的字符串参数
  # 将提供的字符串拆分，返回一维对象，有二级返回列表，没二级返回向量
  #
  lev1 <- unlist(strsplit(str, split = split1))
  if(is.na(split2)) {
    return(lev1)
  } else {
    list.lev2 <- lapply(lev1, function(group) {
      group <- unlist(strsplit(group, split = split2))
      if(all(sapply(as.double(group), is.double))) {
        return(as.double(group))
      } else {
        return(group)
      }
    })
    return(list.lev2)
  }
}
GroupInfo <- function(group, sample) {
  # l.a$group的基本信息
  # 
  group <- as.numeric(Str2Vector(group, split2 = NA))
  sample <- Str2Vector(sample, split2 = NA)
  if(length(group) != length(sample)) stop("Wrong group info @ GroupInfo")
  group.length <- cumsum(group)
  total.col <- sum(group)
  list.group <- mapply(
    function(group.sub, group.length.sub) {
      # 真·分组信息，返回分组包含分组开始列和结束列的列表
      # 
      group.info <- rep(F, total.col)
      start <- group.length.sub - group.sub + 1
      stop <- group.length.sub
      group.info[start:stop] <- T
      attr(group.info, "type") <- "ori"
      return(group.info)
    }, 
    group.sub = group, group.length.sub = group.length, SIMPLIFY = F
  )
  names(list.group) <- sample
  return(list.group)
}
AddGroup <- function(meta, add.sample) {
  # 在l.a$group中增加组合分组
  # 
  add.sample <- Str2Vector(add.sample, split2 = NA)
  add.sample.name <- sub("~.*$", "", add.sample)
  add.sample.body <- sub("^.*~", "", add.sample)
  list.add.sample <- lapply(
    add.sample.body, function(asb) {
      asb <- Str2Vector(asb, split1 = "\\+", split2 = NA)
      all(asb %in% names(meta)) || stop("Wrong sample info @ AddGroup")
      add.group <- apply(data.frame(meta[asb]), MARGIN = 1, any)
      attr(add.group, "type") <- "add"
      return(add.group)
    }
  )
  names(list.add.sample) <- add.sample.name
  meta <- c(meta, list.add.sample)
}
AddQC <- function(meta, qc) {
  # 在l.a$group中添加qc信息
  # 
  len <- length(meta[[1]])
  meta <- lapply(
    meta, function(sub.meta) {
      type <- attr(sub.meta, "type")
      out <- c(sub.meta, rep(F, qc))
      attr(out, "type") <- type
      return(out)
    }
  )
  meta[['QC']] <- c(rep(F, len), rep(T, qc))
  attr(meta[['QC']], "type") <- "qc"
  return(meta)
}
BuildProjectFramework <- function() {
  # 2.0版新加函数，根据参数信息建立项目框架，主要是创建目录
  #
  sapply(
    c(l.a$output, l.a$out.temp), function(dir) {
      if(!dir.exists(dir)) dir.create(dir)
    }
  )
  sapply(
    c('results', 'others'), function(dir) {
      dir.t <- paste(l.a$output, dir, sep = '/')
      if(!dir.exists(dir.t)) dir.create(dir.t)
      if (dir == 'results' & !is.na(pol)) {
        dir.p <- paste(dir.t, pol, sep = '/')
        if(!dir.exists(dir.p)) dir.create(dir.p)
      } # 只有result文件夹里需要有POS和NEG
    }
  )
  # 以上为创建零层和第一层文件夹
  if("stat" %in% l.a$task) {  # 创建统计结果文件夹
    dir <- paste0(
      l.a$output, '/results/', ifelse(is.na(pol), '', paste0(pol, '/')),
      'Statistical Analysis'
    )
    if(!dir.exists(dir)) dir.create(dir)
    temp <- l.a$subDir
    if (length(l.a$group[sapply(l.a$group, attr, "type") == "ori"]) > 2)
    {temp <- c(l.a$subDir, 'TOTAL')
    }
    if(!is.na(l.a$qc)) temp <- c(temp, 'TOTAL with QC')
    sapply(
      temp, function(subDir) {
        subDir <- paste(dir, subDir, sep = '/')
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
  }
  if("heatmap" %in% l.a$task) {  # 创建层次聚类结果文件夹
    dir <- paste0(
      l.a$output, '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'Hierarchical Clustering Analysis'
    )
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      c(l.a$subDir), function(subDir) {
        subDir <- paste(dir, subDir, sep = '/')
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
  }
  if("zscore" %in% l.a$task) {  # 创建层次聚类结果文件夹
    dir <- paste0(
      l.a$output, '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'Zscore'
    )
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      c(l.a$subDir), function(subDir) {
        subDir <- paste(dir, subDir, sep = '/')
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
  }
  if("corrplot" %in% l.a$task) {  # 创建层次聚类结果文件夹
    dir <- paste0(
      l.a$output, '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'cor_heatmap'
    )
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      c(l.a$subDir), function(subDir) {
        subDir <- paste(dir, subDir, sep = '/')
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
  }
  if("anova" %in% l.a$task){
    dir <- paste0(
      l.a$output, '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'ANOVA'
    )
    if(!dir.exists(dir)) dir.create(dir)
  }
  
  if("forecast" %in% l.a$task) {  # 创建模型预测结果文件夹
    dir <- paste0(
      l.a$output, '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'Prediction Analysis'
    )
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      c(l.a$subDir), function(subDir) {
        subDir <- paste(dir, subDir, sep = '/')
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
  }
  if("lipidbubble" %in% l.a$task) {  # 创建脂质气泡图结果文件夹
    dir <- paste0(
      l.a$output, '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'Bubble Plot'
    )
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      c(l.a$subDir), function(subDir) {
        subDir <- paste(dir, subDir, sep = '/')
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
  }
  if("kegg" %in% l.a$task) {  # 创建KEGG分析结果文件夹
    dir <- paste0(
      l.a$output, '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'KEGG Analysis'
    )
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      l.a$subDir, function(subDir) {
        subDir <- paste(dir, subDir, sep = '/')
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
    dir <- paste0(l.a$out.temp, '/KEGG.in',
                  ifelse(is.na(pol), '', pol))
    if(!dir.exists(dir)) dir.create(dir)
    dir <- paste0(l.a$out.temp, '/KEGG.out',
                  ifelse(is.na(pol), '', pol))
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      l.a$subDir, function(subDir) {
        fn <- paste0(dir, '/', subDir, '.txt')
        if(!file.exists(fn)) file.create(fn, showWarnings = F)
      }
    )
  }
  if("pathway" %in% l.a$task) {
    dir <- paste0(
      l.a$output, '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'Pathway Analysis'
    )
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      c(l.a$subDir), function(subDir) {
        subDir <- paste(dir, subDir, sep = '/')
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
    dir <- paste0(l.a$out.temp, '/MA.out',
                  ifelse(is.na(pol), '', pol))
    if(!dir.exists(dir)) dir.create(dir)
    sapply(
      l.a$subDir, function(subDir) {
        subDir <- paste0(dir, '/', subDir)
        if(!dir.exists(subDir)) dir.create(subDir)
      }
    )
  }
}
ReadInput <- function() {  # 原始数据输入
  require(openxlsx)
  files <- list.files(l.a$input)
  if (l.a$type == 'GC') {
    input <- grep('compare.*\\.csv$|similarity.*\\.csv', files, value = T)
  } else {
    input <- grep(
      paste0(".*(", pol, ')?.*(\\.csv|\\.xlsx)'), 
      files, ignore.case = T, value = T
    )
  }
  if(length(input) > 1) {
    cat(
      noquote(
        'Multiple compare chart detected, input number to specify:\n'
      )
    )
    print(input)
    select <- as.numeric(readline('Select one file: '))
    input <- input[select]
  }
  subfix <- sub("^.*\\.(.*?)$", "\\1", input)
  if (subfix == "csv") {
    df.ori <- read.csv(
      paste(l.a$input, input, sep = '/'), 
      stringsAsFactors = F, check.names = F, na.strings = c("", "NA")
    )
  } else if (subfix == "xlsx") {
    df.ori <- openxlsx::read.xlsx(
      paste(l.a$input, input, sep = '/'), 
      na.strings = c("NA", "N/A", "N/F", "")
    )
  }
  {  # 清理LC物质名称的蛋疼格式，用的时候小心点 #TEST
    # 两个\n
    df.ori[grepl("^\\n\\n", df.ori)] <- sub(
      "^\\n\\n", "", df.ori[grepl("^\\n\\n", "", df.ori)]
    )
    # (drug)
    df.ori[grepl("^Acadesine .*rug\\)$", df.ori)] <- sub(
      " .*?rug\\)$", 
      " \\(drug\\)", 
      df.ori[grepl("^Acadesine .*rug\\)$", df.ori)]
    )
    # Dehydroascorbic acid  (Oxidized vitamin C)两个空格
    df.ori[grepl("^Dehydroascorbic acid\\s+(Oxidized vitamin C)$", df.ori)] <- sub(
      "\\s+", "\\s", 
      df.ori[grepl("^Dehydroascorbic acid\\s+(Oxidized vitamin C)$", df.ori)]
    )
    # Uridine diphosphate glucose (UDP-D-Glucose)
    df.ori[grepl("^Uridine diphosphate glucose.*DP-D-Glucose\\)$", df.ori)] <- sub(
      "^Uridine diphosphate glucose.*?DP-D-Glucose\\)$", 
      "^Uridine diphosphate glucose \\(UDP_D_Glucose\\)", 
      df.ori[grepl("^Uridine diphosphate glucose.*DP-D-Glucose\\)$", df.ori)]
    )
    # 3,4-Dihydroxybenzoate (Protocatechuic acid)
    df.ori[grepl("^3,4-Dihydroxybenzoate.*Protocatechuic acid.*$", df.ori)] <- sub(
      ".*", 
      "3,4-Dihydroxybenzoate \\(Protocatechuic acid\\)", 
      df.ori[grepl("^3,4-Dihydroxybenzoate.*Protocatechuic acid.*$", df.ori)]
    )
    # Acadesine (drug)
    df.ori[grepl("^Acadesine.*rug\\)$", df.ori)] <- sub(
      ".*", 
      "Acadesine \\(drug\\)", 
      df.ori[grepl("^Acadesine.*rug\\)$", df.ori)]
    )
    # 
    # df.ori[grepl("", df.ori)] <- sub(
    #     , 
    #     , 
    #     
    # )
  }  #TEST
  output.file <- paste0(
    l.a$output, '/results/', ifelse(is.na(pol), '', paste0(pol, '/')), 
    'Statistical Analysis/', 
    ifelse(is.na(pol), "compare similarity.", paste0(pol, "-ori.")), 
    subfix
  )
  # 把compare表放到Statistical Analysis文件夹下
  if (subfix == "csv") {
    write.csv(df.ori, file = output.file, row.names = F)
  } else if (subfix == "xlsx") {
    write.xlsx(df.ori, file = output.file)
  }
  # file.copy(
  #     from = paste0(l.a$input, input), 
  #     to = paste0(
  #         l.a$output, '/results/', ifelse(is.na(pol), '', paste0(pol, '/')), 
  #         'Statistical Analysis/', 
  #         ifelse(is.na(pol), "compare similarity.", paste0(pol, "-ori.")), 
  #         subfix
  #     )
  # )
  return(df.ori)
}
# DeDup <- function(df.in) {
#     df.omit <- df.in[!is.na(df.in$`MS2 name`), ]
#     ms2.names <- unique(df.omit$`MS2 name`)
#     l.dedup <- lapply(
#         ms2.names, function(ms2.name) {
#             df.sub <- df.omit[df.omit$`MS2 name` == ms2.name, ]
#             df.sub <- df.sub[df.sub$`MS2 score` == max(df.sub$`MS2 score`), ]
#             if(nrow(df.sub) != 1) {
#                 line.sum <- apply(
#                     df.sub[11:length(df.sub)], MARGIN = 1, sum
#                 )
#                 df.sub <- df.sub[which(line.sum == max(line.sum)), ]
#             }
#             return(df.sub)
#         }
#     )
#     df.dedup <- do.call(rbind, l.dedup)
#     df.out <- rbind(df.dedup, df.in[is.na(df.in$`MS2 name`), ])
#     return(df.out[order(df.out$id), ])
# }
DataReshape <- function(df.in) {
  # 进行数据整形
  # Value: 
  #     df.area：整理好的area信息
  #     
  if(l.a$type == "GC") {
    #    df.area <- df.in[
    #      -1, 
    #     c(1:l.a$desc, seq(l.a$desc + 4, ncol(df.in), by = 4))
    #    ]
    #    colnames(df.area)[(l.a$desc + 1):ncol(df.area)] <- names(df.in)[
    #      seq((l.a$desc + 1), ncol(df.in), by = 4)
    #    ]
    df.area <- df.in
    if(!is.numeric(df.area$Similarity)) {
      df.area$Similarity <- as.numeric(df.area$Similarity)
    }
    df.area$Similarity <- round(df.area$Similarity)
    df.area$Similarity[is.na(df.area$Similarity)] <- 0
  } else {
    df.area <- df.in
  }
  colnames(df.area)[1:l.a$desc] <- l.a$desc.names
  if(!"id" %in% colnames(df.area)) {
    df.area <- cbind(id = seq_len(nrow(df.area)), df.area)
    l.a$desc <<- l.a$desc + 1
  }
  # 变量重命名
  if(!is.na(l.a$is)) {
    # 如果有内标，把内标行放到数据框最下
    is.row <- which(df.area$id == l.a$is)
    df.area <- df.area[
      c(1:(is.row - 1), (is.row + 1):nrow(df.area), is.row), 
    ]
  }
  # id号放到数据框最左
  df.area <- cbind(id = df.area$id, df.area[, colnames(df.area) != "id"])
  if(!is.na(l.a$qc)) {  # 如果有QC，把QC组放到数据框最右
    num.qc <- grep('QC', colnames(df.area), ignore.case = T)
    if(length(num.qc) != l.a$qc) stop("Need to Specify QC @ DataReshape")
    df.area <- df.area[, c((1:ncol(df.area))[-num.qc], num.qc)]
  }
  df.area[, -(1:l.a$desc)] <- sapply(df.area[, -(1:l.a$desc)], as.numeric)
  # 转换数字类型的变量格式为数字
  df.description <<- as.data.frame(df.area[, 1:l.a$desc])
  m.area <- as.matrix(df.area[, -(1:l.a$desc)])
  m.area[which(m.area == 0, arr.ind = T)] <- NA
  # 把0替换为NA
  df.area <- as.data.frame(m.area)
  # 拆分描述部分和峰面积部分
}
RawFilter <- function(df.in) {
  # 对实验组数据中的离群点进行筛选，包括四分位距法和rsd法（即CV法）
  # 四分位距法：设超出四分位点1.5倍四分位间距的点为离群点，删除该离群点。此方
  #     法目前仍有疑虑，当实验组数据特别离散时，可能会掩盖某些确实存在的生物学
  #     现象
  # rsd法：认为QC组rsd大于0.3的检测不稳定，删除该物质的所有检测数据
  # 注意内标行不参与过滤
  #
  if(!is.na(l.a$is)){  # 内标所在行不参与过滤
    IS.serial <- which(df.description$id == l.a$is)
    line.IS <- df.in[IS.serial, ]
    df.in <- df.in[-IS.serial, ]
  }
  if(l.a$filter == 'sfw') {  # 四分位距法部分
    if(is.na(l.a$qc)) {
      df.exp <- df.in
    } else {
      df.exp <- df.in[, 1:(ncol(df.in) - l.a$qc)]
    }
    # 提取需要过滤的实验组数据
    list.exp <- sapply(1:nrow(df.exp), function(i, df) {
      # 按行循环
      #
      line <- df[i, ]
      quan <- quantile(line, probs = c(0.25,0.75), na.rm = T)
      # 求该物质所有实验组数据的四分位点
      line.max <- quan[[2]] + 1.5 * (quan[[2]] - quan[[1]])
      line.min <- quan[[1]] - 1.5 * (quan[[2]] - quan[[1]])
      line[line > line.max | line < line.min] <- NA
      # 计算四分位距并进行过滤
      return(line)
    }, df = df.exp, simplify = F)
    df.exp <- do.call(rbind, list.exp)
    if (!is.na(l.a$qc)) {
      df.in <- cbind(df.exp, df.in[, -(1:(ncol(df.in) - l.a$qc))])
    } else {
      df.in <- df.exp
    }
  } else if(l.a$filter == 'rsd') {  # rsd法部分
    stopifnot(!is.na(l.a$qc))  # 检测参数，rsd法qc不能为NA
    df.qc <- df.in[, (ncol(df.in) - l.a$qc + 1):ncol(df.in)]
    # 提取需要过滤的QC组数据
    df.qc$rsd <- apply(df.qc, MARGIN = 1, function(line) {
      # 按行循环
      #
      line.mean <- mean(line, na.rm = T)
      line.sd <- sd(line, na.rm = T)
      rsd <- line.sd / line.mean
      # 计算该物质的QC组rsd
    })
    df.description <<- df.description[df.qc$rsd < 0.3, ]
    df.in <- df.in[df.qc$rsd < 0.3, ]  # 卡值取0.3
  } else if(l.a$filter == 'none') { #不过滤
    df.in <- df.in
  }
  if (!is.na(l.a$is)) {  # 内标所在行不参与过滤
    df.in <- rbind(df.in, line.IS)  # 内标所在行不参与过滤
  } else {
    df.in <- df.in
  }
}
RecodeNA <- function(df.in) {
  # 判断某物质是否值得保留，若保留则填充实验组数据中的NA
  # 保留条件：1、鉴定组数大于总组数的50%
  #           2、任意一个样本的重复实验中：鉴定数大于重复次数的50%
  # 用实验组area最小值的二分之一填充
  #
  keep <- apply(
    df.in, MARGIN = 1, function(line) {
      if(is.na(l.a$qc)) {
        line.data <- line
        line.qc <- NULL
      } else {
        line.data <- line[1:(ncol(df.in) - l.a$qc)]
        line.qc <- line[(ncol(df.in) - l.a$qc + 1):ncol(df.in)]
      }
      if (
        !is.null(line.qc) & 
        sum(!is.na(line.qc)) < 0.5 * length(line.qc)
      ) {
        return(F)
      }
      temp <- sapply(
        l.a$group[sapply(l.a$group, attr, "type") == "ori"], 
        function(group) {
          line.sub <- line[group]
          if(sum(!is.na(line.sub)) >= 0.5 * length(line.sub)) {
            return(T)
          } else {
            return(F)
          }
        }
      )
      ifelse(any(temp), T, F)
    }
  )
  df.in <- df.in[keep, ]
  df.description <<- df.description[keep, ]
  # 删除不要保留的行，df.in/df.description在此处发生变化
  if(is.na(l.a$qc)) {
    df.data <- df.in
  } else {
    df.data <- df.in[, 1:(ncol(df.in) - l.a$qc)]
  }
  exp.min <- min(apply(df.data, MARGIN = 1, min, na.rm = T))
  # 求实验组area最小值
  df.in <- sapply(df.in, function(column) {
    # 填充最小值的二分之一
    column[is.na(column)] <- exp.min / 2
    return(column)
  })
  # 填充NA
  return(df.in)
}
Normalize <- function(df.in) {
  # 对实验组和QC组进行进一步归一化，方法有内标法和面积和法
  # 内标法：认为加入的内标在各组数据间相同，利用内标进行归一化
  #     注意如果用内标法，则l.a$is不能为NA
  # 面积和法：认为检测物质总量一定，得出的结果转换为
  #     各个物质相对物质总量的相对值
  # 
  # Args:
  #     df.in：data frame，总数据框
  #     norm：string，标明归一化方法，'neibiao' or 'mianji'
  #     is：numeric，内标的id号
  #
  if(l.a$norm == 'neibiao') {  # 内标法
    IS.serial <- which(df.description$id == l.a$is)
    df.in <- sapply(
      as.data.frame(df.in), function(column) {
        column <- column / column[IS.serial]  # 内标归一
      }
    )
    df.in <- df.in[-IS.serial, ]  
    df.description <<- df.description[-IS.serial, ]
    # 校正后删除内标行
  } else if(l.a$norm == 'mianji') {  # 面积法
    if(!is.na(l.a$is)) {  # 校正前如果有内标，删除内标行
      IS.serial <- which(df.description$id == l.a$is)
      df.in <- df.in[-IS.serial, ]
      df.description <<- df.description[-IS.serial, ]
    }
    df.in <- sapply(as.data.frame(df.in), function(column) {
      # 利用面积矫正
      #
      sum.col <- sum(column)
      column <- column / sum.col
    })
  } else {
    return(df.in)  # 杜指导的脚本里就返回空值
  }
  return(df.in)
}
# 数据输出 ----
OutputMean <- function(df.in) {
  # 输出原始数据和MEAN表到客户报告文件夹
  # 
  if(any(sapply(l.a$group, length) != ncol(m.final))) {
    stop("Wrong data frame length @ OutputMean")
  }
  list.data <- lapply(
    1:length(l.a$group), function(i) {
      # 分组计算mean
      # 
      df.sub <- as.data.frame(df.in[, l.a$group[[i]], drop = F])
      df.sub$temp <- apply(df.sub, MARGIN = 1, mean)
      names(df.sub)[ncol(df.sub)] <- paste('Mean', names(l.a$group)[i])
      if(attr(l.a$group[[i]], "type") == "add") {
        df.sub <- df.sub[, ncol(df.sub), drop = F]
      }
      # 修改mean列命名
      return(df.sub)
    }
  )
  df.mean <- do.call(cbind, list.data)
  df.mean <- cbind(df.description, df.mean)
  # 合并输出列表
  wb <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
  addWorksheet(
    wb, sheetName = 'file explanation', gridLines = T, tabColour = '#4F81BD' #设置表格的外观
  )
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    valign = 'center', halign = 'center'
  )
  desc1 <- sapply(colnames(df.mean)[1:l.a$desc], Explanation4Desc)
  desc2 <- sapply(
    colnames(df.mean)[(ncol(df.description) + 1):ncol(df.mean)], 
    function(title) {
      if(grepl('^Mean ', title)) {
        paste0(sub('^Mean ', '', title), '样本重复实验的相对定量值均值')
      } else {
        paste0(title, '样品相对定量值')
      }
    }
  )
  df.desc <- data.frame(
    '表头名称' = colnames(df.mean), '描述' = c(desc1, desc2)
  )
  writeDataTable(
    wb, sheet = 1, df.desc, 
    headerStyle = hs, tableStyle = 'TableStyleMedium2', withFilter = T, startRow = 1, startCol = 1, colName = TRUE
  )
  setColWidths(wb, sheet = 1, cols = 1, widths = 20)
  setColWidths(wb, sheet = 1, cols = 2, widths = 80)
  setRowHeights(wb, sheet = 1, rows = 1, 24)
  setRowHeights(wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), 16)
  addStyle(
    wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), cols = 1:2, 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  addWorksheet(
    wb, sheetName = 'Mean', gridLines = T, tabColour = '#4F81BD'
  )
  writeDataTable(
    wb, sheet = 2, df.mean, headerStyle = hs, 
    tableStyle = 'TableStyleMedium2', withFilter = T, startRow = 1, startCol = 1, colName = TRUE
  )
  setColWidths(wb, sheet = 2, cols = 1:ncol(df.mean), widths = 18)
  setRowHeights(wb, sheet = 2, rows = 1, 24)
  setRowHeights(wb, sheet = 2, rows = 2:(nrow(df.mean) + 1), 16)
  addStyle(
    wb, sheet = 2, cols = grep('Mean ', colnames(df.mean)), 
    rows = 1:(nrow(df.mean) + 1), gridExpand = T, stack = T, 
    style = createStyle(fontColour = '#e41a1c')
  )
  addStyle(
    wb, sheet = 2, rows = 2:(nrow(df.mean) + 1), cols = 1:ncol(df.mean), 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  saveWorkbook(
    wb, overwrite = T, 
    file = paste0(
      l.a$output, '/results/', ifelse(is.na(pol), '', paste0(pol, '/')), 
      'Statistical Analysis/', ifelse(is.na(pol), '', paste0(pol, '-')), 
      'Mean.xlsx'
    )
  )
  return(df.mean)
}
Explanation4Desc <- function(entry) {
  # 为df.description提供表头信息
  # 调用效果比如sapply(colnames(df.description, Explanation4Desc))这样
  # 
  entry <- tolower(entry)  # 不区分大小写
  exp <- switch(
    entry, 
    "id" = "该物质在本次定性分析中的唯一数据编号", 
    "peak" = "定性分析得到的物质名称", 
    "ms2 name" = "二级质谱定性匹配分析得到的物质名称", 
    "compound name" = "质谱定性匹配分析得到的物质名称", 
    "similarity" = paste0(
      '定性分析中该物质与标准库中物质的匹配程度，取值[0, 1000], ', 
      '越高说明定性出的物质越准确'
    ),
    "adduct type" = '加合离子类型',
    "molecular weight" =  '检测离子的分子量(为实际检测值)',
    "monoisotopic molecular weight" = '匹配物质的单同位素精确分子量(为理论值)',
    "formula" = '定性匹配分析得到物质的分子式或预测化学式（unknown）',
    "ms2 formula" = '二级定性匹配分析得到物质的分子式',
    "ms2 score" = '二级质谱定性匹配的打分，取值[0, 100]，越大越好', 
    "score" = '定性匹配的打分，取值[0, 1]，越大越好', 
    "ms1 name" = "一级定性离子匹配的名称，参考HMDB数据库或KEGG数据库",
    "ppm" = "质荷比定性匹配的精确度，越接近0越好",
    "ms1 ppm" = "一级质荷比定性的精确度，取值(-10, 10)，越接近0越好", 
    "ms2 ppm" = "二级定性匹配物质精确分子量与单同位素精确分子量之间的精确度，越接近0越好",
    "rt" = '该物质的色谱保留时间', 
    "rt [min]" = '该物质的色谱保留时间',
    "rt (min)" = '该物质的色谱保留时间',
    "mz" = "检测离子的质荷比", 
    "m/z" = "物质特征离子的质荷比",
    "reference m/z" = '参考化合物的质荷',
    "mass" = "物质特征离子的质荷比",
    "ms2 kegg id" = "二级匹配化合物KEGG(京都基因与基因组百科全书)数据库的ID",
    "ms1 kegg id" = "一级匹配物质对应KEGG(京都基因与基因组百科全书)数据库的ID",
    "kegg id" = "KEGG(京都基因与基因组百科全书)数据库的ID",
    "kegg" = "KEGG(京都基因与基因组百科全书)数据库的ID",
    "count" = '该物质在所有实验组中检测到的次数', 
    "type" = paste0(
      '匹配模式，包括二级正/反向匹配（MS2 forward/reverse）', 
      '一级平均分子量/单同位素分子量匹配', 
      '和无匹配等多种情况'
    ),
    # lipid
    "lipid name" = "脂质名字",
    "lipidion" = "定性匹配分析得到的脂质脂质离子",
    "lipidgroup" = "确定的脂类。脂类是根据总脂肪酸碳数和不饱和数(总和组成)分组显示的脂类。如PC(16:0/18:2)+H脂类组:PC(34:2)+H",
    "class" = "脂质分类",
    "fattyacid" = "脂肪酸链",
    "calcmz" = "脂质离子的理论质合比",
    "obsmz" = "峰值的测量质合比",
    "ionformula" = "脂质离子的组成化学式",
    "mean oo" = "该物质在该组对比内的一个实验组的相对定量值均值", 
    "mean xx" = "该物质在该组对比内的另一个实验组的相对定量值均值", 
    "vip" = "该物质在该组对比的OPLS-DA模型得到的变量投影重要度", 
    "p-value" = paste0(
      "该物质在该组对比的t检验得到的P值，", 
      "P值 = 假设是正确但是被拒绝的概率 = 阴性结果个数/结果总个数，", 
      "是对与样本数据的一个检验概率"
    ), 
    "q-value" = paste0(
      "假设检验统计量（P值）经多重假设检验校正之后的结果，", 
      "Q值 = 被拒绝但却是正确的概率", 
      " = 假阳性结果个数/推测为阳性结果的个数，", 
      "是对统计检验得到的推论的一种检验概率，对P值的再统计"
    ), 
    "fold change" = "该物质在该组对比两组实验间的倍数关系", 
    "log_foldchange" = "FOLD CHANGE取以2为底的对数", 
    "暂无说明"
  )
  return(exp)
}
OutputSIMCA <- function(df.in) {
  # 输出供SIMCA输入的数据
  # 
  rownames(df.in) <- df.description$id
  df.out <- as.data.frame(t(df.in))
  l.class <- l.a$group[sapply(l.a$group, attr, "type") != "add"]
  class.id <- rep(names(l.class), sapply(l.class, sum))
  df.out$`classid` <- class.id
  df.out <- df.out[c(ncol(df.out), 1:(ncol(df.out) - 1))]
  if(!is.na(l.a$qc)) {
    width <- nrow(df.out)
    df.data <- df.out[1:(width - l.a$qc), ]
    df.out <- rbind(df.data, df.out[(width - l.a$qc + 1):width, ])
  }
  # 把QC放在最下
  write.csv(
    df.out, file = paste0(
      l.a$out.temp, '/', 
      ifelse(is.na(pol), '', paste0(pol, "-")), 
      "SIMCA.csv"
    )
  )
  return(df.out)
}
# 统计分析 ----
MVDA <- function(df.in,
                 statisticalMethod = l.a$statisticalMethod) {
  # 计算并输出每组对比的统计分析
  #
  df.data <- df.in[, -1]
  class.id <- MakeClassID(
    l.a$group[sapply(l.a$group, attr, "type") != "add"], 
    use.global.level = T
  )
  mvda.model <- list()
  mvdaPlotAttr <- list()
  if(!is.na(l.a$qc)) {  # 有QC才做TOTAL with QC
    df.model <- df.data
    pca <- opls(
      x = df.model, 
      log10L = l.a$transf[1], scaleC = l.a$scaling[1], 
      crossvalI = min(nrow(df.model), 7), 
      info.txtC = 'none', fig.pdfC = 'none'
    )
    if(pca@summaryDF$pre < 3) {  # 主成分不少于3，方便做三维图
      pca <- opls(
        x = df.model, predI = 3, 
        log10L = l.a$transf[1], scaleC = l.a$scaling[1], 
        crossvalI = min(nrow(df.model), 7), 
        info.txtC = 'none', fig.pdfC = 'none'
      )
    }
    
    picDescription <- new("picDescription")
    
    # 判断记录置信区间以外的点 
    df.p <- as.data.frame(pca@scoreMN[, 1:2])
    colnames(df.p) <- c('x', 'y')
    n <- nrow(df.p)
    hfn <- 2*(n-1)*(n^2-1)/(n^2*(n-2))*qf(l.a$ci, 2, (n-2))
    rr <- (df.p$x)^2/(var(df.p$x)*hfn) + (df.p$y)^2/(var(df.p$y)*hfn)
    picDescription@pcaOutOfCi <- length(rr[rr > 1])
    mvdaPlotAttr <- c(mvdaPlotAttr, "TOTAL with QC" = picDescription)
    rm(df.p,n,hfn,rr,picDescription)
    
    mvda.model <- c(mvda.model, "TOTAL with QC" = pca)
    ModelScorePlot(  # 2D得分图
      model = pca, class.id = class.id, 
      output = paste0(
        l.a$output, '/results/', 
        ifelse(is.na(pol), '', paste0(pol, '/')), 
        'Statistical Analysis/TOTAL with QC'
      )
    )
    if(l.a$pca3D) {
      ModelScorePlot3DNew(  # 3D得分图
        model = pca, class.id = class.id, 
        output = paste0(
          l.a$output, '/results/', 
          ifelse(is.na(pol), '', paste0(pol, '/')), 
          'Statistical Analysis/TOTAL with QC'
        )
      )
    }
  }
  if(sum(sapply(l.a$group, attr, "type") == "ori") > 2) {
    if(!is.na(l.a$qc)) {
      df.model <- df.data[!l.a$group$QC, ]
    } else {
      df.model <- df.data
    }
    class.id <- MakeClassID(
      l.a$group[sapply(l.a$group, attr, "type") == "ori"], 
      use.global.level = T
    )
    pca <- opls(
      x = df.model, 
      log10L = l.a$transf[1], scaleC = l.a$scaling[1], 
      crossvalI = min(nrow(df.model), 7), 
      info.txtC = 'none', fig.pdfC = 'none'
    )
    if(pca@summaryDF$pre < 3) {  # 为三维图加主成分
      pca <- opls(
        x = df.model, predI = 3, 
        log10L = l.a$transf[1], scaleC = l.a$scaling[1], 
        crossvalI = min(nrow(df.model), 7), 
        info.txtC = 'none', fig.pdfC = 'none'
      )
    }
    picDescription <- new("picDescription")
    # 判断记录置信区间以外的点 
    df.p <- as.data.frame(pca@scoreMN[, 1:2])
    colnames(df.p) <- c('x', 'y')
    n <- nrow(df.p)
    hfn <- 2*(n-1)*(n^2-1)/(n^2*(n-2))*qf(l.a$ci, 2, (n-2))
    rr <- (df.p$x)^2/(var(df.p$x)*hfn) + (df.p$y)^2/(var(df.p$y)*hfn)
    picDescription@pcaOutOfCi <- length(rr[rr > 1])
    mvdaPlotAttr <- c(mvdaPlotAttr, "TOTAL" = picDescription)
    rm(df.p,n,hfn,rr,picDescription)
    
    mvda.model <- c(mvda.model, "TOTAL" = pca)
    ModelScorePlot(  # 2D得分图
      model = pca, class.id = class.id, 
      output = paste0(
        l.a$output, '/results/', 
        ifelse(is.na(pol), '', paste0(pol, '/')), 
        'Statistical Analysis/TOTAL'
      )
    )
    if(l.a$pca3D) {
      ModelScorePlot3DNew(  # 3D得分图
        model = pca, class.id = class.id, 
        output = paste0(
          l.a$output, '/results/', 
          ifelse(is.na(pol), '', paste0(pol, '/')), 
          'Statistical Analysis/TOTAL'
        )
      )
    }
  }
  list.compare <- lapply(
    1:length(l.a$list.compare), function(i) {
      pair <- l.a$list.compare[[i]]
      df.a <- df.in[l.a$group[[pair[1]]], -1]
      df.b <- df.in[l.a$group[[pair[2]]], -1]
      row.ab <- nrow(df.a) + nrow(df.b)
      class.id.sub <- MakeClassID(l.a$group[pair], use.global.level = T)
      mean.a <- sapply(df.a, mean)
      mean.b <- sapply(df.b, mean)
      pca <- opls(
        x = rbind(df.a, df.b), predI = 2,
        log10L = l.a$transf[1], scaleC = l.a$scaling[1], 
        crossvalI = min(row.ab, 7), 
        info.txtC = 'none', fig.pdfC = 'none'
      )
      if(pca@summaryDF$pre < 2) {
        pca <- opls(
          x = rbind(df.a, df.b), predI = 2, 
          log10L = l.a$transf[1], scaleC = l.a$scaling[1], 
          crossvalI = min(row.ab, 7), 
          info.txtC = 'none', fig.pdfC = 'none'
        )
      }
      picDescription <- new("picDescription")
      # 判断记录置信区间以外的点 
      df.p <- as.data.frame(pca@scoreMN[, 1:2])
      colnames(df.p) <- c('x', 'y')
      n <- nrow(df.p)
      hfn <- 2*(n-1)*(n^2-1)/(n^2*(n-2))*qf(l.a$ci, 2, (n-2))
      rr <- (df.p$x)^2/(var(df.p$x)*hfn) + (df.p$y)^2/(var(df.p$y)*hfn)
      picDescription@pcaOutOfCi <- length(rr[rr > 1])
      rm(df.p,n,hfn,rr)
      
      ModelScorePlot(  # PCA得分图
        model = pca, class.id = class.id.sub, 
        output = paste0(
          l.a$output, '/results/', 
          ifelse(is.na(pol), '', paste0(pol, '/')), 
          'Statistical Analysis/', l.a$subDir[i]
        )
      )
      
      ModelLoadingPlot( #PCA loading plot
        model = pca,
        output = paste0(
          l.a$output, '/results/', 
          ifelse(is.na(pol), '', paste0(pol, '/')), 
          'Statistical Analysis/', l.a$subDir[i]))
      if (l.a$pls.da) {  # PLS-DA
        set.seed(123)  # 固定随机数种子
        pls.da <- opls(
          x = rbind(df.a, df.b), y = as.character(class.id.sub),
          predI = 2,
          permI = 200, log10L = l.a$transf[2], 
          scaleC = l.a$scaling[2], crossvalI = min(row.ab, 7), 
          info.txtC = 'none', fig.pdfC = 'none'
        )
        # 判断记录置信区间以外的点 
        df.p <- as.data.frame(pls.da@scoreMN[, 1:2])
        colnames(df.p) <- c('x', 'y')
        # 判断是否有点位于置信区间外
        n <- nrow(df.p)
        hfn <- 2*(n-1)*(n^2-1)/(n^2*(n-2))*qf(l.a$ci, 2, (n-2))
        rr <- (df.p$x)^2/(var(df.p$x)*hfn) + (df.p$y)^2/(var(df.p$y)*hfn)
        picDescription@plsDaOutOfCi <- length(rr[rr > 1])
        # 判断是否有效区分两组样本
        min2min <- tapply(df.p$x,pls.da@suppLs$y,min)
        max2max <- tapply(df.p$x,pls.da@suppLs$y,max)
        picDescription@plsDaClassify[1] <- all(min2min[1]*max2max[1] > 0,
                                               min2min[2]*max2max[2] > 0,
                                               any(min2min[1]<0,min2min[2]<0),
                                               any(max2max[1]>0,max2max[2]>0))
        picDescription@plsDaClassify[2] <- any(min2min[1] > max2max[2],
                                               min2min[2] > max2max[1])
        # 判断置换检验是否通过
        df.q <- as.data.frame(pls.da@suppLs$permMN[,c(2,3,7)])
        df.q <- df.q[df.q$sim > 0,]
        picDescription@plsDaPermutate[1] <- df.q$`R2Y(cum)`[df.q$sim ==1]
        picDescription@plsDaPermutate[2] <- df.q$`Q2(cum)`[df.q$sim ==1]
        picDescription@plsDaPermutate[3] <- all(df.q$`Q2(cum)`[df.q$sim ==1] > df.q$`Q2(cum)`[df.q$sim < 1])
        lmY <- df.q$`Q2(cum)` - df.q$`Q2(cum)`[df.q$sim ==1]
        lmX <- df.q$sim - 1
        lmQ2 <-lm(lmY ~ lmX - 1)
        picDescription@plsDaPermutate[4] <- df.q$`Q2(cum)`[df.q$sim ==1] <= lmQ2$coefficients
        picDescription@plsDaPermutate[5] <- lmQ2$coefficients > 0
        
        rm(df.p,n,hfn,rr,lmQ2,min2min,max2max,df.q,lmY,lmX)
        
        VIP <- getVipVn(pls.da)
        ModelScorePlot(  # OPLS-DA得分图
          model = pls.da, class.id = class.id.sub, 
          output = paste0(
            l.a$output, '/results/', 
            ifelse(is.na(pol), '', paste0(pol, '/')), 
            'Statistical Analysis/', 
            l.a$subDir[i]
          )
        )
        PermutationPlot(  # 置换检验图
          model = pls.da, 
          output = paste0(
            l.a$output, '/results/', 
            ifelse(is.na(pol), '', paste0(pol, '/')), 
            'Statistical Analysis/', 
            l.a$subDir[i]
          )
        )
      }
      set.seed(123)  # 固定随机数种子
      opls.da <- opls(
        x = rbind(df.a, df.b), y = as.character(class.id.sub), 
        predI = 1, ortho = 1, permI = 200, 
        log10L = l.a$transf[2], scaleC = l.a$scaling[2], 
        crossvalI = min(row.ab, 7), 
        info.txtC = 'none', fig.pdfC = 'none'
      )
      # 判断记录置信区间以外的点 
      df.p <- cbind(
        as.data.frame(opls.da@scoreMN[, 1]), 
        as.data.frame(opls.da@orthoScoreMN[,1])
      )
      
      colnames(df.p) <- c('x', 'y')
      # 判断是否有点位于置信区间外
      n <- nrow(df.p)
      hfn <- 2*(n-1)*(n^2-1)/(n^2*(n-2))*qf(l.a$ci, 2, (n-2))
      rr <- (df.p$x)^2/(var(df.p$x)*hfn) + (df.p$y)^2/(var(df.p$y)*hfn)
      picDescription@oplsDaOutOfCi <- length(rr[rr > 1])
      # 判断是否有效区分两组样本
      min2min <- tapply(df.p$x,opls.da@suppLs$y,min)
      max2max <- tapply(df.p$x,opls.da@suppLs$y,max)
      picDescription@oplsDaClassify[1] <- all(min2min[1]*max2max[1] > 0,
                                              min2min[2]*max2max[2] > 0,
                                              any(min2min[1]<0,min2min[2]<0),
                                              any(max2max[1]>0,max2max[2]>0))
      picDescription@oplsDaClassify[2] <- any(min2min[1] > max2max[2],
                                              min2min[2] > max2max[1])
      # 判断置换检验是否通过
      df.q <- as.data.frame(opls.da@suppLs$permMN[,c(2,3,7)])
      df.q <- df.q[df.q$sim > 0,]
      picDescription@oplsDaPermutate[1] <- df.q$`R2Y(cum)`[df.q$sim ==1]
      picDescription@oplsDaPermutate[2] <- df.q$`Q2(cum)`[df.q$sim ==1]
      picDescription@oplsDaPermutate[3] <- all(df.q$`Q2(cum)`[df.q$sim ==1] > df.q$`Q2(cum)`[df.q$sim < 1])
      lmY <- df.q$`Q2(cum)` - df.q$`Q2(cum)`[df.q$sim ==1]
      lmX <- df.q$sim - 1
      lmQ2 <-lm(lmY ~ lmX - 1)
      picDescription@oplsDaPermutate[4] <- df.q$`Q2(cum)`[df.q$sim ==1] <= lmQ2$coefficients
      picDescription@oplsDaPermutate[5] <- lmQ2$coefficients > 0
      
      rm(df.p,n,hfn,rr,lmQ2,min2min,max2max,df.q)
      
      VIP <- getVipVn(opls.da)
      ModelScorePlot(  # OPLS-DA得分图
        model = opls.da, class.id = class.id.sub, 
        output = paste0(
          l.a$output, '/results/', 
          ifelse(is.na(pol), '', paste0(pol, '/')), 
          'Statistical Analysis/', 
          l.a$subDir[i]
        )
      )
      PermutationPlot(  # 置换检验图
        model = opls.da, 
        output = paste0(
          l.a$output, '/results/', 
          ifelse(is.na(pol), '', paste0(pol, '/')), 
          'Statistical Analysis/', 
          l.a$subDir[i]
        )
      )
      # select statistical method
      # 0: t test
      # 1: Wilcoxon rank test
      if (statisticalMethod == 0){ 
        `P-VALUE` <- sapply(
          1:ncol(df.a), function(k) {
            # 先检验方差齐性，再计算p-value
            #
            col.a <- df.a[, k]
            col.b <- df.b[, k]
            var <- var.test(col.a, col.b)$p.value > 0.05
            # f检验，方差齐性
            p <- t.test(col.a, col.b, var.equal = var)$p.value  # t检验
          }
        )
      } else if (statisticalMethod == 1){
        `P-VALUE` <- sapply(
          1:ncol(df.a), function(k) {
            col.a <- df.a[, k]
            col.b <- df.b[, k]
            p <- wilcox.test(col.a, col.b)$p.value  # 秩和检验
          }
        )
      }
      
      q <- fdrtool(  # {fdrtool}
        `P-VALUE`, statistic = 'pvalue', plot = F, verbose = F
      )$qval
      fc <- mean.a / mean.b  # fold change
      log.fc <- log2(fc)  # log2 fold change
      df.combine <- data.frame(
        mean.a, mean.b, VIP, `P-VALUE`, q, fc, log.fc
      )
      names(df.combine) <- c(
        paste('MEAN', names(l.a$group)[pair[1]]), 
        paste('MEAN', names(l.a$group)[pair[2]]), 
        'VIP', 'P-VALUE', 'Q-VALUE', 'FOLD CHANGE', 'LOG_FOLDCHANGE'
      )
      df.out <- cbind(df.description, df.combine)
      if (l.a$type == 'GC') {
        df.out <- df.out[order(-df.out$Similarity), ]
      }
      VolcanoPlot(  # 火山图
        df.out, output = paste0(
          l.a$output, '/results/',
          ifelse(is.na(pol), '', paste0(pol, '/')),
          'Statistical Analysis/',
          l.a$subDir[i]
        )
      )
      l.return <- if (l.a$pls.da) {
        list(df.out, picDescription ,pca, pls.da, opls.da)
      } else {
        list(df.out, picDescription , pca, opls.da)
      }
    }
  )
  list.compare.result <- lapply(
    1:length(list.compare), function(i) {list.compare[[i]][[1]]}
  )
  mvda.model.each <- c(
    lapply(
      1:length(list.compare), function(i) {list.compare[[i]][[3]]}
    ), 
    lapply(
      1:length(list.compare), function(i) {list.compare[[i]][[4]]}
    ),
    if (l.a$pls.da) {
      lapply(
        1:length(list.compare), function(i) {list.compare[[i]][[5]]}
      )
    } else {
      NULL
    }
  )
  mvdaPlotAttr.each <- c(lapply(
    1:length(list.compare), function(i) {list.compare[[i]][[2]]}
  ))
  names(mvdaPlotAttr.each) <- l.a$subDir
  mvdaPlotAttr <- c(mvdaPlotAttr,mvdaPlotAttr.each)
  assign('mvdaPlotAttr', mvdaPlotAttr, envir = .GlobalEnv)
  if (l.a$pls.da) {
    names(mvda.model.each) <- rep(l.a$subDir, 3)
  }else{
    names(mvda.model.each) <- rep(l.a$subDir, 2)
  }
  mvda.model <- c(mvda.model,mvda.model.each)
  assign('mvda.model', mvda.model, envir = .GlobalEnv)
  return(list.compare.result)
}
MakeClassID <- function(l.group, use.global.level = F) {
  # 根据提供的l.a$group的子集，返回对应子集的class.id
  # 返回值为因子，因子水平由use.global.level参数决定
  # 
  out <- rep(names(l.group), sapply(l.group, sum))
  if(use.global.level) {
    f.out <- factor(out, levels = names(l.a$group))
  } else {
    f.out <- factor(out, levels = names(l.group))
  }
  return(f.out)
}
ModelScorePlot <- function(model, class.id, output) {
  # PCA/OPLS-DA模型得分图
  # 
  if(model@typeC == 'PCA' | model@typeC == 'PLS-DA') {
    df.p <- as.data.frame(model@scoreMN[, 1:2])
  } else if(model@typeC == 'OPLS-DA') {
    df.p <- cbind(
      as.data.frame(model@scoreMN[, 1]), 
      as.data.frame(model@orthoScoreMN[,1])
    )
  }
  #comp. <- model@modelDF[1]
  #comp.xy <- paste(rownames(comp.), " [", 100 * comp.$R2X, "%]", sep = "")
  #names(comp.xy) <- rownames(comp.)
  
  colnames(df.p) <- c('x', 'y')
  n <- nrow(df.p)
  hfn <- 2*(n-1)*(n^2-1)/(n^2*(n-2))*qf(l.a$ci, 2, (n-2))
  rv <- seq(0, 2*pi, length.out = 100)
  df.ell <- data.frame(  # 置信区间数据
    x = sqrt(var(df.p$x)*hfn)*cos(rv), 
    y = sqrt(var(df.p$y)*hfn)*sin(rv)
  )
  p <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_polygon(data = df.ell, aes(x, y), color = 'grey', fill = NA)
  if(l.a$ellipse) {
    p <- p + stat_ellipse(
      data = df.p, geom = 'polygon', level = l.a$ci, 
      aes(x, y, fill = class.id, color = class.id), alpha = I(0.1)
    ) +
      stat_ellipse(
        data = df.p, geom = 'blank', level = l.a$ci, 
        aes(-x, -y, fill = class.id)
      )
  }
  p <- p + geom_point(
    data = df.p, size = 5, aes(x, y, shape = class.id, color = class.id,fill = class.id)
  ) +
    geom_blank(
      data = df.p, aes(-x, -y, shape = class.id, color = class.id,fill = class.id)
    ) +
    scale_shape_manual(  # 形状参数
      values = rep_len(
        #c("\u2605","\u25C4","\u25BC","\u25B2"), 
        c(21:25),
        length(levels(class.id))
      )[sort(unique(as.numeric(class.id)))]
    ) +
    scale_color_manual(  # 颜色参数
      values = rep_len(
        l.a$palette, length(levels(class.id))
      )[sort(unique(as.numeric(class.id)))]
    ) +
    scale_fill_manual(
      values = rep_len(
        l.a$palette, length(levels(class.id))
      )[sort(unique(as.numeric(class.id)))]
    ) +
    ggtitle(paste(model@typeC, "score plot")) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(size = 0.3),
      aspect.ratio = 1,
      panel.grid.minor = element_line(size = 0.3),
      plot.caption = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 8),
      plot.title = element_text(hjust = 0.5),
      legend.key.size = unit(0.3, "cm"),
      legend.title = element_blank(), 
      legend.key = element_blank(), 
      panel.border = element_rect(color = 'black', size = 0.8)
    ) +
    if(model@typeC == 'PCA') {
      labs(x = "PC[1]",
           y = "PC[2]",
           caption = paste("R2X[1] = ", model@modelDF["p1", "R2X"], "  ", "R2X[2]=", model@modelDF["p2", "R2X"], "  ", "Ellipse: Hotelling's T2(95%)", sep = "")) 
    } else if(model@typeC == 'OPLS-DA') {
      labs(x = 't[1]',
           y = 'to[1]',
           caption = paste("R2X[1] = ", model@modelDF["p1", "R2X"], "   ", "R2Xo[2]=", model@modelDF["o1", "R2X"], "  ", "Ellipse: Hotelling's T2(95%)", sep = ""))
    } else if(model@typeC == 'PLS-DA') {
      labs(x = 't[1]',
           y = 't[2]',
           caption = paste("R2X[1] = ", model@modelDF["p1", "R2X"], "   ", "R2X[2]=", model@modelDF["p2", "R2X"], "  ", "Ellipse: Hotelling's T2(95%)", sep = ""))
    }
  ggsave(  # 输出位置
    paste0(output, '/', model@typeC, ' score plot.jpg'), p, 
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/', model@typeC, ' score plot.pdf'), p,
    device = cairo_pdf,
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
}

## loading plot
ModelLoadingPlot <- function(model, output) {
  if(model@typeC == 'PCA') {
    loading.df <- as.data.frame(
      model@loadingMN
    )
  }
  pcaLoading.plot <- ggplot(loading.df, aes(x = p1, y = p2)) +
    geom_point(fill = "green", shape = 21, color = "steelblue")+
    ggtitle("PCA loading plot")+
    xlab("p1")+
    labs(caption = paste("R2X[1] = ", model@modelDF["p1", "R2X"], "   ", "R2X[2]=", model@modelDF["p2", "R2X"], "  ", "Ellipse: Hotelling's T2(95%)", sep = "")) +
    ylab("p2")+                
    theme_bw() +
    theme(
      panel.grid.major = element_line(size = 0.3),
      aspect.ratio = 1,
      panel.grid.minor = element_line(size = 0.3),
      plot.caption = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 8), #x axis label (element_text; inherits from axis.title)
      plot.title = element_text(hjust = 0.5),
      panel.border = element_rect(color = 'black', size = 0.8)
    )
  ggsave(  # 输出位置
    paste0(output, '/', model@typeC, ' loading plot.jpg'),pcaLoading.plot, 
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/', model@typeC, ' loading plot.pdf'), pcaLoading.plot,
    device = cairo_pdf,
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
}
# ModelScorePlot3D <- function(model, class.id, output) {
#     # PCA3D得分图
#     # 
#     require(pca3d)
#     df.data <- as.data.frame(model@scoreMN[, 1:3])
#     colnames(df.data) <- c("PC1", "PC2", "PC3")
#     shape <- rep_len(
#         c('sphere', 'tetrahedron', 'cube'), length(levels(class.id))
#     )
#     names(shape) <- levels(class.id)
#     pca3d::pca3d(
#         m.data, axes.color = 'black', new = F, radius = 1.5, group = class.id,
#         shape = shape, legend = c('topright'), show.scale = T,
#         show.plane = F, show.ellipses = l.a$ellipse, show.centroids = F,
#         ellipse.ci = l.a$ci, palette = l.a$palette
#     )
#     makeMoviePCA(dir = output, fps = 30)
#     rgl.close()
#     file.rename(
#         paste0(output, '/movie.gif'), paste0(output, '/PCA score plot 3D.gif')
#     )
# }
ModelScorePlot3DNew <- function(model, class.id, output) {
  # 新的3D PCA得分图~~撒花~
  # class.id参数必须是因子，不然会出现图例的排序问题
  # 
  require(plotly)
  require(htmlwidgets)
  df.data <- as.data.frame(model@scoreMN[, 1:3])
  df.data <- cbind(class.id, df.data)
  colnames(df.data) <- c("class.id", "PC1", "PC2", "PC3")
  axis.x <- list(
    range = c(- max(abs(df.data$PC1)), max(abs(df.data$PC1)))
  )
  axis.y <- list(
    range = c(- max(abs(df.data$PC2)), max(abs(df.data$PC2)))
  )
  axis.z <- list(
    range = c(- max(abs(df.data$PC3)), max(abs(df.data$PC3)))
  )
  p <- plot_ly(
    df.data, x = ~PC1, y = ~PC2, z = ~PC3, 
    type = "scatter3d", mode = "markers", 
    color = ~class.id, 
    colors = rep_len(l.a$palette, length.out = length(levels(class.id))), 
    symbol = ~class.id, 
    symbols = rep_len(
      #c("\u2605","\u25C4","\u25BC","\u25B2"),
      c(16, 15, 17, 18),
      length.out = length(levels(class.id))
    ), 
    text = ~row.names(df.data), 
    marker = list(size = 4, opacity = 0.6)
  ) %>% 
    layout(
      legend = list(
        bgcolor = "#E2E2E2",
        bordercolor = "#FFFFFF"
      ), 
      scene = list(
        xaxis = axis.x, yaxis = axis.y, zaxis = axis.z
      )
    )
  old.wd <- getwd()
  setwd(output)
  on.exit(setwd(old.wd))
  htmlwidgets::saveWidget(p, file = "PCA score plot 3D.html")
}
PermutationPlot <- function(model, output) {
  # 置换检验图
  # 
  df.p <- as.data.frame(model@suppLs$permMN)
  df.p <- df.p[df.p$sim > 0,]
  fit.r2 <- lm(
    I(`R2Y(cum)` - df.p$`R2Y(cum)`[1]) ~ I(sim - 1) + 0, data = df.p
  )
  fit.q2 <- lm(
    I(`Q2(cum)` - df.p$`Q2(cum)`[1]) ~ I(sim - 1) + 0, data = df.p
  )
  df.lm <- data.frame(
    x1 = 0, x2 = 1, 
    y1.r = df.p$`R2Y(cum)`[1] - as.numeric(fit.r2$coefficients), 
    y2.r = df.p$`R2Y(cum)`[1], 
    y1.q = df.p$`Q2(cum)`[1] - as.numeric(fit.q2$coefficients), 
    y2.q = df.p$`Q2(cum)`[1]
  )
  df.melt <- reshape2::melt(
    df.p[-1, c('R2Y(cum)', 'Q2(cum)', 'sim')], id.vars = 'sim'
  )
  p <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_segment(
      aes(x = x1, y = y1.r, xend = x2, yend = y2.r), data = df.lm, 
      size = 0.8, color = 'gray30', linetype = 'longdash'
    ) +
    geom_segment(
      aes(x = x1, y = y1.q, xend = x2, yend = y2.q), data = df.lm, 
      size = 0.8, color = 'gray30', linetype = 'longdash'
    ) +
    geom_point(
      data = df.melt, aes(sim, value, shape = variable, color = variable), 
      size = 3
    ) +
    geom_point(
      data = df.p[1, ], aes(sim, `R2Y(cum)`), 
      shape = 16, color = "green", size = 3, show.legend = F
    ) +
    geom_point(
      data = df.p[1, ], aes(sim, `Q2(cum)`), 
      shape = 15, color = "blue", size = 3, show.legend = F
    ) +
    scale_color_manual(values = c("green", "blue")) +
    scale_shape_manual(values = c(16, 15)) +
    theme_bw() +
    theme(
      legend.title = element_blank(), 
      legend.key = element_blank(), 
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = 'black', size = 0.8)
    ) +
    labs(
      title = substitute(
        paste(
          'Intercepts: ', R^2, 'Y(cum) = ', a, ', ', 
          Q^2, '(cum) = ', b, sep = ''
        ), 
        list(
          a = paste0('(0, ', round(df.lm$y1.r, digits = 2), ')'), 
          b = paste0('(0, ', round(df.lm$y1.q, digits = 2), ')')
        )
      ), 
      x = 'Correlation Coefficient', 
      y = substitute(paste(R^2, 'Y(cum) and ', Q^2, '(cum)', sep = ''))
    )
  ggsave(
    paste0(output, '/', model@typeC, ' permutation plot.jpg'), p, 
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/', model@typeC, ' permutation plot.pdf'), p, 
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
}
VolcanoPlot <- function(df.deg, output) {
  # 火山图
  #
  df.deg$deg <- 'not significant'
  df.deg$deg[eval(l.a$deg.exp) & df.deg$`FOLD CHANGE` > 1] <- 'up-regulated'
  df.deg$deg[eval(l.a$deg.exp) & df.deg$`FOLD CHANGE` < 1] <- 'down-regulated'
  df.deg <- df.deg[, c('VIP', 'P-VALUE', 'LOG_FOLDCHANGE', 'deg')]
  colnames(df.deg) <- c('vip', 'p', 'fc', 'deg')
  df.deg$deg <- factor(
    df.deg$deg, 
    levels = c('down-regulated', 'not significant', 'up-regulated')
  )
  p <- ggplot() +
    geom_blank(data = df.deg, aes(-fc, -log10(p), size = vip)) +
    geom_hline(
      yintercept = c(-log10(l.a$deg.p)), 
      color = 'grey50', linetype = 'dashed'
    ) +
    geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
    geom_point(
      data = df.deg[df.deg$deg == "not significant", ], 
      aes(fc, -log10(p), size = vip, color = deg)
    ) +
    geom_point(
      data = df.deg[df.deg$deg != "not significant", ], 
      aes(fc, -log10(p), size = vip, color = deg)
    ) +
    scale_color_manual(
      name = 'Status', values = c(
        'down-regulated' = '#619cffa0', 'not significant' = '#b3b3b350', 
        'up-regulated' = '#f8766da0'
      ), guide = guide_legend(order = 1, override.aes = list(size = 3))
    ) +
    scale_size_continuous(
      name = 'VIP', guide = guide_legend(order = 2), range = c(2, 6), 
      breaks = c(min(df.deg$vip), max(df.deg$vip)), 
      labels = c('0.0', round(max(df.deg$vip), 1))
    ) +
    theme_bw() +
    theme(
      legend.key = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = 'black', size = 1)
    ) +
    labs(
      x = expression(paste(log[2], ' Fold Change')), 
      y = expression(paste(-log[10], ' ', italic('P'), '-value'))
    )
  ggsave(
    paste0(output, '/volcano plot.jpg'), p, 
    width = 10, height = 8, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/volcano plot.pdf'), p, 
    width = 10, height = 8, units = 'in', dpi = 600
  )
}
OutputModelInfo <- function(list.in) {
  # 输出统计模型信息表
  # 
  require(ropls)
  require(openxlsx)
  # 以下为创建模型信息表，按列创建列表对象
  list.out <- list()
  list.out$`Model` <- paste0(  # 第一列
    "Model ", 
    sprintf(paste0("%0", nchar(length(list.in)), "d"), 1:length(list.in))
  )
  list.out$`Type` <- sapply(mvda.model, slot, "typeC")
  {  # A，第二列
    SwitchTemp <- function(m, type) {
      switch(
        type, PCA = m@summaryDF$pre, 
        `PLS-DA` = paste(
          m@summaryDF$pre, m@summaryDF$ort, 0, sep = "+"
        ), 
        `OPLS-DA` = paste(
          m@summaryDF$pre, m@summaryDF$ort, 0, sep = "+"
        )
      )
    }
    list.out$`A` <- mapply(
      SwitchTemp, m = mvda.model, type = list.out$`Type`
    )
  }
  {  # N
    SwitchTemp <- function(m, type) {
      switch(
        type, PCA = m@descriptionMC['samples', 1], 
        `PLS-DA` = m@descriptionMC['samples', 1], 
        `OPLS-DA` = m@descriptionMC['samples', 1]
      )
    }
    list.out$`N` <- mapply(
      SwitchTemp, m = mvda.model, type = list.out$`Type`
    )
  }
  {  # R2X
    SwitchTemp <- function(m, type) {
      switch(
        type, PCA = m@summaryDF$`R2X(cum)`, 
        `PLS-DA` = m@summaryDF$`R2X(cum)`, 
        `OPLS-DA` = m@summaryDF$`R2X(cum)`
      )
    }
    list.out$`R2X(cum)` <- mapply(
      SwitchTemp, m = mvda.model, type = list.out$`Type`
    )
  }
  {  # R2Y
    SwitchTemp <- function(m, type) {
      switch(
        type, PCA = "", 
        `PLS-DA` = m@summaryDF$`R2Y(cum)`, 
        `OPLS-DA` = m@summaryDF$`R2Y(cum)`
      )
    }
    list.out$`R2Y(cum)` <- mapply(
      SwitchTemp, m = mvda.model, type = list.out$`Type`
    )
  }
  {  # Q2
    SwitchTemp <- function(m, type) {
      switch(
        type, PCA = "", 
        `PLS-DA` = m@summaryDF$`Q2(cum)`, 
        `OPLS-DA` = m@summaryDF$`Q2(cum)`
      )
    }
    list.out$`Q2(cum)` <- mapply(
      SwitchTemp, m = mvda.model, type = list.out$`Type`
    )
  }
  list.out$`Title` <- names(mvda.model)  # 最后一列
  df.out <- as.data.frame(list.out, check.names = F)
  {  # 创建模型信息表，写第一个sheet
    hs <- createStyle(
      fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
      valign = 'center', halign = 'center'
    )
    wb <- createWorkbook(creator = 'Sensichip')
    modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
    addWorksheet(
      wb, sheetName = "file explanation", 
      gridLines = T, tabColour = "#4F81BD"
    )
    df.desc <- data.frame(
      "表头名称" = colnames(df.out), 
      "描述" = sapply(colnames(df.out), Explanation2ModelInfo)
    )
    writeDataTable(  # 表格式
      wb, sheet = 1, df.desc, headerStyle = hs, 
      tableStyle = 'TableStyleMedium2', withFilter = T
    )
    setColWidths(wb, sheet = 1, cols = 1, 30)  # 第一列列宽
    setColWidths(wb, sheet = 1, cols = 2, 100)  # 第二列列宽
    setRowHeights(wb, sheet = 1, rows = 1, 24)  # 表头行高
    setRowHeights(wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), 16)  # 表行高
    addStyle(
      wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), cols = 1:ncol(df.desc), 
      gridExpand = T, stack = T, 
      style = createStyle(valign = 'center')
    )
  }
  {  # 写第二个sheet
    addWorksheet(
      wb, sheetName = "Models Information", 
      gridLines = T, tabColour = '#4F81BD'
    )
    writeDataTable(
      wb, sheet = 2, df.out, headerStyle = hs, 
      tableStyle = 'TableStyleMedium2', withFilter = T
    )
    setColWidths(wb, sheet = 2, cols = 1:ncol(df.out), widths = 20)
    setRowHeights(wb, sheet = 2, rows = 1, 24)
    setRowHeights(wb, sheet = 2, rows = 2:(nrow(df.out) + 1), 16)
    addStyle(
      wb, sheet = 2, rows = 2:(nrow(df.out) + 1), cols = 1:ncol(df.out), 
      gridExpand = T, stack = T, 
      style = createStyle(valign = 'center', halign = "left")
    )
  }
  saveWorkbook(
    wb, overwrite = T, 
    file.path(
      l.a$output, 'results', ifelse(is.na(pol), '', paste0(pol, '/')), 
      'Statistical Analysis', 
      paste0(
        ifelse(is.na(pol), "", paste0(pol, "-")), 
        'Models Information.xlsx'
      )
    )
  )
  return(df.out)
}
Explanation2ModelInfo <- function(entry) {
  exp <- switch(
    entry, 
    "Model" = "多元统计分析模型编号", 
    "Type" = "模型类型，包括PCA，OPLS-DA或其他类型",
    "A" = "模型主成分个数，OPLS-DA模型的三个数字分别为预测主成分，正交主成分", 
    "N" = "模型的观测个数（此处即为样本数）", 
    "R2X(cum)" = "代表模型对X变量的累计解释度", 
    "R2Y(cum)" = "代表模型对Y变量的累计解释度", 
    "Q2(cum)" = "模型的可预测性", 
    "Title" = "该模型对应的数据对象"
  )
  return(exp)
}
OutputDegPeak <- function(list.in) {
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    valign = 'center', halign = 'center'
  )
  wb1 <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb1, fontSize = 8, fontName = 'Arial')
  addWorksheet(
    wb1, sheetName = 'file explanation', 
    gridLines = T, tabColour = '#4F81BD'
  )
  wb2 <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb2, fontSize = 8, fontName = 'Arial')
  addWorksheet(
    wb2, sheetName = 'file explanation', 
    gridLines = T, tabColour = '#4F81BD'
  )
  desc1 <- colnames(list.in[[1]])
  desc1[grepl("^MEAN ", desc1)] <- c("MEAN OO", "MEAN XX")
  desc2 <- sapply(desc1, Explanation4Desc)
  df.desc <- data.frame("表头名称" = desc1, "描述" = desc2)
  writeDataTable(
    wb1, sheet = 1, df.desc, 
    headerStyle = hs, tableStyle = 'TableStyleMedium2', withFilter = T
  )
  setColWidths(wb1, sheet = 1, cols = 1, widths = '30')
  setColWidths(wb1, sheet = 1, cols = 2, widths = '180')
  setRowHeights(wb1, sheet = 1, rows = 1, 24)
  setRowHeights(wb1, sheet = 1, rows = 2:(nrow(df.desc) + 1), 16)
  addStyle(
    wb1, sheet = 1, rows = 2:(nrow(df.desc) + 1), cols = 1:nrow(df.desc), 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  filter <- paste0(
    "本次差异代谢物筛选条件为VIP大于", l.a$deg.vip, 
    "且P-value小于", l.a$deg.p
  )
  if(l.a$deg.named) {
    filter <- paste0(filter, "，同时物质有确定的定性名称")
  }
  df.desc <- rbind(df.desc, c('', filter))
  writeDataTable(
    wb2, sheet = 1, df.desc, 
    headerStyle = hs, tableStyle = 'TableStyleMedium2', withFilter = T
  )
  setColWidths(wb2, sheet = 1, cols = 1, widths = '30')
  setColWidths(wb2, sheet = 1, cols = 2, widths = '180')
  setRowHeights(wb1, sheet = 1, rows = 1, 24)
  setRowHeights(wb2, sheet = 1, rows = 2:nrow(df.desc), 16)
  setRowHeights(wb2, sheet = 1, rows = nrow(df.desc) + 1, 40)
  addStyle(
    wb2, sheet = 1, rows = 2:(nrow(df.desc) + 1), cols = 1:nrow(df.desc), 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  addStyle(
    wb2, sheet = 1, rows = nrow(df.desc) + 1, cols = 2, stack = T, 
    style = createStyle(fontSize = 12, fontColour = 'red')
  )
  sapply(
    1:length(l.a$list.compare), function(i) {
      df.deg <- list.in[[i]]
      addWorksheet(
        wb1, sheetName = l.a$subDir[i], 
        gridLines = T, tabColour = '#4F81BD'
      )
      writeDataTable(
        wb1, sheet = i + 1, df.deg, headerStyle = hs, 
        tableStyle = 'TableStyleMedium2', withFilter = T
      )
      setColWidths(wb1, sheet = i + 1, cols = 1:ncol(df.deg), 30)
      setRowHeights(wb1, sheet = i + 1, rows = 1, 24)
      setRowHeights(wb1, sheet = i + 1, rows = 2:(nrow(df.deg) + 1), 16)
      addStyle(
        wb1, sheet = i + 1, rows = 2:(nrow(df.deg) + 1), 
        cols = 1:ncol(df.deg), gridExpand = T, stack = T, 
        style = createStyle(valign = 'center')
      )
      df.deg <- df.deg[eval(l.a$deg.exp), ]
      if (l.a$type == "GC" | l.a$type == "lipid" ) df.deg <- df.deg else 
        df.deg <- subset.data.frame(df.deg,
                                    subset = df.deg$`MS2 name`!="NA" | df.deg$`MS1 name`!="NA",
                                    select = c(colnames(df.deg))
        )
      df.group <- l.a$list.compare[[i]]
      m.data <- cbind(m.final[,l.a$group[df.group[1]][[1]]],m.final[,l.a$group[df.group[2]][[1]]])
      row.names(m.data) <- df.description$id
      df.deg <- cbind(df.deg,m.data[as.character(df.deg$id),])
      addWorksheet(
        wb2, sheetName = l.a$subDir[i], 
        gridLines = T, tabColour = '#4F81BD'
      )
      writeDataTable(
        wb2, sheet = i + 1, df.deg, headerStyle = hs, 
        tableStyle = 'TableStyleMedium2', withFilter = T
      )
      setColWidths(wb2, sheet = i + 1, cols = 1:ncol(df.deg), 30)
      setRowHeights(wb2, sheet = i + 1, rows = 1, 24)
      setRowHeights(wb2, sheet = i + 1, rows = 2:(nrow(df.deg) + 1), 16)
      addStyle(
        wb2, sheet = i + 1, style = createStyle(fontColour = 'red'), 
        cols = 1:ncol(df.deg), gridExpand = T, stack = T, 
        rows = which(df.deg$`FOLD CHANGE` > 1) + 1
      )
      addStyle(
        wb2, sheet = i + 1, style = createStyle(fontColour = 'blue'), 
        cols = 1:ncol(df.deg), gridExpand = T, stack = T, 
        rows = which(df.deg$`FOLD CHANGE` < 1) + 1
      )
      addStyle(
        wb2, sheet = i + 1, rows = 2:(nrow(df.deg) + 1), 
        cols = 1:ncol(df.deg), gridExpand = T, stack = T, 
        style = createStyle(valign = 'center')
      )
    }
  )
  saveWorkbook(
    wb1, overwrite = T, 
    paste0(
      l.a$output, '/results/', 
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'Statistical Analysis/', 
      ifelse(is.na(pol), '', paste0(pol, '-')),
      'Statistical Analysis Results.xlsx'
    )
  )
  saveWorkbook(
    wb2, overwrite = T, 
    paste0(
      l.a$output, 
      '/results/',
      ifelse(is.na(pol), '', paste0(pol, '/')),
      'Statistical Analysis/', 
      ifelse(is.na(pol), '', paste0(pol, '-')), 
      'Differentially Expressed Metabolites.xlsx'
    )
  )
}
SaveWorkspace <- function() {
  file <- paste0(
    l.a$out.temp, "/", ifelse(is.na(pol), "", paste0(pol, ".")), "l.a.RData"
  )
  save.image(file = file)
}
ReportAssistant <- function() {
  # 生成一个帮助写报告的小文档
  # 
  file <- paste0(
    l.a$out.temp, "/", ifelse(is.na(pol), "", paste0(pol, ".")), "l.a.RData"
  )
  load(file)  # 载入存档
  group <- l.a$group[sapply(l.a$group, attr, "type") != "qc"]  # 去掉QC
  group.length <- sapply(group, sum)
  str <- paste0(
    "XXXX样本，", 
    "共", length(group), "组，", 
    "分别为", ShitConverter(names(group), "组，")
  )
  temp.str <- ifelse(
    var(group.length) == 0, 
    yes = paste0("均为", sum(group[[1]]), "例，"), 
    no = paste0(
      "分别为", ShitConverter(group.length, "例，")
    )
  )
  str <- paste0(
    str, "各组对应的生物学重复", temp.str, "进行基于", 
    switch(l.a$type, GC = "GC-TOFMS", LC = "LC-QTOFMS", QE = "QE"), 
    "的代谢组学分析，共计", 
    sum(apply(data.frame(group), 1, any)), "例样本。"
  )
  output <- paste0(
    l.a$out.temp, "/", ifelse(is.na(pol), "", paste0(pol, "-")), 
    "ReportAssist.txt"
  )
  write(str, output)
  str <- paste0(  # 样品信息表
    "\n详细分组\t", 
    sub("\t$", "\n", ShitConverter(names(group), "组\t")), 
    "样本数量\t", 
    sub("\t$", "\n", ShitConverter(sapply(group, sum), "例\t")), 
    collapse = ""
  )
  write(str, output, append = T)
  temp.str <- sapply(
    l.a$list.compare, function(pair) {
      paste0(names(group)[pair[1]], "组对", names(group)[pair[2]], "组")
    }
  )
  str <- paste0(
    "\n共", length(l.a$list.compare), "组对比，", 
    ifelse(length(l.a$list.compare) == 1, "即", "分别为"), 
    sub("，$", "。", ShitConverter(temp.str, "，"))
  )
  write(str, output, append = T)
}
ShitConverter <- function(vec, appendix) {
  # 在vec的每一个元素后面加上appendix，
  # 合并为一个字符串，并用“和”替换掉最后一个逗号
  # 
  temp.str <- paste(c(vec, ""), collapse = appendix)
  temp.str <- sub("^(.*)，(.+?)$", "\\1和\\2", temp.str)
  return(temp.str)
}
# 层次聚类分析 ----
Heatmap <- function() {
  # 热图函数，大修了一波
  # 
  require(dplyr)
  require(pheatmap)
  require(openxlsx)
  wb <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    valign = 'center', halign = 'center'
  )
  df.final <- cbind(df.description, m.final)
  if (l.a$type == "GC") df.final <- df.final[order(-df.final$Similarity), ]
  lapply(
    1:length(l.a$list.compare), function(i) {
      df.deg <- list.compare.result[[i]]
      pair <- l.a$list.compare[[i]]
      df.meta <- dplyr::left_join(
        df.final, df.deg[, grepl("^id$|^MEAN", names(df.deg))], 
        by = "id"
      )
      names(df.meta)[
        grepl("Peak|MS2.name|(MS2 name)|(compound name)|(lipid name)", names(df.meta))
      ] <- "met"
      names(df.meta)[
        grepl("Similarity|MS2.score|(MS2 score)|^score$", names(df.meta))
      ] <- "score"
      df.meta <- df.meta[eval(l.a$deg.exp), ]
      #df.meta <- df.meta[order(-df.meta$score), ]
      if(l.a$heatmap.ignore) {  # 根据是否去重，处理样本名称
        df.meta <- HeatmapFilter(df.meta)
      } else {
        dup <- duplicated(df.meta$met) | 
          duplicated(df.meta$met, fromLast = T)
        df.meta$met[dup] <- paste(
          df.meta$met, df.meta$id
        )[dup]
      }
      if(is.null(df.meta)) {
        df.data <- data.frame(X1 = "No Differentially Expressed Metabolites")
      } else {
        list.data <- lapply(
          pair, function(j) {
            df.sub <- df.meta[, l.a$desc + which(l.a$group[[j]])]
            if(nrow(df.sub) == 1) {
              return(df.sub)
            } else {
              x <- pheatmap::pheatmap(
                df.sub, scale = "row", silent = T
              )$tree_col$order
              return(df.sub[, x])
            }
          }
        )
        df.data <- do.call(cbind, list.data)
        df.data <- data.frame(
          df.meta[, grep("^id+|^met+", colnames(df.meta)), drop = F], 
          df.data
        )
        if ("met" %in% colnames(df.data)) {
          rownames(df.data) <- df.data$met
          df.heat <- df.data[, 3:ncol(df.data)]
        } else {
          rownames(df.data) <- df.data$id
          df.heat <- df.data[, -1]
        }
        label.col <- c(
          rep(names(l.a$group)[pair], sapply(l.a$group, sum)[pair])
        )
        fn.jpg <- paste0(
          l.a$output, '/results/', 
          ifelse(is.na(pol), '', paste0(pol, '/')), 
          'Hierarchical Clustering Analysis/',
          l.a$subDir[i], '/heatmap.jpg'
        )
        fn.pdf <- paste0(
          l.a$output, '/results/', 
          ifelse(is.na(pol), '', paste0(pol, '/')), 
          'Hierarchical Clustering Analysis/', 
          l.a$subDir[i], '/heatmap.pdf'
        )
        cluster.rows <- ifelse(nrow(df.data) == 1, F, T)
        pheatmap::pheatmap(
          df.data[, 3:ncol(df.data)], 
          cluster_cols = F, cluster_rows = cluster.rows, 
          color = colorRampPalette(c("blue","white", "red"))(1000),
          scale = 'row', fontsize = 5, 
          cellwidth = 12, cellheight = 24, border_color ="gray", 
          labels_col = label.col, filename = fn.jpg
        )
        hm <- pheatmap::pheatmap(
          df.data[, 3:ncol(df.data)], 
          cluster_cols = F, cluster_rows = cluster.rows,
          color = colorRampPalette(c("blue","white", "red"))(1000),
          scale = 'row', fontsize = 5,
          cellwidth = 12, cellheight = 24, border_color = "gray",
          labels_col = label.col, filename = fn.pdf
        )
        if(!all(is.na(hm$tree_row))) df.data <- df.data[hm$tree_row$order, ]
      }
      addWorksheet(
        wb, sheetName = l.a$subDir[i],
        gridLines = T, tabColour = '#4F81BD'
      )
      writeDataTable(
        wb, sheet = i, df.data, headerStyle = hs,
        tableStyle = 'TableStyleMedium2', withFilter = T
      )
      setColWidths(wb, sheet = i, cols = 1:ncol(df.data), 'auto')
      setRowHeights(wb, sheet = i, rows = 1, 24)
      setRowHeights(wb, sheet = i, rows = 2:(nrow(df.data) + 1), 16)
      addStyle(
        wb, sheet = i, rows = 2:(nrow(df.data) + 1),
        cols = 1:ncol(df.data), gridExpand = T, stack = T,
        style = createStyle(valign = 'center')
      )
    }
  )
  saveWorkbook(
    wb, overwrite = T,
    paste0(
      l.a$output, '/results/', 
      ifelse(is.na(pol), '', paste0(pol, '/')), 
      'Hierarchical Clustering Analysis/', 
      ifelse(is.na(pol), '', paste0(pol, '-')), 
      'hierarchical clustering data matrix.xlsx'
    )
  )
}
HeatmapFilter <- function(df.meta) {
  df.meta <- df.meta[
    !is.na(df.meta$met) & !grepl("Analyte|unknown", df.meta$met), 
  ]
  list.keep <- lapply(
    unique(df.meta$met), function(met) {
      df.sub <- df.meta[df.meta$met == met, ]
      if(nrow(df.sub) > 1) {
        df.sub <- df.sub[df.sub$score == max(df.sub$score), ]
      }
      if(nrow(df.sub) > 1) {
        mean.area <- apply(
          df.sub[, grepl("^MEAN", names(df.sub))], 1, mean
        )
        df.sub <- df.sub[which.max(mean.area), ]
      }
      #browser()
      if(nrow(df.sub) != 1) {
        stop("Wrong dedup procedure @ HeatmapFilter")
      }
      return(df.sub)
    }
  )
  df.meta <- do.call(rbind, list.keep)
  return(df.meta)
}

# 脂质气泡图绘制 ----
lipidBubble <- function(dataFrame = list.compare.result, taskLine = l.a$subDir, outputBefore = l.a$output, limit = 500) {
  if (l.a$type != "lipid") {
    stop("该选项仅对LC脂质组有效！")
  }
  require(ggplot2)
  # 装备存储容器
  outputAll <- paste0(outputBefore, "/results/", pol, "/Bubble Plot/")
  wb <- createWorkbook(creator = "Sensichip")
  modifyBaseFont(wb, fontSize = 10, fontName = "Arial")
  hs <- createStyle(
    fontSize = 12, fontColour = "white", fgFill = "#4F81BD",
    valign = "center", halign = "center"
  )
  # 分析、分装
  for (i in 1:length(taskLine)) {
    # 整理数据
    data4draw <- dataFrame[[i]]
    data4draw <- data4draw[, c("id", "lipid name", "P-VALUE", "FOLD CHANGE")]
    data4draw <- na.omit(data4draw)
    data4draw$`Lipid species` <- gsub("\\(.*\\)$", "", data4draw$`lipid name`)
    data4draw$`Relative diffrence(%)` <- (data4draw$`FOLD CHANGE` - 1) * 100
    dataOut <- data4draw[, c("id", "lipid name", "Lipid species", "P-VALUE", "Relative diffrence(%)")]
    # 单写csv 没有样式
    # write.csv(dataOut,paste0(outputAll,taskLine[i],'/Bubble plot.csv'),row.names = F)
    # 存储整理好的数据
    addWorksheet(
      wb,
      sheetName = taskLine[i],
      gridLines = T, tabColour = "#4F81BD"
    )
    writeDataTable(
      wb,
      sheet = i, dataOut, headerStyle = hs,
      tableStyle = "TableStyleMedium2", withFilter = T
    )
    setColWidths(wb, sheet = i, cols = 1, widths = "10")
    setColWidths(wb, sheet = i, cols = 2, widths = "30")
    setColWidths(wb, sheet = i, cols = 3, widths = "30")
    setColWidths(wb, sheet = i, cols = 4, widths = "30")
    setColWidths(wb, sheet = i, cols = 5, widths = "30")
    setRowHeights(wb, sheet = i, rows = 1, 28)
    setRowHeights(wb, sheet = i, rows = 2:(nrow(dataOut) + 1), 20)
    addStyle(
      wb,
      sheet = i,
      style = createStyle(
        halign = "left", valign = "center", wrapText = T
      ), gridExpand = T, stack = T,
      cols = 1:ncol(dataOut),
      rows = 1:(nrow(dataOut) + 1)
    )
    # 绘制图形、存储图形
    data4draw$`label` <- paste0(as.integer(data4draw$`Relative diffrence(%)`), "%")
    data4draw$`Relative diffrence(%)`[data4draw$`Relative diffrence(%)` < -limit] <- (data4draw$`Relative diffrence(%)`[data4draw$`Relative diffrence(%)` < -limit] + limit) / 10 - limit
    data4draw$`Relative diffrence(%)`[data4draw$`Relative diffrence(%)` > limit] <- (data4draw$`Relative diffrence(%)`[data4draw$`Relative diffrence(%)` > limit] - limit) / 10 + limit
    data4draw$`Singificance` <- -log10(data4draw$`P-VALUE`)
    p <- ggplot(data = data4draw, aes(x = `Relative diffrence(%)`, y = `Lipid species`)) +
      geom_point(data = data4draw[data4draw$`P-VALUE` >= 0.05, ], aes(size = `Singificance`, color = `Lipid species`), pch = 19, color = "grey50", alpha = 0.7) +
      geom_point(data = data4draw[data4draw$`P-VALUE` < 0.05, ], aes(size = `Singificance`, color = `Lipid species`), pch = 19, alpha = 0.7) +
      geom_vline(xintercept = 0, lty = 2, color = "grey50") +
      # geom_text(data = data4draw[abs(data4draw$`Relative diffrence(%)`) > limit, ], aes(label = `label`)) +
      geom_rug(sides = "b", alpha = 0.7) +
      scale_size_continuous(range = c(1, 10), breaks = c(min(data4draw$`Singificance`), max(data4draw$`Singificance`))) +
      labs(
        x = expression("Relative diffrence(%)"),
        y = expression("Lipid species")
      ) +
      # facet_wrap(~face,nrow=1,scales="free_x")+
      scale_x_continuous(breaks = seq(-limit, limit, 100)) +
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(colour = "#000000", size = 15),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA)
      )
    ggsave(
      paste0(outputAll, taskLine[i], "/Bubble plot.jpg"), p,
      width = 12, height = 7.5, units = "in", dpi = 600
    )
    ggsave(
      paste0(outputAll, taskLine[i], "/Bubble plot.pdf"), p,
      width = 12, height = 7.5, units = "in", dpi = 600
    )
  }
  # 盖上盖子
  saveWorkbook(wb,
               overwrite = T,
               paste0(outputAll, "Bubble plot.xlsx")
  )
}
# KEGG分析 ----
compareMeetMap <- function(list.compare.result = list.compare.result,
                           df.map = df.map
){
  # mainly inherit from old script
  # accept the compare list and map dataframe
  # make a list comtaining a raw kegg in of each compare
  # the raw kegg in has MS2.name, LogFc, deg(boolean) and KEGG Id columns
  exp.cpd <- expression(
    if(l.a$type == "GC") {
      !grepl("Analyte|unknown", df.map$Peak)
    } else if(l.a$type == "LC") {
      !is.na(df.map$`MS2 name`)
    }
  )
  
  cpd.col <- switch(l.a$type, GC = "Peak", LC = "MS2 name")
  list.met <- lapply(  # 读取所有对比的代谢物信息
    list.compare.result, function(df.deg) {
      df.deg$deg <- 0
      df.deg$deg[eval(l.a$deg.exp)] <- 1
      df.deg <- df.deg[, c(cpd.col, 'LOG_FOLDCHANGE', "deg")]
      df.deg <- df.deg[eval(exp.cpd), ]
      df.deg <- left_join(
        df.deg, df.map[, c(cpd.col, "KEGG COMPOUND ID")], 
        by = cpd.col
      )
      colnames(df.deg) <- c(cpd.col, "LOG_FOLDCHANGE", "deg", "KEGG")
      df.deg <- na.omit(df.deg)
      return(df.deg)
    }
  )
  
  return(list.met)
}

KEGGColPathway <- function(list.met = compareMeetMap(list.compare.result, df.map),
                           keggInFileName = l.a$subDir,
                           keggOrg = l.a$keggOrg,
                           keggDatabaseDir = "D:/share/KEGGDB/",
                           keggInDir = file.path(
                             l.a$out.temp, 
                             paste0("KEGG.in", ifelse(is.na(pol), "", pol))),
                           keggInExcelInputFile = file.path(
                             l.a$out.temp, 
                             paste0("KEGG.in", ifelse(is.na(pol), "", pol)),
                             "keggInForLocal.xlsx"),
                           keggPicOutputDir = file.path(
                             l.a$output, 
                             "results", 
                             ifelse(is.na(pol), "", pol), 
                             "KEGG Analysis"),
                           minHit = l.a$keggHitMin
                           
){
  # mainly inherit from old script
  # one function do all the tricks involve kegg
  # input requires the compare and map result, list.met
  # output the match pics and files
  
  list.kegg <- base::mapply(
    keggPartSingle,
    df.in = list.met,
    keggInFileName = keggInFileName,
    MoreArgs = list(
      keggInDir = keggInDir,
      keggOrg = keggOrg,
      keggDatabaseDir = keggDatabaseDir,
      minHit = minHit
    ), 
    SIMPLIFY = F
  )
  #
  PathwayList(list.kegg)  # 输出通路表
  # write kegg in excel
  writeKeggInputExcel(inputDir = keggInDir, outputDir = keggInDir)
  # draw kegg pics
  drawKeggPics(org = keggOrg,
               inputFile = keggInExcelInputFile,
               keggDatabaseDir = keggDatabaseDir,
               minHit = minHit)
}

keggPartSingle <- function(df.in,
                           keggInFileName,
                           keggInDir,
                           keggOrg,
                           keggDatabaseDir,
                           minHit
){
  # use to accept multiple arguments and pass to sub function
  
  # get sig and all kegg id and logFc
  sigAndAllRawKeggIn <- getSigAndAllRawKeggIn(df.in) 
  # write kegg in txt
  # writeKeggIn(sigAndAllRawKeggIn$sig, 
  #             outputDir = keggInDir, 
  #             fileName = keggInFileName)
  # get sig and all pathway
  df.out.keggp2 <- getSigAndAllKeggPathway(df.sig = sigAndAllRawKeggIn$sig,
                                           df.all = sigAndAllRawKeggIn$all,
                                           org = keggOrg,
                                           keggDatabaseDir = keggDatabaseDir,
                                           minHit = minHit,
                                           comparison = keggInFileName)
  return(df.out.keggp2)
}



getSigAndAllRawKeggIn <- function(df.in){
  # Mainly inherit from old script
  # Accept the raw kegg in
  # Separate and filter the sig and all kegg id with logFc
  # Return the a list contains both
  # Each is formatted as only one logFc column with kegg id as rownames  
  df.in.l <- df.in[, c("KEGG", "LOG_FOLDCHANGE", "deg")]
  colnames(df.in.l) <- c("kegg", "fc", "deg")
  df.in.l <- na.omit(df.in.l)
  uniqueKeggId <- unique(df.in.l$kegg)
  uniqueKeggId <- uniqueKeggId[uniqueKeggId != ""] # remove empty elements
  uniqueKeggId <- uniqueKeggId[uniqueKeggId != " "] # remove empty elements
  list.in <- lapply(
    uniqueKeggId, function(kegg) {
      # 对于多个代谢物对应一个KEGG ID的情况，并不能无脑取第一个
      # 多个代谢物中同时存在差异和非差异代谢物的情况，优先差异代谢物
      # 
      df.sub <- df.in.l[df.in.l$kegg == kegg, ]
      df.sub <- df.sub[order(df.sub$deg, decreasing = T), ]
      df.sub[1, ]
    }
  )
  df.in.l <- do.call(rbind, list.in)
  rownames(df.in.l) <- df.in.l$kegg
  df.all <- df.in.l[, c(-1, -ncol(df.in.l)), drop = F]
  df.sig <- df.in.l[df.in.l$deg == 1, c(-1, -ncol(df.in.l)), drop = F]
  
  sigAndAllRawKeggIn <- list(sig = df.sig, all = df.all)
  return(sigAndAllRawKeggIn)
}

writeKeggIn <- function(rawKeggIn,
                        outputDir = "."                        ,
                        fileName = "keggIn"
){
  # mainly inherit from old script
  # accept raw kegg in which has one logFc column and kegg id as rownames
  # turn logFc to red or blue as up or download regulated symbols
  # write the kegg in files
  
  if (!dir.exists(outputDir)) dir.create(outputDir, recursive = T)
  df.out <- as.numeric(as.matrix(rawKeggIn))
  df.out <- ifelse(df.out > 0, "red", "blue")
  attributes(df.out) <- attributes(as.matrix(rawKeggIn))
  write.table(
    df.out, 
    file = file.path(
      outputDir, paste0(fileName,".txt")
    ), 
    col.names = F, quote = F, sep = "\t"
  )
}


getSigAndAllKeggPathway <- function(df.sig, 
                                    df.all,
                                    org = l.a$keggOrg,
                                    keggDatabaseDir = "D:/share/KEGGDB/",
                                    minHit = l.a$keggHitMin,
                                    comparison = "This comparison",
                                    #type = l.a$type
                                    type = "QE"){
  # mainly inherit from old script
  # accept the df.sig and df.all
  # get the kegg id of each, fetch kegg pathway
  # return a list contains kegg pathway
  if(type == "QE"){
    df.null <- data.frame(
      "KEGG pathway" = character(), "Compound" = character(), 
      stringsAsFactors = F, check.names = F
    )
    if(length(df.sig) == 0){
      print(paste0(
        comparison,
        " has zero DEM found. Please change the deg.p or the keggHitMin."))
      return(df.null)
    }
    
    {  # 确定通路
      
      load(file = paste0(keggDatabaseDir, org,"_KGML/listCg.RData"))
      
      pathwayCompoundMatchSig <- lapply(listCg, function(x){
        query <- df.sig
        query[query %in% x]
      })
      load(file="D:/share/KEGGDB/keggCpdInfo.RData")
      df.ptw <- data.frame(
        pathwayId = names(pathwayCompoundMatchSig),
        Freq = sapply(pathwayCompoundMatchSig, length),
        cpd = sapply(pathwayCompoundMatchSig, function(x) {
          z <- c()
          for(i in x){
            y <- grep(i,keggCpdInfo,value = T)
            z <- c(z,y)
          }
          paste(z, collapse = "; ")
        })
      )
      df.ptw <- df.ptw[df.ptw$Freq >= minHit, ]
      df.ptw <- df.ptw[base::order(df.ptw$Freq, decreasing = T), ]
      
      load(paste0(keggDatabaseDir, org, "_KGML/pathwaysInfo.RData"))
      
      keggOrg.path <- pathwaysInfo$pathAnnotation
      names(keggOrg.path) <- pathwaysInfo$pathNo # TODO sub no need
      
      if(nrow(df.ptw) == 0){
        print(paste0(
          comparison,
          " has zero kegg pathway found. Please change the deg.p or the keggHitMin."))
      }
      df.out.keggp2 <- data.frame(
        "KEGG pathway" = paste0(df.ptw$pathwayId," ",keggOrg.path[df.ptw$pathwayId],"(",df.ptw$Freq,")"),
        "Compound" = df.ptw$cpd, 
        stringsAsFactors = F, check.names = F
      )
      
    }
    return(df.out.keggp2)
  }else{
    if(nrow(df.sig) == 0){
      print(paste0(
        comparison,
        " has zero DEM found. Please change the deg.p or the keggHitMin."))
    }
    
    {  # 确定通路
      
      load(file = paste0(keggDatabaseDir, org,"_KGML/listCg.RData"))
      
      pathwayCompoundMatchSig <- lapply(listCg, function(x){
        query <- rownames(df.sig)
        query[query %in% x]
      })
      
      pathwayCompoundMatchAll <- lapply(listCg, function(x){
        query <- rownames(df.all)
        query[query %in% x]
      })
      
      df.ptw <- data.frame(pathwayId = names(pathwayCompoundMatchSig), 
                           Freq = sapply(pathwayCompoundMatchSig, length))
      df.ptw <- df.ptw[df.ptw$Freq >= minHit, ]
      df.ptw <- df.ptw[base::order(df.ptw$Freq, decreasing = T), ]
      
      load(paste0(keggDatabaseDir, org, "_KGML/pathwaysInfo.RData"))
      keggOrg.path <- pathwaysInfo$pathAnnotation
      names(keggOrg.path) <- pathwaysInfo$pathNo # TODO sub no need
      
      if(nrow(df.ptw) == 0){
        print(paste0(
          comparison,
          " has zero kegg pathway found. Please change the deg.p or the keggHitMin."))
      }
      
      df.out.keggp2 <- data.frame(
        "pathway" = df.ptw$pathwayId,
        "description" = keggOrg.path[df.ptw$pathwayId], 
        stringsAsFactors = F, check.names = F
      )
      
      list.temp <- lapply(
        df.ptw$pathwayId, function(pathwayId){
          cpdsSig <- unlist(pathwayCompoundMatchSig[pathwayId])
          cpdsAll <- unlist(pathwayCompoundMatchAll[pathwayId])
          return(
            list(
              length(cpdsSig), base::paste(cpdsSig, collapse = ";"),
              length(cpdsAll), base::paste(cpdsAll, collapse = ";")
            )
          )
        }
      )
      
      df.out.keggp2$`# compounds (dem)` <- sapply(list.temp, `[[`, 1)
      df.out.keggp2$`compounds (dem)`<- sapply(list.temp, `[[`, 2)
      df.out.keggp2$`# compounds (all)` <- sapply(list.temp, `[[`, 3)
      df.out.keggp2$`compounds (all)`<- sapply(list.temp, `[[`, 4)
    }
    return(df.out.keggp2)
  }
  
}

writeKeggInputExcel <- function(inputDir = "./",
                                outputDir = "./")
{
  # read all kegg in files
  # write an excel file fits for kegg local script input
  # the excel file has sheets contain pathway name and color
  # the sheets name should be compares name
  
  keggInFiles <- list.files(inputDir, pattern = "\\.txt$")
  
  wb <- createWorkbook(creator = 'Sensichip')
  for (compare in keggInFiles){
    addWorksheet(
      wb, sheetName = gsub(".txt$", "", compare)
    )
    keggIn <-  read.table(paste0(inputDir, '/', compare),
                          sep = "\t", header = FALSE)
    colnames(keggIn) <- c("KEGG ID", "Colour")
    writeDataTable(
      wb, sheet = gsub(".txt$", "", compare), 
      keggIn
    )
  }
  
  saveWorkbook(
    wb, overwrite = T, 
    file.path(
      outputDir, "/keggInForLocal.xlsx"
    )
  )
}



checkKeggOrgExistence <- function(org, dir = "D:/share/KEGGDB/"){
  # check if certain org kegg database exists
  if(dir.exists(paste0(dir, org, "_KGML"))){
    print(paste0(org, " database exists"))
    return(TRUE)
  } else {
    print(paste0(org, " database doesn't exist"))
    print("Plase use downloadKeggDatbase() to download")
    return(FALSE)
  }
  
}

downloadKeggDatabase <- function(org, 
                                 dir = "D:/share/KEGGDB/",
                                 downloadPath = T,
                                 downloadGene = F){
  #
  #
  ncharOrg <- nchar(org)
  url <- "http://rest.kegg.jp/"# no need to change at most time
  dir.create(paste0(dir, org,"_KGML"), showWarnings = FALSE)
  
  # function: to get all the pathways of one org 
  if (downloadPath){
    pathwayList <- function(org){
      resp <- GET(paste0(url,"list/pathway/",org))
      Raw <- strsplit(httr::content(resp,as = "text"),"\n")
      pathways <- unlist(lapply(Raw, function(x){strsplit(x, "\t")}))
      pathNo <- pathways[1:length(pathways) %% 2 == 1]
      pathAnnotation <- pathways[1:length(pathways) %% 2 == 0]
      pathwaysInfo <- list(pathNo = pathNo, pathAnnotation = pathAnnotation)
      return(pathwaysInfo)
    }
    pathwaysInfo <- pathwayList(org)
    pathwaysInfo$pathNo <- sub("path:", "", pathwaysInfo$pathNo)
    save(pathwaysInfo, file = paste0(dir, org, "_KGML/pathwaysInfo.RData"), compress = TRUE)
    
    listPath <- pathwaysInfo$pathNo
    print(paste("Total", length(listPath), "pathways found"))
    
    # function: to get all KGMLs of one org
    downloadKGML <- function(Fx){
      fileDir <- paste0(dir, org,"_KGML/",Fx,".xml")
      if (!file.exists(fileDir)){
        resp <- GET(paste0(url,"get/",Fx,"/kgml"))
        html <- httr::content(resp, as = "parsed",encoding = "UTF-8")
        write_xml(html,fileDir)
        print(paste(fileDir, "downloaded"))
      } else {
        print(paste(fileDir, "exists. Download skipped."))
      }
    }
    lapply(listPath,downloadKGML)
    # function: to get all maps of one org
    downloadMap <- function(Fx){
      fileDir <- paste0(dir, org,"_KGML/",Fx,"_map.RData")
      if (!file.exists(fileDir)){
        imgHtml <- GET(paste0(url,"get/",Fx,"/image"))
        imgArray <- httr::content(imgHtml, as = "parsed")
        save(imgArray, file = fileDir , compress = TRUE)
        print(paste(fileDir, "downloaded"))
      } else {
        print(paste(fileDir, "exists. Download skipped."))
      }
    }
    lapply(listPath,downloadMap)
    
    # function: to get cpd info
    downloadCpdInfo <- function(){
      resp <- GET(paste0("http://rest.kegg.jp/list/cpd"))
      Raw <- strsplit(httr::content(resp,as = "text"),"\n")
      rawV <- unlist(Raw)
      rawV <- gsub("\t"," ",rawV)
      rawF <- gsub(";.*","",rawV)
      keggCpdInfo <- rawF
      save(keggCpdInfo,file="D:/share/KEGGDB/keggCpdInfo.RData")
    }
    
    listCgDir <- paste0(dir, org,"_KGML/listCg.RData")
    if (!file.exists(listCgDir)){
      # function: to get all the genes of one pathway
      pathwayGene <- function(path){
        resp <- GET(paste0(url,"link/",org,"/",path))
        Raw <- strsplit(httr::content(resp,as = "text"),"\n")
        geneNo <- lapply(Raw,function(x){substr(x,(13+ncharOrg*2),nchar(x))})
        return(unlist(geneNo))
      }
      # function: to get all the cpds of one pathway
      pathwayCpd <- function(path){
        resp <- GET(paste0(url,"link/cpd/map",substr(path,(ncharOrg+1),(ncharOrg+5))))
        Raw <- strsplit(httr::content(resp,as = "text"),"\n")
        cpdNo <- lapply(Raw,function(x){substr(x,19,nchar(x))})
        return(unlist(cpdNo))
      }
      # list the genes and cpds of each pathway
      listGene <- lapply(listPath,pathwayGene)
      names(listGene) <- listPath
      listCpd <- lapply(listPath,pathwayCpd)
      names(listCpd) <- listPath
      listCg <- lapply(listPath,function(x){c(listGene[[x]],listCpd[[x]])})
      names(listCg) <- listPath
      save(listCg,file =listCgDir, compress = TRUE)
      print(paste(listCgDir, "downloaded"))
    } else {
      print(paste(listCgDir, "exists. Download skipped."))
    }
  }
  
  
  # gene name and entry
  if (downloadGene){
    genesDir <- paste0(dir, org, "_KGML/genesInfo.RData")
    if (!file.exists(genesDir)){
      getGeneList <- function(org){
        resp <- GET(paste0(url,"list/",org)) # download a specie whole genes html
        Raw <- strsplit(httr::content(resp,as = "text"),"\n") # get the text
        genes <- unlist(lapply(Raw, function(x){strsplit(x, "\t")})) # seperate the entry and annotation
        
        geneNum <- genes[1:length(pathways) %% 2 == 1] 
        geneNumPure <- unlist(strsplit(geneNum, ":"))[1:length(geneNum) %%2 == 0] 
        
        geneAnnotation <- genes[1:length(pathways) %% 2 == 0]
        geneAnnotationList <- strsplit(geneAnnotation, "; ")
        geneName <- unlist(lapply(geneAnnotationList, function(x){x <- x[1]}))
        geneNames <- strsplit(geneName, ", ")
        # geneDescription <- unlist(lapply(geneAnnotationList, 
        #                                  function(x){
        #                                    ifelse(length(x) == 1,"", x[2])})) # needless for now
        genesInfo <- geneNames
        names(genesInfo) <- geneNumPure
        return(genesInfo)
      }
      genesInfo <- getGeneList(org)
      save(pathwaysInfo, file = paste0(dir, org, "_KGML/genesInfo.RData"), compress = TRUE)
      print(paste(genesDir, "downloaded"))
    } else{
      print(paste(genesDir, "exists. Download skipped."))
    }
  }
  
  print(paste(org, "database download completed"))
}

drawKeggPics <- function(org, 
                         inputFile = "keggInForLocal.xlsx", 
                         keggDatabaseDir = "D:/share/KEGGDB/",
                         outputDir = file.path(l.a$output, "results", ifelse(is.na(pol), "", pol), 
                                               "KEGG Analysis"),
                         minHit = l.a$keggHitMin){
  # load listCg
  load(file = paste0(keggDatabaseDir, org,"_KGML/listCg.RData"))
  # function to draw circle on the map
  polygonCircle <- function(x,y,r,bg,fg= "black",lwd = 1){
    Cx <- x+sin(seq(0,2*pi,by = 0.01))*(r+0.5)
    Cy <- y+cos(seq(0,2*pi,by = 0.01))*(r+0.5)
    polygon(Cx,Cy,col = bg,border = fg,lwd = lwd,density = 1000)
  }
  # above this is ok to run once ,if u have not change org 
  aVsB <- getSheetNames(inputFile)
  for(i in aVsB) {
    dir.create(paste0(outputDir,'/',i), showWarnings = FALSE)
    #browser()
    input <-  read.xlsx(inputFile,
                        sheet = i,
                        colNames = T,
                        rowNames = F)
    queryColor <- as.vector(paste(input[, 1], input[, 2], sep = "%09"))
    query <- input[, 1]
    pathHitR <- lapply(listCg, function(x) {
      queryColor[query %in% x]
    })
    pathHit <- subset(pathHitR,subset = grepl(switch(minHit, "1" ="%09","2" = "%09.*%09"),pathHitR)) # only run map which has more than 2 DEC|DEG
    
    
    #rawUrls <- lapply(pathHit, pasteNv)
    #rUrls <- unlist(rawUrls)
    # colorUrls <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?",
    #                     names(rUrls),
    #                     rUrls)
    #write.csv(colorUrls,file = paste0(i,"/","Urls.csv"))
    if (length(pathHit) != 0){
      for (k in 1:length(pathHit)) {
        j <- pathHit[k] # easy for test
        map <- names(j)
        load(file = paste0(keggDatabaseDir, org,"_KGML/",map,"_map.RData"))
        html <- read_xml(paste0(keggDatabaseDir, org,"_KGML/",map,".xml"))
        entrys <- xml_find_all(html,".//entry")
        
        
        
        a <- dim(imgArray)[[2]]
        b <- dim(imgArray)[[1]]
        CairoPNG(paste0(outputDir,"/", i,"/",map,".png"),
                 width = 2*a, height = 2*b, bg = "transparent",quality = 100)
        #pdf(paste0(i,"/",map,".pdf"),width = a/100, height = b/100, bg = "transparent") 
        #if u need high resolution image ,use pdf
        par(mai = c(0,0,0,0))
        plot(
          0,
          xlim = c(0,a),
          ylim = c(0,b),
          asp = 1,
          type = 'n',
          xlab = NA,
          ylab = NA,
          frame.plot = F,
          axes = F
        )
        rasterImage(imgArray, 0, 0, a, b,bg = "transparent")
        
        addColor <- function(Fx){
          type <- xml_attr(Fx,"type")
          if(type == "rectangle") {
            x<- as.numeric(xml_attr(Fx,"x"))
            y<- as.numeric(xml_attr(Fx,"y"))
            width<-as.numeric(xml_attr(Fx,"width"))
            height<-as.numeric(xml_attr(Fx,"height"))
            polygon(
              c(x - width / 2, x + width / 2, x + width / 2, x - width / 2),
              c(b - y - height / 2, b - y - height / 2, b - y + height /
                  2, b - y + height / 2),
              col = Fcol,
              border = Fcol,
              lwd = 1,
              density = 5
            )
          } else if (type == "circle") {
            x<- as.numeric(xml_attr(Fx,"x"))
            y<- as.numeric(xml_attr(Fx,"y"))
            width<-as.numeric(xml_attr(Fx,"width"))
            height<-as.numeric(xml_attr(Fx,"height"))
            polygonCircle(x,b-y,max(width,height)/2,bg = Fcol,lwd = 1) 
          } else if (type == "line") {
            coords <- as.numeric(unlist(strsplit(xml_attr(Fx,"coords"),",")))
            xCoords <- coords[seq(1,length(coords),2)]
            yCoords <- coords[seq(2,length(coords),2)]
            lines(xCoords,b-yCoords,col = Fcol,lwd = 2)
          }
        }
        
        addColors <- function(Fx){
          name <- unlist(strsplit(Fx,"%09"))[1]
          Fcol <<- unlist(strsplit(Fx,"%09"))[2]
          entry <- subset(entrys,grepl(paste0(':',name,'$'),xml_attr(entrys,"name"))|grepl(paste0(':',name,''),xml_attr(entrys,"name")))
          graphic <- xml_find_all(entry,".//graphics")
          lapply(graphic,addColor)
        }
        k <- unlist(j)
        lapply(k,addColors)
        
        dev.off()
      }
      print(paste0(i, " kegg pics done"))
    } else {
      print(paste0(i, " kegg pics passed"))
    }
  }
}

# search and color pathway
# function to paste elements in vector 
pasteNv <- function(x,y = "",sep = "/"){
  y <- names(x)
  for(i in x){y <- paste(y,i,sep = sep)}
  return(y)
}
# function to draw circle on the map
polygonCircle <- function(x,y,r,bg,fg= "black",lwd = 1){
  Cx <- x+sin(seq(0,2*pi,by = 0.01))*(r+0.5)
  Cy <- y+cos(seq(0,2*pi,by = 0.01))*(r+0.5)
  polygon(Cx,Cy,col = bg,border = fg,lwd = lwd,density = 1000)
}

MakeMAInput <- function(list.in) {
  dir <- file.path(
    l.a$out.temp, 
    paste0("MA.in", ifelse(is.na(pol), "", pol))
  )
  if(!dir.exists(dir)) dir.create(dir, recursive = T)
  list.out <- lapply(
    1:length(list.in), function(i) {
      df.temp <- list.in[[i]]
      df.temp <- df.temp[df.temp$deg == 1, "KEGG", drop = F]
      write.table(
        df.temp, 
        file = file.path(dir, paste0(l.a$subDir[i], ".txt")), 
        quote = F, row.names = F, col.names = F
      )
      return(df.temp)
    }
  )
  return(NULL)
}

MetaboliteMapping <- function() {
  require(dplyr)
  {  # 创建mapping表
    load(
      file = switch(
        l.a$type, GC = paste0(l.a$input,"df.fiehn.RData"), LC = paste0(l.a$input,"ZhuMetLib.RData")
      )
    )
    if (l.a$type == "GC") {
      df.map <- left_join(
        df.description[, c("id", "Peak"), drop = F], 
        df.fiehn, 
        by = "Peak"
      )
    } else if (l.a$type == "LC") {
      df.map <- left_join(
        df.description[, c("id", "MS2 name"), drop = F], 
        df.zhu, 
        by = "MS2 name"
      )
    }
    exp.cpd <- expression(
      if(l.a$type == "GC") {
        !grepl("Analyte|unknown", df.map$Peak)
      } else if(l.a$type == "LC") {
        !is.na(df.map$`MS2 name`)
      }
    )
    df.map <- df.map[
      eval(exp.cpd), 
      colnames(df.map) != "cpd.bi"
    ]
    colnames(df.map)[colnames(df.map) == "cpd.KEGG"] <- "KEGG COMPOUND name"
  }
  {  # 创建workbook对象
    wb <- createWorkbook(creator = 'Sensichip')
    modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
    hs <- createStyle(
      fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
      halign = 'center', valign = 'center'
    )
    addWorksheet(
      wb, sheetName = 'file explanation', gridLines = T, 
      tabColour = '#4F81BD'
    )
    df.desc <- data.frame(
      "表头名称" = colnames(df.map), 
      "描述" = sapply(colnames(df.map), Explanation2KEGGMetMapping)
    )
    writeDataTable(  # 表格式
      wb, sheet = 1, df.desc, headerStyle = hs, 
      tableStyle = 'TableStyleMedium2', withFilter = T
    )
    setColWidths(wb, sheet = 1, cols = 1, 30)  # 第一列列宽
    setColWidths(wb, sheet = 1, cols = 2, 100)  # 第二列列宽
    setRowHeights(wb, sheet = 1, rows = 1, 24)  # 表头行高
    setRowHeights(wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), 16)  # 表行高
    addStyle(
      wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), cols = 1:ncol(df.desc), 
      gridExpand = T, stack = T, 
      style = createStyle(valign = 'center')
    )
  }
  {  # 写mapping表
    addWorksheet(
      wb, sheetName = "Metabolite Mapping", 
      gridLines = T, tabColour = '#4F81BD'
    )
    writeDataTable(
      wb, sheet = 2, df.map, headerStyle = hs, 
      tableStyle = 'TableStyleMedium2', withFilter = T
    )
    setColWidths(
      wb, sheet = 2, cols = 1:(ncol(df.map) - 2), 18
    )
    setColWidths(wb, sheet = 2, cols = ncol(df.map) - 1, 80)
    setColWidths(wb, sheet = 2, cols = ncol(df.map), 160)
    setRowHeights(wb, sheet = 2, rows = 1, 24)
    addStyle(
      wb, sheet = 2, 
      style = createStyle(
        halign = 'left', valign = 'center', wrapText = T
      ), 
      gridExpand = T, stack = T, 
      cols = 1:ncol(df.map), 
      rows = 1:(nrow(df.map) + 1)
    )
    freezePane(wb, sheet = 2, firstRow = T)
  }
  saveWorkbook(
    wb, overwrite = T, 
    file.path(
      l.a$output, 'results', ifelse(is.na(pol), '', paste0(pol, '/')),
      paste0(ifelse(is.na(pol), '', paste0(pol, '-')),'Metabolite Mapping.xlsx')
    )
  )
  return(df.map)
}
PathwayList <- function(list.in) {
  # 
  list.in <- lapply(  # 统一输出表格的列名
    list.in, `colnames<-`, 
    c(
      "Pathway", "Description", 
      "# compounds (dem)", "Compounds (dem)", 
      "# compounds (all)", "Compounds (all)"
    )
  )
  wb <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    halign = 'center', valign = 'center'
  )
  addWorksheet(
    wb, sheetName = 'file explanation', gridLines = T, tabColour = '#4F81BD'
  )
  df.desc <- data.frame(
    "表头名称" = colnames(list.in[[1]]), 
    "描述" = sapply(colnames(list.in[[1]]), Explanation2KEGGPathway)
  )
  writeDataTable(  # 表格式
    wb, sheet = 1, df.desc, headerStyle = hs, 
    tableStyle = 'TableStyleMedium2', withFilter = T
  )
  setColWidths(wb, sheet = 1, cols = 1, 30)  # 第一列列宽
  setColWidths(wb, sheet = 1, cols = 2, 100)  # 第二列列宽
  setRowHeights(wb, sheet = 1, rows = 1, 24)  # 表头行高
  setRowHeights(wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), 16)  # 表行高
  addStyle(
    wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), cols = 1:ncol(df.desc), 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  lapply(
    1:length(l.a$list.compare), function(i) {
      addWorksheet(
        wb, sheetName = l.a$subDir[i], 
        gridLines = T, tabColour = '#4F81BD'
      )
      writeDataTable(
        wb, sheet = i + 1, list.in[[i]], headerStyle = hs, 
        tableStyle = 'TableStyleMedium2', withFilter = T
      )
      setColWidths(wb, sheet = i + 1, cols = 1, 18)
      setColWidths(wb, sheet = i + 1, cols = 2, 40)
      setColWidths(wb, sheet = i + 1, cols = 3, 18)
      setColWidths(wb, sheet = i + 1, cols = 4, 80)
      setColWidths(wb, sheet = i + 1, cols = 5, 18)
      setColWidths(wb, sheet = i + 1, cols = 6, 80)
      setRowHeights(wb, sheet = i + 1, rows = 1, 24)
      addStyle(
        wb, sheet = i + 1, 
        style = createStyle(
          halign = 'left', valign = 'center', wrapText = T
        ), 
        gridExpand = T, stack = T, 
        cols = 1:ncol(list.in[[i]]), 
        rows = 1:(nrow(list.in[[i]]) + 1)
      )
      freezePane(wb, sheet = i + 1, firstRow = T)
    }
  )
  saveWorkbook(
    wb, overwrite = T, 
    file.path(
      l.a$output, '/results/', ifelse(is.na(pol), '', paste0(pol, '/')), 
      'KEGG Analysis/',paste0(ifelse(is.na(pol), '', paste0(pol, '-')),'KEGG Pathway.xlsx')
    )
  )
}
Explanation2KEGGMetMapping <- function(entry) {
  # 
  # entry <- tolower(entry)  # 不区分大小写
  exp <- switch(
    entry, 
    "id" = "该物质在本次定性分析中的唯一数据编号", 
    "Peak" = "代谢物名称", 
    "KEGG COMPOUND name" = "该代谢物在KEGG COMPOUND数据库中对应的名称", 
    "KEGG COMPOUND ID" = "该物质在KEGG COMPOUND数据库中的索引", 
    "Formula" = "该物质的化学组成", 
    "Exact Mass" = "该物质的精确分子量", 
    "NIKKAJI" = "该物质在NIKKAJI数据库中的索引", 
    "ChEBI" = "该物质在ChEBI数据库中的索引", 
    "PubChem" = "该物质在PubChem数据库中的编号", 
    "CAS" = "该物质的CAS号", 
    "HMDB" = "该物质在HMDB数据库中的索引", 
    "Kingdom" = "该物质在HMDB数据库中的一级分类", 
    "Super Class" = "该物质在HMDB数据库中的二级分类", 
    "Class" = "该物质在HMDB数据库中的三级分类", 
    "Sub Class" = "该物质在HMDB数据库中的四级分类", 
    "KEGG Link" = "该物质的KEGG COMPOUND数据库链接", 
    "Pathway" = "该物质映射的KEGG PATHWAY", 
    # 2017/10/30 New Added
    "MS2 name" = "二级质谱定性匹配分析得到的物质名称",
    "labid" = "该物质在本地数据库的编号",
    "mz" = "物质特征离子的质荷比",
    "formula" = "该物质对应的分子式",
    "METLIN" = "该物质在METLIN数据库中的索引",
    "Super.Class" = "该物质在HMDB数据库中的二级分类",
    "Class.y" = "该物质在HMDB数据库中的三级分类中文名称",
    "Sub.Class" = "该物质在HMDB数据库中的四级分类",
    "pathway" = "该物质对应的KEGG数据库代谢通路",
    
    "暂无说明"
  )
  return(exp)
}
Explanation2KEGGPathway <- function(entry) {
  exp <- switch(
    entry, 
    "Pathway" = "该通路的KEGG PATHWAY数据库ID", 
    "Description" = "该通路名称", 
    "# compounds (dem)" = "该通路内差异代谢物的数量", 
    "Compounds (dem)" = "该通路内差异代谢物的KEGG COMPOUND数据库ID", 
    "# compounds (all)" = "该通路内检测到的所有代谢物的数量", 
    "Compounds (all)" = "该通路内检测到的所有代谢物的KEGG COMPOUND数据库ID", 
    "暂无说明"
  )
  return(exp)
}

KEGGPart1 <- function(list.met,
                      keggInFileName = l.a$subDir,
                      keggDatabaseDir = "D:/share/KEGGDB/",
                      keggInDir = file.path(
                        l.a$out.temp, 
                        paste0("KEGG.in", ifelse(is.na(pol), "", pol)))
) {
  # output kegg in files only
  
  for (i in 1:length(list.met)){
    # get sig and all kegg id and logFc
    sigAndAllRawKeggIn <- getSigAndAllRawKeggIn(as.data.frame(list.met[i])) 
    # write kegg in txt
    writeKeggIn(sigAndAllRawKeggIn$sig, 
                outputDir = keggInDir, 
                fileName = keggInFileName[i])
  }
}

KEGGPart2 <- function() {
  require(openxlsx)
  
  list.cpd <- lapply(
    1:length(l.a$list.compare), function(i) {
      cpd <- readLines(
        paste0(
          l.a$out.temp, '/KEGG.out',
          ifelse(is.na(pol), '', pol), '/', 
          l.a$subDir[i], '.txt'
        ), warn = F
      )
      cpd.title <- c(grep('\\(\\d+\\)$', cpd), length(cpd) + 1)
      list.temp <- lapply(
        1:length(cpd.title[1:length(cpd.title) - 1]), function(j) {
          if(length(cpd) == 0) {
            return(
              data.frame(
                'KEGG Pathway' = '', 'Compound' = '', 
                check.names = F
              )
            )
          } else {
            sub.cpd <- cpd[
              (cpd.title[j] + 1):(cpd.title[j + 1] - 1)
            ]
            sub.cpd <- sub.cpd[sub.cpd != '']
            data.frame(
              'KEGG Pathway' = cpd[cpd.title[j]], 
              'Compound' = paste(sub.cpd, collapse = ';'), 
              check.names = F
            )
          }
        }
      )
      data.frame(
        do.call(rbind.data.frame, list.temp), 
        stringsAsFactors = F, check.names = F
      )
    }
  )
  wb <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    valign = 'center', halign = 'center'
  )
  lapply(
    1:length(list.cpd), function(i) {
      addWorksheet(
        wb, sheetName = l.a$subDir[i], 
        gridLines = T, tabColour = '#4F81BD'
      )
      writeDataTable(
        wb, sheet = i, list.cpd[[i]], headerStyle = hs, 
        tableStyle = 'TableStyleMedium2', withFilter = T
      )
      setColWidths(wb, sheet = i, cols = 1, widths = '50')
      setColWidths(wb, sheet = i, cols = 2, widths = '100')
      setRowHeights(wb, sheet = i, rows = 1, 24)
      addStyle(
        wb, sheet = i, 
        style = createStyle(
          halign = 'left', valign = 'center', wrapText = T
        ), gridExpand = T, stack = T, 
        cols = 1:ncol(list.cpd[[i]]), 
        rows = 1:(nrow(list.cpd[[i]]) + 1)
      )
    }
  )
  saveWorkbook(
    wb, overwrite = T, 
    paste0(
      l.a$output, '/results/', ifelse(is.na(pol), '', paste0(pol, '/')), 
      'KEGG Analysis/',paste0(ifelse(is.na(pol), '', paste0(pol, '-')),'KEGG Pathway.xlsx')
    )
  )
}


KEGGPart3 <- function(){
  keggInDir = file.path(l.a$out.temp,
                        paste0("KEGG.in", ifelse(is.na(pol), "", pol)))
  wb <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    valign = 'center', halign = 'center'
  )
  keggout <- KeggCheck()
  lapply(l.a$subDir,function(i){
    if(i %in% l.a$subDir[keggout]) {
      df.keggPM <- data.frame(x1 = "no kegg pathway found")
    } else {
      df.sig <- read.table(paste0(keggInDir,"/",i,".txt"),sep = "\t",header = F)
      query <- unique(df.sig[,1])
      df.keggPM <- getSigAndAllKeggPathway(df.sig = query)
    }
    addWorksheet(
      wb, sheetName = i, 
      gridLines = T, tabColour = '#4F81BD'
    )
    writeDataTable(
      wb, sheet = i,df.keggPM, headerStyle = hs, 
      tableStyle = 'TableStyleMedium2', withFilter = T
    )
    setColWidths(wb, sheet = i, cols = 1, widths = '50')
    setColWidths(wb, sheet= i, cols = 2, widths = '100')
    setRowHeights(wb, sheet = i, rows = 1, 24)
    addStyle(
      wb, sheet = i, 
      style = createStyle(
        halign = 'left', valign = 'center', wrapText = T
      ), gridExpand = T, stack = T, 
      cols = 1:ncol(df.keggPM), 
      rows = 1:(nrow(df.keggPM) + 1)
    )
  })
  saveWorkbook(
    wb, overwrite = T, 
    paste0(
      l.a$output, '/results/', ifelse(is.na(pol), '', paste0(pol, '/')), 
      'KEGG Analysis/',paste0(ifelse(is.na(pol), '', paste0(pol, '-')),'KEGG Pathway.xlsx')
    )
  )
  writeKeggInputExcel(keggInDir, keggInDir)
  drawKeggPics(
    org = l.a$keggOrg,
    inputFile = paste0(keggInDir,"/keggInForLocal.xlsx"),
    outputDir = file.path(l.a$output, "results", ifelse(is.na(pol), "", pol), "KEGG Analysis/")
  )
}
# MA通路分析 ----
biotreeMA <- function(inputFile = "",
                      outputDir = "",
                      dataBaseDir = paste0(l.a$input,"localMAdataBase/"),
                      fileCmpd = paste0(l.a$input,"localMAdataBase/compound_db.rds"),
                      fileSyn = paste0(l.a$input,"localMAdataBase/syn_nms.rds"),
                      queryBy = "name", # 'kegg' 'hmdb' 'name' 为输入查询的Id类型
                      maOrg = l.a$maOrg)
  
  
{ # 使用代码选择模式物种
  
  # given a list of compound names or ids, find matched name or ids from selected databases
  CrossReferencing <- function(q.type, hmdb=T, pubchem=T, chebi=F, kegg=T, metlin=F){
    
    # record the filter for 8 major databases
    return.cols <<- c(hmdb, pubchem, chebi, kegg, metlin);
    
    # record all the data
    if(!exists("name.map")){
      name.map <<- list();
    }
    
    # distribute job
    dataSet$q.type <<- q.type;
    MetaboliteMappingExact(q.type);
    
    # do some sanity check
    todo.inx <-which(is.na(name.map$hit.inx));
    
  }
  # Mapping from different metabolite IDs
  # For compound names to other id, can do exact or approximate match
  # For other IDs, except HMDB ID, all other may return multiple /non-unique hits
  # multiple hits or non-unique hits will all users to manually select
  MetaboliteMappingExact<-function(q.type){
    qvec <- dataSet$cmpd;
    
    # variables to record results
    hit.inx = vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
    match.values = vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
    match.state = vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0 
    
    cmpd.db <- readRDS(fileCmpd);
    if(q.type == "hmdb"){
      hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
      match.values <- cmpd.db$name[hit.inx];
      match.state[!is.na(hit.inx)] <- 1;
    }else if(q.type == "pubchem"){
      hit.inx <- match(tolower(qvec), tolower(cmpd.db$pubchem));
      match.values <- cmpd.db$name[hit.inx];
      match.state[!is.na(hit.inx)] <- 1;
    }else if(q.type == "chebi"){
      hit.inx <- match(tolower(qvec), tolower(cmpd.db$chebi));
      match.values <- cmpd.db$name[hit.inx];
      match.state[!is.na(hit.inx)] <- 1;
    }else if(q.type == "metlin"){
      hit.inx <- match(tolower(qvec), tolower(cmpd.db$metlin));
      match.values <- cmpd.db$name[hit.inx];
      match.state[!is.na(hit.inx)] <- 1;
    }else if(q.type == "kegg"){
      hit.inx <- match(tolower(qvec), tolower(cmpd.db$kegg));
      hit.inx2 <- match(tolower(qvec), rev(tolower(cmpd.db$kegg)));
      
      # unique hits
      nonuniq.hits <- hit.inx + hit.inx2 != nrow(cmpd.db) + 1;
      hit.inx[nonuniq.hits] <- NA;
      match.values <- cmpd.db$name[hit.inx];
      match.state[!is.na(hit.inx)] <- 1;
      
    }else if(q.type == "name"){
      # first find exact match to the common compound names
      hit.inx <- match(tolower(qvec), tolower(cmpd.db$name));
      match.values <- cmpd.db$name[hit.inx];
      match.state[!is.na(hit.inx)] <- 1;
      
      # then try to find exact match to synanyms for the remaining unmatched query names one by one
      syn.db <- readRDS(fileSyn)
      syns.list <-  syn.db$syns.list;
      todo.inx <-which(is.na(hit.inx));
      if(length(todo.inx) > 0){
        for(i in 1:length(syns.list)){
          syns <-  syns.list[[i]];
          hitInx <- match(tolower(qvec[todo.inx]), tolower(syns));
          
          hitPos <- which(!is.na(hitInx));
          if(length(hitPos)>0){
            # record matched ones
            orig.inx<-todo.inx[hitPos];
            hit.inx[orig.inx] <- i;                  
            # match.values[orig.inx] <- syns[hitInx[hitPos]];  # show matched synnames
            match.values[orig.inx] <- cmpd.db$name[i];    # show common name
            match.state[orig.inx] <- 1;
            
            # update unmatched list
            todo.inx<-todo.inx[is.na(hitInx)];
          }
          if(length(todo.inx) == 0) break;
        }
      }
    }else{
      print(paste("Unknown compound ID type:", q.type));
      # guess a mix of kegg and hmdb ids
      hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
      hit.inx2 <- match(tolower(qvec), tolower(cmpd.db$kegg));
      nohmdbInx <- is.na(hit.inx);
      hit.inx[nohmdbInx]<-hit.inx2[nohmdbInx]
    }
    # empty memory
    gc();
    
    name.map$hit.inx <<- hit.inx;
    name.map$hit.values <<- match.values;
    name.map$match.state <<- match.state;
  }
  
  GetMappingResultTable<-function(){
    
    qvec <- dataSet$cmpd;
    if(is.null(qvec)){
      return();
    }
    
    
    hit.inx<-name.map$hit.inx;
    hit.values<-name.map$hit.values;
    match.state<-name.map$match.state;
    
    # contruct the result table with cells wrapped in html tags
    # the unmatched will be highlighted in different background
    csv.res<-matrix("", nrow=length(qvec), ncol=8);
    colnames(csv.res)<-c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "METLIN", "Comment");
    cmpd.db <- readRDS(fileCmpd);
    for (i in 1:length(qvec)){
      hit <-cmpd.db[hit.inx[i], ,drop=F];
      
      csv.res[i, ]<-c(qvec[i],
                      ifelse(match.state[i]==0, "NA", hit.values[i]),
                      ifelse(match.state[i]==0, "NA", hit$hmdb_id),
                      ifelse(match.state[i]==0, "NA", hit$pubchem_id),
                      ifelse(match.state[i]==0, "NA", hit$chebi_id),
                      ifelse(match.state[i]==0, "NA", hit$kegg_id),
                      ifelse(match.state[i]==0, "NA", hit$metlin_id),
                      match.state[i]);
    }
    # return only columns user selected
    
    # add query and match columns at the the beginning, and 'Detail' at the end
    return.cols <- c(TRUE, TRUE, return.cols, TRUE);
    #browser()
    csv.res <- csv.res[,return.cols, drop=F];
    
    # store the value for report
    dataSet$map.table <<- csv.res;
    write.csv(csv.res, file=paste0(outputDir,"/name_map.csv"), row.names=F);
  }
  GetFinalNameMap<-function(){
    
    hit.inx<-name.map$hit.inx;
    hit.values<-name.map$hit.values;
    match.state<-name.map$match.state;
    
    qvec <- dataSet$cmpd;
    nm.mat<-matrix(nrow=length(qvec), ncol=3);
    colnames(nm.mat)<-c("query", "hmdb", "kegg");
    cmpd.db <- readRDS(fileCmpd);
    for (i in 1:length(qvec)){
      hit <-cmpd.db[hit.inx[i], ,drop=F];
      if(match.state[i]==0){
        hmdb.hit <- NA;
        kegg.hit <- NA;
      }else{
        hmdb.hit <- ifelse(nchar(hit.values[i])==0, NA, hit.values[i]);  ####mark1
        kegg.hit <- ifelse(nchar(hit$kegg_id)==0, NA, hit$kegg_id);
      }
      nm.mat[i, ]<-c(qvec[i],hmdb.hit, kegg.hit);
    }
    return(as.data.frame(nm.mat));
  }
  
  # method is "fisher" or "hyperg"
  CalculateOraScore<-function(nodeImp, method){
    
    # make a clean dataSet$cmpd data based on name mapping
    # only valid kegg id will be used
    nm.map <- GetFinalNameMap();
    valid.inx <- !(is.na(nm.map$kegg)| duplicated(nm.map$kegg));
    ora.vec <- nm.map$kegg[valid.inx];
    q.size<-length(ora.vec);
    if(is.na(ora.vec) || q.size==0) {
      nres.mat<-matrix(0, nrow=1, ncol=10);
      colnames(nres.mat)<-c("Total", "Expected", "Hits", "Raw p", "-log(p)", "Holm adjust", "FDR", "Impact","Hits Cpd","Total Cpd");
      nres.mat[1,] <- 0
      rownames(nres.mat) <- "no_result"
      write.csv(nres.mat, file=paste0(outputDir,"/pathway_results.csv"));
      
    }else{
      require(KEGGgraph);
      current.mset <- metpa$mset.list;
      uniq.count <- metpa$uniq.count;
      
      # check if need to be filtered against reference metabolome
      if(use.metabo.filter && exists('metabo.filter.kegg')){
        current.mset <- lapply(current.mset, function(x){x[x %in% metabo.filter.kegg]});
        analSet$ora.filtered.mset <<- current.mset;
        uniq.count <- length(unique(unlist(current.mset, use.names=FALSE)));
      }
      
      hits <- lapply(current.mset, function(x){x[x %in% ora.vec]});
      hit.num<-unlist(lapply(hits, function(x){length(x)}), use.names=FALSE);
      set.size<-length(current.mset);
      set.num<-unlist(lapply(current.mset, length), use.names=FALSE);
      
      # prepare for the result table
      res.mat<-matrix(0, nrow=set.size, ncol=10);
      rownames(res.mat)<-names(current.mset);
      colnames(res.mat)<-c("Total", "Expected", "Hits", "Raw p", "-log(p)", "Holm adjust", "FDR", "Impact","Hits Cpd","Total Cpd");
      
      if(nodeImp == "rbc"){
        imp.list <- metpa$rbc;
      }else{
        imp.list <- metpa$dgr;
      }
      
      res.mat[,1]<-set.num;
      res.mat[,2]<-q.size*(set.num/uniq.count);
      res.mat[,3]<-hit.num;
      
      
      if(method == "fisher"){
        res.mat[,4]<-GetFisherPvalue(hit.num, q.size, set.num, uniq.count)
      }else{
        res.mat[,4]<-phyper(hit.num-1, set.num, uniq.count-set.num, q.size, lower.tail=F)
      }
      res.mat[,5] <- -log(res.mat[,4]);
      
      # adjust for multiple testing problems
      res.mat[,6] <- p.adjust(res.mat[,4], "holm");
      res.mat[,7] <- p.adjust(res.mat[,4], "fdr");
      
      # calculate the sum of importance
      res.mat[,8] <- mapply(function(x, y){sum(x[y])}, imp.list, hits);
      
      res.mat <- res.mat[hit.num>0,,drop = F];
      ord.inx<-order(res.mat[,4], res.mat[,8]);
      analSet$ora.mat <<- signif(res.mat[ord.inx,,drop = F],5);
      #browser()
      analSet$ora.hits <<- hits;
      analSet$node.imp <<- nodeImp;
      
      save.mat <- analSet$ora.mat;
      rownames(save.mat) <- GetORA.pathNames();
      save.mat[,9] <- GetORA.demNames()
      save.mat[,10] <- GetORA.totalNames()
      if(nrow(save.mat)>=1){
        write.csv(save.mat, file=paste0(outputDir,"/pathway_results.csv"));
      }else{
        nres.mat<-matrix(0, nrow=1, ncol=10);
        colnames(nres.mat)<-c("Total", "Expected", "Hits", "Raw p", "-log(p)", "Holm adjust", "FDR", "Impact","Hits Cpd","Total Cpd");
        nres.mat[1,] <- 0
        rownames(nres.mat) <- "no_result"
        write.csv(nres.mat, file=paste0(outputDir,"/pathway_results.csv"));
      }
    }
  }
  
  GetORA.pathNames<-function(){
    hit.inx <- match(rownames(analSet$ora.mat), metpa$path.ids);
    return(names(metpa$path.ids)[hit.inx]);
  }
  GetORA.demNames<-function(){
    hit.ins <- match(rownames(analSet$ora.mat),names(analSet$ora.hits))
    return(sapply(analSet$ora.hits[hit.ins],function(x)paste0(names(x)," cpd:",x,collapse = "; ")))
  }
  GetORA.totalNames<-function(){
    hit.ins <- match(rownames(analSet$ora.mat),names(metpa$mset.list))
    return(sapply(metpa$mset.list[hit.ins],function(x)paste0(names(x)," cpd:",x,collapse = "; ")))
  }
  
  
  
  dataSet <- list()
  dataSet$cmpd <- read.csv(inputFile,header = F,sep = "\t",stringsAsFactors = F)[,1,drop = T]
  CrossReferencing(q.type = queryBy)
  GetMappingResultTable()
  GetFinalNameMap()
  use.metabo.filter <- F
  analSet <- list()
  load(paste0(dataBaseDir,maOrg,".rda"))
  CalculateOraScore(nodeImp = "rbc", method = "hyperg")
  
}
MetaboAnalyst <- function() {
  require(openxlsx)
  require(treemap)
  require(ggplot2)
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    valign = 'center', halign = 'center'
  )
  wb1 <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb1, fontSize = 8, fontName = 'Arial')
  addWorksheet(
    wb1, sheetName = 'file explanation', 
    gridLines = T, tabColour = '#4F81BD'
  )
  df.disc1 <- data.frame(
    '表头名称' = c(
      'Pathway', 'Total', 'Hits', 'Raw p', '-ln(p)', 
      'Holm adjust', 'FDR', 'Impact','Hits Cpd','Total Cpd'
    ), 
    '描述' = c(
      '代谢通路名称', '该通路中的代谢物个数', 
      '差异代谢物命中该通路的个数', 
      '代谢通路富集分析的P值', 'P值取以e为底的负对数（负自然底对数）', 
      '经Holm-Bonferroni方法进行多重假设检验校正后的P值', 
      paste0(
        '经错误发现率（false discovery rate, FDR）方法', 
        '进行多重假设检验校正后的P值'
      ), 
      '代谢通路拓扑分析的Impact值',
      '命中该通路的差异代谢物名称及KEGG ID',
      '该通路包含的全部代谢物名称及KEGG ID'
    )
  )
  writeDataTable(
    wb1, sheet = 1, df.disc1, 
    headerStyle = hs, tableStyle = 'TableStyleMedium2', withFilter = T
  )
  setColWidths(wb1, sheet = 1, cols = 1:2, 'auto')
  setRowHeights(wb1, sheet = 1, rows = 1, 24)
  setRowHeights(wb1, sheet = 1, rows = 2:(nrow(df.disc1) + 1), 16)
  addStyle(
    wb1, sheet = 1, rows = 2:(nrow(df.disc1) + 1), cols = 1:ncol(df.disc1), 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  wb2 <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb2, fontSize = 8, fontName = 'Arial')
  addWorksheet(
    wb2, sheetName = 'file explanation', 
    gridLines = T, tabColour = '#4F81BD'
  )
  df.disc2 <- data.frame(
    '表头名称' = c(
      'Query', 'Match', 'HMDB', 'PubChem', 'KEGG', 'Comment'
    ), 
    '描述' = c(
      '输入物质名称或编号', '匹配物质名称', '该物质在HDMB数据库中的id', 
      '该物质在PubChem数据库中的id', 
      '该物质在KEGG数据库中的id', 
      '匹配情况：0表示无匹配，1表示精确匹配，2表示模糊匹配'
    )
  )
  writeDataTable(
    wb2, sheet = 1, df.disc2, 
    headerStyle = hs, tableStyle = 'TableStyleMedium2', withFilter = T
  )
  setColWidths(wb2, sheet = 1, cols = 1:2, widths = 'auto')
  setRowHeights(wb2, sheet = 1, rows = 1, 24)
  setRowHeights(wb2, sheet = 1, rows = 2:(nrow(df.disc2) + 1), 16)
  addStyle(
    wb2, sheet = 1, rows = 2:(nrow(df.disc2) + 1), cols = 1:ncol(df.disc2), 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  lapply(
    1:length(l.a$list.compare), function(i) {
      dir <- paste0(
        l.a$out.temp, '/MA.out', 
        ifelse(is.na(pol), '', pol), 
        '/', l.a$subDir[i]
      )
      files <- list.files(dir)
      zip.file <- grep('\\.zip$', files, value = T)
      if(length(zip.file) == 1) {
        unzip(zipfile = paste(dir, zip.file, sep = '/'), exdir = dir)
      }else{
        dirIn <- paste0(
          l.a$out.temp, '/MA.in', 
          ifelse(is.na(pol), '', pol), 
          '/', l.a$subDir[i],".txt"
        )
        cmpd <- read.csv(dirIn,header = F,sep = "\t",stringsAsFactors = F)[,1,drop = T]
        biotreeMA(inputFile = dirIn,
                  outputDir = dir)
      }
      df.ma <- read.csv(
        paste0(dir, '/pathway_results.csv'), 
        stringsAsFactors = F, check.names = F
      )#[,1:9]
      j <- which(names(df.ma) == '-log(p)')
      k <- which(names(df.ma) == 'Impact')
      l <- which(names(df.ma) == '')
      colnames(df.ma)[c(j, k, l)] <- c('p', 'impact', 'pathway')
      if(length(df.ma$`p`)==0){
        df.ma <- rbind(df.ma,"no_result")
      }
      if(length(unique(df.ma$`p`)) >= 2 & length(unique(df.ma$`impact`)) >= 2){
        if(l.a$ma.treemap) {
          df.ma$impact.mod <- df.ma$impact + 0.001
          jpeg(
            filename = paste0(
              l.a$output, '/results/', 
              ifelse(is.na(pol), '', paste0(pol, '/')), 
              'Pathway Analysis/', 
              l.a$subDir[i], '/treemap.jpg'
            )
            , res = 600, width = 6000, height = 6000
          )
          treemap::treemap(
            df.ma, index = 'pathway', vSize = 'impact.mod', 
            vColor = 'p', title = '', 
            title.legend = expression(
              paste(-ln, ' ', italic('P'), '-value')
            ), 
            type = 'value', palette = '-RdYlGn', 
            range = c(0, 3), mapping = c(0, 1.5, 3)
          )
          dev.off()
          pdf(
            paste0(
              l.a$output, '/results/', 
              ifelse(is.na(pol), '', paste0(pol, '/')), 
              'Pathway Analysis/', 
              l.a$subDir[i], '/treemap.pdf'
            )
          )
          treemap::treemap(
            df.ma, index = 'pathway', vSize = 'impact.mod', 
            vColor = 'p', title = '', 
            title.legend = expression(
              paste(-ln, ' ', italic('P'), '-value')
            ), 
            type = 'value', palette = '-RdYlGn', 
            range = c(0, 3), mapping = c(0, 1.5, 3)
          )
          dev.off()
        }
        if(l.a$ma.bubble) {
          p <- DrawBubblePlot(df.ma)
          ggsave(
            paste0(
              l.a$output, '/results/', 
              ifelse(is.na(pol), '', paste0(pol, '/')), 
              'Pathway Analysis/', 
              l.a$subDir[i], '/Bubble Plot.pdf'
            ), p, 
            width = 9.6, height = 7.2, units = 'in', dpi = 600
          )
          ggsave(
            paste0(
              l.a$output, '/results/', 
              ifelse(is.na(pol), '', paste0(pol, '/')), 
              'Pathway Analysis/', 
              l.a$subDir[i], '/Bubble Plot.jpg'
            ), p, 
            width = 9.6, height = 7.2, units = 'in', dpi = 600
          )
        }
      }else{
        writeLines("no_result",
                   paste0(
                     l.a$output, '/results/', 
                     ifelse(is.na(pol), '', paste0(pol, '/')), 
                     'Pathway Analysis/', 
                     l.a$subDir[i], '/no_result.txt'
                   )
        )
        
      }
      df.ma <- df.ma[
        , colnames(df.ma) != 'Expected' & 
          colnames(df.ma) != 'impact.mod'
      ]
      colnames(df.ma)[1:8] <- c(
        'Pathway', 'Total', 'Hits', 'Raw p', '-ln(p)', 
        'Holm adjust', 'FDR', 'Impact'
      )
      addWorksheet(
        wb1, sheetName = l.a$subDir[i], 
        gridLines = T, tabColour = '#4F81BD'
      )
      writeDataTable(
        wb1, sheet = i + 1, df.ma, headerStyle = hs, 
        tableStyle = 'TableStyleMedium2', withFilter = T
      )
      setColWidths(
        wb1, sheet = i + 1, cols = 1:ncol(df.ma), widths = 'auto'
      )
      setRowHeights(wb1, sheet = i + 1, rows = 1, 24)
      setRowHeights(
        wb1, sheet = i + 1, rows = 2:(nrow(df.ma) + 1), 16
      )
      addStyle(
        wb1, sheet = i + 1, rows = 2:(nrow(df.ma) + 1), 
        cols = 1:ncol(df.ma), 
        gridExpand = T, stack = T, 
        style = createStyle(valign = 'center')
      )
      df.map <- read.csv(
        paste0(dir, '/name_map.csv'), 
        stringsAsFactors = F, check.names = F
      )
      if(nrow(df.map) == 0){df.map <- rbind(df.map,"no_result")}
      df.map <- df.map[df.map$Query != 'KEGG COMPOUND id', ]
      addWorksheet(
        wb2, sheetName = l.a$subDir[i], 
        gridLines = T, tabColour = '#4F81BD'
      )
      writeDataTable(
        wb2, sheet = i + 1, df.map, headerStyle = hs, 
        tableStyle = 'TableStyleMedium2', withFilter = T
      )
      setColWidths(
        wb2, sheet = i + 1, cols = 1:ncol(df.map), widths = 'auto'
      )
      setRowHeights(wb2, sheet = i + 1, rows = 1, 24)
      setRowHeights(
        wb2, sheet = i + 1, rows = 2:(nrow(df.map) + 1), 16
      )
      addStyle(
        wb2, sheet = i + 1, rows = 2:(nrow(df.map) + 1), 
        cols = 1:ncol(df.map), 
        gridExpand = T, stack = T, 
        style = createStyle(valign = 'center')
      )
      
      # on.exit(
      #     file.remove(  # 删除文档以节省硬盘空间
      #         file.path(dir, files)[!grepl("\\.zip$", files)]
      #     )
      #)
      
    }
  )
  saveWorkbook(
    wb1, overwrite = T, 
    paste0(
      l.a$output, '/results/', 
      ifelse(is.na(pol), '', paste0(pol, '/')), 
      'Pathway Analysis/', 
      ifelse(is.na(pol), '', paste0(pol, '-')), 
      'Pathway Analysis.xlsx'
    )
  )
  saveWorkbook(
    wb2, overwrite = T, 
    paste0(
      l.a$output, '/results/', 
      ifelse(is.na(pol), '', paste0(pol, '/')), 
      'Pathway Analysis/', 
      ifelse(is.na(pol), '', paste0(pol, '-')), 
      'Metabolite Mapping.xlsx'
    )
  )
}
DrawBubblePlot <- function(df.ma) {
  require(ggplot2)
  require(ggrepel)
  if(length(which(df.ma$impact > 0)) >= 5) {
    filter <- expression(
      df.ma$impact >= sort(df.ma$impact, decreasing = T)[5]
    )
  } else {
    filter <- expression(df.ma$impact >= 0)
  }
  p <- ggplot() +
    geom_point(
      data = df.ma, shape = 21, color = 'black', stroke = 1, 
      mapping = aes(x = impact, y = p, size = impact, fill = p)
    ) +
    geom_label_repel(
      data = df.ma[eval(filter) | df.ma$`Raw p` < 0.05, ], 
      mapping = aes(impact, p, fill = p, label = pathway), 
      color = 'black', size = 3, 
      label.padding = unit(0.2, 'lines'), 
      point.padding = unit(0.5, 'lines'), 
      min.segment.length = unit(0.1, "lines"), 
      segment.color = 'grey50', segment.size = 1, 
      show.legend = F
    ) +
    theme_bw() +
    theme(
      #panel.grid.major =element_blank(),
      panel.grid.major = element_line(
        color = '#3f00f6', linetype = 'dashed'
      ), 
      panel.grid.minor = element_blank()
    ) +
    scale_size_continuous(
      name = 'Impact', range = c(8, 16), 
      breaks = c(min(df.ma$impact), max(df.ma$impact)), 
      labels = c(min(df.ma$impact), round(max(df.ma$impact), 3)), 
      guide = guide_legend(title.hjust = 0.5)
    ) + 
    scale_fill_gradient(
      name = expression(paste(-ln, ' ', italic('P'), '-value')), 
      low = '#ffeda0', high = '#f03b20', 
      guide = guide_colorbar(
        title.position = 'top', title.hjust = 0.5, 
        direction = 'horizontal'
      )
    ) +
    labs(
      x = 'Impact', 
      y = expression(paste(-ln, ' ', italic('P'), '-value'))
    )
}
MetaboAnalyst2KEGGCol <- function(chr.pol) {
  # 从MetaboAnalyst的分析结果直接制作KEGGCol文件
  # 使用前确保MetaboAnalyst的分析结果.zip文件已存放至指定的文件夹
  # 输出对应的KEGGCol.txt
  # 函数返回列表，所有元素全是NULL
  # KEGG的后续分析应该能套用KEGGPart2()函数
  # 
  lapply(
    chr.pol, function(pol) {
      #load(
      #    file.path(
      #        l.a$out.temp, 
      #       paste0(ifelse(is.na(pol), '', paste0(pol, ".")), "l.a.RData")   
      #   ))
      lapply(
        1:length(l.a$list.compare), function(i) {
          dir <- file.path(
            l.a$out.temp, paste0("MA.out", ifelse(is.na(pol), '', pol)), 
            l.a$subDir[i]
          )
          files <- list.files(dir)
          zip.file <- grep('\\.zip$', files, value = T)
          if(length(zip.file) == 1){
            unzip(zipfile = file.path(dir, zip.file), exdir = dir)
          }
          df.map2 <- read.csv(
            file.path(dir, "name_map.csv"), na.strings = c("NA", "")
          )
          df.pathway <- read.csv(
            file.path(dir, "pathway_results.csv"), na.strings = c("NA", "")
          )[1,1]
          if(df.pathway != "no_result") {
            #files <- list.files(dir)
            # file.remove(  # 删除文档以节省硬盘空间
            #     file.path(dir, files)[!grepl("\\.zip$", files)]
            # )
            df.map2 <- df.map2[!is.na(df.map2$KEGG), c("KEGG", "Query")]
            #筛选出差异物质
            df.meta <- ExtractionDiffMet(list.compare.result)[[i]]
            #names(df.meta)[
            #  grepl("Peak|(MS2 name)|(compound name)", names(df.meta))
            #] <- "met"
            # names(df.meta)[
            # grepl("Similarity|MS2.score|(MS2 score)|^score$", names(df.meta))
            #] <- "score"
            #删除差异物质重复物质
            df.meta <- HeatmapFilter(df.meta)
            df.meta$met <- sub("^\n\n|\\s$", "", df.meta$met)  
            # 去掉开头双换行符和/或结尾单空格
            df.temp <- df.meta[
              !is.na(df.meta$met), c("met", "LOG_FOLDCHANGE")
            ]
            df.map2$col <- sapply(
              df.map2$Query, function(metabolite) {
                fc <- df.temp$LOG_FOLDCHANGE[
                  df.temp$met == metabolite
                ]
                if(length(fc) > 1){stop(paste0("发现重名物质",metabolite,"，请检查重名！"))}
                ifelse(fc > 0, "red", "blue")
              }
            )
            df.map2 <- df.map2[, -2]
            
            write.table(
              df.map2, 
              file = file.path(
                l.a$out.temp, paste0("KEGG.in", ifelse(is.na(pol), '', pol)), 
                paste0(l.a$subDir[i], ".txt")
              ), quote = F, sep = '\t', row.names = F, col.names = F
            )
          } else {
            writeLines("no_result",
                       paste0(
                         l.a$output, '/results/', 
                         ifelse(is.na(pol), '', paste0(pol, '/')), 
                         'KEGG Analysis/', 
                         l.a$subDir[i], '/no_result.txt'
                       )
            )
          }
        }
      )
    }
  )
}
## 返回通路分析没有富集到通路的对比序号
KeggCheck <- function() {
  for (i in 1:length(l.a$list.compare)) {
    dir <- file.path(
      l.a$out.temp, paste0("MA.out", ifelse(is.na(pol), '', pol)), 
      l.a$subDir[i]
    )
    df.pathway <- read.csv(
      file.path(dir, "pathway_results.csv"), na.strings = c("NA", "")
    )[1,1]
    if(df.pathway == "no_result") {
      l.a$kegg.out[length(l.a$kegg.out)+1] <- i
    }
  }
  return(l.a$kegg.out)
}
# 转录组关联 ----
MetaboAnalyst2MetScape <- function(chr.pol) {
  # 从MetaboAnalyst的结果中提取LC项目差异代谢物Mapping信息，回到compare表中
  # 提取P值和fold change，输出表格到temp/MetScape文件夹
  # 
  require(dplyr)
  lapply(
    chr.pol, function(pol) {
      load(file.path(
        l.a$out.temp, 
        paste0(ifelse(is.na(pol), '', paste0(pol, ".")), "l.a.RData")   
      ))
      dir.out <- file.path(
        l.a$out.temp, paste0("MetScape", ifelse(is.na(pol), '', pol))
      )
      dir.create(dir.out, showWarnings = F, recursive = T)
      file.cys <- file.path(dir.out, paste0(pol, ".cys"))
      if(!file.exists(file.cys)) file.create(file.cys)
      lapply(
        1:length(l.a$list.compare), function(i) {
          dir.in <- file.path(
            l.a$out.temp, paste0("MA.out", ifelse(is.na(pol), '', pol)), 
            l.a$subDir[i]
          )
          files <- list.files(dir.in)
          zip.file <- grep('\\.zip$', files, value = T)
          length(zip.file) == 1 || return(NULL)
          unzip(zipfile = file.path(dir.in, zip.file), exdir = dir.in)
          df.map <- read.csv(
            file.path(dir.in, "name_map.csv"), 
            na.strings = c("NA", ""), stringsAsFactors = F
          )
          file.remove(  # 删除文档以节省硬盘空间
            file.path(dir.in, files)[!grepl("\\.zip$", files)]
          )
          df.map <- df.map[!is.na(df.map$KEGG), c("Query", "KEGG")]
          df.meta <- list.compare.result[[i]]
          names(df.meta)[
            grepl("Peak|(MS2 name)|(compound name)", names(df.meta))
          ] <- "met"
          df.temp <- df.meta[
            !is.na(dfmeta$`MS2 name`), 
            c("MS2 name", "P-VALUE", "LOG_FOLDCHANGE")
          ]
          names(df.temp)[names(df.temp) == "MS2 name"] <- "Query"
          df.map <- dplyr::left_join(df.map, df.temp, by = "Query")
          df.map <- df.map[, -1]
          write.table(
            df.map, file = file.path(
              dir.out, paste0(l.a$subDir[i], ".txt")
            ), quote = F, sep = '\t', row.names = F, col.names = T
          )
          dir.result <- file.path(
            l.a$output, "results", 
            ifelse(is.na(pol), '', pol), 
            "Network Analysis", 
            l.a$subDir[i]
          )
          dir.create(dir.result, recursive = T, showWarnings = F)
          f <- file.path(
            dir.result, 
            paste0(pol, "_", l.a$subDir[i], "_network.csv")
          )
          if(!file.exists(f)) file.create(f)
        }
      )
    }
  )
}
# ROC曲线 ----
ROC4GC <- function() {
  # 给GC数据做ROC曲线
  # 需要读取SIMCA.csv和差异代谢物筛选表
  # 对每个差异代谢物进行单独的ROC曲线绘制，使用{plotROC}包
  # 在结果目录下创建ROC Curve文件夹，为每组对比创建ROC曲线
  # 函数本身没有返回值
  # 
  require(reshape2)
  require(ggplot2)
  require(plotROC)
  load(file.path(l.a$out.temp, "l.a.RData"))
  df.roc <- df.out.SIMCA
  lapply(
    1:length(l.a$list.compare), function(i) {
      df.deg <- list.compare.result[[i]]
      df.deg <- df.deg[eval(l.a$deg.exp), c("id", "Peak")]
      df.deg <- df.deg[!is.na(df.deg$`Peak`), ]
      df.deg <- df.deg[!grepl("unknown", df.deg$`Peak`), ]
      # 临时加的，去掉unknown
      pair <- l.a$list.compare[[i]]
      df.sub <- df.roc[
        mapply(`|`, l.a$group[[pair[1]]], l.a$group[[pair[2]]]), 
        c("classid", df.deg$id)
      ]
      colnames(df.sub) <- c("classid", df.deg$`Peak`)
      df.sub[, -1] <- scale(df.sub[, -1])
      df.sub$D <- as.numeric(
        factor(
          df.sub$classid, 
          levels = names(l.a$group)[c(pair[2], pair[1])]
        )
      ) - 1
      df.sub <- reshape2::melt(
        df.sub, id.vars = c('classid', 'D')
      )
      df.sub$variable <- as.character(df.sub$variable)
      dir <- file.path(
        l.a$output, "results","ROC Curve", l.a$subDir[i]
      )
      if(!dir.exists(dir)) dir.create(dir, recursive = T)
      list.aucs <- lapply(
        unique(df.sub$variable), function(m) {
          p <- ggplot(df.sub[df.sub$variable == m,],
                      aes(d = D, m = value)) + geom_roc(
                        labelround = 2,
                        hjust = 0.4,
                        vjust = 0.9,
                        linealpha = 0.6
                      ) + geom_rocci(labels = F) +
            style_roc(guide = TRUE)
          auc <- calc_auc(p)$AUC
          if(auc >= 0.5) {
            p <- p +
              annotate(
                'label', x = 0.9, y = 0.1, 
                fill = '#e31a1c', color = 'white', 
                label = paste0(
                  'AUC = ', round(auc, digits = 2)
                )
              )
          } else {
            p <- p +
              annotate(
                'label', x = 0.1, y = 0.9, 
                fill = '#1f78b4', color = 'white', 
                label = paste0(
                  'AUC = ', round(auc, digits = 2)
                )
              )
          }
          m.sub <- gsub(
            "/|:|\\*|\\?|\"|<|>|\\|", "", m, perl = T
          )
          ggsave(
            file.path(dir, paste0(m.sub, '.jpg')),
            p, dpi = 600, 
            height = 7, width = 7, units = "in"
          )
          ggsave(
            file.path(dir, paste0(m.sub, '.pdf')),
            p, dpi = 600, 
            height = 7, width = 7, units = "in"
          )
          return(list(m, m.sub, auc))
        }
      )
      df.aucs <- data.frame(
        "Peak" = sapply(list.aucs, `[[`, 1, USE.NAMES = F), 
        "pic name" = sapply(list.aucs, `[[`, 2, USE.NAMES = F), 
        "auc" = sapply(list.aucs, `[[`, 3, USE.NAMES = F), 
        check.names = F
      )
      df.aucs <- df.aucs[order(df.aucs$`pic name`), ]
      write.csv(
        df.aucs, file = file.path(dir, "auc.csv"), 
        row.names = F
      )
    }
  )
}

ROC4LC <- function() {
  # 给LC数据做ROC曲线
  # 需要读取正负离子模式下的SIMCA.csv和差异代谢物筛选表
  # 对每个差异代谢物进行单独的ROC曲线绘制，使用{plotROC}包
  # 分别在正负离子模式的结果目录下创建ROC Curve文件夹，为每组对比创建ROC曲线
  # 函数本身没有返回值
  # 
  require(reshape2)
  require(ggplot2)
  require(plotROC)
  lapply(
    c("POS", "NEG"), function(pol) {
      load(file.path(l.a$out.temp, paste0(pol, ".l.a.RData")))
      df.roc <- df.out.SIMCA
      lapply(
        1:length(l.a$list.compare), function(i) {
          df.deg <- list.compare.result[[i]]
          df.deg <- df.deg[eval(l.a$deg.exp), c("id", "MS2 name")]
          df.deg <- df.deg[!is.na(df.deg$`MS2 name`), ]
          df.deg <- df.deg[!grepl("unknown", df.deg$`MS2 name`), ]
          # 临时加的，去掉unknown
          pair <- l.a$list.compare[[i]]
          df.sub <- df.roc[
            mapply(`|`, l.a$group[[pair[1]]], l.a$group[[pair[2]]]), 
            c("classid", df.deg$id)
          ]
          colnames(df.sub) <- c("classid", df.deg$`MS2 name`)
          df.sub[, -1] <- scale(df.sub[, -1])
          df.sub$D <- as.numeric(
            factor(
              df.sub$classid, 
              levels = names(l.a$group)[c(pair[2], pair[1])]
            )
          ) - 1
          df.sub <- reshape2::melt(
            df.sub, id.vars = c('classid', 'D')
          )
          df.sub$variable <- as.character(df.sub$variable)
          dir <- file.path(
            l.a$output, "results", pol, "ROC Curve", l.a$subDir[i]
          )
          if(!dir.exists(dir)) dir.create(dir, recursive = T)
          list.aucs <- lapply(
            unique(df.sub$variable), function(m) {
              p <- ggplot(
                df.sub[df.sub$variable == m, ], 
                aes(d = D, m = value)
              ) + geom_roc(
                labelround = 2,
                hjust = 0.4,
                vjust = 0.9,
                linealpha = 0.6
              ) + geom_rocci(labels = F) +
                style_roc(guide = TRUE)
              
              auc <- calc_auc(p)$AUC
              if(auc >= 0.5) {
                p <- p +
                  annotate(
                    'label', x = 0.9, y = 0.1, 
                    fill = '#e31a1c', color = 'white', 
                    label = paste0(
                      'AUC = ', round(auc, digits = 2)
                    )
                  )
              } else {
                p <- p +
                  annotate(
                    'label', x = 0.1, y = 0.9, 
                    fill = '#1f78b4', color = 'white', 
                    label = paste0(
                      'AUC = ', round(auc, digits = 2)
                    )
                  )
              }
              m.sub <- gsub(
                "/|:|\\*|\\?|\"|<|>|\\|", "", m, perl = T
              )
              ggsave(
                file.path(dir, paste0(m.sub, '.jpg')),
                p, dpi = 600, 
                height = 7, width = 7, units = "in"
              )
              ggsave(
                file.path(dir, paste0(m.sub, '.pdf')),
                p, dpi = 600, 
                height = 7, width = 7, units = "in"
              )
              return(list(m, m.sub, auc))
            }
          )
          df.aucs <- data.frame(
            "MS2 name" = sapply(list.aucs, `[[`, 1, USE.NAMES = F), 
            "pic name" = sapply(list.aucs, `[[`, 2, USE.NAMES = F), 
            "auc" = sapply(list.aucs, `[[`, 3, USE.NAMES = F), 
            check.names = F
          )
          df.aucs <- df.aucs[order(df.aucs$`pic name`), ]
          write.csv(
            df.aucs, file = file.path(dir, "auc.csv"), 
            row.names = F
          )
        }
      )
    }
  )
}

# ANOVA ----
oneWayANOVA <- function(includedSampleNames = ifelse(is.na(l.a$qc), 
                                                     paste(names(l.a$group), collapse = " "),
                                                     paste(
                                                       names(l.a$group)[1:length(names(l.a$group)) - 1],
                                                       collapse = " "
                                                     )
),
dfData = df.out.SIMCA[, -1],
dfDescription = df.description,
outputDir =  paste0(
  l.a$output, '/results/',
  ifelse(is.na(pol), '', paste0(pol, '/')),
  'ANOVA/'
)
)
{ # input the sample name that wanted to be included in the one way anova analysis in a single string
  # default is all the sample name from l.a$group except the qc
  # i.e.
  # oneWayANOVA("A B C")
  if(!dir.exists(outputDir)) dir.create(path = outputDir,recursive = TRUE)
  l.a$anova[length(l.a$anova)+1] <- includedSampleNames
  # sample name transform
  sampleNames <- unlist(strsplit(includedSampleNames, split = " "))
  
  # get included sample data
  includedSampleBoolean <- l.a$group[sampleNames]
  includedAllSampleBoolean <- Reduce(`|`, includedSampleBoolean)
  includedSampleData <- dfData[includedAllSampleBoolean, ]
  
  # sample group information
  sampleGroupInfo <- NULL
  for (sampleName in sampleNames){
    sampleGroupInfo[l.a$group[[sampleName]]] <- sampleName
  }
  sampleGroupInfo <- sampleGroupInfo[includedAllSampleBoolean]
  if (NA %in% sampleGroupInfo){stop("NA exisits in sample group information")}
  
  # anova
  pOfOneWayANOVA  <- NULL
  for (i in 1:ncol(includedSampleData)){
    singleTestData <- data.frame(data = includedSampleData[,i], sampleGroupInfo = sampleGroupInfo)
    anovaResult <- performOneWayANOVA(singleTestData)
    pOfOneWayANOVA  <- c(pOfOneWayANOVA , anovaResult)
  }
  qOfOneWayANOVA <- p.adjust(pOfOneWayANOVA, method = "BH", n = length(pOfOneWayANOVA))
  
  # make matrix
  resultDataframe <- cbind(dfDescription, pOfOneWayANOVA, qOfOneWayANOVA, t(includedSampleData))
  colnames(resultDataframe) <- c(colnames(dfDescription), "ANOVA P-VALUE", "Q-VALUE", rownames(includedSampleData))
  # file explanation
  desc1 <- sapply(colnames(resultDataframe)[1:(ncol(resultDataframe) - nrow(includedSampleData))],
                  Explanation4Desc)
  desc2 <- sapply(
    colnames(resultDataframe)[(ncol(resultDataframe) - nrow(includedSampleData) + 1):ncol(resultDataframe)],
    function(title) {
      if(grepl('^Mean ', title)) {
        paste0(sub('^Mean ', '', title), '样本重复实验的相对定量值均值')
      } else {
        paste0(title, '样品相对定量值')
      }
    }
  )
  df.desc <- data.frame(
    '表头名称' = colnames(resultDataframe), '描述' = c(desc1, desc2)
    ,stringsAsFactors = FALSE)
  
  #输出满足筛选条件的数据
  resultDataframe <- subset.data.frame(resultDataframe,
                                       subset = resultDataframe$`MS2 name`!="NA" | resultDataframe$`MS1 name`!="NA",
                                       select = c(colnames(resultDataframe))
  )
  
  # write output xlsx
  writeOutputXlsx(
    resultDataframe[which(resultDataframe$`ANOVA P-VALUE` < 0.05),],
    
    df.desc,
    subOutputDir = paste0(outputDir, includedSampleNames, "/"),
    "ANOVA.xlsx",
    secondSheetName = ifelse(nchar(includedSampleNames) > 25,
                             paste0(substr(includedSampleNames, 1, 25), "..."),
                             includedSampleNames)
  )
  
  # heatmap
  resultDataframeSig <- resultDataframe[which(resultDataframe$`ANOVA P-VALUE` < 0.05),]
  resultDataframeSig <- resultDataframeSig[which(!is.na(resultDataframeSig$`MS2 name`)),]
  #resultDataframeSig <- na.omit(resultDataframeSig)
  #write.table(resultDataframe,file = paste0(outputDir, includedSampleNames, "/","ANOVA.txt"))
  # jpg
  pheatmap::pheatmap(
    resultDataframeSig[,(l.a$desc+3):ncol(resultDataframeSig)], 
    cluster_cols = T, cluster_rows = T, 
    scale = 'row', fontsize = 5, 
    color = colorRampPalette(c("blue","black", "yellow"))(1000),
    cellwidth = 16, cellheight = 5, border_color = NA, 
    labels_row = resultDataframeSig$`MS2 name`, 
    filename = paste0(outputDir, includedSampleNames, "/", "ANOVA heatmap.jpg")
  )
  # pdf
  pheatmap::pheatmap(
    resultDataframeSig[, (l.a$desc+3):ncol(resultDataframeSig)], 
    cluster_cols = F, cluster_rows = T, 
    scale = 'row', fontsize = 5,
    color = colorRampPalette(c("blue","black", "yellow"))(1000),
    cellwidth = 16, cellheight = 5, border_color = NA, 
    labels_row = resultDataframeSig$`MS2 name`, 
    filename = paste0(outputDir, includedSampleNames, "/", "ANOVA heatmap.pdf")
  )
  
}

performOneWayANOVA <- function(df.in){
  # performOneWayANOVA function will accept a data frame with two columns, which are data and sample group information.
  # it will first test the homogeneity of variance by bartlett.test(),
  # then pass the result to oneway.test() var.queal parameter.
  # One way ANOVA will be test by oneway.test()
  # it will return the p value.
  
  # homogeneity test of variance
  isEqual <- ifelse(bartlett.test(data~sampleGroupInfo, df.in)$p.value < 0.05,
                    FALSE, TRUE)
  
  # one way anova
  pOfOneWayANOVA <- oneway.test(data~sampleGroupInfo, df.in, var.equal = isEqual)$p.value
  
  return(pOfOneWayANOVA)
}

writeOutputXlsx <- function(dataFrame, 
                            df.desc, 
                            subOutputDir,
                            outputFileName,
                            secondSheetName = "Sheet2"){
  
  wb <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
  
  # write first sheet
  addWorksheet(
    wb, sheetName = 'file explanation', gridLines = T, tabColour = '#4F81BD'
  )
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    valign = 'center', halign = 'center'
  )
  writeDataTable(
    wb, sheet = 1, df.desc, 
    headerStyle = hs, tableStyle = 'TableStyleMedium2', withFilter = T
  )
  # first sheet style
  setColWidths(wb, sheet = 1, cols = 1, widths = 20)
  setColWidths(wb, sheet = 1, cols = 2, widths = max(nchar(df.desc[,2])) * 2)
  setRowHeights(wb, sheet = 1, rows = 1, 24)
  setRowHeights(wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), 16)
  addStyle(
    wb, sheet = 1, rows = 2:(nrow(df.desc) + 1), cols = 1:2, 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  
  # write second sheet
  addWorksheet(
    wb, sheetName = secondSheetName, gridLines = T, tabColour = '#4F81BD'
  )
  writeDataTable(
    wb, sheet = 2, dataFrame, headerStyle = hs, 
    tableStyle = 'TableStyleMedium2', withFilter = T
  )
  # second sheet style
  setColWidths(wb, sheet = 2, cols = 1:ncol(dataFrame), widths = 18)
  setRowHeights(wb, sheet = 2, rows = 1, 24)
  setRowHeights(wb, sheet = 2, rows = 2:(nrow(dataFrame) + 1), 16)
  addStyle(
    wb, sheet = 2, rows = 2:(nrow(dataFrame) + 1), cols = 1:ncol(dataFrame), 
    gridExpand = T, stack = T, 
    style = createStyle(valign = 'center')
  )
  
  # save xlsx
  if(!dir.exists(subOutputDir)) dir.create(subOutputDir)
  saveWorkbook(
    wb, overwrite = T, 
    file = paste0(subOutputDir, outputFileName)
  )
}
#grou_by
#为长数据添加分组
group_by <- function(group_samples = group_samples,Data = Data){
  for(i in 1:length(group_samples)){
    Data$Group[Data$variable %in% group_samples[[i]]] <- names(group_samples)[i]
  }
  return(Data)
}

# Z-score图----
zscore <- function(list.in = list.compare.result) {
  wb <- createWorkbook(creator = "Sensichip")
  modifyBaseFont(wb, fontSize = 8, fontName = "Arial")
  hs <- createStyle(
    fontSize = 10, fontColour = "white", fgFill = "#4F81BD",
    valign = "center", halign = "center"
  )
  main.list <- ExtractionDiffMet(list.compare.result)
  lapply(
    1:length(l.a$list.compare), function(i) {
      df.plot <- main.list[[i]]
      df.plot <- RepnameFilter(df.plot)
      if (is.null(df.plot)) {
        de_df <- data.frame(X1 = "No Differentially Expressed Metabolites")
        writeLines(
          "no_result",
          paste0(
            l.a$output, "/results/",
            ifelse(is.na(pol), "", paste0(pol, "/")),
            "Zscore/",
            l.a$subDir[i], "/no_result.txt"
          )
        )
      } else {
        rank_data <- arrange(df.plot, desc(met))
        de_df <- scale(t(rank_data[, -c(1:(l.a$desc + 7))]))
        de_df <- cbind.data.frame(`met` = rank_data[, grepl("met", colnames(rank_data))], t(de_df))
        if (nrow(de_df) > 50) {
          de_df <- de_df[c(1:50), ]
        } else {
          de_df <- de_df
        }
        ZscoreData <- melt(de_df, value.name = "Zscore")
        group_samples <- lapply(1:length(l.a$group), function(j) {
          group_list <- colnames(df.area)[l.a$group[[j]]]
        })
        names(group_samples) <- names(l.a$group)
        ZscoreData <- group_by(group_samples[l.a$list.compare[[i]]], ZscoreData)
        zscoreplot <- ggplot(data = ZscoreData, aes(x = Zscore, y = met, color = Group)) +
          geom_point(shape = 21, stroke = 1.2, size = 2, alpha = 0.6) +
          xlab("z-score") +
          theme(axis.title.y = element_blank())
        ggsave(
          filename = paste0(
            l.a$output, "/results/",
            ifelse(is.na(pol), "", paste0(pol, "/")),
            "Zscore/", l.a$subDir[i], "/", l.a$subDir[i], ".png"
          ),
          plot = zscoreplot,
          device = "png",
          height = 20,
          width = 16,
          units = "cm"
        )
        ggsave(
          filename = paste0(
            l.a$output, "/results/",
            ifelse(is.na(pol), "", paste0(pol, "/")),
            "Zscore/", l.a$subDir[i], "/", l.a$subDir[i], ".pdf"
          ),
          plot = zscoreplot,
          device = "pdf",
          height = 20,
          width = 16,
          units = "cm"
        )
      }
      addWorksheet(
        wb,
        sheetName = l.a$subDir[i],
        gridLines = T, tabColour = "#4F81BD"
      )
      writeDataTable(
        wb,
        sheet = i, de_df, headerStyle = hs,
        tableStyle = "TableStyleMedium2", withFilter = T
      )
      setColWidths(wb, sheet = i, cols = 1:ncol(de_df), "auto")
      setRowHeights(wb, sheet = i, rows = 1, 24)
      setRowHeights(wb, sheet = i, rows = 2:(nrow(de_df) + 1), 16)
      addStyle(
        wb,
        sheet = i, rows = 2:(nrow(de_df) + 1),
        cols = 1:ncol(de_df), gridExpand = T, stack = T,
        style = createStyle(valign = "center")
      )
    }
  )
  saveWorkbook(
    wb,
    overwrite = T,
    paste0(
      l.a$output, "/results/",
      ifelse(is.na(pol), "", paste0(pol, "/")),
      "Zscore/",
      ifelse(is.na(pol), "", paste0(pol, "-")),
      "Zscore data matrix.xlsx"
    )
  )
}
RepnameFilter <- function(df.meta) {
  df.meta <- df.meta[
    !is.na(df.meta$met) & !grepl("Analyte|unknown", df.meta$met), 
  ]
  list.keep <- lapply(
    unique(df.meta$met), function(met) {
      df.sub <- df.meta[df.meta$met == met, ]
      if(nrow(df.sub) > 1) {
        df.sub <- df.sub[df.sub$score== max(df.sub$score), ]
      }
      if(nrow(df.sub) > 1) {
        mean.area <- apply(
          df.sub[, grepl("^MEAN", names(df.sub))], 1, mean
        )
        df.sub <- df.sub[which.max(mean.area), ]
      }
      #browser()
      if(nrow(df.sub) != 1) {
        stop("Wrong dedup procedure @ HeatmapFilter")
      }
      return(df.sub)
    }
  )
  df.meta <- do.call(rbind, list.keep)
  return(df.meta)
}
###QC质控####
containQC <- function(data = data){
  dir.create(paste0(
    l.a$output, '/results/', 
    ifelse(is.na(pol), '', paste0(pol, '/')), 
    'QC/')
  )
  tiff(filename = paste0(
    l.a$output, '/results/', 
    ifelse(is.na(pol), '', paste0(pol, '/')), 
    'QC/', "Correlation coefficient graph",".tiff"
  )
  )
  psych::pairs.panels(data[,grep("QC",colnames(data))],method = "pearson", pch=19,
                      hist.col = "#00AFBB",
                      digits= 4,
                      density = TRUE,  
                      ellipses = TRUE )
  dev.off()
}

### Extraction different compound----
ExtractionDiffMet <- function(list.in = list.compare.result){
  diffMet.list <- lapply(
    1:length(l.a$list.compare), function(i) {
      df.deg <- list.in[[i]]
      df.deg <- df.deg[eval(l.a$deg.exp),]
      names(df.deg)[
        grepl("Peak|MS2.name|(MS2 name)|(compound name)|(lipid name)", names(df.deg))
      ] <- "met"
      names(df.deg)[
        grepl("Similarity|MS2.score|(MS2 score)|^score$", names(df.deg))
      ] <- "score"
      df.meta <- df.deg[which(df.deg$met!="NA"),]
      df.group <- l.a$list.compare[[i]]
      m.data <- cbind.data.frame(m.final[,l.a$group[df.group[1]][[1]]],m.final[,l.a$group[df.group[2]][[1]]])
      row.names(m.data) <- df.description$id
      diffMet <- cbind(df.meta, m.data[as.character(df.meta$id), ])
    }
  )
  return(diffMet.list)
}
### Extraction different peak---
ExtractionDiffPeak <- function(list.in = list.compare.result){
  DiffPeak.list <- lapply(
    1:length(l.a$list.compare), function(i) {
      df.deg <- list.in[[i]]
      df.deg <- df.deg[eval(l.a$deg.exp),]
      names(df.deg)[
        grepl("Peak|MS2.name|(MS2 name)|(compound name)|(lipid name)", names(df.deg))
      ] <- "met"
      names(df.deg)[
        grepl("Similarity|MS2.score|(MS2 score)|^score$", names(df.deg))
      ] <- "score"
      if (l.a$type == "QC" | l.a$type == "lipid") {
        df.meta <- subset.data.frame(df.deg,
                                     subset = df.deg$met!="NA", select = c(colnames(df.deg)))
      } else {
        df.meta <- subset.data.frame(df.deg,
                                     subset = df.deg$met!="NA" | df.deg$`MS1 name`!="NA",
                                     select = c(colnames(df.deg))
        )}
      df.group <- l.a$list.compare[[i]]
      df.final <- cbind(df.description$id, m.final)
      m.data <- cbind(m.final[,l.a$group[df.group[1]][[1]]],m.final[,l.a$group[df.group[2]][[1]]])
      DiffPeak <- cbind(df.meta,m.data[df.meta$id,])
    }
  )
  return(DiffPeak.list)
}
### 维恩图----
Venn <- function(comparSN = 1:length(list.compare.result)) {
  library(VennDiagram)
  dir.create(paste0(
    l.a$output, "/results/",
    ifelse(is.na(pol), "", paste0(pol, "/")),
    "Venn/"
  ))
  if (length(list.compare.result) > 5) {
    cat("Venn No drawing, Allow maximum comparison to exceed the limit!")
  } else {
    diffMet.list <- ExtractionDiffMet(list.compare.result)
    diffMet.list <- diffMet.list[comparSN]
    venn.list <- lapply(1:length(diffMet.list), function(i) {
      unique(diffMet.list[[i]]$met)
    })
    names(venn.list) <- c(l.a$subDir[comparSN])
    fill.color <- c("deeppink2", "chartreuse", "cadetblue", "brown3", "goldenrod2")

    venn.diagram(venn.list,
      filename = paste0(
        l.a$output, "/results/",
        ifelse(is.na(pol), "", paste0(pol, "/")),
        "Venn/", "Venn.tiff"
      ),
      cex = 1.4,
      col = "black",
      resolution = 300,
      imagetype = "tiff",
      alpha = 0.50, lwd = 1.2, cat.cex = 1.4,
      fill = fill.color[1:length(diffMet.list)],
      margin = 0.15
    )
    venn.diagram(venn.list,
      filename = paste0(
        l.a$output, "/results/",
        ifelse(is.na(pol), "", paste0(pol, "/")),
        "Venn/", "Venn.png"
      ),
      cex = 0.3,
      col = "black",
      imagetype = "png",
      alpha = 0.50,
      lwd = 0.3, cat.cex = 0.3,
      height = 960,
      width = 960,
      fill = fill.color[1:length(diffMet.list)],
      margin = 0.15
    )
  }
}
### 差异代谢物相关系数热图----
cor_heatmap <- function() {
  library(corrplot)
  wb <- createWorkbook(creator = 'Sensichip')
  modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
  hs <- createStyle(
    fontSize = 10, fontColour = 'white', fgFill = '#4F81BD', 
    valign = 'center', halign = 'center'
  )
  main.list <- ExtractionDiffMet(list.compare.result)
  lapply(1:length(list.compare.result), function(i){
    Sub.df <- main.list[[i]]
    Sub.df <- RepnameFilter(Sub.df)
    row.names(Sub.df) <- Sub.df$met
    if(is.null(Sub.df)) {
      Sub.df <- data.frame(X1 = "No Differentially Expressed Metabolites")
      writeLines("no_result",
                 paste0(
                   l.a$output, '/results/', 
                   ifelse(is.na(pol), '', paste0(pol, '/')), 
                   'cor_heatmap/', 
                   l.a$subDir[i], '/no_result.txt'
                 ))
    } else {
      Sub.df1 <- Sub.df [, -c(1:(l.a$desc+7))]
      cor_mat = cor(t(Sub.df1) ,method = c("pearson"))
      
      fn.jpg <- paste0(
        l.a$output, '/results/', 
        ifelse(is.na(pol), '', paste0(pol, '/')), 
        'cor_heatmap/',l.a$subDir[i], "/",
        l.a$subDir[i], '.png'
      )
      
      fn.pdf <- paste0(
        l.a$output, '/results/', 
        ifelse(is.na(pol), '', paste0(pol, '/')), 
        'cor_heatmap/', l.a$subDir[i], "/",
        l.a$subDir[i], '.pdf'
      )
      
      pdf(file = fn.pdf,
          height =unit(min(24, log10(ncol(cor_mat)) * 12), "cm"),
          width = unit(min(24, log10(ncol(cor_mat)) * 12), "cm"),
          pointsize = min(25, 25 - log10(ncol(cor_mat)))
      )
      
      corrplot(cor_mat, type = "lower", order = "hclust", tl.col = "black", tl.srt = 45, tl.cex = 0.5)
      dev.off()
      
      png(file = fn.jpg,
          height = min(2400, log10(ncol(cor_mat)) * 960),
          width = min(2400, log10(ncol(cor_mat)) * 960),
          pointsize = min(25, 25 - log10(ncol(cor_mat)))
      )
      corrplot(cor_mat, type = "lower", order = "hclust", tl.col = "black", tl.srt = 45, tl.cex = 0.5)
      dev.off()
    }
    addWorksheet(
      wb, sheetName = l.a$subDir[i],
      gridLines = T, tabColour = '#4F81BD'
    )
    
    writeDataTable(
      wb, sheet = i, Sub.df, headerStyle = hs,
      tableStyle = 'TableStyleMedium2', withFilter = T
    )
    setColWidths(wb, sheet = i, cols = 1:ncol(Sub.df), 'auto')
    setRowHeights(wb, sheet = i, rows = 1, 24)
    setRowHeights(wb, sheet = i, rows = 2:(nrow(Sub.df) + 1), 16)
    addStyle(
      wb, sheet = i, rows = 2:(nrow(Sub.df) + 1),
      cols = 1:ncol(Sub.df), gridExpand = T, stack = T,
      style = createStyle(valign = 'center')
    )
  }
  )
  saveWorkbook(
    wb, overwrite = T,
    paste0(
      l.a$output, '/results/', 
      ifelse(is.na(pol), '', paste0(pol, '/')), 
      'cor_heatmap/',
      ifelse(is.na(pol), '', paste0(pol, '-')), 
      'cor_heatmap data matrix.xlsx'
    )
  )
}
# MAIN 主函数----
# 基础分析 ----
biotree <- function(
  # 基础参数
  type = 'GC',  # "GC" or "LC" or "QE" "lipid"
  qc, # integer or NA
  is,  # integer or NA
  group,  # string, eg: '10 10 10'
  sample,  # string, eg: 'exp1 exp2 exp3'
  add.sample = NA, # 正确的写法为："合并后的分组名称~分组1名称+分组2名称.."
  compare,  # string, eg: '1:2 1:3 2:3'
  vsSymbol = '-',# 对比符号 例如：A-B ；支持自定义，例如' VS ',效果为：A VS B
  desc = NA,  # integer or NA
  desc.names = NA,  # string or NA
  # 统计模型参数
  filter = 'sfw',  # string, 'rsd' or 'sfw'
  norm = 'neibiao',  # string, 'neibiao' or 'mianji'
  task = 'basic stat heatmap kegg pathway', # string
  ci = 0.95, # double
  pca3D = T, 
  pls.da = T, # 是否做PLS-DA
  deg.vip = 1, # double
  deg.p = 0.05,  # double
  deg.named = F, 
  keggHitMin = 2, #kegg图表最少命中数 可选 1 或者 2
  palette ='preset2', #画图颜色参数
  ellipse = F, # 分组椭圆
  scaling = 'Ctr UV', # string,pca和PLS-DA/OPLS-DA标准化化方式 可选'None, Ctr, Par, UV'，四取二可重复
  transf = 'F F', # string,pca和PLS-DA/OPLS-DA是否log10转换 eg: 'F T'，填'T/F'或'TRUE/FALSE', 需加引号
  # 聚类分析参数
  heatmap.ignore = T, # 是否去除热图的未知物
  # KEGG分析参数
  keggOrg = "hsa",  # KEGG分析的物种
  maOrg = "hsa", # MA分析的物种
  kegg.deg.only = T,  # 是否在KEGG通路图上只显示差异代谢物
  kegg.pathview = F,  # 要不要用pathview包做KEGG图
  kegg.db = tempdir(),   # 如果用pathview，那么KEGG数据库的本地位置
  kegg.ori = T,  # 要不要做原始形态的KEGG图
  kegg.in = T,  # 是否生成KEGG.in文件夹里的小东西供debug
  # 通路分析参数
  ma.treemap = F, # 通路分析选择
  ma.bubble = T, 
  # 其他参数
  input = paste0("D:/sensichip/"),  # string
  output = paste0('D:/sensichip','/report'),  # string
  out.temp = paste0("D:/sensichip",'/temp'),
  statisticalMethod = 0, # 统计检验方法 0: t test. 1: Wilcoxon rank test 秩和检验
  anova = anova
) {
  # 主函数，做全套~
  # 参数：
  #     qc：质控组个数，本次实验有多少组QC（即混样组），没有就设NA
  #     is：内标序号，内标物的Analyte序号，没有就设NA
  #     group：分组信息，包含分组信息的字符串，每组实验个数按顺序用空格隔开
  #     sample：每组名称，由于分析执行单和compare表中的组别名称可能会有差异，
  #         以分析执行单为准
  #     add.sample：正确的写法为："合并后的分组名称~分组1名称+分组2名称.."
  #     compare：差异比较策略，包含group间两两比较的信息，用冒号连接，空格隔
  #         开，分子在前，分母在后
  #     desc：TODO
  #     desc.names：TODO
  #     filter：过滤法，
  #         'rsd': relative standard deviation, or CV
  #         'sfw': interquartile range, 四分位距
  #         注意：使用rsd法时，qc不能为NA
  #     norm：归一化方法，
  #         'none',不归一化
  #         'neibiao'：利用内标进行归一化
  #         'mianji'：利用TIC进行归一化（存疑）
  #     task：重要，标明需要做哪些分析的参数，必有basic和stat，用空格隔开，
  #         可选的内容包括：basic，stat，heatmap，kegg，pathway，roc和
  #         metscape
  #     type：重要，用了GC还是LC平台
  #     ci：置信区间，默认0.95（即95%）
  #     deg.p：TODO
  #     deg.vip：TODO
  #     deg.named：TODO
  #     input：输入文件
  #     output：输出文件位置
  #     out.temp：SIMCA等中间文件位置
  #     statisticalMethod: 统计检验方法选择
  l.a <<- ArgsVarify( #  确认和整理参数
    type, qc, is, group, sample, add.sample, compare, desc, desc.names, 
    filter, norm, task, ci, pca3D, pls.da, deg.vip, deg.p, deg.named, 
    palette, ellipse, heatmap.ignore, 
    keggOrg, maOrg,kegg.deg.only, kegg.db, kegg.pathview, kegg.ori, kegg.in, 
    ma.treemap, ma.bubble, 
    scaling, transf, input, output, out.temp ,keggHitMin,
    statisticalMethod, vsSymbol ,anova
  )
  BuildProjectFramework()  # 创建项目框架
  df.ori <<- ReadInput()  # 读取原始数据
  if (is.na(pol)) {
    df.area <<- DataReshape(df.ori)  # 原始数据整形
    df.filtered <<- RawFilter(df.area)  # 过滤离群点/不稳定点
    m.recoded <<- RecodeNA(df.area)  # 重编码缺失值
  } else {
    df.area <<- DataReshape(df.ori)
    m.recoded <<- RecodeNA(df.area)
  }
  m.final <<- Normalize(m.recoded)  # 归一化
  #containQC <<- containQC(m.final)  #QC相关系数图
  df.out.mean <<- OutputMean(m.final)  # 输出MEAN表
  df.out.SIMCA <<- OutputSIMCA(m.final)  # SIMCA分析前置格式整理部分
  list.compare.result <<- MVDA(df.out.SIMCA)  # 单变量统计分析
  #list.compare.result.Peak <<- ExtractionDiffPeak(list.compare.result) #选择带有ms1,或ms2 name的peak,合并峰面积输出
  OutputModelInfo(mvda.model)  # 输出模型信息表
  OutputDegPeak(list.compare.result)  # 输出统计表和DEG-PEAK表
#  if (l.a$type == "QE" | l.a$type == "LC") {
    if ("pathway" %in% l.a$task) {
      dir <- file.path(l.a$out.temp,
                       paste0("MA.in", ifelse(is.na(pol), "", pol)))
      if (!dir.exists(dir))
        dir.create(dir, recursive = T)  # 输出MA用的输入文件
      lapply(1:length(l.a$list.compare), function(i) {
        df.deg <- list.compare.result[[i]]
        if (l.a$type == "lipid") {
        inMa <- unique(df.deg[eval(l.a$deg.exp), "lipid name"])
        } else inMa <- unique(df.deg[eval(l.a$deg.exp), "MS2 name"])
        inMa <- gsub("[ ]*$","",inMa)
        inMa <- inMa[!is.na(inMa)]
        if (length(inMa) == 0) {
          inMa <- "No Differentially Expressed Metabolites"
        }
        write.table(
          inMa,
          file = file.path(dir, paste0(l.a$subDir[[i]], ".txt")),
          col.names = F,row.names = F,
          quote = F,
          sep = "\t"
        )
        
      })
    }
#  } else {
    # 老的代码，GC、LC maping到本地数据库fiehn、zhuMetlab在做通路分析和KEGG
    # 通路分析和KEEG的接受文件输入为KEGG ID    
    #    df.map <<- MetaboliteMapping()  # 输出代谢物映射表
    #    if (any("kegg" %in% l.a$task, "pathway" %in% l.a$task)){
    #      list.met <- compareMeetMap(list.compare.result, df.map)
    #      if ("pathway" %in% l.a$task){
    #        MakeMAInput(list.met)  # 输出MA用的输入表格
    #      }
    #      if ("kegg" %in% l.a$task){
    #        KEGGPart1(list.met)
    #      }
    #    }
    
  #}
  
  SaveWorkspace()  # 保存工作空间，给别的脚本用
  ReportAssistant()  # 撰写部分报告
}
# 高级分析 ----
{
  {  # GC部分
    # Heatmap()
    # KEGGColPathway()
    # MetaboAnalyst()
  }
  {  # LC部分
    # pol <<- "POS"
    # load(file.path(l.a$out.temp, paste0(pol, ".l.a.RData")))
    # Heatmap()
    # KEGGColPathway()
    # MetaboAnalyst()
    # pol <<- "NEG"
    # load(file.path(l.a$out.temp, paste0(pol, ".l.a.RData")))
    # Heatmap()
    # KEGGColPathway()
    # MetaboAnalyst()
    
    # pol <<- "POS"
    # load(file.path(l.a$out.temp, paste0(pol, ".l.a.RData")))
    # MetaboAnalyst()
    # pol <<- "NEG"
    # load(file.path(l.a$out.temp, paste0(pol, ".l.a.RData")))
    # MetaboAnalyst()
    # QE
    # MetaboAnalyst()
    # MetaboAnalyst2KEGGCol("POS")
    # pol <<- "POS"
    # load(file.path(l.a$out.temp, paste0(pol, ".l.a.RData")))
    # KEGGPart2()
    
    # pol <<- "NEG"
    # load(file.path(l.a$out.temp, paste0(pol, ".l.a.RData")))
    # KEGGPart2()
    
    # ROC4LC()
  }
  {  # 关联分析
    # MetaboAnalyst2MetScape()  # LC专用，制作MetScape需要的文件和框架
  }
}

# 以下为MA中可供使用的模式动物的'代码'查询表
# =====================================================================
# Mammals
# Homo sapiens (human) [80] "hsa"
# Mus musculus (mouse) [82] "mmu"
# Rattus norvegicus (rat) [81] "rno"
# Bos taurus (cow) [81] "bta"
# Birds
# Gallus gallus (chicken) [78] "gga"
# Fish
# Danio rerio (zebrafish) [81] "dre"
# Insects
# Drosophila melanogaster (fruit fly) [79] "dme"
# Nematodes
# Caenorhabditis elegans (nematode) [78] "cel"
# Fungi
# Saccharomyces cerevisiae (yeast) [65] "sce"
# Plants
# Oryza sativa japonica (Japanese rice) [83] "osa"
# Arabidopsis thaliana (thale cress) [87] "ath"
# Parasites
# Schistosoma mansoni [69] "smm"
# Plasmodium falciparum 3D7 (Malaria) [47] "pfa"
# Trypanosoma brucei [54] "tbr"
# Prokaryotes
# Escherichia coli K-12 MG1655 [87] "eco"
# Bacillus subtilis [80] "bsu"
# Pseudomonas putida KT2440 [89] "ppu"
# Staphylococcus aureus N315 (MRSA/VSSA) [73] "sau"
# Thermotoga maritima [57] "tma"
# Synechococcus elongatus PCC7942 [75] "syf"
# Mesorhizobium loti [86] "mlo"

# 使用示例 ----
# 
# kegg单独绘图--- drawKeggPics(org = 'mmu',inputFile = "~/keggInForLocal.xlsx",outputDir = normalizePath("~"))

# 特别申明 分组名称禁用"-"符号 建议使用"_"代替

# GC--------------------------------------------------------------------------------------
# pol <- NA
# biotree(  # GC
#   type = 'GC',
#   qc = 12, is = 762,
#   group = '6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 5 5 5 3',
#   sample = 'Group1_1 Group1_2 Group1_3 Group1_4 Group2_1 Group2_2 Group2_3 Group2_4 Group3_1 Group3_2 Group3_3 Group3_4 Group4_1 Group4_2 Group4_3 Group4_4 Group5_1 Group5_2 Group5_3 Group5_4',
#   compare = "2:1 3:1 4:1 6:5 7:5 8:5 10:9 11:9 12:9 14:13 15:13 16:13 18:17 19:17 20:17",
#   task = "basic stat heatmap kegg pathway",
#   norm = "mianji",
#   palette = "presetWQ",
#   keggOrg = "mmu",
#   keggHitMin = 2,
#   maOrg = "mmu",
#   kegg.deg.only = T,
#   kegg.ori = F,
#   kegg.pathview = F
# )
# Heatmap()
# KEGGColPathway()
# MetaboAnalyst()

# QE--------------------------------------------------------------------------------------
# pol <- "POS"
# biotree(  # POS
#   type = "QE", qc = NA, is = NA, group = "6 6",
#   desc = 4, desc.names = c("MS2 name", "rt", "mz","MS2 score"),
#   sample = "wdd_T wdd_B", compare = "2:1",
#   keggOrg = "osa",
#   maOrg = "osa",
#   kegg.deg.only = T,
#   keggHitMin = 2,
#   deg.p = 0.05,
#   task = "basic stat heatmap kegg pathway",
#   norm = "mianji"
# )
#Heatmap()
#MetaboAnalyst()
#MetaboAnalyst2KEGGCol(pol) 
#KEGGPart3() #新增模块

# LC--------------------------------------------------------------------------------------
# pol <- "POS"
# biotree(
#   # NEG
#   type = "LC", #QE 请填写 QE , 不能填写为 LC
#   qc = 5,
#   is = NA,
#   group = "9 9 9",
#   # desc = 4, desc.names = c("MS2 name","MS2 score","mz", "rt"),
#   sample = "Treatment Control Model",
#   compare = "3:2 1:2 1:3",
#   task = "basic stat heatmap kegg pathway",
#   filter = "rsd",
#   norm = "mianji",
#   keggOrg = "mmu",
#   maOrg = "mmu", 
#   kegg.deg.only = T
# )

# Heatmap()
# KEGGColPathway()
# MetaboAnalyst() 

# 报告系统 ----
if(0){# 报告系统（暂时只有GC、LC及QE分析可用）----
  save(l.a,pol,file = "~/autoReportTemp.RData")
  repeat{
    newMessage <- readline(
      cat("请更新'",l.a$input,"forAutoReport/headMessage.xlsx'文件中的项目基本信息，完成后输入'#'：\n", sep = ""))
    if(newMessage == "#"){break()}
  }
  source(paste0(l.a$input,"/forAutoReport/oneStepGetFinalHtmlForBiotree.R"))
  file.remove("~/autoReportTemp.RData")
}
# 使用记录 ----

# LC--------------------------------------------------------------------------------------
pol <- "POS" 
biotree(
  type = 'QE',
  qc = 8,
  is = NA,
  group = "51 14 30",
  sample = "PKU MHP Control",
  #add.sample = "F~B+C+D G~C+D",
  compare = "1:3 2:3 1:2",
  task = "basic stat heatmap zscore corrplot kegg pathway", # basic stat heatmap zscore corrplot venn kegg pathway roc anova
  vsSymbol = ' - ',
  desc = 12,  # integer or NA
  desc.names = c('id','MS2 name','Formula','Monoisotopic Molecular Weight','MS2 ppm','MS2 score',
                 'MS1 name','MS1 ppm','KEGG','Molecular Weight','m/z','RT [min]'
),  # string or NA
  filter = "rsd",
  norm = "mianji",
  pls.da = T,
  heatmap.ignore = T,
  palette = "preset2",
  #ellipse = T, 
  scaling = 'UV UV',
  transf = 'F T',
  keggOrg = "hsa",
  maOrg = "hsa",
  kegg.deg.only = T,
  keggHitMin = 2,
  kegg.ori = F,
  ma.treemap = F, # 通路分析马赛克图
  ma.bubble = T, 
  input = paste0("D:/sensichip/"),  # string
  output = paste0('D:/sensichip','/report'),  # string
  out.temp = paste0("D:/sensichip",'/temp'),
  statisticalMethod = 0,
  anova = NA
)
### LC
#KEGGColPathway()
#MetaboAnalyst()
#Heatmap()
#zscore()
#ROC4LC()
### QE
Heatmap()
zscore()
cor_heatmap()
MetaboAnalyst()
MetaboAnalyst2KEGGCol("POS") 
KEGGPart3()
#lipidBubble()
#Venn()
#oneWayANOVA("3QP1 3QP2 3QP3 3QP4")
#ROC4LC()

