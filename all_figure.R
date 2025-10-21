###FIG1B,C 
human <- readRDS("human.rds") #### OR MOUSE
library("scplotter")
CellDimPlot(human,
            group_by = "celltype", ####split_by
            reduction = "UMAP")

###FS1A,B,C 
CellDimPlot(human, 
            group_by = "SubCellType", 
            reduction = "UMAP",
            theme = "theme_blank")



###FIG1H,I  
gene=c("CD8A","CD3D","NKG7","GNLY","GZMB",
       "CSF3R","S100A12","S100A8","CD14",
       "CD79A","IGKC","MS4A1",
       "PPBP","PF4",
       "HBB","HBA1","HBA2"
)


DotPlot(scRNA_harmony, features = gene) +
  scale_color_gradient(low = "lightblue", high = "red",limits = c(-1, 3)) +  # Custom color gradient
  #scale_size(range = c(2, 8)) +  # Adjust the size of the dots
  theme_minimal() + theme_bw() +  # Use a minimal theme
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45, color = "black", size = 10),  # Adjust axis text color and size
    axis.text.y.left = element_text(color = "black", size = 10),  # Adjust y-axis text color and size
    panel.grid.major = element_blank(),  # Remove major grid lines
    legend.text = element_text(size = 6),  # 调整图例文字字体大小
    legend.title = element_text(size = 6),  # 调整图例标题字体大小
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "right"  # Position the legend on the right
  ) +
  labs(x = "", y = "", title = "")  # Customize axis labels

###############FIG1D


human_data <- readRDS("human.rds")
human_avg <- AverageExpression(human_data,
                               group.by = "cell_type",
                               assays = "RNA")$RNA

hum.data <- readRDS("mouse.rds")
hum_avg <- AverageExpression(hum.data,
                             group.by = "cell_type",
                             assays = "RNA")$RNA

common_genes <- intersect(rownames(human_avg), rownames(hum_avg))

correlation_matrix <- cor(human_avg[common_genes,],
                          hum_avg[common_genes,],
                          method = "spearman")

library(pheatmap)
pheatmap(correlation_matrix,
         display_numbers = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Cell Type Correlation")


######################FIG1E,F

meta=human@meta.data ###OR MOUSE

library(ggplot2)
library(dplyr)

ggplot(meta, aes(x = age, fill = group1)) +
  geom_histogram(binwidth = 5, position = "stack") +
  labs(title = "donor_age",
       x = "Age",
       y = "count") +
  scale_fill_manual(values=sample_color) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white"))


######################FIG1G


#
library(Seurat)
library(ggplot2)
library(dplyr)


##############
sub=human
library(readxl)
sample_info <- read_excel("sample.xlsx")

sample_info$sample=gsub("_","",sample_info$sample)


head(sample_info)


sample_info <- sample_info[, c("sample", "age", "sex","group2")]
rownames(sample_info) <- sample_info$sample

# 将 `age` 和 `sex` 信息添加到 Seurat 对象的 `meta.data` 中
# 将 age 和 sex 信息转换为向量，并添加到 Seurat 对象的 meta.data 中
sub@meta.data$age <- as.vector(sample_info$age[match(sub$orig.ident, sample_info$sample)])
scRNA_harmony@meta.data$sex <- as.vector(sample_info$sex[match(scRNA_harmony$orig.ident, sample_info$sample)])


Idents(sub)=sub@meta.data$orig.ident

sub <- subset(sub, idents = setdiff(levels(Idents(sub)), c("F016","F017","F024")))



Idents(sub)=sub@meta.data$age

####age group for FS1E
group1= subset(sub, subset = age < 35) 
group2 <- subset(sub, subset = age >= 35 & age < 50)
group3= subset(sub, subset = age >= 50) 

#####
exp<- AverageExpression(sub)$RNA #####OR AGE GROUP

exp <- cor(exp, method= "spearman")
pheatmap::pheatmap(exp,cluster_rows = T,cluster_cols = T,)






###################F2B


rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)  
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)


library(future)
plan("multiprocess", workers =12)
options(future.globals.maxSize = 2000 * 1024^2)


scRNA_harmony <- readRDS("neu.rds")

type <- unique(scRNA_harmony$type)


common_types <- type[type %in% unique(scRNA_harmony@meta.data$type[scRNA_harmony@meta.data$group1 == "TWC"]) & 
                       type %in% unique(scRNA_harmony@meta.data$type[scRNA_harmony@meta.data$group1 == "LS"])]


r.deg <- data.frame()


result_list <- list()


for (i in 1:length(common_types)) {
  
  Idents(scRNA_harmony) <- "type"
  
  
  deg <- FindMarkers(scRNA_harmony, ident.1 = "TWC", ident.2 = "LS",
                     group.by = "group1", subset.ident = common_types[i])
  
  
  deg$gene <- rownames(deg)
  
  
  rownames(deg) <- 1:nrow(deg)
  
  
  deg$celltype <- common_types[i]
  deg$unm <- i - 1
  write.csv(deg, file = paste0(common_types[i], 'TWC_LS.csv'))
  
  result_list[[i]] <- deg
}


final_result <- do.call(rbind, result_list)


write.csv(final_result,file="r.deg_TWC_LS.csv",quote=F)

r.deg1=final_result

table(r.deg$unm) 
#############################################################################

r.deg1 <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0)


r.deg1$threshold <- as.factor(ifelse(r.deg1$avg_log2FC > 0 , 'Up', 'Down'))
dim(r.deg1)

r.deg1$adj_p_signi <- as.factor(ifelse(abs(r.deg1$avg_log2FC )> 1  , 'Highly', 'Lowly'))
r.deg1$thr_signi <- paste0(r.deg1$threshold, "_", r.deg1$adj_p_signi)
r.deg1$unm %<>% as.vector(.) %>% as.numeric(.)


##top_up_label <- r.deg1 %>% 
##subset(., threshold%in%"Up") %>% 
##group_by(unm) %>% 
##top_n(n = 5, wt = avg_log2FC) %>% 
##as.data.frame()

top_up_label <- r.deg1 %>% 
  subset(., threshold%in%"Up") %>% 
  group_by(unm) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  as.data.frame()


#top_down_label <- r.deg1 %>% 
#subset(., threshold %in% "Down") %>% 
#group_by(unm) %>% 
#top_n(n = -5, wt = avg_log2FC) %>% 
# as.data.frame()

top_down_label <- r.deg1 %>% 
  subset(., threshold %in% "Down") %>% 
  group_by(unm) %>% 
  top_n(n = -5, wt = avg_log2FC) %>% 
  as.data.frame()

top_label <- rbind(top_up_label,top_down_label)
top_label$thr_signi %<>% 
  factor(., levels = c("Up_Highly","Down_Highly","Up_Lowly","Down_Lowly"))


background_position <- r.deg1 %>%
  dplyr:: group_by(unm) %>%
  dplyr:: summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
  as.data.frame()
## `summarise()` ungrouping output (override with `.groups` argument)
background_position$unm %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$unm - 0.4
background_position$end <- background_position$unm + 0.4

### 
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$unm - 0.5
cluster_bar_position$end <- cluster_bar_position$unm + 0.5
cluster_bar_position$unm %<>% 
  factor(., levels = c(0:max(as.vector(.))))

## 设置填充颜色
cols_thr_signi <- c("Up_Highly" = "#d7301f",
                    "Down_Highly" = "#225ea8",
                    "Up_Lowly" = "black",
                    "Down_Lowly" = "black")
cols_cluster <- c("0" = "#35978f",
                  "1" = "#8dd3c7",
                  "2" = "#ffffb3",
                  "3" = "#b3de69",
                  "4" = "#fb8072",
                  "5" = "#80b1d3",
                  "6" = "#fdb462",
                  "7" = "#bebada")



p= ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                            ymax = Max),
            fill = "#525252", alpha = 0.1) + ###添加灰色背景色
  geom_jitter(data = r.deg1, aes(x =unm, y = avg_log2FC, colour = thr_signi),
              size = 1,position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi) +
  scale_x_continuous(limits = c(-0.5, max(r.deg1$unm) + 0.5),
                     breaks = seq(0, max(r.deg1$unm), 1),
                     label = seq(0, max(r.deg1$unm),1)) + #修改坐标轴显示刻度
  
  # 根据top_label标注基因名
  geom_text_repel(data = top_label, aes(x =unm, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
  
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,
                                             ymax = 0.4, fill = unm), color = "black", alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  labs(x = "Cluster", y = "average log2FC") +
  theme(
    axis.text.y = element_text(colour = 'black', size = 14),
    axis.text.x = element_text(colour = 'black', size = 14, vjust = 25))+
  
  theme_bw()


plot1 <- p + theme(panel.grid.minor = element_blank(), ##去除网格线
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 14),
                   axis.text.x = element_text(colour = 'black', size = 14, vjust = 28.5), #调整x轴坐标,vjust的值按照最终结果稍加调整
                   #panel.border = element_blank(), ## 去掉坐标轴,也就是框里的线
                   axis.ticks.x = element_blank(), ## 去掉的坐标刻度线
                   axis.line.y = element_line(colour = "black")) #添加y轴坐标轴

plot1
ggsave(filename = "TWC_LS.pdf", plot = plot1, width = 5.5, height = 3)


###################F2C

################GO


ggplot(data = GO,
       aes(x = Count, y = Description)) +  ######aes(x = Count, y = reorder(pathway, Count))) +
  geom_point(aes(size = Count, color = -log10(pvalue))) + #气泡大小及颜色设置
  scale_color_distiller(palette = "Spectral",direction = -1) +
  labs(x = "Gene Number",
       y = "",
       title = "Dotplot of Enriched GO Pathways",
       size = "gene number") + #图例名
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))

###################F2E


# 加载所需的包
library(Seurat)
library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 



library(reshape2)

pB2_df <- table(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$type) %>% melt()
##pB2_df <- table(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$Tissue) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")

unique(pB2_df$Cluster)


##调整顺序pB2_df$Sample <- factor(pB2_df$Sample,levels = c("Tumor","Normal"))
sample_color <- col_vector[1:100] 
table(scRNA_harmony@meta.data[["orig.ident"]])
##table(scRNA_harmony@meta.data[["Tissue"]])
pB2 <- ggplot(data = pB2_df, aes(x = Cluster, y = Number, fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB2

pB3 <- ggplot(data = pB2_df, aes(x =Number, y = Cluster , fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB3

pB4 <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8, position="fill") +
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "") +
  theme(axis.text.y = element_text(size = 10, colour = "black")) +
  theme(axis.text.x = element_text(size = 10, colour = "black")) +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.text = element_text(size = 8),  # 调整图例文字字体大小
    legend.title = element_text(size = 8),  # 调整图例标题字体大小
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),  # 调整图例键大小
    legend.spacing.x = unit(0.5, 'cm')  # 调整图例之间的间距
  )


###################F2F



library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)


in.dat<-sub@meta.data
in.dat


R_oe<-calTissueDist(in.dat,
                    byPatient=F,
                    colname.cluster="type",
                    colname.patient="orig.ident",
                    colname.tissue="group1",
                    method="chisq",
                    min.rowSum=0)
R_oe


##############################

library(ComplexHeatmap)
library(circlize)

# Define the function to categorize R_oe values
categorize_roe <- function(value) {
  if (value > 1.5) {
    return("+++")
  } else if (value >  1) {
    return("++")
  } else if (value >= 0.5 && value <= 1) {
    return("+")
  } else if (value > 0 && value < 0.5) {
    return("+/−")
  } else if (value == 0) {
    return("−")
  }
}




# Apply the function to categorize R_oe
R_oe_categorized <- apply(as.matrix(R_oe), c(1, 2), categorize_roe)

# Define the color mapping for the categories
color_mapping <- colorRamp2(
  c(0, 0.5, 1, 2), 
  c("blue", "green", "yellow", "orange")
)

# Create the heatmap
Heatmap(as.matrix(R_oe),
        col = color_mapping,
        show_heatmap_legend=T,
        cluster_rows=T,
        cluster_columns=T,
        row_names_side='right',
        show_column_names=T,
        show_row_names=T,
        
        row_names_gp=gpar(fontsize=10),
        column_names_gp=gpar(fontsize=10),
        heatmap_legend_param=list(
          title="R_o/e",
          at=seq(0, 2, by=0.5),
          labels=seq(0, 2, by=0.5)),
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(R_oe_categorized[i, j], x, y, gp=gpar(fontsize=10))
        }
)



###################F2G

doi: 10.1016/j.exger.2013.11.019. Epub 2013 Dec 4. PMID: 24316038.


###################F2H

doi: 10.1038/s41586-021-04345-x. Epub 2021 Dec 22. PMID: 34937051; PMCID: PMC8828466.


###################F2I

library(limma)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd("/home/xiao_yu/old/merge_all_data/doublefinder/1zhikong/sub_b/score/")
# 选择一个配色方案，比如Set1
my_colors <- brewer.pal(n = 10, name = "Set3")
display.brewer.all()
#########AddModuleScore
upgene=read.table("induce.txt", header=F, sep="\t", check.names=F)

genes_Aucell=upgene$V1

genes_Aucell=genes_Aucell[genes_Aucell%in%rownames(sub)]


library(AUCell)
cells_rankings<-AUCell_buildRankings(sub@assays$RNA@counts,nCores=10,plotStats=TRUE,splitByBlocks=TRUE)##基因排序
cells_AUC<-AUCell_calcAUC(genes_Aucell,cells_rankings,aucMaxRank=nrow(cells_rankings)*0.05)##计算AUC值。
cells_assignment<-AUCell_exploreThresholds(cells_AUC,plotHist=TRUE,nCores=1,assign=TRUE)##挑选阈值
aucs<-getAUC(cells_AUC)
sub$AUC<-aucs[1,]

FeatureStatPlot(sub,features = "AUC", reduction = "UMAP")

############################################
library(CytoTRACE)
library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)

sce=sub
exp1 <- as.matrix(sce@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 3,]
results <- CytoTRACE(exp1,ncores = 10)
phenot <- sce$type
phenot <- as.character(phenot)
names(phenot) <- rownames(sce@meta.data)
emb <- sce@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb)
CytoTRACE=read.table("CytoTRACE_plot_table.txt", header=T, sep="\t", check.names=F, row.names=1)
identical(rownames(CytoTRACE),rownames(sub@meta.data))
sub$CytoTRACE=CytoTRACE$CytoTRACE
library("SCP",lib="/home/xiao_yu/rpackage")

FeatureDimPlot( srt = sub, 
                features = c("CytoTRACE"), 
                reduction = "UMAP"#, theme_use = "theme_blank"
)


#######################################F2J


library(stringr)
library(Seurat)
library(patchwork)
library(SummarizedExperiment)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(SCENIC)
install.packages(qs)
library(qs)
library(BiocParallel)
library(ggrepel)
library(tidyverse)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggnetwork)



register(MulticoreParam(workers = 20, progressbar = TRUE))
loom=open_loom("out_SCENIC.loom")
sce <- sub
regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons = regulonsToGeneLists(regulons_incidMat) 
class(regulons)

regulonAUC = get_regulons_AUC(loom,column.attr.name='RegulonsAUC')

regulonAucThresholds = get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings = get_embeddings(loom)
embeddings

sub_regulonAUC = regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
dim(sub_regulonAUC)

identical(colnames(sub_regulonAUC), colnames(sce))

cellTypes = data.frame(row.names = colnames(sce),
                       celltype = sce$scissor)

selectedResolution = "celltype"
cellsPerGroup = split(rownames(cellTypes),cellTypes[,selectedResolution]) 

sub_regulonAUC = sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]
dim(sub_regulonAUC)

selectedResolution = "celltype"
rss = calcRSS(AUC=getAUC(sub_regulonAUC),
              cellAnnotation=cellTypes[colnames(sub_regulonAUC),selectedResolution]) 
rss=na.omit(rss)

saveRDS(rss,file = "rss.rds")

pdf('RSS_plot.pdf',width = 10, height = 20)
rssPlot = plotRSS(rss,
                  labelsToDiscard = NULL, 
                  zThreshold = 1, 
                  cluster_columns = FALSE, 
                  order_rows = T,
                  thr = 0.01, 
                  varName = "cellType",
                  col.low = '#330066',
                  col.mid = '#66CC66',
                  col.high= '#FFCC33',
                  revCol = F,
                  verbose = TRUE)

rssPlot$plot
dev.off()



#########
for (i in unique(sce$celltype)) {
  pdf(paste0('RSS_plot_',i,'.pdf'),width = 8, height = 8)
  print(plotRSS_oneSet(rss, setName = i) )
  dev.off()
}


############################################ FS2A,B


library(dplyr)
library(Seurat)
library(monocle)

data <- as(as.matrix(sub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data =sub@meta.data)
fData <- data.frame(gene_short_name = row.names(sub), row.names = row.names(sub))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=20, relative_expr = TRUE)
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
saveRDS(mycds,file = "mycds.rds")



# PLOT
plot_cell_trajectory(mycds, color_by = "Pseudotime")|
  plot_cell_trajectory(mycds, color_by = "State")|
  plot_cell_trajectory(mycds, color_by = "group")


##############################################FS2C-E


library(monocle)
library(ggplot2)


pseudotime <- pData(mycds)$Pseudotime
group3 <- pData(mycds)$group  
score <- read.table("/home/docxiao/old/data2/CytoTRACE_plot_table.txt", header=T, sep="\t", check.names=F)  # 
score=score$CytoTRACE
# 
plot_data <- data.frame(
  Pseudotime = pseudotime,
  Score = score,
  Group = group3
)



library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

sample_color <- col_vector[1:100] 

ggplot(plot_data, aes(x = Pseudotime, y = Score, color = Group)) +
  geom_line() +
  scale_color_manual(values = sample_color) +  # 为每个 group3 赋予不同颜色
  labs(x = "Pseudotime", y = "Score") +
  theme_classic() + facet_wrap(~ Group, nrow = 1) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )



# 
pseudotime <- pData(mycds)$Pseudotime
group <- pData(mycds)$type

# 
plot_data <- data.frame(
  Pseudotime = pseudotime,
  Group = group
)

#  ggplot2 
ggplot(plot_data, aes(x = Pseudotime)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "#7FC97F") + 
  labs(x = "Pseudotime", y = "Cell Count") +
  theme_classic() +
  facet_wrap(~ Group, nrow = 1) +  # 
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )



##################
mycds=detectGenes(mycds ,min_expr = 0.1)
head(fData(mycds))
mycds_expressed_genes=row.names(subset(fData(mycds), num_cells_expressed>=10))
mycds_filtered=mycds[mycds_expressed_genes,]
cds_subset=mycds_filtered[gene,]



################## FS2F
plot_genes_in_pseudotime(cds_subset[gene]], 
                         color_by = 'seurat_clusters',ncol = 5)


################## FS2G
plot_pseudotime_heatmap(cds_subset[c(NFKB,HLA,IFIT,S100,CXCL,corgene),],
                        num_clusters = 7,
                        cores = 1,
                        show_rownames = T,
                        cluster_rows = T)

################## F3P

library(ggsankey)
library(tidyverse)
merged_data$Species <- ifelse(merged_data$group1 %in% c("Plain", "Plateau"), "Mouse", "Human")

sankey_data=merged_data@meta.data[,c(10,12,18,20,21)]
sankey_data=meta.data


sankey_data$hcell <- ifelse(sankey_data$group == "Human", sankey_data$celltype, NA)
sankey_data$mcell=ifelse(sankey_data$group == "Mouse", sankey_data$celltype,
                         NA)

sankey_data$group2 <- ifelse(is.na(sankey_data$group2), 
                             "hBlood",
                             paste0("m", sankey_data$group2))


hcell <- na.omit(unique(sankey_data$hcell))
hcell <- unique(sankey_data$hcell)[!is.na(unique(sankey_data$hcell))]
clust=unique(sankey_data$cluster)
mcell <- unique(sankey_data$mcell)[!is.na(unique(sankey_data$mcell))]

sankey_data <- sankey_data %>%  make_long(hcell, cluster, mcell)
sankey_data$node <- factor(sankey_data$node) #,                            levels = rev(c(input, h1[h1 != "Residual"],                                      h2[h2 != "Residual"],                                      h3[h3 != "Residual"],                                      h4[h4 != "Residual" & h4 != "SUMOylation"],                                      h5,h6,outcome)))
sankey_data$next_node <- factor(sankey_data$next_node) #,                            levels = rev(c(h1[h1 != "Residual"],                                          h2[h2 != "Residual"],                                          h3[h3 != "Residual"],                                          h4[h4 != "Residual" & h4 != "SUMOylation"],                                          h5,h6,outcome)))
sankey_data$x <- factor(sankey_data$x) #,                         levels = c("input", "h1", "h2", "h3", "h4", "h5", "h6", "outcome"),                        labels = c("Inputs", "H1", "H2", "H3", "H4", "H5", "H6", "Outcome"))
#sankey_data$node <- factor(sankey_data$node,                            levels = rev(c(input, h1[h1 != "Residual"],                                      h2[h2 != "Residual"],                                      h3[h3 != "Residual"],                                      h4[h4 != "Residual" & h4 != "SUMOylation"],                                      h5,h6,outcome)))
#sankey_data$next_node <- factor(sankey_data$next_node,                            levels = rev(c(h1[h1 != "Residual"],                                          h2[h2 != "Residual"],                                          h3[h3 != "Residual"],                                          h4[h4 != "Residual" & h4 != "SUMOylation"],                                          h5,h6,outcome)))
#sankey_data$x <- factor(sankey_data$x,                         levels = c("input", "h1", "h2", "h3", "h4", "h5", "h6", "outcome"),                        labels = c("Inputs", "H1", "H2", "H3", "H4", "H5", "H6", "Outcome"))
reds <- my_colors

colors <- c(setNames(reds, hcell), setNames(reds, clust), setNames(reds, mcell))

# 
sankey_data_filtered <- sankey_data %>%
  filter(!(is.na(next_x) & is.na(next_node)))


y_range <- range(table(sankey_data_filtered$x))


# 4. 绘图
ggplot(sankey_data,
       aes(x = x,
           next_x = next_x,
           node = node,
           next_node = next_node,
           fill = factor(node))) +  
  geom_sankey(width = 0.15,
              smooth = 5,
              space = 15,
              na.rm = T,
              position = "identity",
              flow.alpha = 0.4,
              node.color = "transparent") +
  geom_sankey_text(aes(label = node),
                   width = 0.2,
                   space = 15,
                   position = "identity",
                   size = 3,
                   color = "black",
                   hjust = 0.1) +
  theme_sankey(base_size = 12) +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.08),
        plot.margin = unit(c(0.1, 0, 0.1, 0), "cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size =12),
        legend.position = "none")


#############################


#############################F3P

library(ggsankey)
library(tidyverse)
merged_data$Species <- ifelse(merged_data$group1 %in% c("Plain", "Plateau"), "Mouse", "Human")

sankey_data=merged_data@meta.data[,c(10,12,18,20,21)]
sankey_data=meta.data


sankey_data$hcell <- ifelse(sankey_data$group == "Human", sankey_data$celltype, NA)
sankey_data$mcell=ifelse(sankey_data$group == "Mouse", sankey_data$celltype,
                         NA)

sankey_data$group2 <- ifelse(is.na(sankey_data$group2), 
                             "hBlood",
                             paste0("m", sankey_data$group2))


hcell <- na.omit(unique(sankey_data$hcell))
hcell <- unique(sankey_data$hcell)[!is.na(unique(sankey_data$hcell))]
clust=unique(sankey_data$cluster)
mcell <- unique(sankey_data$mcell)[!is.na(unique(sankey_data$mcell))]

sankey_data <- sankey_data %>%  make_long(hcell, cluster, mcell)
sankey_data$node <- factor(sankey_data$node) #,                            levels = rev(c(input, h1[h1 != "Residual"],                                      h2[h2 != "Residual"],                                      h3[h3 != "Residual"],                                      h4[h4 != "Residual" & h4 != "SUMOylation"],                                      h5,h6,outcome)))
sankey_data$next_node <- factor(sankey_data$next_node) #,                            levels = rev(c(h1[h1 != "Residual"],                                          h2[h2 != "Residual"],                                          h3[h3 != "Residual"],                                          h4[h4 != "Residual" & h4 != "SUMOylation"],                                          h5,h6,outcome)))
sankey_data$x <- factor(sankey_data$x) #,                         levels = c("input", "h1", "h2", "h3", "h4", "h5", "h6", "outcome"),                        labels = c("Inputs", "H1", "H2", "H3", "H4", "H5", "H6", "Outcome"))
#sankey_data$node <- factor(sankey_data$node,                            levels = rev(c(input, h1[h1 != "Residual"],                                      h2[h2 != "Residual"],                                      h3[h3 != "Residual"],                                      h4[h4 != "Residual" & h4 != "SUMOylation"],                                      h5,h6,outcome)))
#sankey_data$next_node <- factor(sankey_data$next_node,                            levels = rev(c(h1[h1 != "Residual"],                                          h2[h2 != "Residual"],                                          h3[h3 != "Residual"],                                          h4[h4 != "Residual" & h4 != "SUMOylation"],                                          h5,h6,outcome)))
#sankey_data$x <- factor(sankey_data$x,                         levels = c("input", "h1", "h2", "h3", "h4", "h5", "h6", "outcome"),                        labels = c("Inputs", "H1", "H2", "H3", "H4", "H5", "H6", "Outcome"))
reds <- my_colors

colors <- c(setNames(reds, hcell), setNames(reds, clust), setNames(reds, mcell))

# 
sankey_data_filtered <- sankey_data %>%
  filter(!(is.na(next_x) & is.na(next_node)))


y_range <- range(table(sankey_data_filtered$x))


# 4. 绘图
ggplot(sankey_data,
       aes(x = x,
           next_x = next_x,
           node = node,
           next_node = next_node,
           fill = factor(node))) +  
  geom_sankey(width = 0.15,
              smooth = 5,
              space = 15,
              na.rm = T,
              position = "identity",
              flow.alpha = 0.4,
              node.color = "transparent") +
  geom_sankey_text(aes(label = node),
                   width = 0.2,
                   space = 15,
                   position = "identity",
                   size = 3,
                   color = "black",
                   hjust = 0.1) +
  theme_sankey(base_size = 12) +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.08),
        plot.margin = unit(c(0.1, 0, 0.1, 0), "cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size =12),
        legend.position = "none")


##################################### F4A

import scanpy as sc
import hotspot
import numpy as np


adata = sc.read_h5ad("anndata_object.h5ad")

# 1. top5000
sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor='seurat')
adata = adata[:, adata.var.highly_variable]


if 'counts' not in adata.layers:
  adata.layers['counts'] = adata.X.copy()


##################
#
sc.pp.filter_genes(adata, min_cells=10)  # 
sc.pp.filter_cells(adata, min_genes=200)  #


# 
zero_genes = np.sum(adata.X > 0, axis=0).A1 == 0  # 


print(adata.var_names[zero_genes])


#############
hs = hotspot.Hotspot(
  adata,
  layer_key="counts",
  model='danb',
  latent_obsm_key="X_umap",  
  umi_counts_obs_key="nCount_RNA"
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)

hs_results = hs.compute_autocorrelations()

hs_genes = hs_results.loc[hs_results.FDR < 0.05].index

local_correlations = hs.compute_local_correlations(hs_genes, jobs=20)

modules = hs.create_modules(min_gene_threshold=30, core_only=False, fdr_threshold=0.05)

module_scores = hs.calculate_module_scores()

##########

module_scores.to_csv('Hotspot_module_scores_filtered_p0.05_gene30.csv')
modules.to_csv('modules_filtered_p0.05_gene30.csv')
hs.results.to_csv('hs.results_filtered_p0.05_gene30.csv')
local_correlations.to_csv('local_correlations_filtered_p0.05_gene30.csv')

import matplotlib.pyplot as plt
# 首先导入matplotlib
import matplotlib.pyplot as plt

# 然后绘图并保存
fig = hs.plot_local_correlations()
plt.savefig('correlation_modules.png',
            dpi=300,              # 
            bbox_inches='tight',  # 
            format='png',         # 
            transparent=False)    # 

plt.savefig('correlation_modules.pdf',
            dpi=300,              
            bbox_inches='tight',  
            format='pdf',         
            transparent=False)    


####################################F4B
library(viridis)
library(patchwork)
module_scores=read.csv("Hotspot_module_scores_filtered_p0.05_gene30.csv",header = T,row.names = 1)
colnames(module_scores) <- str_c('Module',1:ncol(module_scores),sep = "_")

# 
module_scores_scaled <- scale(module_scores)
module_scores_scaled <- as.data.frame(module_scores_scaled)


sub <- AddMetaData(sub,metadata = module_scores)
plt <- list()
for (i in 1:ncol(module_scores)) {  
  plt[[i]] <- FeaturePlot(sub,
                          features = str_c('Module_',i))+    
    scale_color_viridis(option = 'H')+
    theme(aspect.ratio = 1,
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size = 7),
          panel.background = element_rect(fill = NA,colour = 'black')
    )
}

######################umap

plt[[ncol(module_scores)+1]] <- UMAPPlot(sub)+
  scale_color_d3('category20')+
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 7),
        panel.background = element_rect(fill = NA,colour = 'black'))
patchwork::wrap_plots(plt,ncol = 5)


#################### F4D

library(igraph)
par(mfrow = c(2,5), mar = c(1,1,1,1))

for (i in sort(unique(modules$Module))) {  
  
  current_data <- data_heatmap[subset(modules,Module == i)$X[1:10],
                               subset(modules,Module == i)$X[1:10]]
  
  
  g1 <- graph.adjacency(
    as.matrix(current_data), 
    mode = "undirected",
    weighted = TRUE, 
    diag = FALSE
  )
  
  
  edge_weights <- E(g1)$weight
  
  
  edge_widths <- scales::rescale(abs(edge_weights), to = c(1, 10))
  
  
  layoutCircle <- layout.circle(g1)
  
  
  plot(g1,
       edge.width = edge_widths,  # 
       edge.color = adjustcolor(col = row_group_color[i], alpha.f = 0.25),
       edge.alpha = .5,
       vertex.color = row_group_color[i],
       vertex.label = as.character(colnames(current_data)),
       vertex.label.dist = 2.8,
       vertex.label.degree = -pi/2,
       vertex.label.color = "black",
       vertex.label.family = "Helvetica",
       vertex.label.font = 3,
       vertex.label.cex = 1,
       vertex.frame.color = "black",
       layout = jitter(layoutCircle),
       vertex.size = 20,
       main = str_c('Module ',i),
       cex.main = 1)
}



legend("bottomright", xpd = TRUE,
       title = "Weight",
       legend = c("Strong", "Medium", "Weak"),
       lwd = c(8, 4, 1),
       col = adjustcolor(row_group_color[i], alpha.f = 0.25))

################################ F4E-F

library(WGCNA)

Idents(sub)=sub@meta.data$group1
sub1<- subset(sub, idents = c("TWC","LS"))

data <- as.matrix(GetAssayData(sub1,layer = 'data',assay = 'RNA')[modules$X,])
WGCNA_data <- t(data)
all(rownames(data) == modules$X)
moduleID <- modules$Module
nGenes  <-  ncol(WGCNA_data)
nSamples  <-  nrow(WGCNA_data)

design <- model.matrix(~0+as.factor(sub1$group2))
colnames(design) <- c("35-49","50-70","<35")


MES0 <- moduleEigengenes(WGCNA_data,moduleID)$eigengenes
MEs = orderMEs(MES0)
moduleTraitCor <- cor(MEs,design,use = "p") 
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
saveRDS(moduleTraitPvalue,file = "WGCNA_p.rds")

textMatrix  <-  str_c(signif(moduleTraitCor,2),"(",
                      signif(moduleTraitPvalue,1),")",sep = '')


addSignificance <- function(pvalue) {
  stars <- ""
  if (pvalue < 0.001) stars <- "***"
  else if (pvalue < 0.01) stars <- "**"
  else if (pvalue < 0.05) stars <- "*"
  return(stars)
}


textMatrix <- matrix(NA, nrow=nrow(moduleTraitCor), ncol=ncol(moduleTraitCor))
for(i in 1:nrow(moduleTraitCor)) {
  for(j in 1:ncol(moduleTraitCor)) {
    stars <- addSignificance(moduleTraitPvalue[i,j])
    textMatrix[i,j] <- paste0(signif(moduleTraitCor[i,j], 2), stars)
  }
}



dim(textMatrix) <- dim(moduleTraitCor)
saveRDS(moduleTraitCor,file = "WGCNA_COR.rds")

par(mar = c(4, 4, 3, 2))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = T,
               cex.text = 0.8,
               #zlim = c(-.1,.1),
               main = paste("Module-trait relationships")) 




module_cols <- paste0("Module_", 1:10)


result_matrix <- matrix(nrow = length(unique(sub1$orig.ident)), 
                        ncol = length(module_cols))
rownames(result_matrix) <- unique(sub1$orig.ident)
colnames(result_matrix) <- module_cols


for(sample in unique(sub1$orig.ident)) {
  for(module in module_cols) {
    result_matrix[sample, module] <- mean(sub1@meta.data[sub1$orig.ident == sample, module])
  }
}


result <- data.frame(
  Sample = rownames(result_matrix),
  Top_Module = apply(result_matrix, 1, function(x) colnames(result_matrix)[which.max(x)]),
  Max_Score = apply(result_matrix, 1, max)
)
saveRDS(result,file = "module_result.rds")
saveRDS(result_matrix,file = "module_result_matrix.rds")


print(result)


library(pheatmap)


pheatmap(result_matrix,
         #main = "Module Scores by Sample Groups",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = TRUE,
         number_format = "%.2f")

############################################### F5




library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)
library(dplyr)


cellchat.TWC<- subset(combined_data, idents = c("TWC"))
cellchat.LS<- subset(combined_data, idents = c("LS"))

data.input <- GetAssayData(cellchat.TWC, assay = "RNA", slot = "data")
identity <- subset(cellchat.TWC@meta.data, select = "celltype3")
cellchat.TWC <- createCellChat(object = data.input, meta = identity,  group.by = "celltype3")

data.input <- GetAssayData(cellchat.LS, assay = "RNA", slot = "data")
identity <- subset(cellchat.LS@meta.data, select = "celltype3")
cellchat.LS <- createCellChat(object = data.input, meta = identity,  group.by = "celltype3")


#############

####
CellChatDB <- CellChatDB.human
##
showDatabaseCategory(CellChatDB)

##
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)


unique(CellChatDB$interaction$annotation)
# use Secreted Signaling for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- CellChatDB
cellchat.TWC@DB <- CellChatDB.use # set the used database in the object



cellchat.TWC <- subsetData(cellchat.TWC)

# 
cellchat.TWC <- identifyOverExpressedGenes(cellchat.TWC)
# 
cellchat.TWC <- identifyOverExpressedInteractions(cellchat.TWC)
#cellchat <- projectData(cellchat, PPI.human)
cellchat.TWC <- computeCommunProb(cellchat.TWC, raw.use = TRUE,population.size=F)
###
cellchat.TWC <- filterCommunication(cellchat.TWC, min.cells = 10)
saveRDS(cellchat.TWC,file = "cellchat.TWC.rds")

############################### 2

cellchat.LS@DB <- CellChatDB.use # set the used database in the object



cellchat.LS <- subsetData(cellchat.LS)
future::plan("multiprocess", workers = 20)
# 
cellchat.LS <- identifyOverExpressedGenes(cellchat.LS)
# 
cellchat.LS <- identifyOverExpressedInteractions(cellchat.LS)
#cellchat <- projectData(cellchat, PPI.human)
cellchat.LS <- computeCommunProb(cellchat.LS, raw.use = TRUE,population.size=F)
###
cellchat.LS <- filterCommunication(cellchat.LS, min.cells = 10)
saveRDS(cellchat.LS,file = "cellchat.LS.rds")

############################# 3


df.net <- subsetCommunication(cellchat.o1)

cellchat.TWC=computeCommunProbPathway(cellchat.TWC)
cellchat.TWC <- aggregateNet(cellchat.TWC)
groupSize <- as.numeric(table(cellchat.TWC@idents))
netVisual_circle(cellchat.TWC@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

cellchat.LS=computeCommunProbPathway(cellchat.LS)
cellchat.LS <- aggregateNet(cellchat.LS)




object.list <- list(LS = cellchat.LS,TWC = cellchat.TWC)
saveRDS(object.list,file = "object.list.rds")

cellchat<-mergeCellChat(object.list,add.names=names(object.list),cell.prefix=TRUE)
saveRDS(cellchat,file = "cellchat.rds")

gg1<-compareInteractions(cellchat,show.legend=F,group=c(1,2))
gg2<-compareInteractions(cellchat,show.legend=F,group=c(1,2),measure="weight")
gg1+gg2

par(mfrow=c(1,2),xpd=TRUE)
netVisual_diffInteraction(cellchat,weight.scale=T,comparison=c(1,2))
netVisual_diffInteraction(cellchat,weight.scale=T,measure="weight",comparison=c(1,2))

###########################################  

gg1<-netVisual_heatmap(cellchat,comparison=c(1,2))
gg2<-netVisual_heatmap(cellchat,measure="weight",comparison=c(1,2))
gg1+gg2


vln_data %>% ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+ 
  facet_grid(vln_data $ gene~.,scales = "free_y")+ 
  scale_fill_manual(values = col)+ 
  scale_x_discrete("")+scale_y_continuous("")+ 
  theme_bw()+ 
  theme( 
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    legend.position = "none" )



weight.max<-getMaxWeight(object.list,attribute=c("idents","count"))
par(mfrow=c(1,2),xpd=TRUE)
for(i in 1:length(object.list)){
  netVisual_circle(object.list[[i]]@net$count,weight.scale=T,label.edge=F,edge.weight.max=weight.max[2],edge.width.max=12,title.name=paste0("Number of interactions - ", names(object.list)[i]))
}


###################
num.link<-sapply(object.list2,function(x){rowSums(x@net$count)+colSums(x@net$count)-diag(x@net$count)})
weight.MinMax<-c(min(num.link),max(num.link))# control the dot size in the different datasets
gg<-list()
for(i in 1:length(object.list2)){
  gg[[i]]<-netAnalysis_signalingRole_scatter(object.list2[[i]],title=names(object.list2)[i],weight.MinMax=weight.MinMax)
}
patchwork::wrap_plots(plots=gg)



cellchat.TWC <- netAnalysis_computeCentrality(cellchat.TWC, slot.name = "netP")
cellchat.LS <- netAnalysis_computeCentrality(cellchat.LS, slot.name = "netP")


object.list2 <- list(LS = cellchat.LS,TWC = cellchat.TWC)
saveRDS(object.list2,file = "object.list2.rds")

cellchat2<-mergeCellChat(object.list2,add.names=names(object.list),cell.prefix=TRUE)
######################



######################
install.packages("uwot")
library(uwot)

###################
cellchat2<-computeNetSimilarityPairwise(cellchat2,type="functional",comparison=c(1,2))
#cellchat2<-netEmbedding(cellchat2,type="structural",comparison=c(1,2,3))
cellchat2<- netEmbedding(cellchat2, umap.method = 'uwot',type = "functional",comparison=c(1,2))

trace('netClustering', edit = T, where = asNamespace("CellChat"))  ###多线程问题
cellchat2<-netClustering(cellchat2,type="functional",comparison=c(1,2),do.parallel = F)
netVisual_embeddingPairwise(cellchat2,type="functional",label.size=3.5,comparison=c(1,2))


########################
cellchat2<-computeNetSimilarityPairwise(cellchat2,type="structural",comparison=c(2,3))
#cellchat2<-netEmbedding(cellchat2,type="structural",comparison=c(1,2,3))
cellchat2<- netEmbedding(cellchat2, umap.method = 'uwot',type = "structural",comparison=c(2,3))

trace('netClustering', edit = T, where = asNamespace("CellChat"))  ###多线程问题
cellchat2<-netClustering(cellchat2,type="structural",comparison=c(2,3),do.parallel = F)
netVisual_embeddingPairwise(cellchat2,type="structural",label.size=2.5,comparison=c(2,3),)

saveRDS(cellchat2,file = "cellchat2.rds")


###################

rankSimilarity(cellchat2,type="structural",comparison=c(1,2))

###################
gg1<-rankNet(cellchat2,mode="comparison",stacked=T,do.stat=T,comparison=c(1,2))

gg2<-rankNet(cellchat,mode="comparison",stacked=F,do.stat=TRUE,comparison=c(1,2))
gg1+gg2

###################

library(ComplexHeatmap)
i=1
pathway.union<-union(object.list2[[i]]@netP$pathways,object.list2[[i+1]]@netP$pathways)

##outgoing
ht1=netAnalysis_signalingRole_heatmap(object.list2[[i]],pattern="outgoing",signaling=pathway.union,title=names(object.list2)[i],width=5,height=6)
ht2=netAnalysis_signalingRole_heatmap(object.list2[[i+1]],pattern="outgoing",signaling=pathway.union,title=names(object.list2)[i+1],width=5,height=6)
draw(ht1+ht2,ht_gap=unit(0.5,"cm"))





cellchat2@meta$datasets=factor(cellchat2@meta$datasets,levels=c("TWC","LS"))#set factor level
plotGeneExpression(cellchat2,signaling=c("APRIL"),split.by="datasets",colors.ggplot=T)

plotGeneExpression(cellchat2,signaling=c("TGFb","RESISTIN","VISFATIN"),split.by="datasets",colors.ggplot=T)



plotGeneExpression(cellchat2, signaling =c("CXCL","CCL","IL1","IL2"), type = 'dot', 
                   ,colors.ggplot=T,split.by="datasets") #

genes <- c("CXCL8","CXCR1","CXCR2","CXCL1","CCL5","CCR1","IL1B","IL1R1","IL1R2","IL1RAP","IL18",
           "IL18R1","IL18RAP","IL1A","IL7","IL7R","IL2RG")

combined_data$group3 <- paste(combined_data$celltype3, combined_data$group2, sep = "_")
Idents(combined_data)=combined_data@meta.data$group3
genes <- grep("^(TGF)", rownames(combined_data@assays$RNA@data), value = TRUE)


aaa=DotPlot(combined_data, features = genes) +
  scale_color_gradient(low = "lightblue", high = "red") +  # Custom color gradient  #,limits = c(-1, 3)
  #scale_size(range = c(2, 8)) +  # Adjust the size of the dots
  theme_minimal() + theme_bw() +  # Use a minimal theme
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45, color = "black", size = 10),  # Adjust axis text color and size
    axis.text.y.left = element_text(color = "black", size = 10),  # Adjust y-axis text color and size
    panel.grid.major = element_blank(),  # Remove major grid lines
    legend.text = element_text(size = 6),  # 
    legend.title = element_text(size = 6),  # 
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "right"  # Position the legend on the right
  ) +
  labs(x = "", y = "", title = "")#+
#coord_flip()  # Customize axis labels


cell_types <- unique(gsub("_(TWC|LS)$", "", unique(combined_data@meta.data$group6)))

# 
ordered_levels <- unlist(lapply(cell_types, function(x) {
  paste0(x, c("_TWC", "_LS"))
}))

# 
combined_data@meta.data$group6 <- factor(combined_data@meta.data$group6, 
                                         levels = ordered_levels)

library(dplyr)


library(Seurat)

# 
combined_data@meta.data$group4 <- ifelse(grepl("Neu", combined_data@meta.data$celltype3), "Neu", combined_data@meta.data$celltype3)
combined_data$group5 <- paste(combined_data$group4, combined_data$group2, sep = "_")
combined_data$group6 <- paste(combined_data$celltype3, combined_data$group1, sep = "_")

Idents(combined_data)=combined_data@meta.data$group6
# 
aaa = DotPlot(combined_data, features = genes) +
  scale_color_gradient(low = "lightblue", high = "red") +
  theme_minimal() + theme_bw() +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45, color = "black", size = 10),
    axis.text.y.left = element_text(color = "black", size = 10),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "", y = "", title = "") +
  scale_y_discrete(limits = rev(ordered_levels))  # 

gene2=c("CXCL8","CXCL1","CXCR1","CXCR2","CXCR5")

library(gridExtra)
library(ggplot2)
library(patchwork)

aaa=DotPlot(combined_data, features = genes) +
  scale_color_gradient(low = "lightblue", high = "red") +
  theme_minimal() + theme_bw() +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45, color = "black", size = 10),
    axis.text.y = element_blank(),     # 
    axis.ticks.y = element_blank(),    # 
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "", y = "", title = "") +
  scale_y_discrete(limits = rev(ordered_levels))

mytheme<- theme(
  axis.text.y = element_blank(),
  axis.line.y = element_blank(),
  axis.title = element_blank(),
  #axis.text.x = element_text(size = 14),
  #legend.text = element_text(size = 14),
  #legend.title = element_text(size = 15),
  plot.margin = margin(t = 5.5, # 
                       r= 5.5, # 
                       b= 5.5, # 
                       l= 0) #
)
aaa<- aaa + mytheme
aaa

######################
#Idents(combined_data)=combined_data@meta.data$group1
#combined_data<- subset(combined_data, idents = c("TWC","LS"))

library(Seurat)
library(readxl)

sample_info <- read_excel("sample.xlsx")
sample_info$sample=gsub("_","",sample_info$sample)
combined_data@meta.data$group2 <- as.vector(sample_info$group7[match(combined_data$orig.ident, sample_info$sample)])


avg_exp <- AverageExpression(combined_data, 
                             features =genes,
                             group.by = c("celltype3", "group1"),
                             slot = "data")$RNA


library(tidyr)
library(dplyr)
library(ggplot2)


plot_data <- avg_exp %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  gather(key = "sample", value = "expression", -gene) %>%
  
  separate(sample, into = c("celltype", "condition"), sep = "_")
plot_data$group <- paste(plot_data$celltype, plot_data$condition, sep = "_")

unique(plot_data$condition)
plot_data$condition <- factor(plot_data$condition, 
                              levels = c("TWC","LS"))
plot_data$group <- factor(plot_data$group, 
                          levels = ordered_levels)

pp <- ggplot(data = plot_data,
             aes(x = 0.5, y = group)) +
  geom_tile(aes(fill = condition)) +
  scale_fill_manual(values = c('#27527a','#8ac3d0',"#470d60")) +
  scale_y_discrete(limits = rev)  # 

pp<- pp +
  scale_x_continuous(expand = c(0,0)) +
  xlim(-20,1) + #
  coord_fixed(ratio = 2) #
pp

mytheme2<- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.title = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 6),
  plot.margin = margin(t = 5.5,
                       r= 0, #
                       b= 5.5,
                       l= 5.5))
pp1<- pp + theme_bw() + mytheme2
pp1

###########
# 
plot_data$group2 <- ifelse(!duplicated(plot_data$celltype), plot_data$celltype, NA)

# 
library(dplyr)
data <- data %>%
  group_by(celltype) %>%
  mutate(group2 = if_else(row_number() == 1, celltype, NA_character_)) %>%
  ungroup()
############


pp2 <- pp1 +
  geom_text(data = plot_data,
            aes(x = -0.1,           
                y = group,
                label = group2),
            size = 2.5,
            hjust = 1)

pp2



pp2+ aaa + plot_layout(guides = 'collect')



######################F6E
library('GSEABase')
library(GSVA)
library(msigdbr)
library(Seurat)

set.seed(123)

exp=AverageExpression(sub,group.by = "group3") 
exp=exp[["RNA"]]

counts2=exp
#counts2=exp

GSVA_hall <- gsva(expr=as.matrix(counts2), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, # 
                  kcdf="Gaussian", #
                  parallel.sz=6) # 



GSVA_hall_df <- as.data.frame(GSVA_hall)
GSVA_hall_df$pathway <- rownames(GSVA_hall)
GSVA_hall_df=GSVA_hall_df[,-15]
GSVA_hall_df=GSVA_hall_df[,c("Cord_blood","GZ<35","GZ_35-49","GZ>70",
                             "SH<35","SH_35-49","SH_50-70","SH>70",
                             "LS<35","LS_35-49","LS_50-70",
                             "TWC<35","TWC_35-49","TWC_50-70")]

pheatmap::pheatmap(GSVA_hall_df , #
                   cluster_rows = T,#
                   cluster_cols =F,#
                   gaps_col = gaps_col,
                   show_colnames=T,fontsize_row = 8,fontsize_col = 8,
                   scale = "row", #
                   color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))




