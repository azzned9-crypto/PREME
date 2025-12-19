#cat gene_counts.txt |sed 's/\/home\/yangp\/P48\/STAR\/mRNA_LU_//g'|sed 's/Aligned.sortedByCoord.out.bam//g'|sed 's/"//g'|awk '{printf "%s ",$1; for(i=9;i<=104;i++) printf "%s ", $i; print ""}'|awk '{for(i=1;i<=95;i+=2) printf "%s ",$i;print ""}'|head
#cat gene_counts.txt |sed 's/\/home\/yangp\/P48\/STAR\/mRNA_LU_//g'|sed 's/Aligned.sortedByCoord.out.bam//g'|sed 's/"//g'|awk '{printf "%s ",$1; for(i=8;i<=102;i+=2) printf "%s ",$i; print ""}'|sed '1d' > TumorCounts.txt


library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(viridis)

setwd("D:/工作文件/PRAME转录组")
Counts <- read.table(file = "TumorCounts.txt",row.names=1,header=T)
Info <-read.csv(file = "GroupInfo2.csv",header=T)
Info$ID <-paste0("X",Info$ID)

#4种肿瘤
    #edgeR流程
    	#由于没有21号样本的分组信息,因此排除掉21号样本
    	Counts <- subset(Counts, select = -which(names(Counts) == "X21T"))
    	y <- DGEList(counts = Counts, genes = rownames(Counts))
        keep <- rowSums(cpm(y)>1) >= 1 
        y <- y[keep,,keep.lib.sizes=FALSE]
        y <- calcNormFactors(y, method = "TMM")
        y <- estimateCommonDisp(y, verbose = TRUE)
        y <- estimateTagwiseDisp(y, verbose = TRUE)

        #分组
        InfoSort <- factor(Info$REC[match(colnames(Counts), Info$ID)])
        sample_info <- data.frame(sample = colnames(Counts),group=InfoSort)
        design <- model.matrix(~0 +InfoSort, data = sample_info)
        sample_info$group<-factor(sample_info$group)
        colnames(design) <- levels(sample_info$group)

        #差异表达分析
             #使用GLM广义线性模型
             fit <- glmFit(y, design)
             lrt <- glmLRT(fit, contrast = c(-1, 1))#design文件中，NO和YES分别对应文件中的第1列和第2列,则c(-1,1)是 YES vs NO,而c(1,-1)则是NO vs YES
        #提取差异表达基因
            de_genes <- topTags(lrt, n = Inf)
            de_genes <- as.data.frame(de_genes)
            de_genes$padj <- p.adjust(de_genes$PValue, method = "BH")
            #转换ID
    		gene_symbols <- bitr(de_genes$genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
    		#处理多对一映射(保留第一个映射值）
    		gene_symbols <- gene_symbols[!duplicated(gene_symbols$ENSEMBL), ]
    		#将转换结果合并回原始数据框
    		de_genes <- merge(de_genes, gene_symbols, by.x = "genes", by.y = "ENSEMBL", all.x = TRUE)
            write.csv(de_genes,"DEGs.csv",row.names=T)
         #设置显著性阈值
             sig_threshold <- 0.05
             logFC_threshold <- 1
             significant_genes <- de_genes[de_genes$FDR < sig_threshold & abs(de_genes$logFC) > logFC_threshold, ]
             write.csv(significant_genes,"SIGN.csv",row.names=T)

             #选出上调基因
             UP <- subset(significant_genes,logFC>0)

    #上调基因的KEGG分析
        KEGG_database <- 'hsa'
        gene <- bitr(UP$genes,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = GO_database)
        KEGG<-enrichKEGG(gene$ENTREZID,organism = KEGG_database,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    #提取KEGG结果
      	 keggResults <-KEGG@result
    #筛选结果,p小于0.05,并且按Count取前10个作为展示
      	 Pkegg <-subset(keggResults,pvalue<0.05)
      	 Pkegg <-Pkegg[,c(1,4,10,14)]
      	 Pkegg$logP <- -log10(Pkegg$pvalue)
    	Pkegg <- Pkegg %>% arrange(Count) %>%   slice_head(n = 10) 
    	#将Description转换为有序因子
    	Pkegg$Description <- factor(Pkegg$Description, levels = Pkegg$Description)

    #绘图
    	#绘制柱状图,并使用scale_fill_viridis()着色
    	ggplot(Pkegg, aes(x = Description, y = Count, fill = logP)) +geom_bar(stat = "identity", width = 0.8) + scale_fill_viridis( name = "-log10(P-value)") +coord_flip() + theme_bw() + theme(panel.grid = element_blank()) +labs( title = "KEGG Pathway Analysis",x = "", y = "Gene Count")

#只有腺癌28例(REC 18 vs NON 10)
    #载入文件
        Counts <- read.table(file = "LUADCounts.txt",row.names=1,header=T)
        Info <-read.csv(file = "LUADInfo.csv",header=T)
        Info$ID <-paste0("X",Info$ID)
    #edgR
        y <- DGEList(counts = Counts, genes = rownames(Counts))
        keep <- rowSums(cpm(y)>1) >= 1 
        y <- y[keep,,keep.lib.sizes=FALSE]
        y <- calcNormFactors(y, method = "TMM")
        y <- estimateCommonDisp(y, verbose = TRUE)
        y <- estimateTagwiseDisp(y, verbose = TRUE)
        InfoSort <- factor(Info$REC[match(colnames(Counts), Info$ID)])
        sample_info <- data.frame(sample = colnames(Counts),group=InfoSort)
        design <- model.matrix(~0 +InfoSort, data = sample_info)
        sample_info$group<-factor(sample_info$group)
        colnames(design) <- levels(sample_info$group)
        fit <- glmFit(y, design)
        lrt <- glmLRT(fit, contrast = c(-1, 1))
    #提取差异表达基因
        de_genes <- topTags(lrt, n = Inf)
        de_genes <- as.data.frame(de_genes)
        de_genes$padj <- p.adjust(de_genes$PValue, method = "BH")
    #转换ID
        gene_symbols <- bitr(de_genes$genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
    #处理多对一映射(保留第一个映射值）
        gene_symbols <- gene_symbols[!duplicated(gene_symbols$ENSEMBL), ]
    #将转换结果合并回原始数据框
        de_genes <- merge(de_genes, gene_symbols, by.x = "genes", by.y = "ENSEMBL", all.x = TRUE)
        write.csv(de_genes,"DEGs.csv",row.names=T)
    #设置显著性阈值
        sig_threshold <- 0.05
        logFC_threshold <- 1
        significant_genes <- de_genes[de_genes$FDR < sig_threshold & abs(de_genes$logFC) > logFC_threshold, ]
        write.csv(significant_genes,"SIGN.csv",row.names=T)