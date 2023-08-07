require(rrvgo)
require(GOxploreR)
require(clusterProfiler)

setwd("~/Dropbox/Barcelona/P1/03_BadiRate/Marine_vs_non/GD_model/res_processing/FDR_set/01_GO_Enrichment/Further_analyses/")

PangenomeExtra <- read.delim("~/Dropbox/Barcelona/P1/Pangenome_annotation/Pangenome_COG_Desc_Names_KEGG_mostFreq.txt")
head(PangenomeExtra)

# pangenome <- read.delim("../../../../../../../Pangenome_annotation/Pangenome_GO_most_freqIDs_0.1.txt", h=F)
# head(pangenome)
# colnames(pangenome) <- c("goid", "OG")
# pangenome$goid <- sub(" ","", pangenome$goid)
# category <- go2ont(unique(pangenome$goid))
# MF <- category$go_id[category$Ontology=="MF"]
# pangenomeFMF <- pangenome[pangenome$goid%in%MF,]
# head(pangenomeFMF)
# write.table(pangenomeFMF, file="~/Dropbox/Barcelona/P1/Pangenome_annotation/Pangenome_GO_MF_most_freqIDs_0.1.txt", quote=F, sep="\t", row.names = F)

pangenomeFBP <- read.delim("~/Dropbox/Barcelona/P1/Pangenome_annotation/Pangenome_GO_BP_most_freqIDs_0.1_Terms.txt")
head(pangenomeFBP)

PangenomeFMF <- read.delim("~/Dropbox/Barcelona/P1/Pangenome_annotation/Pangenome_GO_MF_most_freqIDs_0.1.txt")
head(PangenomeFMF)

# Expanded OGs
 expanded <- read.table("../../Rates_significant_FDR_LR_expanded.txt")
 expanded <- row.names(expanded)
 
 expNetRates <- read.table("../../Net_NetDiv_Rates_significant_FDR_LR_expanded_ordered.txt", h = T)
 head(expNetRates)
 dim(expNetRates)
 
 expExtraRates <- cbind.data.frame(expNetRates, PangenomeExtra[row.names(expNetRates),])
 head(expExtraRates)
 # # Repeat enrichment on filtered annotation
 # BKG_ids <- pangenomeFBP[,1:2]
 # colnames(BKG_ids) <- c("go_id","gene_id")
 # head(BKG_ids)
 # BKGterms <- go2term(unique(BKG_ids$go_id))
 # head(BKGterms)
 # res <- enricher(gene = expanded, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
 #                 TERM2GENE = BKG_ids, TERM2NAME = BKGterms, minGSSize = 10, maxGSSize = 1000)
 # res
 # head(res@result, n=30)[,1:6]
 # res@result[1:117,1:6]
 # 
 # # write.table(res@result, "Enrichment_Expanded_GO_BP_Freq0.1.txt", sep="\t", quote=F)
 # # write.table(res@result[which(res@result$p.adjust<0.05),], "Enrichment_Expanded_significant_GO_BP_F0.1.txt", sep="\t", quote=F)
 
 expEnriched <- read.delim("Enrichment_Expanded_significant_GO_BP_F0.1.txt")
 head(expEnriched)
 expEnriched[,1:6]
 dim(expEnriched)
 
 # Compute semantic similarity
  expEnrichSimMat <- calculateSimMatrix(expEnriched$ID, orgdb = "org.Ce.eg.db", ont="BP", method="Wang") #method = c("Resnik", "Lin", "Rel", "Jiang", "Wang")
  scores <- setNames(-log10(expEnriched$p.adjust), expEnriched$ID)

 # Scatter Plot
 library(scico)
 rates <- read.table("../../Net_NetDiv_Rates_significant_FDR_LR_expanded_ordered.txt", h = T)
 head(rates)      
 m <- expEnrichSimMat 
 #df <- res@result[row.names(enrichmentBP_leaf_simMat),]
 res <- expEnriched
 df <- res[row.names(m),]
 head(df)
 #dim(df)
 ogcomposition <- apply(df, 1, function(x){strsplit(x[8], split = "/")[[1]]})
 df <- cbind.data.frame(df, "MedianRate"=unlist(lapply(ogcomposition, function(x){median(rates[x,])})))
 df <- cbind.data.frame(df, "MeanRate"=unlist(lapply(ogcomposition, function(x){mean(rates[x,])})))
 head(df)
 mds <- cmdscale(as.dist(1-m))   
 
 df <- cbind.data.frame(df, mds)
 colnames(df)[9:10] <- c("MDS1", "MDS2")
 head(df)
 
 fastest25 <- df[order(log(df$MeanRate), decreasing = T),][1:25,1]
 enriched25 <- df[order(-log(df$p.adjust), decreasing = T),][!df$ID%in%fastest25,][1:25,1] # 25 most enriched, not included in the previous set
 
 p <- ggplot(df, aes(x= MDS1, y=MDS2)) +
      geom_point(aes(size=-log(p.adjust), color = log(MeanRate)), stroke = 0) +
      scale_colour_scico(alpha=0.7,
                          palette = "lajolla",
                          end = 0.8,
                          limits =c(-5,5),
                          oob = scales::squish) +
      scale_size(range=c(3, 13)) +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            axis.text.y = element_text(size=1),
            legend.key.size = unit(0.1, 'cm'),
            legend.title = element_text(size=2),
            legend.text = element_text(size=2)
             ) +
      ggrepel::geom_text_repel(aes(label = Description, col=NULL), seed = 10,
                                #data = subset(df, -MeanRate > sort(-df$MeanRate, decreasing = F)[length(df$MeanRate) - 20]), #label xx fastest
                                data = df[c(fastest25, enriched25),],
                                box.padding = grid::unit(1, "lines"),
                                size = 4,
                                segment.linetype = 2,
                                segment.curvature = -0.1,
                                segment.ncp = 3,
                                segment.angle = 10,
                                segment.inflect = T,
                                nudge_x = 0.1,
                                max.overlaps = 50,
                                alpha = 0.7) 
 p  
 
 
 pdf("Enriched_BP_GO_F0.1_Expanded_ScatterPlot_MeanRates_Label25Fastest&25Enriched.pdf", h=7, w=14)
 p
 dev.off()  
 
 # Get only terms related to metabolism
  expenrichMet <- expEnriched[grep("metabol", x=expEnriched$Description),]
  tail(expenrichMet)[,1:7]
 # Reduce redundancy using semantic similarity
   expEnrichMetSimMat <- calculateSimMatrix(expenrichMet$ID, orgdb = "org.Ce.eg.db", ont="BP", method="Wang")
   scores <- setNames(-log10(expenrichMet$p.adjust), expenrichMet$ID)
   expEnrichMetSimMatRed <- reduceSimMatrix(expEnrichMetSimMat, scores, threshold = 0.7, orgdb = "org.Ce.eg.db")
 
   pdf("Enriched_Metabolism_GO_Terms_F0.1_Expanded_OGs_Treemap_w_scores_Threshold0.8.pdf", width = 12, height = 8)
    treemapPlot(expEnrichMetSimMatRed, size = "score")
   dev.off()
 
 # DotPlot
   df <- expEnriched[grep("metabol", expEnriched$Description),c(2,6)]
   p1 <- ggplot(df, aes(x = -log(p.adjust), y = reorder(Description,-log(p.adjust)))) +
     geom_point(aes(color=-log(p.adjust)), size=4) +
     scale_colour_continuous(type = "viridis"
                             #end = 0.8,
                             #limits =c(-5,5),
                             #oob = scales::squish
     ) +
     theme_light()  +
     theme(panel.grid.major.x = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text.y = element_text(size=5),
           legend.key.size = unit(0.3, 'cm'),
           legend.title = element_text(size=4),
           legend.text = element_text(size=4))
   
   p1     
   pdf(file="Enriched_Metabolism_GO_Terms_F0.1_Expanded_OGs_DotPlot.pdf", h=3, w=6)
   p1
   dev.off()
   
   # Get only terms related to transport
   expenrichTrans <- expEnriched[grep("transport", x=expEnriched$Description),]
   tail(expenrichTrans)[,1:7]
   # Reduce redundancy using semantic similarity
   expEnrichTransSimMat <- calculateSimMatrix(expenrichTrans$ID, orgdb = "org.Ce.eg.db", ont="BP", method="Wang")
   scores <- setNames(-log10(expenrichTrans$p.adjust), expenrichTrans$ID)
   expEnrichTransSimMatRed <- reduceSimMatrix(expEnrichTransSimMat, scores, threshold = 0.5, orgdb = "org.Ce.eg.db")
   
   pdf("Enriched_Transport_GO_Terms_F0.1_Expanded_OGs_Treemap_w_scores_Threshold0.5.pdf", width = 8, height = 6)
   treemapPlot(expEnrichTransSimMatRed, size = "score")
   dev.off()
   
   # DotPlot
   df <- expEnriched[grep("transport", expEnriched$Description),c(2,6)]
   p1 <- ggplot(df, aes(x = -log(p.adjust), y = reorder(Description,-log(p.adjust)))) +
     geom_point(aes(color=-log(p.adjust)), size=3) +
     scale_colour_continuous(type = "viridis"
                             #end = 0.8,
                             #limits =c(-5,5),
                             #oob = scales::squish
     ) +
     theme_light()  +
     theme(panel.grid.major.x = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text.y = element_text(size=5),
           legend.key.size = unit(0.3, 'cm'),
           legend.title = element_text(size=4),
           legend.text = element_text(size=4))
   
   p1     
   pdf(file="Enriched_Transport_GO_Terms_F0.1_Expanded_OGs_DotPlot.pdf", h=3, w=6)
   p1
   dev.off()

  # OGs annotated to xx Term
  OGsIntAbs <- unlist(str_split(expEnriched[expEnriched$Description=="intestinal absorption",]$geneID, "/"))
  OGsRegTmTrans <- unlist(str_split(expEnriched[expEnriched$Description=="regulation of transmembrane transport",]$geneID, "/"))
  OGsRegIonTmTrans <- unlist(str_split(expEnriched[expEnriched$Description=="regulation of ion transmembrane transport",]$geneID, "/"))
  OGsAnionTrans <- unlist(str_split(expEnriched[expEnriched$Description=="anion transport",]$geneID, "/"))
  OGsRenalWatHom <- unlist(str_split(expEnriched[expEnriched$Description=="renal water homeostasis",]$geneID, "/"))
  OGsCHOMet <- unlist(str_split(expEnriched[expEnriched$Description=="carbohydrate metabolic process",]$geneID, "/"))
  OGsCHOHom <- unlist(str_split(expEnriched[expEnriched$Description=="carbohydrate homeostasis",]$geneID, "/"))
  OGsCellRMetalIon <- unlist(str_split(expEnriched[expEnriched$Description=="cellular response to metal ion",]$geneID, "/"))
  
  #OGs annotated to a word
  OGsTransport <- unique(unlist(str_split(expEnriched[grep("transport",expEnriched$Description),]$geneID, "/")))
  
  # Get Info
  PangenomeExtra[OGsIntAbs,3:6]
  PangenomeExtra[OGsRegTmTrans,3:6]
  PangenomeExtra[OGsRegIonTmTrans,c(2,3,4,5,6,8,9,10)]
  PangenomeExtra[OGsTransport,c(3,5,6)]
  PangenomeExtra[OGsAnionTrans, 3:6]
  PangenomeExtra[OGsRenalWatHom, 3:6]
  PangenomeExtra[OGsCHOMet, 3:6]
  PangenomeExtra[expanded[grep("Secondary",PangenomeExtra[expanded,3])], 3:6]
  
  write.table(PangenomeExtra[OGsRegTmTrans,c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/Reg_transmembrane_transport.txt", quote=F, sep="\t")
  write.table(PangenomeExtra[OGsRegIonTmTrans,c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/Reg_Ion_transmembrane_transport.txt", quote=F, sep="\t")
  write.table(PangenomeExtra[OGsTransport,c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/transport.txt", quote=F, sep="\t")
  
  write.table(PangenomeExtra[OGsCHOMet,c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/CarbohydrateMetabolicProcess.txt", quote=F, sep="\t")
  
  write.table(PangenomeExtra[expanded[grep("Secondary",PangenomeExtra[expanded,3])],c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/SecondayMetabolism.txt", quote=F, sep="\t")
  
  write.table(PangenomeExtra[OGsCHOHom,c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/CarbohydrateHomeostasis.txt", quote=F, sep="\t")
  write.table(PangenomeExtra[OGsCHOHom,c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/CarbohydrateHomeostasis.txt", quote=F, sep="\t")

  ###########################################


# Contracting OGs
contracted <- read.table("../../Net_NetDiv_Rates_significant_FDR_LR_contracted_ordered.txt")
contracted <- row.names(contracted)
length(contracted)

# Enrich
  BKG_ids <- pangenomeFBP[,1:2]
  colnames(BKG_ids) <- c("go_id","gene_id")
  head(BKG_ids)
  BKGterms <- go2term(unique(BKG_ids$go_id))
  head(BKGterms)
  res <- enricher(gene = contracted, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  TERM2GENE = BKG_ids, TERM2NAME = BKGterms, minGSSize = 10, maxGSSize = 2000)
  head(res@result, n=30)[,1:6]
  
  write.table(res@result, file="Enrichment_Contracted_OG_BP_Freq0.1.txt", quote=F, sep="\t")
  
  
  #
  conExtra <- PangenomeExtra[row.names(contracted),]
  conExtra <- cbind.data.frame(contracted, conExtra)
  head(conExtra)  
  write.table(conExtra,"Contracted_ExtraAnnotations.txt", quote=F, sep="\t", row.names = F)
  
###########
  # Parallel expansions transitions
  
  parExp <- read.table("~/Dropbox/Barcelona/P1/03_BadiRate/Marine_vs_non/GD_model/res_processing/FDR_set/04_Reconstructions/Node_analyses/Expansions/Convergently_expanded_OGs_N173_170_147_DROS1.txt")
  parExp <- parExp$V1
  
  # Enrich COG categories
  PangenomeCOG <- read.delim("~/Dropbox/Barcelona/P1/Pangenome_annotation/COG/Pangenome_COG_IDs_ExpandedMultiple_Terms_most_freq.txt")
  BKG_ids <- PangenomeCOG[,2:1]
  colnames(BKG_ids) <- c("go_id","gene_id")
  head(BKG_ids)
  BKGterms <- unique(PangenomeCOG[,2:3])
  BKGterms
  res <- enricher(gene = parExp, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  TERM2GENE = BKG_ids, TERM2NAME = BKGterms, minGSSize = 5, maxGSSize = 15000)
  res
  head(res@result, n=30)[,1:6]
  
  parExpEnriched <- res@result[1:120,]
  
  write.table(res@result, "Enrichment_ParallelExpanded_COG.txt", sep="\t", quote=F)
  #write.table(res@result[which(res@result$p.adjust<0.05),], "Enrichment_ParallelExpanded_significant_GO_BP_F0.1.txt", sep="\t", quote=F)
  
  # Repeat enrichment on filtered annotation
  BKG_ids <- pangenomeFBP[,1:2]
  colnames(BKG_ids) <- c("go_id","gene_id")
  head(BKG_ids)
  BKGterms <- go2term(unique(BKG_ids$go_id))
  head(BKGterms)
  res <- enricher(gene = parExp, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  TERM2GENE = BKG_ids, TERM2NAME = BKGterms, minGSSize = 10, maxGSSize = 500)
  res
  head(res@result, n=30)[,1:6]
  
  parExpEnriched <- res@result[1:120,]

  # write.table(res@result, "Enrichment_ParallelExpanded_GO_BP_Freq0.1.txt", sep="\t", quote=F)
  # write.table(res@result[which(res@result$p.adjust<0.05),], "Enrichment_ParallelExpanded_significant_GO_BP_F0.1.txt", sep="\t", quote=F)
  
  parExpEnriched <- read.delim("~/Dropbox/Barcelona/P1/03_BadiRate/Marine_vs_non/GD_model/res_processing/FDR_set/01_GO_Enrichment/Further_analyses/Enrichment_ParallelExpanded_significant_GO_BP_F0.1.txt")
  head(parExpEnriched)
  
  # Plots
  # Reduce redundancy using semantic similarity
  parExpEnrichSimMat <- calculateSimMatrix(parExpEnriched$ID, orgdb = "org.Ce.eg.db", ont="BP", method="Wang") #method = c("Resnik", "Lin", "Rel", "Jiang", "Wang")
  scores <- setNames(-log10(parExpEnriched$p.adjust), parExpEnriched$ID)
  parExpEnrichSimMatRed <- reduceSimMatrix(parExpEnrichSimMat, scores, threshold = 0.9, orgdb = "org.Ce.eg.db")
  
  pdf(file="Enrichment_ParallelExpanded_Treemap_0.9.pdf", h=15, w=20)
  treemapPlot(parExpEnrichSimMatRed, size = "score")
  dev.off()
  
  # Scatter Plot
  rates <- read.table("../../Net_NetDiv_Rates_significant_FDR_LR_expanded_ordered.txt", h = T)
  head(rates)      
  rates <- rates[parExp,,drop=F]
  m <- parExpEnrichSimMat #Calculate similarity matrix as above
  #df <- res@result[row.names(enrichmentBP_leaf_simMat),]
  res <- parExpEnriched
  df <- res[row.names(m),]
  head(df)
  #dim(df)
  ogcomposition <- apply(df, 1, function(x){strsplit(x[8], split = "/")[[1]]})
  df <- cbind.data.frame(df[, c(1:4,6,9)], "MedianRate"=unlist(lapply(ogcomposition, function(x){median(rates[x,])})))
  df <- cbind.data.frame(df, "MeanRate"=unlist(lapply(ogcomposition, function(x){mean(rates[x,])})))
  head(df)
  mds <- cmdscale(as.dist(1-m))   
  
  df <- cbind.data.frame(df, mds)
  colnames(df)[9:10] <- c("MDS1", "MDS2")
  head(df)
  
  # Get 50 labels, 25 fastest, 25 highest pval
  fastest25 <- df[order(log(df$MeanRate), decreasing = T),][1:25,1]
  enriched25 <- df[order(-log(df$p.adjust), decreasing = T),][!df$ID%in%fastest25,][1:25,1] # 25 most enriched, not included in the previous set
  
  p <- ggplot(df, aes(x= MDS1, y=MDS2)) +
    geom_point(aes(size=-log(p.adjust), color = log(MeanRate)), stroke = 0) +
    scale_colour_scico(alpha=0.7,
                       palette = "lajolla",
                       end = 0.8,
                       limits =c(-5,5),
                       oob = scales::squish) +
    scale_size(range=c(4, 14)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          axis.text.y = element_text(size=1),
          legend.key.size = unit(0.1, 'cm'),
          legend.title = element_text(size=2),
          legend.text = element_text(size=2)
    ) +
    ggrepel::geom_text_repel(aes(label = Description, col=NULL), seed = 10,
                             #data = subset(df, -MeanRate > sort(-df$MeanRate, decreasing = F)[length(df$MeanRate) - 20]), #label xx fastest
                             data = df[c(fastest25, enriched25),],
                             box.padding = grid::unit(1, "lines"),
                             size = 4,
                             segment.linetype = 2,
                             segment.curvature = -0.1,
                             segment.ncp = 3,
                             segment.angle = 10,
                             segment.inflect = T,
                             nudge_x = 0.1,
                             max.overlaps = 50,
                             alpha = 0.7) 
  p  
  
  
  pdf("Enriched_BP_GO_F0.1_ParallelExpanded_ScatterPlot_MeanRates_Label25Fastest&25Enriched.pdf", h=7, w=14)
  p
  dev.off()  
  
  
  # OGs annotated to word or term
  OGsWaterHomeo <- unique(unlist(str_split(parExpEnriched[grep("water homeostasis",parExpEnriched$Description),]$geneID, "/")))
  OGsHyperOsmRes <- unique(unlist(str_split(parExpEnriched[grep("hyperosmotic",parExpEnriched$Description),]$geneID, "/")))
  OGsStarvation <- unique(unlist(str_split(parExpEnriched[grep("starvation",parExpEnriched$Description),]$geneID, "/")))
  OGsCellRespMetIon <- unique(unlist(str_split(parExpEnriched[grep("cellular response to metal ion",parExpEnriched$Description),]$geneID, "/")))
  
  # Get Info
  PangenomeExtra[OGsWaterHomeo,c(2,3,4,5,6,8,9,10)]
  PangenomeExtra[OGsHyperOsmRes,c(2,3,4,5,6,8,9,10)]
  PangenomeExtra[OGsStarvation,c(2,3,4,5,6,8,9,10)]
  
  OGsCellRespMetIon
  
  write.table(PangenomeExtra[OGsStarvation,c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/ParalelExp_CellularResponseToStarvation.txt", quote=F, sep="\t")
  write.table(PangenomeExtra[OGsCellRespMetIon,c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/ParalelExp_CellularResponseToMetalIon.txt", quote=F, sep="\t")
  write.table(PangenomeExtra[unique(unlist(str_split(parExpEnriched[grep("regulation of insulin secretion",parExpEnriched$Description),]$geneID, "/"))),
                             c(2,3,4,5,6,8,9,10)],
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/ParalelExp_RegulationInsulinSecretion.txt", quote=F, sep="\t")

  parExpAnnotExtra <- PangenomeExtra[parExp,c(2,3,4,5,6,8,9,10)]
  parExpAnnotExtra <- parExpAnnotExtra[order(row.names(parExpAnnotExtra)),]
  write.table(parExpAnnotExtra,
              "~/Dropbox/Barcelona/P1/12_FunctionalExploration/ParalelExp_All.txt", quote=F, sep="\t")
  parExpAnnotExtra    
  
  
  # Which annotated to Transport
  
  parExpAnnotExtra[grep("ransp", parExpAnnotExtra$Description),]
  
 ## Innovations
  
  iN148 <- read.table("../../04_Reconstructions/Node_analyses/Innovations/N148_Innovations.txt")$V1
  iN174 <- read.table("../../04_Reconstructions/Node_analyses/Innovations/N174_Innovations.txt")$V1
  iN197 <- read.table("../../04_Reconstructions/Node_analyses/Innovations/N197_Innovations.txt")$V1
  iN180 <- read.table("../../04_Reconstructions/Node_analyses/Innovations/N180_Innovations.txt")$V1
  
  OGs_N148 <- read.table("../../04_Reconstructions/Node_analyses/OGs_Node148.txt")$V1
  OGs_N174 <- read.table("../../04_Reconstructions/Node_analyses/OGs_Node174.txt")$V1
  OGs_N197 <- read.table("../../04_Reconstructions/Node_analyses/OGs_Node197.txt")$V1
  OGs_N180 <- read.table("../../04_Reconstructions/Node_analyses/OGs_Node180.txt")$V1
  
  #Enrich. Background is against individual pseudogenome? I think better against pangenome
  enrich <- function(background, og_list, type){
    BKG_ids <- background[,1:2]
    colnames(BKG_ids)<-c("term_id", "gene_id")
    if (type=="go") { 
      BKGterms <- go2term(unique(BKG_ids$term_id))
    } else if (type=="kegg") { 
      BKGterms <- ko2name(unique(BKG_ids$term_id))
    } else {stop("Invalid database")}
    colnames(BKGterms) <- c("term_id", "Term")  
    res_test <- enricher(gene = og_list, pAdjustMethod = "BH",
                         pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                         TERM2GENE = BKG_ids, TERM2NAME = BKGterms)
    return(res_test)
  }
  
  goEnrich_iN148 <- enrich(pangenomeFBP, iN148, "go")
  goEnrich_iN174 <- enrich(pangenomeFBP, iN174, "go")
  goEnrich_iN197 <- enrich(pangenomeFBP, iN197, "go")
  goEnrich_iN180 <- enrich(pangenomeFBP, iN180, "go")

  PangenomeExtra[iN148,]  
  PangenomeExtra[iN174,]
  PangenomeExtra[iN197,]
  PangenomeExtra[iN180,]
  
  #Enrich all innovations
  goEnrich_iAll <- enrich(pangenomeFBP, c(iN148,iN174,iN180, iN197), "go")
  
  #Plot all GO terms (without enrichment)
  
  iN180_go <- unique(pangenomeFBP[pangenomeFBP$OG%in%iN180, c(1,3)])
  iN148_go <- unique(pangenomeFBP[pangenomeFBP$OG%in%iN148, c(1,3)])
  iN174_go <- unique(pangenomeFBP[pangenomeFBP$OG%in%iN174, c(1,3)])
  iN197_go <- unique(pangenomeFBP[pangenomeFBP$OG%in%iN197, c(1,3)])
  
  iN148_go_p <- prioritizedGOTerms(lst = iN148_go$GO_id, sp = F, domain = "BP", organism = "Caenorhabditis elegans")$HF
  iN174_go_p <- prioritizedGOTerms(lst = iN174_go$GO_id, sp = F, domain = "BP", organism = "Caenorhabditis elegans")$HF
  iN180_go_p <- prioritizedGOTerms(lst = iN180_go$GO_id, sp = F, domain = "BP", organism = "Caenorhabditis elegans")$HF
  iN197_go_p <- prioritizedGOTerms(lst = iN197_go$GO_id, sp = F, domain = "BP", organism = "Caenorhabditis elegans")$HF
  
  
  matBP <- GOSemSim::godata("org.Ce.eg.db", ont = "BP")
  simMat_iN148 <- calculateSimMatrix(iN148_go_p, method = "Wang", orgdb="org.Ce.eg.db")
  simMat_iN174 <- calculateSimMatrix(iN174_go_p, method = "Wang", orgdb="org.Ce.eg.db")
  simMat_iN180 <- calculateSimMatrix(iN180_go_p, method = "Wang", orgdb="org.Ce.eg.db")
  simMat_iN197 <- calculateSimMatrix(iN197_go_p, method = "Wang", orgdb="org.Ce.eg.db")
  
  pdf(file="Treemap_iN180_allGO_reduced_withEnrichmentScores.pdf", h=12, w=16)
  treemapPlot(reduceSimMatrix(simMat_iN148,
                              threshold = 0.9,
                              orgdb="org.Ce.eg.db",
                              #scores = setNames(-log(goEnrich_iN180@result[row.names(simMat_iN180_p),6]),
                               #                 nm = row.names(simMat_iN180_p))
                              
  )
  )
  dev.off()  
  
  