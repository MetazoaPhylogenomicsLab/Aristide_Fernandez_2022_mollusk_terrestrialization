#### Performs enrichment (overrepresentation test) on set of gene families, reduces and plots ####
### It is exemplified here with one of the datasets and it is shown how to reproduce the analyses and figures in the manuscript.
#### Author: Leandro Aristide 2022 ####

library(clusterProfiler)
library(rrvgo)
library(GOxploreR)
library(scico)

#### Create convenience function to perform  enrichment ####
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
  
#### GO #####
    # Load annotations
    pangenome <- read.delim("Pangenome_GO_unique_IDs_and_terms.txt")
    head(pangenome)
  # Keep only Biological Processes Terms
    category <- go2ont(unique(pangenome$goid))
    BP <- category$go_id[category$Ontology=="BP"]
    pangenome_BP <- pangenome[pangenome$goid%in%BP,]
    
  # Load OGs to test
    test <- read.table("FastExpanding.txt", h = F) # e.g. Fast Expanding OGs; Figure 3a in manuscript.
    head(test)
    test <- test$V1
    length(test)

  # Enrich
    res <- enrich(pangenome_BP, test, "go") #Background annotations would change according to the specific enrichment (e.g. pangenome or pseudogenome)
    head(res@result)
    
  # Keep enriched terms (p adjusted < 0.05)
    enriched <- res@result[res@result$p.adjust < 0.05,]
    nrow(enriched)
  
  # As enriched terms are too many for plotting, we use a clustering method to reduce redundancy and summarize
    
    # First we used a more stringent threshold (p adjusted < 0.01)
      enriched <- res@result[res@result$p.adjust < 0.01,]
      nrow(enriched)
    # Compute similarity matrix
      matBP <- GOSemSim::godata("org.Ce.eg.db", ont = "BP")
      enrichment_simMat <- calculateSimMatrix(enriched$ID, orgdb = "org.Ce.eg.db", ont = "BP", method = "Wang", semdata = matBP)
      scores <- setNames(-log(enriched$p.adjust[match(row.names(enrichment_simMat), enriched$ID)]), 
                         row.names(enrichment_simMat))
    # Reduce by clustering
      enrichedBP_clustered <- reduceSimMatrix(enrichment_simMat, threshold = 0.7, orgdb = "org.Ce.eg.db", scores = scores)
      enrichedBP_clustered <- unique(enrichedBP_clustered$parent)

  ## Now plot ###
      m <- enrichment_simMat[enrichedBP_clustered, enrichedBP_clustered]
      df <- res[enrichedBP_clustered,]
      head(df)
  # Scatterplot points proportional to median rate. Calculate median rates for each term.
      rates <- read.table("FastExpanding_NetRates.txt", h = T)
      head(rates)      
      ogcomposition <- apply(df, 1, function(x){strsplit(x[8], split = "/")[[1]]})
      df <- cbind.data.frame(df[, c(1:4,6,9)], "MedianRate" = unlist(lapply(ogcomposition, function(x){median(rates[x, ])})))
      df <- cbind.data.frame(df, "MeanRate" = unlist(lapply(ogcomposition, function(x){mean(rates[x,])})))
      head(df)
  # Multidimensional scaling on similarity matrix
    mds <- cmdscale(as.dist(1-m))   
    df <- cbind.data.frame(df, mds)
    colnames(df)[9:10] <- c("MDS1", "MDS2")
    head(df)
    
  # Create plot
    p <- ggplot(df, aes(x = MDS1, y = MDS2)) +
      geom_point(aes(size = -log(p.adjust), color = log(MeanRate))) +
      scale_colour_scico(alpha = 0.7,
                         palette = "lajolla",
                         end = 0.8,
                         limits = c(-5,5),
                         oob = scales::squish) +
      scale_size(range = c(4, 14)) +
      theme_minimal() +
      theme(panel.grid.major = element_blank()) +
      ggrepel::geom_label_repel(aes(label = Description, col = NULL),
                                data = subset(df, MeanRate > sort(df$MeanRate, decreasing = F)[length(df$MeanRate) - 25]), #label 25 fastest
                                box.padding = grid::unit(2, "lines"),
                                size = 4,
                                segment.linetype = 2,
                                segment.curvature = -0.1,
                                segment.ncp = 3,
                                segment.angle = 20,
                                segment.inflect = T,
                                nudge_x = 0.1,
                                max.overlaps = 50,
                                alpha = 0.7)
    p
    
    pdf("FastExpanding_Reduced_EnrichedTerms_ScatterPlot.pdf", h=6, w=12)
    p
    dev.off()
    
#### KEGG pathways ####
  # Annotations
    pangenome_kegg <- read.delim("Pangenome_KeggPathways_BRITE_uniqueIDs.txt", h = T)
    head(pangenome_kegg)
    pangenome_kegg <- pangenome_kegg[,c(2,1)]
    colnames(pangenome_kegg) <- c("kegg_id", "gene_id")
    
  # Enrich
    res_kegg <- enrich(background = pangenome_kegg, og_list = test, type = "kegg")
    
  # Circular DotPlot
    library(dplyr)
    library(forcats)
    library(paletteer)
    library(stringr)
    
    # Prepare dataframe for plotting, add rates as above
      dfk <- res_kegg[res_kegg$p.adjust < 0.05 & res_kegg$qvalue < 0.2, ]      
      ogcomposition <- apply(dfk, 1, function(x){strsplit(x[8], split = "/")[[1]]})
      dfk <- cbind.data.frame(dfk[, c(1:4,6,9)], "MedianRate" = unlist(lapply(ogcomposition, function(x){median(rates[x, ])})))
      dfk <- cbind.data.frame(dfk, "MeanRate" = unlist(lapply(ogcomposition, function(x){mean(rates[x, ])})))
      dfk <- cbind.data.frame(dfk, "GeneRatioNum" = apply(dfk, 1, function(x){eval(parse(text = x[3]))}))
      head(dfk)
    
    # Add BRITE categories
      brite <- read.delim("Pangenome_KeggPathways_BRITE_uniqueIDs.txt", h=T)
      brite <- unique(brite[brite$KEGG_ID%in%dfk$ID, 2:5])
      row.names(brite) <- brite$KEGG_ID
      dfk <- cbind.data.frame(dfk, brite[dfk$ID,3:4])
      head(dfk)    
      # Assign empty, remove Human Diseases categories, and non-relevant/redundant pathways
      dfk[is.na(dfk$BRITE_1),]$BRITE_1 <- rep("NoCat", times = 5)
      dfk[dfk$BRITE_2== "",]$BRITE_2 <- rep("NoCat", times = 5)
      dfk <- dfk[!dfk$BRITE_1 == "Human Diseases",]
      dfk <- dfk[!dfk$Description == "Mitophagy - yeast",]
      dfk <- dfk[!dfk$Description == "Hippo signaling pathway - multiple species",]
      dfk <- dfk[!dfk$Description == "Longevity regulating pathway - multiple species",]
      dfk <- dfk[!dfk$Description == "Longevity regulating pathway - worm",]
      dfk <- dfk[!dfk$Description == "Circadian rhythm - plant",]
    
    # Define labels angle
      dfk2 <- cbind.data.frame(dfk, "Tr_p.adjust"=-log(dfk$p.adjust),"IDnum" = 1:nrow(dfk)) # Truncate scale to fit plot (see below)
      head(dfk2)
      angle <-  90 - 360 * (dfk2$IDnum - 0.5)/nrow(dfk2)
    # Establish ordering variables
      label_data <- dfk2[order(dfk2$BRITE_1, dfk2$BRITE_2, log(dfk2$MeanRate)),] 
      label_data$IDnum <- 1:nrow(dfk2)
      label_data$hjust <- ifelse( angle < -90, 1, 0)
      label_data$angle <- ifelse(angle < -90, angle+180, angle) #flip angle to make readable
      dfk2 <- dfk2[order(dfk2$BRITE_1, dfk2$BRITE_2, log(dfk2$MeanRate)),] #if ordering by a continuous variable, do it directly in ggplot call
    
    # Define background colors according to BRITE category
      freqs <- as.data.frame(table(dfk2$BRITE_1))
      rects <- data.frame("xstart" = c(1, 7, 18, 19, 49, 54) - 0.5,
                          "xend" = c(7, 18, 19 , 49, 54, 101) - 0.5,
                          col = freqs$Var1)
    # Truncate range of p.values for plotting
      range(dfk2$Tr_p.adjust)
      dfk2 %>% mutate(dfk2, Tr_p.adjust=pmin(Tr_p.adjust, 8)) -> dfk2
      range(dfk2$Tr_p.adjust)
    # Colors for labels
      B12 <- paste0(dfk2$BRITE_1,";", dfk2$BRITE_2)
      B12 <- as.data.frame(table(B12))
      B12 <- cbind.data.frame(str_split(B12$B12,";",simplify = T), B12$Freq)
      colnames(B12) <- c("B1","B2", "Freq")
      pals <- palettes_dynamic_names[palettes_dynamic_names$package == "cartography",]
      pals <- pals$palette[1:n_distinct(B12$B1)]
      B2_freq <- as.data.frame(unlist(lapply(split(B12,B12$B1), nrow)))
      cols_i <- NULL
      for (i in 1:nrow(B2_freq)){
        cols_i  <- c(cols_i,
                     setNames(paletteer_dynamic(paste0("cartography::",pals[i]), n=B2_freq[i,1]), nm = B12[B12[,1]==row.names(B2_freq)[i],2]))
      }
      colors <- cols_i[dfk2$BRITE_2]
    
    #
    p2 <- ggplot(dfk2, aes(y = Tr_p.adjust, x = fct_inorder(Description, NA))) + 
          geom_point(aes(size = GeneRatioNum, color = log(MeanRate))) + 
          labs(y = "Kegg Pathway", x = "Adjusted p-value (-log)", fill = "p.adjust") +
          scale_size(range=c(2, 6)) +
          scale_color_scico(alpha = 0.8,
                            palette = "lajolla",
                            end = 0.8,
                            limits = c(-4, 4),
                            oob = scales::squish
                            ) + 
          theme_minimal() +
          theme(panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_line(size = 0.05, color = "gray50"),
                axis.text.x = element_blank(),
                axis.title = element_blank()
                ) +
      geom_rect(data = rects, inherit.aes = F,
                aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), alpha = 0.4) +
      scale_fill_discrete(name = "BRITE Category") +
      ylim(0, 9) + 
      coord_polar(start = 0, clip = "off") +
      geom_text(data = label_data,
                aes(x= IDnum, y=9, label= Description, hjust=hjust),
                color=colors,
                fontface="bold",
                alpha=0.6,
                size=2.5,
                angle= label_data$angle,
                inherit.aes = FALSE)
    
    p2 
    pdf("FastExpanding_KEGG_EnrichedTerms_CircularDotPlot.pdf", h=10, w=10)
    p2
    dev.off()  
    
    
    
    
