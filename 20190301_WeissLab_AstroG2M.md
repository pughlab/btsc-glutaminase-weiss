
----
# Weiss Lab Single Cell Collaboration


L.Richards  
March 2019  

---  

Collaborators: Ian Restall, Weiss Lab, University of Calgary  
Analysis Dir: /mnt/work1/users/pughlab/projects/Weiss_Astro_scRNAseq  



----
# 1.0 Introduction
---
Cultures resistant to glutamine inhibitor ...
- express higher levels of astrocyte markers
- express higher levels of glutamate transporters, EAAT1/2

**Anlaysis Plan:**
- score cells for astrocyte signature
- score cells for oligodendrocyte signature
- score cells for stem signature
- score cells for cell cycle phase
- plot frequency of all these cells across samples (dot plot)
- plot expression of EAAT1/2 for +/- of each category
- plot the different scores against eachother
- do not need to normalize the cell lines all together, as a first pass, perform the scoring on each sample separately

**Objective:**   
Compare all these gene signature scores to see if any phenotypes correlate with in vitro response to a gluatminase inhibitor. 

---
# 2.0 Cohort
---

We will be using all the BT lines from the SU2C scRNAseq sequencing effort:
> BT127_L  
> BT14_L  
> BT48_L  
> BT67_L  
> BT73_L  
> BT84_L  
> BT89_L  
> BT94_L  

I have previously QC'ed and optimized clustering for these samples.

---
# 3.0 Format gene signatures
---

Combine [Ian's signatures](https://docs.google.com/spreadsheets/d/1GIrYQaDiJKW5Bb1-J5ZiPuc-xreA6lU2WIr3sdU0ogM/edit?usp=sharing) with cell cycle signatures from Regev Lab.


```R
cc.sig <- paste0(dat.dir, "regev_lab_cell_cycle_genes.txt")
cc.genes <- readLines(con = cc.sig)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
```


```R
gene.sigs <- read.csv("~/Desktop/Samwise/projects/Weiss_Astro_scRNAseq/data/WeissLab_GeneSignatures_Dec2018.csv",
                     stringsAsFactors = F)
```


```R
gene.list <- list(gene.sigs[,1][gene.sigs[,1] != ""],
                  gene.sigs[,2][gene.sigs[,2] != ""],
                  gene.sigs[,3][gene.sigs[,3] != ""],
                  gene.sigs[,4][gene.sigs[,4] != ""],
                  gene.sigs[,5][gene.sigs[,5] != ""],
                  gene.sigs[,6][gene.sigs[,6] != ""],
                  gene.sigs[,7][gene.sigs[,7] != ""],
                  gene.sigs[,8][gene.sigs[,8] != ""],
                  gene.sigs[,9][gene.sigs[,9] != ""],
                  s.genes,
                  g2m.genes
                 )



names(gene.list) <- c(colnames(gene.sigs), "s.genes", "g2m.genes")
str(gene.list)

gene.list

save(gene.list, file = "~/Desktop/Samwise/projects/Weiss_Astro_scRNAseq/data/Weiss_GeneSignatures.RData")
```

    List of 11
     $ Suva_Astro        : chr [1:65] "APOE" "SPARCL1" "ALDOC" "CLU" ...
     $ Suva_Astro_noEAATs: chr [1:63] "APOE" "SPARCL1" "ALDOC" "CLU" ...
     $ Suva_Stemness     : chr [1:63] "SOX4" "CCND2" "SOX11" "RBM6" ...
     $ Suva_G1S          : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
     $ Suva_G2M          : chr [1:55] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
     $ Suva_Oligo        : chr [1:35] "OLIG1" "SNX22" "GPR17" "DLL3" ...
     $ GO_Astro_Diff     : chr [1:39] "TSPAN2" "CDK6" "LAMC3" "CNTF" ...
     $ Housekeeping      : chr [1:98] "ACTB" "B2M" "HNRPLL" "HPRT" ...
     $ Genes_of_Interest : chr [1:3] "EAAT12" "SLC1A2" "SLC1A3"
     $ s.genes           : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
     $ g2m.genes         : chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...



```R
gene.list[[1]][ ! gene.list[[1]] %in% rownames(eb1S@data)]
"SLC1A1" %in% rownames(eb1S@data)
```






TRUE


---
# 4.0 Score and visualize signatures
---

Results from this section previously sent to Ian in January.


```R
suppressMessages(library(Seurat))
suppressMessages(library(RColorBrewer))

plot.cor <- function(x, y, title, xtitle, ytitle){
    
    p.cor <- cor(x, y,method = "pearson")
    a <- cor.test(x, y)
    
    plot(x, 
         y,
         cex = 01,
         xlab = xtitle,
         ylab = ytitle,
         main = title,
         col = "black"
     )
    
    abline(lm(y~x), 
       col="red",
       lty = 1,
       lwd = 2
      )
    
    legend("top",
      c(paste0("r = ", round(p.cor, 2), ", p =", round(a$p.value, 3))),
       bty = "n",
      cex = 1
      )
}
```


```R
setwd("~/pughlab/projects/Weiss_Astro_scRNAseq/scoring/")


dat.dir <- "~/pughlab/projects/Weiss_Astro_scRNAseq/data/"
files <- list.files(dat.dir, pattern = "_L_res.")
files

samples <- gsub('.{14}$', '', files)
samples

gene.sig.file <- "~/pughlab/projects/Weiss_Astro_scRNAseq/data/Weiss_GeneSignatures.RData"
```


<ol class=list-inline>
	<li>'BT127_L_res.0.1.RData'</li>
	<li>'BT147_L_res.0.4.RData'</li>
	<li>'BT48_L_res.0.1.RData'</li>
	<li>'BT67_L_res.0.2.RData'</li>
	<li>'BT73_L_res.0.1.RData'</li>
	<li>'BT84_L_res.0.1.RData'</li>
	<li>'BT89_L_res.0.3.RData'</li>
	<li>'BT94_L_res.0.2.RData'</li>
</ol>




<ol class=list-inline>
	<li>'BT127_L'</li>
	<li>'BT147_L'</li>
	<li>'BT48_L'</li>
	<li>'BT67_L'</li>
	<li>'BT73_L'</li>
	<li>'BT84_L'</li>
	<li>'BT89_L'</li>
	<li>'BT94_L'</li>
</ol>




```R
for (i in 7:8){
    
    print("")
    print("")
    print("###################")
    print(samples[i])
    print("###################")
    print("")
    print("")
    
    load.file <- paste0(dat.dir, files[i])
    print(load.file)
    load(load.file)
    colnames(eb1S@meta.data)[grep("res", colnames(eb1S@meta.data))] <- "Cluster"

    print("Filter Gene Signatures (Only include genes found in data matrix)")
    print("Below are genes NOT found in data:")
    load(gene.sig.file)

        for (j in 1:length(gene.list)){
    
            print("")
            print(names(gene.list)[j])
            aa <- gene.list[[j]][ ! gene.list[[j]] %in% rownames(eb1S@data)]
            print(aa)
            print("")  
            gene.list[[j]] <- gene.list[[j]][gene.list[[j]] %in% rownames(eb1S@data)]
    
        }
    
    print("Final Gene List Lengths for Scoring")
    print(str(gene.list))

    print("")
    print("")
    print("Scoring Cell Cycle Phase")
    print("")
    print("")

    eb1S <- CellCycleScoring(object = eb1S, g2m.genes = gene.list$g2m.genes, s.genes = gene.list$s.genes)
    
    print("Scoring Gene Signatures")
    print("")
    print("")

    eb1S <- AddModuleScore(eb1S,
                           genes.list = gene.list[1:9],
                           n.bin = 25,
                           seed.use = 1,
                           ctrl.size = 100,
                           random.seed = 123,
                           enrich.name = names(gene.list)[1:9]
                      )


    MKI67 <- data.frame(eb1S@data[rownames(eb1S@data) == "MKI67", ])
    colnames(MKI67) <- "MKI67"

    #SLC1A2 <- data.frame(eb1S@data[rownames(eb1S@data) == "SLC1A2", ])
    #colnames(SLC1A2) <- "SLC1A2"

    SLC1A2 <- data.frame(colnames(eb1S@data), 0)
    rownames(SLC1A2) <- colnames(eb1S@data)
    SLC1A2[,1] <- NULL
    colnames(SLC1A2) <- "SLC1A2"
    SLC1A2

    SLC1A3 <- data.frame(eb1S@data[rownames(eb1S@data) == "SLC1A3", ])
    colnames(SLC1A3) <- "SLC1A3"

    dat <- cbind(eb1S@dr$tsne@cell.embeddings, eb1S@meta.data, MKI67, SLC1A2, SLC1A3 )

    #print cluster tSNE for sample
    #colour by cluster

    plot.title <- paste(samples[i], "     ", nrow(dat), "cells")
    #print(plot.title)

    sample_tSNE <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=Cluster)) +
                   geom_point(alpha = 0.7, size = 2) +
                   labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
                   scale_colour_brewer(palette = "Dark2") +
                   theme_bw() +
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())

    cell.cycle_tSNE <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=Phase)) +
                   geom_point(alpha = 0.7, size = 2,) +
                   labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
                   scale_color_manual(values=c("#999999", "darkred", "#E69F00")) +
                   theme_bw() +
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())
    
    SLC1A2 <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=SLC1A2)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
         scale_colour_gradientn(colours = c("grey", "darkblue")) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())

    SLC1A3 <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=SLC1A3)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
         scale_colour_gradientn(colours = c("grey", "darkblue")) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())

    MKI67 <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=MKI67)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
         scale_colour_gradientn(colours = c("grey", "darkblue")) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())


p <- dat$Suva_Astro1
upper <- max(abs(max(p)), abs(min(p)))
lower <- upper * -1

Astro <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=Suva_Astro1)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
       
         scale_colour_gradientn(colours = rev(brewer.pal(10, "Spectral")), limits = c(lower, upper)) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())



p <- dat$Suva_Astro_noEAATs2
upper <- max(abs(max(p)), abs(min(p)))
lower <- upper * -1

Astro_no <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=Suva_Astro_noEAATs2)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
       
         scale_colour_gradientn(colours = rev(brewer.pal(10, "Spectral")), limits = c(lower, upper)) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())


p <- dat$Suva_Oligo6
upper <- max(abs(max(p)), abs(min(p)))
lower <- upper * -1

Oligo <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=Suva_Oligo6)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
       
         scale_colour_gradientn(colours = rev(brewer.pal(10, "Spectral")), limits = c(lower, upper)) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())

p <- dat$Suva_Stemness3
upper <- max(abs(max(p)), abs(min(p)))
lower <- upper * -1

Stemness <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=Suva_Stemness3)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
       
         scale_colour_gradientn(colours = rev(brewer.pal(10, "Spectral")), limits = c(lower, upper)) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())


p <- dat$GO_Astro_Diff7
upper <- max(abs(max(p)), abs(min(p)))
lower <- upper * -1

Astro_diff <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=GO_Astro_Diff7)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
       
         scale_colour_gradientn(colours = rev(brewer.pal(10, "Spectral")), limits = c(lower, upper)) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())


p <- dat$Genes_of_Interest9
upper <- max(abs(max(p)), abs(min(p)))
lower <- upper * -1

GenesI <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=Genes_of_Interest9)) + 
        geom_point(alpha = 0.6, size = 1.5) +  
       
         scale_colour_gradientn(colours = rev(brewer.pal(10, "Spectral")), limits = c(lower, upper)) +
        labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
        theme_bw() + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())



    


print("")
print("")
print("Saving data")

new.dat.name <- paste0(dat.dir, samples[i], "_WeissScoring.Rdata")
print(new.dat.name)
save(dat, file = new.dat.name)

pdf.name <- paste0(samples[i], "_Scoring_WeissCollab.pdf") 
print("")
print("")
print("Plotting")
print(pdf.name)

pdf(pdf.name)

print(sample_tSNE)
print(cell.cycle_tSNE)
print(SLC1A2)
print(SLC1A3)
print(MKI67)
print(Astro)
print(Astro_no)
print(Astro_diff)
print(Stemness)
print(Oligo)
print(GenesI)



plot.cor(dat$SLC1A2, dat$SLC1A3, plot.title, "SLC1A2 Expression", "SLC1A3 Expression")
plot.cor(dat$SLC1A2, dat$MKI67, plot.title, "SLC1A2 Expression", "MKI67 Expression")
plot.cor(dat$SLC1A3, dat$MKI67, plot.title, "SLC1A3 Expression", "MKI67 Expression")

plot.cor(dat$Suva_Astro1, 
         dat$G2M.Score, 
         plot.title, 
         "Suva Astro Signature Score", 
         "G2M Signature Score"
        )
plot.cor(dat$Suva_Astro1, 
         dat$Suva_Stemness3, 
         plot.title, 
         "Suva Astro Signature Score", 
         "Suva Stemness Signature Score"
        )

plot.cor(dat$Suva_Astro1, 
         dat$Suva_Stemness3, 
         plot.title, 
         "Suva Astro Signature Score", 
         "Suva Stemness Signature Score"
        )


plot.cor(dat$Suva_Astro1, 
         dat$Suva_Stemness3, 
         plot.title, 
         "Suva Astro Signature Score", 
         "Suva Stemness Signature Score"
        )

plot.cor(dat$G2M.Score, 
         dat$Suva_Stemness3, 
         plot.title, 
         "G2M Signature Score", 
         "Suva Stemness Signature Score"
        )

plot.cor(dat$G2M.Score, 
         dat$Genes_of_Interest9, 
         plot.title, 
         "G2M Signature Score", 
         "SLC1A2/3 Score"
        )


plot.cor(dat$Suva_Stemness3, 
         dat$Genes_of_Interest9, 
         plot.title, 
         "Suva Stemness Signature Score", 
         "SLC1A2/3 Score"
        )

plot.cor(dat$Suva_Astro1, 
         dat$Genes_of_Interest9, 
         plot.title, 
         "Suva Astrocyte Signature Score", 
         "SLC1A2/3 Score"
        )


dev.off()
 
    
    
#}


```

    [1] ""
    [1] ""
    [1] "###################"
    [1] "BT94_L"
    [1] "###################"
    [1] ""
    [1] ""
    [1] "~/Desktop/Samwise/projects/Weiss_Astro_scRNAseq/data/BT94_L_res.0.2.RData"
    [1] "Filter Gene Signatures (Only include genes found in data matrix)"
    [1] "Below are genes NOT found in data:"
    [1] ""
    [1] "Suva_Astro"
    [1] "GABRB1" "DOK5"   "SLC1A2" "LIX1"   "NOG"    "HSPB8" 
    [1] ""
    [1] ""
    [1] "Suva_Astro_noEAATs"
    [1] "GABRB1" "DOK5"   "LIX1"   "NOG"    "HSPB8" 
    [1] ""
    [1] ""
    [1] "Suva_Stemness"
     [1] "SPDYE7P"  "SPDYE1"   "NCRUPAR"  "NELL2"    "CCL5"     "EVI2A"   
     [7] "LYZ"      "POU5F1"   "MBOAT1"   "LOC90834" "SPDYE5"  
    [1] ""
    [1] ""
    [1] "Suva_G1S"
    [1] "MLF1IP"
    [1] ""
    [1] ""
    [1] "Suva_G2M"
    character(0)
    [1] ""
    [1] ""
    [1] "Suva_Oligo"
    [1] "NEU4"  "OMG"   "LIMS2"
    [1] ""
    [1] ""
    [1] "GO_Astro_Diff"
     [1] "TSPAN2" "LAMC3"  "CNTF"   "DRD1"   "LIF"    "S100A8" "S100A9" "TAL1"  
     [9] "GCM1"   "PPAP2B"
    [1] ""
    [1] ""
    [1] "Housekeeping"
    [1] "HNRPLL"  "HPRT"    "PRPS1L1" "PRPS1L3" "RPL10L"  "RPL3L"   "RPS6KA6"
    [1] ""
    [1] ""
    [1] "Genes_of_Interest"
    [1] "EAAT12" "SLC1A2"
    [1] ""
    [1] ""
    [1] "s.genes"
    [1] "MLF1IP"
    [1] ""
    [1] ""
    [1] "g2m.genes"
    character(0)
    [1] ""
    [1] "Final Gene List Lengths for Scoring"
    List of 11
     $ Suva_Astro        : chr [1:59] "APOE" "SPARCL1" "ALDOC" "CLU" ...
     $ Suva_Astro_noEAATs: chr [1:58] "APOE" "SPARCL1" "ALDOC" "CLU" ...
     $ Suva_Stemness     : chr [1:52] "SOX4" "CCND2" "SOX11" "RBM6" ...
     $ Suva_G1S          : chr [1:42] "MCM5" "PCNA" "TYMS" "FEN1" ...
     $ Suva_G2M          : chr [1:55] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
     $ Suva_Oligo        : chr [1:32] "OLIG1" "SNX22" "GPR17" "DLL3" ...
     $ GO_Astro_Diff     : chr [1:29] "CDK6" "ABL1" "GFAP" "DLL1" ...
     $ Housekeeping      : chr [1:91] "ACTB" "B2M" "PSMB2" "PSMB4" ...
     $ Genes_of_Interest : chr "SLC1A3"
     $ s.genes           : chr [1:42] "MCM5" "PCNA" "TYMS" "FEN1" ...
     $ g2m.genes         : chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
    NULL
    [1] ""
    [1] ""
    [1] "Scoring Cell Cycle Phase"
    [1] ""
    [1] ""
    [1] "Scoring Gene Signatures"
    [1] ""
    [1] ""



    Error in as.data.frame.default(x[[i]], optional = TRUE): cannot coerce class "structure("dgCMatrix", package = "Matrix")" to a data.frame
    Traceback:


    1. data.frame(eb1S@data[rownames(eb1S@data) == "SLC1A2", ])

    2. as.data.frame(x[[i]], optional = TRUE)

    3. as.data.frame.default(x[[i]], optional = TRUE)

    4. stop(gettextf("cannot coerce class \"%s\" to a data.frame", deparse(class(x))), 
     .     domain = NA)


---
# 5.0 Classify cells for gene signatures 
---

Ian is interested in correlating the proportion of cycling, astrocytic cells to drug response. 

Anlaysis:  
(1) Classify cells as astrocytic  
(2) Classify cells as cycling  
(3) Identify double positive cells, classified as being both astrocytes and G2M phase.   




```R
suppressMessages(library(Seurat))
```


<ol class=list-inline>
	<li>'BT127_L_res.0.1.RData'</li>
	<li>'BT147_L_res.0.4.RData'</li>
	<li>'BT48_L_res.0.1.RData'</li>
	<li>'BT67_L_res.0.2.RData'</li>
	<li>'BT73_L_res.0.1.RData'</li>
	<li>'BT84_L_res.0.1.RData'</li>
	<li>'BT89_L_res.0.3.RData'</li>
	<li>'BT94_L_res.0.2.RData'</li>
</ol>




<ol class=list-inline>
	<li>'BT127_L'</li>
	<li>'BT147_L'</li>
	<li>'BT48_L'</li>
	<li>'BT67_L'</li>
	<li>'BT73_L'</li>
	<li>'BT84_L'</li>
	<li>'BT89_L'</li>
	<li>'BT94_L'</li>
</ol>




```R
dat.dir <- "~/pughlab/projects/Weiss_Astro_scRNAseq/data/"
files <- list.files(dat.dir, pattern = "_L_res.")
files <- files[-c(1,5)]

samples <- gsub('.{14}$', '', files)
samples

i <- 1
load.file <- paste0(dat.dir, files[i])
print(load.file)
load(load.file)
colnames(eb1S@meta.data)[grep("res", colnames(eb1S@meta.data))] <- "Cluster"

meta <- eb1S@meta.data
rownames(meta) <- paste(meta$orig.ident,rownames(meta), sep = "_")
head(meta)

raw.data <- eb1S@raw.data
colnames(raw.data) <- paste(meta$orig.ident, colnames(eb1S@raw.data), sep = "_")
head(raw.data)

merged.meta <- meta

#make a seurat object

merged <- CreateSeuratObject(raw.data, 
                             min.cells = 0,
                             min.genes = 0, 
                             normalization.method = NULL,
                             do.scale = FALSE, 
                             do.center = FALSE,
                             display.progress = TRUE
                            )


#now append all the other samples to this object created above

for (i in 2:length(files)){
    
load.file <- paste0(dat.dir, files[i])
print(load.file)
load(load.file)
colnames(eb1S@meta.data)[grep("res", colnames(eb1S@meta.data))] <- "Cluster"

meta <- eb1S@meta.data
rownames(meta) <- paste(meta$orig.ident,rownames(meta), sep = "_")
merged.meta <- rbind(merged.meta, meta)
dim(merged.meta)

raw.data <- eb1S@raw.data


merged <- AddSamples(object = merged, 
                     new.data = raw.data, 
                     do.normalize = FALSE,
                     add.cell.id = as.character(unique(eb1S@meta.data$orig.ident)),
                    )
    
}

table(merged@meta.data$orig.ident)

#add the meta.data 
merged <- AddMetaData(merged, metadata = merged.meta)

#normalize and scale the data
merged <-   NormalizeData(merged, 
                          assay.type = "RNA",
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000,
                          display.progress = TRUE
                         )
merged <- ScaleData(object = merged, vars.to.regress = c('percent.mito', "nUMI"))


```


```R
save(merged, file = "BTSC_Weiss_merged.RData")
```


```R
#now score all the data with the gene signatures

gene.sig.file <- "~/pughlab/projects/Weiss_Astro_scRNAseq/data/Weiss_GeneSignatures.RData"
load(gene.sig.file)

        for (j in 1:length(gene.list)){
    
            print("")
            print(names(gene.list)[j])
            aa <- gene.list[[j]][ ! gene.list[[j]] %in% rownames(merged@data)]
            print(aa)
            print("")  
            gene.list[[j]] <- gene.list[[j]][gene.list[[j]] %in% rownames(merged@data)]
    
        }

save(gene.list, file = "Weiss_Filtered_GeneList.RData")
```


```R
    load("Weiss_Filtered_GeneList.RData")

    print("")
    print("")
    print("Scoring Cell Cycle Phase")
    print("")
    print("")

    merged <- CellCycleScoring(object = merged, g2m.genes = gene.list$g2m.genes, s.genes = gene.list$s.genes)
    
    print("Scoring Gene Signatures")
    print("")
    print("")

    merged <- AddModuleScore(merged,
                           genes.list = gene.list[1:9],
                           n.bin = 25,
                           seed.use = 1,
                           ctrl.size = 100,
                           random.seed = 123,
                           enrich.name = names(gene.list)[1:9]
                      )
```


```R
#define random cutoff for scoring my simulating 100 gene lists

for (i in 1:length(gene.list)){
    
    sig.name <- names(gene.list)[i]
    print(sig.name)
    
    sig.length <- length(gene.list[[i]])
    print(sig.length)
    
    print("Generating random lists")
    random.sigs <- lapply(1:100, function(x) sample(rownames(merged@data),size = sig.length))
    
    names(random.sigs) <- paste0("Random", 1:100, "_", sig.name)
        
    file1 <- paste0(sig.name, "_random_lists.Rdata")
    print(file1)
    save(random.sigs, file = file1)
    
    print("Scoring random lists with AddModuleScore")
    
    merged <- AddModuleScore(merged,
                       genes.list = random.sigs,
                       n.bin = 25,
                       seed.use = 1,
                       ctrl.size = 100,
                       random.seed = 123,
                       enrich.name = names(random.sigs)
                      )
    print("")
    print("")
    print("")
    print("")
    
}
```


```R
#subet out random gene sigs
random.scores <- merged@meta.data[grep("Random", colnames(merged@meta.data))]
save(random.scores, file = "ADDMODULESCORE_GENESIG_RANDOM_SCORES.RData")

#remove all the random from the merged meta.data
merged@meta.data <- merged@meta.data[,-grep("Random", colnames(merged@meta.data))]
```


```R
#define 5% cutoff from random sigs

sigs <- c("Suva_Astro",
          "Suva_Astro_noEAATs",
          "Suva_Stemness",
          "Suva_G1S",
          "Suva_G2M"  ,
          "Suva_Oligo" ,
          "GO_Astro_Diff",
          "Housekeeping",
          "Genes_of_Interest"        
          )       


for (sig.name in sigs){
    
    print(sig.name)
    range <- grep(sig.name, colnames(random.scores))
    print(range)
    
    bb <- as.vector(as.matrix(random.scores[,range]))
    print(length(bb))
    
    d <- density(bb)
    cutoff <- qnorm(0.95, mean = mean(bb), sd = sd(bb))
    print(cutoff)
    
    filename <- paste0(sig.name, "AddModuleScore_RandomListCutoff.pdf")
    pdf(filename)

    plot_raw <- plot(d,
                 main = paste0(sig.name," - Random Signatures"),
                 xlab = "AddModuleScore()"
                )
    legend("topright",
       legend = c(paste("Mean:", round(mean(bb),3)),
                          paste("Sd:", round(sd(bb), 3)),
                          paste("Cutoff:", round(cutoff,3))
                     ),
               bty = 'n',
               cex = 1
          )
        abline(v=cutoff, lty = 2, col = "red")

    dev.off()
}

```


```R
save(merged, file = "BTSC_Weiss_merged_scores.RData")

meta <- merged@meta.data
save(meta, file = "BTSC_Weiss_merged_scores_meta.RData")
```

---
# 6.0 Proportion of cells
---

Visualize the proportion of cells classified as Astro, G2M and double pos per sample.  
Use stacked bar charts for this.


```R
#### write out classification columns for Ian 
meta$Suva_Astro_Classification <- meta$Suva_Astro1 >= 0.071 
meta$Suva_Astro_noEAATs_Classification <- meta$Suva_Astro_noEAATs2 >=0.068
meta$Suva_Stemness_Classification <- meta$Suva_Stemness3 >= 0.072
meta$Suva_G1S_Classification <- meta$Suva_G1S4 >= 0.087
meta$Suva_G2M_Classification <- meta$Suva_G2M5 >= 0.075
meta$GO_Astro_Diff_Classification <- meta$GO_Astro_Diff7 >= 0.091
meta$GenesOfInterest_Classification <- meta$Genes_of_Interest9 >= 0.34

#### merge MKI67 column to meta data

MKI67 <- data.frame(merged@scale.data["MKI67", ])
colnames(MKI67) <- "MKI67"
head(MKI67)
meta$MKI67 <- MKI67$MKI67

###remove Suva oligo, housekeeping, s.genes, g2m.genes

meta <- meta[ ,-c(1,2,4,5,14)]
head(meta)


#### write out table for Ian

write.table(meta, 
            file = "Weiss_GeneSig_Scores.txt", 
            quote = FALSE, 
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE
           )

```


```R
library(ggplot2)

setwd("~/Desktop/Samwise/projects/Weiss_Astro_scRNAseq/")
load("BTSC_Weiss_merged_scores_meta.RData")

meta$Astro_Classification <- meta$Suva_Astro1 >= 0.069 
#cutoff from analysisin 5 for suva astro signature
```


```R
meta$Astro_G2M <- meta$Suva_Astro1 >= 0.069 & meta$Phase == "G2M"
```


```R
#plot the proportion of astrocytic cells for each sample

pdf("~/Desktop/AstroG2M_correlations.pdf")

props <- prop.table(table(meta$orig.ident, meta$Astro_Classification),1)
props

barplot(t(props),
        las = 2,
        col = c("#a6611a", "#018571"),
        ylab = "Proportion Cells ",
        main = "Figure 1. \n Classification: Green = Astrocytic, Brown = Not Astrocytic"
       )


props.CC <- prop.table(table(meta$orig.ident, meta$Phase),1)
props.CC

barplot(t(props.CC),
        las = 2,
        col = c("#7fc97f", "#beaed4", "#fdc086"),
        ylab = "Proportion Cells ",
        main = "Figure 2. \n Classification: Green = G1, Purple = G2M , Orange = S"
       )


props.dp <- prop.table(table(meta$orig.ident, meta$Astro_G2M),1)
props.dp

barplot(t(props.dp),
        las = 2,
        col = c("black", "grey"),
        ylab = "Proportion Cells ",
        main = "Figure 3. \n Classification: Black = Double Negative, Grey = Astro+G2M"
 )


dev.off()

```


             
                  FALSE      TRUE
      BT127_L 0.3317175 0.6682825
      BT147_L 0.8990719 0.1009281
      BT48_L  0.8568507 0.1431493
      BT67_L  0.8585691 0.1414309
      BT73_L  0.9417829 0.0582171
      BT84_L  0.4961024 0.5038976
      BT89_L  1.0000000 0.0000000
      BT94_L  0.8528209 0.1471791



             
                      G1        G2M          S
      BT127_L 0.69598338 0.14265928 0.16135734
      BT147_L 0.38283063 0.27726218 0.33990719
      BT48_L  0.23449216 0.27880027 0.48670757
      BT67_L  0.35274542 0.26123128 0.38602329
      BT73_L  0.56519102 0.21043056 0.22437841
      BT84_L  0.81403118 0.05345212 0.13251670
      BT89_L  0.31019037 0.35946249 0.33034714
      BT94_L  0.40392478 0.24448078 0.35159444



             
                    FALSE        TRUE
      BT127_L 0.927977839 0.072022161
      BT147_L 0.988399072 0.011600928
      BT48_L  0.978186776 0.021813224
      BT67_L  0.992512479 0.007487521
      BT73_L  0.996967859 0.003032141
      BT84_L  0.986636971 0.013363029
      BT89_L  1.000000000 0.000000000
      BT94_L  0.986917416 0.013082584


    Warning message:
    “Removed 2 rows containing non-finite values (stat_smooth).”Warning message:
    “Removed 2 rows containing missing values (geom_point).”Warning message:
    “Removed 2 rows containing missing values (geom_text_repel).”



    Warning message:
    “Removed 2 rows containing non-finite values (stat_smooth).”Warning message:
    “Removed 2 rows containing missing values (geom_point).”Warning message:
    “Removed 2 rows containing missing values (geom_text_repel).”



    Warning message:
    “Removed 2 rows containing non-finite values (stat_smooth).”Warning message:
    “Removed 2 rows containing missing values (geom_point).”Warning message:
    “Removed 2 rows containing missing values (geom_text_repel).”




<strong>pdf:</strong> 2



```R
proportions <- data.frame(cbind(props, props.CC, props.dp))
colnames(proportions) <- c("Astro_FALSE", "Astro_TRUE", "G1", "G2M", "S", "DoubleNeg", "Astro_G2M")
proportions$PercentViability <- c(NA, 0.7, 0.82, 0.86, NA, 0.71, 0.03, 0.15)
proportions

```


<table>
<thead><tr><th></th><th scope=col>Astro_FALSE</th><th scope=col>Astro_TRUE</th><th scope=col>G1</th><th scope=col>G2M</th><th scope=col>S</th><th scope=col>DoubleNeg</th><th scope=col>Astro_G2M</th><th scope=col>PercentViability</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L</th><td>0.3317175  </td><td>0.6682825  </td><td>0.6959834  </td><td>0.14265928 </td><td>0.1613573  </td><td>0.9279778  </td><td>0.072022161</td><td>  NA       </td></tr>
	<tr><th scope=row>BT147_L</th><td>0.8990719  </td><td>0.1009281  </td><td>0.3828306  </td><td>0.27726218 </td><td>0.3399072  </td><td>0.9883991  </td><td>0.011600928</td><td>0.70       </td></tr>
	<tr><th scope=row>BT48_L</th><td>0.8568507  </td><td>0.1431493  </td><td>0.2344922  </td><td>0.27880027 </td><td>0.4867076  </td><td>0.9781868  </td><td>0.021813224</td><td>0.82       </td></tr>
	<tr><th scope=row>BT67_L</th><td>0.8585691  </td><td>0.1414309  </td><td>0.3527454  </td><td>0.26123128 </td><td>0.3860233  </td><td>0.9925125  </td><td>0.007487521</td><td>0.86       </td></tr>
	<tr><th scope=row>BT73_L</th><td>0.9417829  </td><td>0.0582171  </td><td>0.5651910  </td><td>0.21043056 </td><td>0.2243784  </td><td>0.9969679  </td><td>0.003032141</td><td>  NA       </td></tr>
	<tr><th scope=row>BT84_L</th><td>0.4961024  </td><td>0.5038976  </td><td>0.8140312  </td><td>0.05345212 </td><td>0.1325167  </td><td>0.9866370  </td><td>0.013363029</td><td>0.71       </td></tr>
	<tr><th scope=row>BT89_L</th><td>1.0000000  </td><td>0.0000000  </td><td>0.3101904  </td><td>0.35946249 </td><td>0.3303471  </td><td>1.0000000  </td><td>0.000000000</td><td>0.03       </td></tr>
	<tr><th scope=row>BT94_L</th><td>0.8528209  </td><td>0.1471791  </td><td>0.4039248  </td><td>0.24448078 </td><td>0.3515944  </td><td>0.9869174  </td><td>0.013082584</td><td>0.15       </td></tr>
</tbody>
</table>



---
# 7.0 Correlate cell types with drug response
---

Use absolute proportion values from the previous section and correlate to drug response values from Ian.  
Drug response metrics: https://docs.google.com/spreadsheets/d/1GIrYQaDiJKW5Bb1-J5ZiPuc-xreA6lU2WIr3sdU0ogM/edit?usp=sharing


```R
library(ggpubr)
```

    Warning message:
    “package ‘ggpubr’ was built under R version 3.4.4”Loading required package: magrittr



```R
a <- ggscatter(proportions, 
          x = "Astro_TRUE", 
          y = "PercentViability", 
          add = "reg.line", 
          add.params = list(color = "red", fill = "darkgrey", cex = 0.1),
          conf.int = TRUE, 
          #cor.coef = TRUE, 
          cor.method = "pearson",
          xlab = "Proportion Astrocytic Cells", 
          ylab = "Percent Viability",
          size = 1.5,
          title = "pearson r=0.4, p=0.43",
          color = "black",
          label = rownames(proportions),
          repel = TRUE,
          font.label = c(8)
         )

b <- ggscatter(proportions, 
          x = "G2M", 
          y = "PercentViability", 
          add = "reg.line", 
          add.params = list(color = "red", fill = "darkgrey", cex = 0.1),
          conf.int = TRUE, 
          #cor.coef = TRUE, 
          cor.method = "pearson",
          xlab = "Proportion G2M Cells", 
          ylab = "Percent Viability",
          size = 1.5,
          title = "pearson r=-0.39, p=0.45",
          color = "black",
          label = rownames(proportions),
          repel = TRUE,
          font.label = c(8)
         )

c <- ggscatter(proportions, 
          x = "Astro_G2M", 
          y = "PercentViability", 
          add = "reg.line", 
          add.params = list(color = "red", fill = "darkgrey", cex = 0.1),
          conf.int = TRUE, 
          #cor.coef = TRUE, 
          cor.method = "pearson",
          xlab = "Proportion G2M Astrocytic Cells", 
          ylab = "Percent Viability",
          size = 1.5,
          title = "pearson r=0.55, p=0.25",
          color = "black",
          label = rownames(proportions),
          repel = TRUE,
          font.label = c(8),
          xlim = c(0, 0.025)
         )





a
b
c
```

    Warning message:
    “Removed 2 rows containing non-finite values (stat_smooth).”Warning message:
    “Removed 2 rows containing missing values (geom_point).”Warning message:
    “Removed 2 rows containing missing values (geom_text_repel).”



    Warning message:
    “Removed 2 rows containing non-finite values (stat_smooth).”Warning message:
    “Removed 2 rows containing missing values (geom_point).”Warning message:
    “Removed 2 rows containing missing values (geom_text_repel).”




![png](output_29_4.png)


    Warning message:
    “Removed 2 rows containing non-finite values (stat_smooth).”Warning message:
    “Removed 2 rows containing missing values (geom_point).”Warning message:
    “Removed 2 rows containing missing values (geom_text_repel).”




![png](output_29_7.png)



![png](output_29_8.png)


---
**End of Analysis**  
March 2019  
L. Richards  

----
