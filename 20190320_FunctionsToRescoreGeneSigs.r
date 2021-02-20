
## write a function to convert gene signature csv to list, check which genes are in the dataset
## genes sigs as columns in csv file

##Input:
##csv = list of gene signatures, columns with headers
##dat.mat = gene expression matrix, gene names must be rownames, form of matrix

##Output:
##cleaned gene lists to only include genes in your dataset
##table of genes in signature, but missing in data
##cleaned gene signatures in list format, saved as RData and csv
## file.prefix = what to name beginning of files

csvConvertToList <- function(csv,
                             dat.mat,
                             file.prefix
                            ){
    
    #read in csv file and data
    #extract rownames
    options(stringsAsFactors = FALSE)
    sigs <- read.csv(csv)
    sigs <- as.list(sigs)
    genes <- rownames(dat.mat)
    
    print("Genes in data matrix:")
    print(length(genes))
    print("")
    print("Number of gene signautes:")
    print(length(sigs))
    
    #clean up gene  signatures
    #remove empty spaces in lists
    #l is list

    removeEmpty <- function(l){
        l[l!= ""]
    }
    sigs <- lapply(sigs, removeEmpty)
   
    ##check how many genes are in data
    ##remove genes not in data frame
    sigs.cleaned <- lapply(sigs, function(x) x[x %in% genes])
       
    genes.removed <- setdiff(unlist(sigs), unlist(sigs.cleaned))
    print("")
    print("The following genes are not in data:")
    print(genes.removed)
           
    ##save files
    save(sigs.cleaned, genes.removed, file = paste0(file.prefix, ".RData") )
        
}

##write a function to score cells and define a cutoff

## INPUT: 
#gene.sigs = object in list format of gene signatures you want to score
#data.matrix = data matrix for classification, must be a seurat object
#num.random = number of random gene lists to evaluate
#results.name = string to append to saved results file

## OUTPUT:
# seruat object with scores
# classification cutoffs for each gene signature
# single classification of each cell
# 


scoreAndClassify_AddModuleScore <- function(gene.sigs, 
                             data,
                             num.random,
                              num.bin,
                             results.name
                                            
                            ) {
  
    start <- Sys.time()
    print("")
    print("Gene signature details:")
    print(Sys.time())
    print(str(gene.sigs))
    

    if(class(data) == "seurat"){
        
        print("")
        print("Expression data is a Seurat Object")
        print("")
        print("Scoring gene signatures with AddModuleScore()")
        print(Sys.time())
        
        scored <- AddModuleScore(data,
                       genes.list = gene.sigs,
                       n.bin = num.bin,
                       seed.use = 1,
                       ctrl.size = num.random,
                       random.seed = 123,
                       enrich.name = names(gene.sigs)     
                      )
        print("")
        print("Scoring gene signatures complete")
        print(Sys.time())
        
        print("")
        print("Generate 100 random gene lists of the same size for each signature")
        print(Sys.time())
        
        size <- as.numeric(lapply(gene.sigs, length)) #extract length of gene signatures
        random.sigs <- replicate(100, lapply(size, function(x) sample(rownames(scored@data), size = x)))
         
        #name the random sigs to match parent signature
        names(random.sigs) <- seq(1:length(random.sigs))
     
        
        for (i in 1:length(gene.sigs)){ 
            n <- seq(from = i, to = length(random.sigs), by = length(gene.sigs))
            names(random.sigs)[n] <- paste0("Random", 1:100, "_", names(gene.sigs)[i])        
        }
            
        print("")
        print("Score random gene signatures")
        print("...this might take a while, go get a coffee....")
        print(length(random.sigs))
        print(Sys.time())
        
        
        scored <- AddModuleScore(scored,
                       genes.list = random.sigs,
                       n.bin = num.bin,
                       seed.use = 1,
                       ctrl.size = num.random,
                       random.seed = 123,
                       enrich.name = names(random.sigs)     
                      )
            
        random.scores <- scored@meta.data[grep("Random", colnames(scored@meta.data))]
        scored@meta.data <- scored@meta.data[,-grep("Random", colnames(scored@meta.data))]
        
        
        print("")
        print("Scoring random signatures complete")
        print(Sys.time())
            
        print("")
        print("Define top 5% cutoff for classification and plot")
        print(Sys.time())
        
        cutoffs <- c()
            
        for (sig.name in names(gene.sigs)){
         
            print(sig.name)
            range <- grep(sig.name, colnames(random.scores))
            #print(range)
            
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
        cutoffs <- c(cutoffs, cutoff)
        
        }
        
        names(cutoffs) <- names(gene.sigs)
            
        print("")
        print("Classify cells")
        print(Sys.time())
        
        meta <- scored@meta.data
            
        for (i in 1:length(gene.sigs)){
        
        subset <- meta[ ,grep(names(gene.sigs[i]), colnames(meta))]
        
        
        col.name <- paste0("Classify_", names(gene.sigs[i]), "_", 
                           round(as.numeric(cutoffs[grep(names(gene.sigs[i]), names(cutoffs))]), 4)
                           )
            
        meta[col.name] <- subset > as.numeric(cutoffs[grep(names(gene.sigs[i]), names(cutoffs))])
        
        }
            
        print("")
        print("Save results and write out csv")
        print(Sys.time())
            
        write.csv(meta, 
                  file = paste0(results.name, ".csv"), 
                  row.names = T
                 )
            
        save(random.sigs,
             gene.sigs,
             scored,
             random.scores,
             meta,
             file = paste0(results.name, ".RData")
            )
            
        
        print("")   
        print("---------------")
        end <- Sys.time()
        print("END OF FUNCTION")
        print(Sys.time())
            
        print("Number of cells:")
        print(ncol(scored@data))
        print("Number of gene signatures:")
        print(length(gene.sigs))
        print("Number of random signatures:")
        print(length(random.sigs))
        print("Total run time:")
        print(end - start)
        print("---------------")

    
      } else { 
        print("Expression matrix is not a Seurat Object")
        }

}




### test run

library(Seurat)
options(stringsAsFactors = FALSE)
load("../BTSC_Weiss_merged.RData")

csvConvertToList("WeissLab_GeneSignatures_Dec2018 .GeneSignatures.csv",
                 merged@data,
                 "WeissSignatures_April2019_TEST"
)

scoreAndClassify_AddModuleScore(sigs.cleaned, 
                                merged, 
                                100, 
                                "Weiss_testRun"
                               )
                                     
