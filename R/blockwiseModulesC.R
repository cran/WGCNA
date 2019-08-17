# Hello,

# I am doing co-expression network with WGCNA on RNA-seq data (70-200 samples). 
# While using blockwiseModules() function, I obtain modules smaller than the module size cut-off. 
# Namely, minModuleSize <= 20, I got modules with just 1, 2, 4, genes. 
# It seems a kind of exception in the function.

# How can this be possible?

# Further, I can say that this happens more frequently when the soft-thresholding power < 4 (unsigned, signed hybrid networks) or power < 9 (signed networks) even if the soft-threshold criterion is satisfied. The occurrence for larger powers is negligible, still some exceptions remain.

# Here the code, the most is default:

net = blockwiseModules(input_train, 
                       power = beta_value, 
                       networkType = NT, 
                       randomSeed = seed, 
                       minModuleSize = 20, 
                       corType = "bicor", 
                       mergeCutHeight = 0.15, 
                       pamStage = TRUE, 
                       pamRespectsDendro = FALSE, 
                       TOMType = "signed", 
                       saveTOMs = FALSE, 
                       maxBlockSize = 5000, 
                       numericLabels = TRUE, 
                       nThreads = 4, 
                       verbose = 3)
