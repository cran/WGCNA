# Function that returns GO terms with best enrichment and highest number of genes.

# So that I don't forget: information about GO categories is contained in the GO.db package

GOenrichmentAnalysis = function(labels, entrezCodes, 
               yeastORFs = NULL,
               organism = "human",
               ontologies = c("BP", "CC", "MF"),
               evidence = "all",
               includeOffspring = TRUE,
               backgroundType = "givenInGO",
               removeDuplicates = TRUE,
               leaveOutLabel = NULL,
               nBestP = 10, pCut = NULL, nBiggest = 0,
               getTermDetails = TRUE,
               verbose = 2, indent = 0 )
{

   sAF = options("stringsAsFactors")
   options(stringsAsFactors = FALSE);
   on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)

   organisms = c("human", "mouse", "rat", "malaria", "yeast", "fly", "bovine", "worm", "canine",
                 "zebrafish", "chicken");
   allEvidence =  c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "ISS", "ISO", "ISA",
                    "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA", "TAS", "NAS", "IC", "ND", "IEA",
                    "NR");
   allOntologies = c("BP", "CC", "MF");

   backgroundTypes = c("allGiven", "allInGO", "givenInGO");

   spaces = indentSpaces(indent);
   orgInd = pmatch(organism, organisms);
   if (is.na(orgInd))
     stop(paste("Unrecognized 'organism' given. Recognized values are ", 
                paste(organisms, collapse = ", ")));

   if (length(evidence)==0)
     stop("At least one valid evidence code must be given in 'evidence'.");
   if (length(ontologies)==0)
     stop("At least one valid ontology code must be given in 'ontology'.");

   if (evidence=="all")
      evidence = allEvidence;

   evidInd = pmatch(evidence, allEvidence);
   if (sum(is.na(evidInd))!=0)
     stop(paste("Unrecognized 'evidence' given. Recognized values are ", 
                paste(allEvidence, collapse = ", ")));

   ontoInd = pmatch(ontologies, allOntologies);
   if (sum(is.na(ontoInd))!=0)
     stop(paste("Unrecognized 'ontologies' given. Recognized values are ", 
                paste(allEvidence, collapse = ", ")));

   backT = pmatch(backgroundType, backgroundTypes);
   if (is.na(backT))
     stop(paste("Unrecognized 'backgroundType' given. Recognized values are ", 
           paste(backgroundTypes,  collapse = ", ")));

   orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
   orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
   reverseMap = c(rep(".egGO2EG", 4), ".sgdGO2ORF", rep(".egGO2EG", 6))

   missingPacks = NULL;
   packageName = paste("org.", orgCodes[orgInd], orgExtensions[orgInd], ".db", sep="");
   if (!eval(parse(text = "require(packageName, character.only = TRUE)")))
     missingPacks = c(missingPacks, packageName);

   if (!eval(parse(text="require(GO.db)")))
     missingPacks = c(missingPacks, "GO.db");

   if (!is.null(missingPacks)) 
     stop(paste("Could not load the requisite package(s)", 
           paste(missingPacks, collapse = ", "), ". Please install the package(s)."))

   if (verbose > 0)
   {
     printFlush(paste(spaces, "GOenrichmentAnalysis: loading annotation data..."));
   }

   if (orgInd==5)
   {
      # Yeast needs special care.
      if (!is.null(yeastORFs))
      {
        entrezCodes = yeastORFs
      } else {
        # Map the entrez IDs to yeast ORFs
        x = eval(parse(text = "org.Sc.sgd:::org.Sc.sgdENTREZID"))
        # x = org.Sc.sgd:::org.Sc.sgdENTREZID
        xx = as.list(x[mapped_genes])
        allORFs = names(xx);
        mappedECs = as.character(sapply(xx, as.character))
        entrez2orf = match(entrezCodes, mappedECs);
        fin = is.finite(entrez2orf);
        newCodes = paste("InvalidCode", c(1:length(entrezCodes)), sep = ".");
        newCodes[fin] = allORFs[entrez2orf[fin]];
        entrezCodes = newCodes;
     }
   } 

   labels = as.matrix(labels);

   nSets = ncol(labels);
   nGivenRaw = nrow(labels);

   if (removeDuplicates)
   {
     # Restrict given entrezCodes such that each code is unique
     ECtab = table(entrezCodes);
     uniqueEC = names(ECtab);
     keepEC = match(uniqueEC, entrezCodes);
     entrezCodes = entrezCodes[keepEC]
     labels = labels[keepEC, , drop = FALSE];
   } else
     keepEC = c(1:nGivenRaw);

   egGO = eval(parse(text = paste(packageName, ":::org.", orgCodes[orgInd], orgExtensions[orgInd], 
                                  "GO", sep = "")));

   if (orgInd==5)
   {
      mapped_genes = as.character(do.call(match.fun("mappedkeys"), list(egGO)));
      encodes2mapped = match(as.character(entrezCodes), mapped_genes);
   } else {
      mapped_genes = as.numeric(as.character(do.call(match.fun("mappedkeys"), list(egGO))));
      encodes2mapped = match(as.numeric(entrezCodes), mapped_genes);
   }
   encMapped = is.finite(encodes2mapped);
   nAllIDsInGO = sum(encMapped);

   mapECodes = entrezCodes[encMapped];
   mapLabels = labels[encMapped, , drop = FALSE];
   nMappedGenes = nrow(mapLabels);

   if (nMappedGenes==0)
     stop(paste("None of the supplied gene identifiers map to the GO database.\n",
                "Please make sure you have specified the correct organism (default is human)."))

   Go2eg = eval(parse(text = paste("AnnotationDbi::as.list(", packageName, ":::org.", orgCodes[orgInd], 
                                           reverseMap[orgInd],")", sep = "")));
   nTerms = length(Go2eg);

   goInfo = as.list(GO.db::GOTERM);
   if (length(goInfo) > 0)
   {
      orgGoNames = names(Go2eg);
      dbGoNames = as.character(sapply(goInfo, GOID));
      dbGoOntologies = as.character(sapply(goInfo, Ontology));
   } else {
      dbGoNames = "";
   }
   goOffSpr = list();
   if (includeOffspring)
   {
     goOffSpr[[1]] = as.list(GOBPOFFSPRING);
     goOffSpr[[2]] = as.list(GOCCOFFSPRING);
     goOffSpr[[3]] = as.list(GOMFOFFSPRING);
   }
   term2info = match(names(Go2eg), names(goInfo));
   termOntologies = dbGoOntologies[term2info];

   if (backT==1)
   {
     nBackgroundGenes = sum(!is.na(entrezCodes));
   } else
   {
     if (backT==2)
     {
       nBackgroundGenes = length(mapped_genes)
     } else
       nBackgroundGenes = nMappedGenes;
   }

   termCodes = vector(mode="list", length = nTerms);
   collectGarbage();
   nExpandLength = 0;
   blockSize = 3000; # For a more efficient concatenating of offspring genes
   nAllInTerm = rep(0, nTerms);
   if (verbose > 0)
   {
      printFlush(paste(spaces, " ..of the", length(entrezCodes), 
                      " Entrez identifiers submitted,", sum(encMapped), 
                      "are mapped in current GO categories."));
      printFlush(paste(spaces, " ..will use", nBackgroundGenes, 
                       "background genes for enrichment calculations."));
      cat(paste(spaces, " ..preparing term lists (this may take a while).."));
      pind = initProgInd();
   }

   for (c in 1:nTerms) if (!is.na(Go2eg[[c]][[1]]))
   {
      te = as.character(names(Go2eg[[c]])); # Term evidence codes
      tc = Go2eg[[c]];
      if (includeOffspring) 
      {
        termOffspring = NULL;
        for (ont in 1:length(goOffSpr))
        {
          term2off = match(names(Go2eg)[c], names(goOffSpr[[ont]]))
          if (!is.na(term2off))
            termOffspring = c(termOffspring, goOffSpr[[ont]][[term2off]]);
        }
        if (length(termOffspring)>0) 
        {
           maxLen = blockSize;
           tex = rep("", maxLen);
           tcx = rep("", maxLen);
           ind = 1;
           len = length(te);
           tex[ ind:len ] = te;
           tcx[ ind:len ] = tc;
           ind = len + 1;
           o2go = match(termOffspring, as.character(names(Go2eg)));
           o2go = o2go[is.finite(o2go)]
           if (length(o2go)>0) for (o in 1:length(o2go)) if (!is.na(Go2eg[[o2go[o]]][[1]]))
           {
             #printFlush(paste("Have offspring for term", c, ": ", names(Go2eg)[c], 
             #           Term(goInfo[[term2info[c]]])));
             newc = Go2eg[[o2go[o]]];
             newe = names(newc);
             newl = length(newe);
             if ((len + newl) > maxLen)
             {
               nExpand = ceiling( (len + newl - maxLen)/blockSize);
               maxLen = maxLen + blockSize * nExpand;
               tex = c(tex, rep("", maxLen - length(tex)));
               tcx = c(tcx, rep("", maxLen - length(tex)));
               nExpandLength = nExpandLength + 1;
             }
             tex[ind:(len + newl)] = newe;
             tcx[ind:(len + newl)] = newc;
             ind = ind + newl;
             len = len + newl;
           }
           te = tex[1:len];
           tc = tcx[1:len];
        }
      }
      use = is.finite(match(te, evidence));
      if (orgInd==5)
      {
         if (backT==2)
         {
            termCodes[[c]] = unique(as.character(tc[use]));
         } else
            termCodes[[c]] = as.character(intersect(tc[use], mapECodes));
      } else {
         if (backT==2)
         {
            termCodes[[c]] = unique(as.character(tc[use]));
         } else
            termCodes[[c]] = as.numeric(as.character(intersect(tc[use], mapECodes)));
      }
      nAllInTerm[c] = length(termCodes[[c]]); 
      if ( (c %%50 ==0) & (verbose > 0)) pind = updateProgInd(c/nTerms, pind);
   }
   if (verbose > 0) 
   {
      pind = updateProgInd(1, pind);
      printFlush("");
   }
   if ((verbose > 5) & (includeOffspring))
      printFlush(paste(spaces, " ..diagnostic for the developer: offspring buffer was expanded",
                       nExpandLength, "times."));

   ftp = function(...) { fisher.test(...)$p.value }

   setResults = list();

   for (set in 1:nSets)
   {
      if (verbose > 0)
        printFlush(paste(spaces, " ..working on label set", set, ".."));
      labelLevels = levels(factor(labels[, set]));
      if (!is.null(leaveOutLabel))
      {
        keep = !(labelLevels %in% as.character(leaveOutLabel));
        if (sum(keep)==0)
          stop("No labels were kept after removing labels that are supposed to be ignored.");
        labelLevels = labelLevels[keep]
      }
      nLabelLevels = length(labelLevels);

      modCodes = list();
      nModCodes = rep(0, nLabelLevels);
      if (backT==1)
      {
        for (ll in 1:nLabelLevels)
        {
           modCodes[[ll]] = entrezCodes[labels[, set]==labelLevels[ll]];
           nModCodes[ll] = length(modCodes[[ll]]);
        }
      } else 
      {
        for (ll in 1:nLabelLevels)
        {
           modCodes[[ll]] = mapECodes[mapLabels[, set]==labelLevels[ll]];
           nModCodes[ll] = length(modCodes[[ll]]);
        }
      }

      countsInTerm = matrix(0, nLabelLevels, nTerms);
      enrichment = matrix(1, nLabelLevels, nTerms);

      for (ll in 1:nLabelLevels)
        countsInTerm[ll, ] = sapply(lapply(termCodes, intersect, modCodes[[ll]]), length)

      nAllInTermMat = matrix(nAllInTerm, nLabelLevels, nTerms, byrow = TRUE);
      nModCodesMat = matrix(nModCodes, nLabelLevels, nTerms);
      tabArr = array(c(countsInTerm, nAllInTermMat - countsInTerm, 
                       nModCodesMat - countsInTerm,
                       nBackgroundGenes - nModCodesMat - nAllInTermMat + countsInTerm),
                     dim = c(nLabelLevels * nTerms, 2, 2));

      if (verbose > 0)
         printFlush(paste(spaces, "   ..calculating enrichments (this may also take a while).."));
      calculate = c(countsInTerm) > 0;
      enrichment[calculate] = apply(tabArr[calculate, , ], 1, ftp, alternative = "g");
    
      dimnames(enrichment) = list (labelLevels, names(Go2eg));
      dimnames(countsInTerm) = list (labelLevels, names(Go2eg));

      bestPTerms = list();
      modSizes = table(labels[ !(labels[, set] %in% leaveOutLabel), set]);

      if (!is.null(pCut) || nBestP > 0)
      {
         printFlush(paste(spaces, "   ..putting together terms with highest enrichment significance.."));
         nConsideredOntologies = length(ontologies)+1;
         for (ont in 1:nConsideredOntologies)
         {
            if (ont==nConsideredOntologies)
            {
               ontTerms = is.finite(match(termOntologies, ontologies))
               bestPTerms[[ont]] = list(ontology = ontologies);
               names(bestPTerms)[ont] = paste(ontologies, collapse = ", ");
            } else {
               ontTerms = termOntologies==ontologies[ont];
               bestPTerms[[ont]] = list(ontology = ontologies[ont]);
               names(bestPTerms)[ont] = ontologies[ont];
            }
            bestPTerms[[ont]]$enrichment = NULL;
            bestPTerms[[ont]]$termInfo = list();
            nOntTerms = sum(ontTerms)
            ontEnr = enrichment[, ontTerms, drop = FALSE];
            order = apply(ontEnr, 1, order);
            for (ll in 1:nLabelLevels)
            {
              if (!is.null(pCut))
              {
                 reportTerms = c(1:nTerms)[ontTerms][ontEnr[ll, ] < pCut];
                 reportTerms = reportTerms[order(ontEnr[ll, ][reportTerms])];
              } else
                 reportTerms =  c(1:nTerms)[ontTerms][order[1:nBestP, ll]];
              nRepTerms = length(reportTerms);
              enrTab = data.frame(module = rep(labelLevels[ll], nRepTerms), 
                                  modSize = rep(modSizes[ll], nRepTerms),
                                  bkgrModSize = rep(nModCodes[ll], nRepTerms), 
                                  rank = c(1:nRepTerms),
                                  enrichmentP = enrichment[ll, reportTerms],
                                  BonferoniP = pmin(rep(1, nRepTerms), 
                                                    enrichment[ll, reportTerms] * nOntTerms),
                                  nModGenesInTerm = countsInTerm[ll, reportTerms],
                                  fracOfBkgrModSize = countsInTerm[ll, reportTerms]/nModCodes[ll],
                                  fracOfBkgrTermSize = 
                                       countsInTerm[ll, reportTerms]/nAllInTerm[reportTerms],
                                  bkgrTermSize = nAllInTerm[reportTerms],
                                  termID = names(Go2eg)[reportTerms],
                                  termOntology =  rep("NA", nRepTerms),
                                  termName = rep("NA", nRepTerms),
                                  termDefinition =  rep("NA", nRepTerms))
              bestPTerms[[ont]]$forModule[[ll]] = list();
              for (rci in 1:length(reportTerms))
              {
                 term = reportTerms[rci];
                 termID = names(Go2eg)[term];
                 dbind = match(termID, dbGoNames);
                 if (is.finite(dbind))
                 {
                   enrTab$termName[rci] = eval(parse(text = "AnnotationDbi:::Term(goInfo[[dbind]])"));
                   enrTab$termDefinition[rci] = 
                         eval(parse(text = "AnnotationDbi:::Definition(goInfo[[dbind]])"));
                   enrTab$termOntology[rci] = 
                         eval(parse(text = "AnnotationDbi:::Ontology(goInfo[[dbind]])"));
                 } 
                 if (getTermDetails)
                 {
                   geneCodes = intersect(modCodes[[ll]], termCodes[[term]])
                   bestPTerms[[ont]]$forModule[[ll]][[rci]] = list(termID = termID,
                           termName = enrTab$termName[rci],
                           enrichmentP = enrTab$enrichmentP[rci],
                           termDefinition = enrTab$termDefinition[rci],
                           termOntology = enrTab$termOntology[rci],
                           geneCodes = geneCodes,
                           genePositions = keepEC[match(geneCodes, entrezCodes)]);
                }
              }
              if (ll==1) 
              {
                 bestPTerms[[ont]]$enrichment = enrTab;
              } else 
                 bestPTerms[[ont]]$enrichment = rbind(bestPTerms[[ont]]$enrichment, enrTab);
           }
         }
      }
           
      biggestTerms = list();
   
      if (nBiggest > 0)
      {
         printFlush(paste(spaces, "   ..putting together terms with largest number of genes in modules.."));
         nConsideredOntologies = length(ontologies)+1;
         for (ont in 1:nConsideredOntologies)
         {
            if (ont==nConsideredOntologies)
            {
               ontTerms = is.finite(match(termOntologies, ontologies))
               biggestTerms[[ont]] = list(ontology = ontologies);
               names(biggestTerms)[ont] = paste(ontologies, collapse = ", ");
            } else {
               ontTerms = termOntologies==ontologies[ont];
               biggestTerms[[ont]] = list(ontology = ontologies[ont]);
               names(biggestTerms)[ont] = ontologies[ont];
            }
            biggestTerms[[ont]]$enrichment = NULL;
            biggestTerms[[ont]]$termInfo = list();
            nOntTerms = sum(ontTerms)
            ontNGenes = countsInTerm[, ontTerms, drop = FALSE];
            order = apply(-ontNGenes, 1, order);
            for (ll in 1:nLabelLevels)
            {
              reportTerms = c(1:nTerms)[ontTerms][order[1:nBiggest, ll]];
              nRepTerms = length(reportTerms);
              enrTab = data.frame(module = rep(labelLevels[ll], nRepTerms),
                                  modSize = rep(modSizes[ll], nRepTerms),
                                  bkgrModSize = rep(nModCodes[ll], nRepTerms), 
                                  rank = c(1:nRepTerms),
                                  enrichmentP = enrichment[ll, reportTerms],
                                  BonferoniP = pmin(rep(1, nRepTerms), 
                                                    enrichment[ll, reportTerms] * nOntTerms),
                                  nModGenesInTerm = countsInTerm[ll, reportTerms],
                                  fracOfModSize = countsInTerm[ll, reportTerms]/nModCodes[ll],
                                  fracOfBkgrTermSize = 
                                       countsInTerm[ll, reportTerms]/nAllInTerm[reportTerms],
                                  bkgrTermSize = nAllInTerm[reportTerms],
                                  termID = names(Go2eg)[reportTerms],
                                  termOntology =  rep("NA", nRepTerms),
                                  termName = rep("NA", nRepTerms),
                                  termDefinition =  rep("NA", nRepTerms))
              biggestTerms[[ont]]$forModule[[ll]] = list();
              for (rci in 1:length(reportTerms))
              {
                 term = reportTerms[rci];
                 termID = names(Go2eg)[term];
                 dbind = match(termID, dbGoNames);
                 if (is.finite(dbind))
                 {
                   enrTab$termName[rci] = eval(parse(text="AnnotationDbi:::Term(goInfo[[dbind]])"));
                   enrTab$termDefinition[rci] =
                       eval(parse(text="AnnotationDbi:::Definition(goInfo[[dbind]])"));
                   enrTab$termOntology[rci] = eval(parse(text="AnnotationDbi:::Ontology(goInfo[[dbind]])"));
                 }
                 if (getTermDetails)
                 {
                   geneCodes = intersect(modCodes[[ll]], termCodes[[term]])
                   biggestTerms[[ont]]$forModule[[ll]][[rci]] = list(termID = termID,
                           termName = enrTab$termName[rci],
                           enrichmentP = enrTab$enrichmentP[rci],
                           termDefinition = enrTab$termDefinition[rci],
                           termOntology = enrTab$termOntology[rci],
                           geneCodes = geneCodes,
                           genePositions = keepEC[match(geneCodes, entrezCodes)]);
                }
              }
              if (ll==1)
              {
                 biggestTerms[[ont]]$enrichment = enrTab;
              } else
                 biggestTerms[[ont]]$enrichment = rbind(biggestTerms[[ont]]$enrichment, enrTab);
           }
         }
      }
      setResults[[set]] = list(countsInTerm = countsInTerm, 
                               enrichmentP = enrichment,
                               bestPTerms = bestPTerms,
                               biggestTerms = biggestTerms);
   }
   
   inGO = rep(FALSE, nGivenRaw);
   inGO[keepEC] = encMapped;
   kept = rep(FALSE, nGivenRaw);
   kept[keepEC] = TRUE;
   if (nSets==1)
   {
      list(keptForAnalysis = kept,
           inGO = inGO,
           countsInTerms = setResults[[1]]$countsInTerm,
           enrichmentP = setResults[[1]]$enrichmentP,
           bestPTerms = setResults[[1]]$bestPTerms,
           biggestTerms = setResults[[1]]$biggestTerms);
   } else {
      list(keptForAnalysis = kept,
           inGO = inGO,
           setResults = setResults);
   }
}

