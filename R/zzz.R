# first and last lib functions

.onAttach = function(libname, pkgname)
{

  x = as.data.frame(installed.packages());
  ourRow = match("WGCNA", x$Package);
  ourVer = as.character(x$Version[ourRow]);

  printFlush("==========================================================================\n*");
  printFlush("*  Package WGCNA version", if (is.finite(ourRow)) ourVer else "unknown", "loaded.\n*")

  if (.useNThreads()==1 && .nProcessorsOnline() > 1)
  {
   printFlush(spaste(
         "*    Important note: It appears that your system supports multi-threading,\n",
         "*    but it is not enabled within WGCNA in R. \n",
         "*    To allow multi-threading within WGCNA with all available cores, use \n",
         "*\n",
         "*          allowWGCNAThreads()\n",
         "*\n",
         "*    within R. Use disableWGCNAThreads() to disable threading if necessary.\n",
         "*    Alternatively, set the following environment variable on your system:\n", 
         "*\n",
         "*          ", .threadAllowVar, "=<number_of_processors>\n",
         "*\n",
         "*    for example \n",
         "*\n",
         "*          ", .threadAllowVar, "=", .nProcessorsOnline(), "\n",
         "*\n",
         "*    To set the environment variable in linux bash shell, type \n",
         "*\n",
         "*           export ", .threadAllowVar, "=", .nProcessorsOnline(), 
         "\n*",
         "\n*     before running R. Other operating systems or shells will", 
         "\n*     have a similar command to achieve the same aim.\n*"));
  }
  printFlush("==========================================================================\n\n");

  impRow = match("impute", x$Package);
  if (is.finite(impRow))
  {
    version = as.character(x$Version[impRow]);
    if (compareVersion(version, "1.12")< 0)
    {
      printFlush(paste("Caution: installed package 'impute' is not the newest version.\n",
            "Older versions can occasionally crash the code or the entire R session.\n",
            "If you already have the newest version available from CRAN, \n",
            "and you still see this warning, please download the impute package \n",
            "from Bioconductor at \n",
            "http://www.bioconductor.org/packages/release/bioc/html/impute.html . \n",
            "If the above link is dead, search for package 'impute' \n",
            "in the Downloads -> Software section of http://www.bioconductor.org .\n"));
    }
  }
}

