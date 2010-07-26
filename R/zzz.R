# first and last lib functions

.onAttach = function(libname, pkgname)
{

  x = as.data.frame(installed.packages());
  ourRow = match("WGCNA", x$Package);
  ourVer = as.character(x$Version[ourRow]);
  printFlush("\nPackage WGCNA version", if (is.finite(ourRow)) ourVer else "unknown", "loaded.\n")
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

