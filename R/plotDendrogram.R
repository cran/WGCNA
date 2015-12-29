#=======================================================================================================
#
# Plot dendrogram
#
#=======================================================================================================

.plotDendrogram = function(tree, labels = NULL, 
                          horiz = FALSE, reverseDirection = FALSE,
                          hang = 0.1, xlab = "", ylab = "", axes = TRUE,
                          cex.labels = 1, ..., adjustRange = FALSE)
{

  hang.gr = hang;
  if (hang < 0) hang.gr = 0.1;
  n = length(tree$order);
  heights = tree$height;
  range = range(heights);
  hang.scaled = hang.gr * (max(heights) - min(heights));
  range[1] = range[1] - hang.scaled;

  indexLim = c(0.5, n+0.5);
  if (adjustRange)
  {
    ctr = mean(indexLim);
    indexLim = ctr + (indexLim - ctr)/1.08;
  }

  nMerge = n-1;
  if (is.null(labels)) labels = tree$labels;
  if (is.null(labels)) labels = rep("", n);
  if (is.na(labels[1])) labels = rep("", n);
  if (is.logical(labels) && labels[1]=="FALSE") labels = rep("", n);

  if (horiz) 
  {
    plot(NA, NA, xlim = if (reverseDirection) range else rev(range), ylim = indexLim,
         axes = axes, yaxt = "none", frame = FALSE, type ="n", xlab = xlab, ylab = ylab, ...);
  } else {
    plot(NA, NA, ylim = if (reverseDirection) rev(range) else range, xlim = indexLim,
         axes = axes, xaxt = "none", frame = FALSE, type ="n", xlab = xlab, ylab = ylab, ...);
  }

  singleton.x = rep(NA, n);
  singleton.x[tree$order] = c(1:n);
  cluster.x = rep(NA, n);

  for (m in 1:nMerge)
  {
     o1 = tree$merge[m, 1]
     o2 = tree$merge[m, 2]
     h = heights[m];
     hh = if (hang>0) h-hang.scaled else range[1];
     h1 = if (o1 < 0) hh else heights[o1];
     h2 = if (o2 < 0) hh else heights[o2];

     x1 = if (o1 < 0) singleton.x[-o1] else cluster.x[o1]
     x2 = if (o2 < 0) singleton.x[-o2] else cluster.x[o2]

     cluster.x[m] = mean(c(x1, x2));

     if (!is.null(labels))
     {
       if (horiz)
       {
          if (o1 < 0) text(h1, x1, spaste(labels[-o1], " "), adj = c(0, 0.5), srt = 0,
                        cex = cex.labels, xpd = TRUE)
          if (o2 < 0) text(h2, x2, spaste(labels[-o2], " "), adj = c(0, 0.5), srt = 0,
                        cex = cex.labels, xpd = TRUE)
       } else {
          if (o1 < 0) text(x1, h1, spaste(labels[-o1], " "), adj = c(1, 0.5), srt = 90,
                        cex = cex.labels, xpd = TRUE)
          if (o2 < 0) text(x2, h2, spaste(labels[-o2], " "), adj = c(1, 0.5), srt = 90,
                        cex = cex.labels, xpd = TRUE)
       }
     }
  
     if (horiz)
     {
       lines(c(h1, h, h, h2), c(x1, x1, x2, x2));
     } else {
       lines(c(x1, x1, x2, x2), c(h1, h, h, h2));
     }
  }
}

