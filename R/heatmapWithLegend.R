# Replacement for the function image.plot

.autoTicks = function(min, max, maxTicks = 6 , tickPos = c(1,2,5))
{
  range = max - min;
  tick0 = range/maxTicks;
  maxTick = max(tickPos);
  # Ticks can only be multiples of tickPos
  mult = 1;
  if (tick0 < maxTick/10)
  {
     while (tick0 < maxTick/10) {tick0 = 10*tick0; mult = mult*10; }
  } else
     while (tick0 >=maxTick ) {tick0 = tick0/10; mult = mult/10; }

  ind = sum(tick0 > tickPos) + 1;
  tickStep = tickPos[ind] / mult;

  lowTick = min/tickStep;
  if (floor(lowTick)!=lowTick) lowTick = lowTick + 1;
  lowTick = floor(lowTick);

  ticks = tickStep * (lowTick:(lowTick + maxTicks+1));
  ticks = ticks[ticks <= max];
  ticks;
}

.plotStandaloneLegend = function(
                            colors,
                            lim,
                            tickLen = 0.7* strwidth("M"),
                            tickGap = 0.3 * strwidth("M"),
                            minBarWidth = strwidth("M"),
                            maxBarWidth = Inf,
                            mar = c(0.5, 0.2, 0.5, 0.1))
{
  par(mar = mar);
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "");
  box = par("usr");
  tickVal = .autoTicks(lim[1], lim[2]);
  maxTickWidth = max(strwidth(tickVal));
  if (maxTickWidth + tickLen + tickGap > box[2]-box[1]-minBarWidth) 
     warning("Some tick labels will be truncated.");
  xMax = max(box[2]-maxTickWidth - tickLen - tickGap, box[1] + minBarWidth);
  if (xMax - box[1] > maxBarWidth) xMax = box[1] + maxBarWidth;
  .plotColorLegend(box[1], xMax,
                   box[3], box[4], 
                   colors = colors,
                   lim = lim,
                   tickLen = tickLen,
                   tickGap = tickGap);
}

.plotColorLegend = function(xmin, xmax, ymin, ymax,
                            colors,
                            tickLen = 0.7* strwidth("M"),
                            tickGap = 0.3 * strwidth("M"),
                            lim, cex.legend = 1)
{
      tickVal = .autoTicks(lim[1], lim[2]);
      tickY = (tickVal - lim[1]) / (lim[2] - lim[1]) * (ymax - ymin) + ymin;
      nTicks = length(tickVal);

      # Ticks:
      for (t in 1:nTicks)
        lines(c(xmax, xmax + tickLen), c(tickY[t], tickY[t]), xpd = TRUE);
      text(rep(xmax + tickLen + tickGap), tickY, tickVal, adj = c(0, 0.5), cex = cex.legend,
           xpd = TRUE);

      # Fill with color:
      nColors = length(colors);
      ybl = (ymax-ymin)/nColors * (0:(nColors-1)) + ymin;
      ytl = (ymax-ymin)/nColors * (1:nColors) + ymin;
      rect(xleft = rep(xmin, nColors), xright = rep(xmax, nColors),
           ybottom = ybl, ytop = ytl, col = colors, border = colors, xpd = TRUE);

      lines(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin), xpd = TRUE );
}


.heatmapWithLegend = function(data, signed, colors, naColor = "grey", zlim = NULL, 
                     reverseRows = TRUE,
                     plotLegend = TRUE,
                     keepLegendSpace = plotLegend,
                     cex.legend = 1, 
                     legendShrink = 0.94,
                     legendSpace = 4 * strwidth("M"),
                     legendWidth = 1 * strwidth("M"),
                     legendGap = 0.7 * strwidth("M"),
                     frame = TRUE,
                     frameTicks = FALSE, tickLen = 0.7* strwidth("M"),
                     ...)
{
  data = as.matrix(data); nCols = ncol(data); nRows = nrow(data);
  if (is.null(zlim)) 
  {
    zlim = range(data, na.rm = TRUE);
    if (signed) zlim = c(-max(abs(zlim)), max(abs(zlim)));
  }

  barplot(1, col = "white", border = "white", axisnames = FALSE,
                  axes = FALSE, ...);

  box = par("usr");
  xminAll = box[1]; 
  xmaxAll = box[2]; 
  yminAll = box[3]; 
  ymaxAll = box[4]; 

  if (!keepLegendSpace && !plotLegend)
  {
     legendSpace = 0;
     legendWidth = 0;
     legendGap = 0;
  }

  ymin = yminAll; 
  ymax = ymaxAll; 
  xmin = xminAll; 
  xmax = xmaxAll - legendSpace;
  if (xmax < xmin) stop("'legendSpace is too large, not enough space for the heatmap."); 

  xStep = (xmax - xmin)/nCols; 
  xLeft = xmin + c(0:(nCols-1)) * xStep;
  xRight = xLeft + xStep; 
  xMid = (xLeft + xRight)/2;

  yStep = (ymax - ymin)/nRows; yBot  = ymin + c(0:(nRows-1)) * yStep;
  yTop  = yBot + yStep; yMid = c(yTop+ yBot)/2;

  if (reverseRows)
  {
    colorMat = numbers2colors(.reverseRows(data), signed, colors = colors, lim = zlim,
                              naColor = naColor)
  } else
    colorMat = numbers2colors(data, signed, colors = colors, lim = zlim, naColor = naColor)

  dim(colorMat) = dim(data);

  for (c in 1:nCols)
  {
    rect(xleft = rep(xLeft[c], nRows), xright = rep(xRight[c], nRows),
         ybottom = yBot, ytop = yTop, col = colorMat[, c], border = colorMat[, c]);
  }
  if (frame) lines( c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin) );

  if (plotLegend)
  {
      # Now plot the legend.
      .plotColorLegend(xmin = xmaxAll - (legendSpace - legendGap),
                       xmax = xmaxAll - (legendSpace - legendGap - legendWidth),
                       ymin = yminAll + (1-legendShrink) * (ymaxAll - yminAll),
                       ymax =  ymaxAll - (1-legendShrink) * (ymaxAll - yminAll),
                       lim = zlim,
                       colors = colors,
                       tickLen = tickLen,
                       cex.legend = cex.legend);
    
  }

  list(xMid = xMid, yMid = if (reverseRows) .reverseVector(yMid) else yMid, 
       box = c(xmin, xmax, ymin, ymax), xLeft = xLeft, xRight = xRight,
       yTop = yTop, yBot = yBot);
  
}

