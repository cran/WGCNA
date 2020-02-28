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
                            ## These dimensions are in inches
                            tickLen = 0.09,
                            tickGap = 0.04,
                            minBarWidth = 0.09,
                            maxBarWidth = Inf,
                            mar = c(0.5, 0.2, 0.5, 0.1))
{
  par(mar = mar);
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "");
  box = par("usr");
  tickVal = .autoTicks(lim[1], lim[2]);
  pin = par("pin");
  xrange = box[2] - box[1];
  tickLen.usr = tickLen/pin[1] * xrange
  tickGap.usr = tickGap/pin[1] * xrange
  minBarWidth.usr = minBarWidth/pin[1] * xrange
  maxBarWidth.usr = maxBarWidth/pin[1] * xrange
  maxTickWidth = max(strwidth(tickVal));
  if (maxTickWidth + tickLen.usr + tickGap.usr > box[2]-box[1]-minBarWidth.usr) 
     warning("Some tick labels will be truncated.");
  xMax = max(box[2]-maxTickWidth - tickLen.usr - tickGap.usr, box[1] + minBarWidth.usr);
  if (xMax - box[1] > maxBarWidth.usr) xMax = box[1] + maxBarWidth.usr;
  .plotColorLegend(box[1], xMax,
                   box[3], box[4], 
                   colors = colors,
                   lim = lim,
                   tickLen.usr = tickLen.usr,
                   tickGap.usr = tickGap.usr);
}

.plotColorLegend = function(xmin, xmax, ymin, ymax,
                            colors,
                            tickLen.usr = 0.7* strwidth("M"),
                            tickGap.usr = 0.3 * strwidth("M"),
                            lim, cex.legend = 1)
{
      tickVal = .autoTicks(lim[1], lim[2]);
      tickY = (tickVal - lim[1]) / (lim[2] - lim[1]) * (ymax - ymin) + ymin;
      nTicks = length(tickVal);

      # Ticks:
      for (t in 1:nTicks)
        lines(c(xmax, xmax + tickLen.usr), c(tickY[t], tickY[t]), xpd = TRUE);
      text(rep(xmax + tickLen.usr + tickGap.usr), tickY, tickVal, adj = c(0, 0.5), cex = cex.legend,
           xpd = TRUE);

      # Fill with color:
      nColors = length(colors);
      ybl = (ymax-ymin)/nColors * (0:(nColors-1)) + ymin;
      ytl = (ymax-ymin)/nColors * (1:nColors) + ymin;
      rect(xleft = rep(xmin, nColors), xright = rep(xmax, nColors),
           ybottom = ybl, ytop = ytl, col = colors, border = colors, xpd = TRUE);

      lines(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin), xpd = TRUE );
}

.boxDimensionsForHeatmapWithLegend = function(
                     data,
                     plotLegend = TRUE,
                     keepLegendSpace = plotLegend,
                     cex.legend = 1,
                     legendShrink = 0.94,
                     ## The following arguments are now in inches
                     legendSpace = 0.5,
                     legendWidth = 0.13,
                     legendGap = 0.09, 
                     startTempPlot = TRUE,
                     plotDevice = "pdf",
                     plotDeviceOptions = list(),
                     width = 7, height = 7,...)
{
  data = as.matrix(data); nCols = ncol(data); nRows = nrow(data);

  if (startTempPlot)
  {
    if (!is.null(plotDevice))
    {
      if (plotDevice == "x11") 
      {
        do.call(match.fun(plotDevice), c(list(width = width, height = height), plotDeviceOptions));
        on.exit(dev.off());
      } else {
        file = tempfile();
        do.call(match.fun(plotDevice), c(list(file = file, width = width, height = height), plotDeviceOptions))
        on.exit({ dev.off(); unlink(file)});
      }
      par(mar = par("mar"));
    }
    barplot(1, col = "white", border = "white", axisnames = FALSE,
                  axes = FALSE, ...);
  }
  pin = par("pin");
  box = par("usr");
  xminAll = box[1];
  xmaxAll = box[2];
  yminAll = box[3];
  ymaxAll = box[4];

  legendSpace.usr = legendSpace/pin[1] * (xmaxAll-xminAll);
  legendWidth.usr = legendWidth/pin[1] * (xmaxAll-xminAll);
  legendGap.usr = legendGap/pin[1] * (xmaxAll-xminAll);

  if (!keepLegendSpace && !plotLegend)
  {
     legendSpace.usr = 0;
     legendWidth.usr = 0;
     legendGap.usr = 0;
  }

  ymin = yminAll;
  ymax = ymaxAll;
  xmin = xminAll;
  xmax = xmaxAll - legendSpace.usr;
  if (xmax < xmin) stop("'legendSpace is too large, not enough space for the heatmap.");
  xStep = (xmax - xmin)/nCols;
  xLeft = xmin + c(0:(nCols-1)) * xStep;
  xRight = xLeft + xStep;
  xMid = (xLeft + xRight)/2;

  yStep = (ymax - ymin)/nRows; yBot  = ymin + c(0:(nRows-1)) * yStep;
  yTop  = yBot + yStep; yMid = c(yTop+ yBot)/2;

  list(xMin = xmin, xMax = xmax, yMin = ymin, yMax = ymax,
       xLeft = xLeft, xRight = xMid, xMid = xMid,
       yTop = yTop, yMid = yMid, yBottom = yBot);
}


.heatmapWithLegend = function(data, signed, 
                     colorMatrix = NULL,
                     colors, naColor = "grey", zlim = NULL, 
                     reverseRows = TRUE,
                     plotLegend = TRUE,
                     keepLegendSpace = plotLegend,
                     cex.legend = 1, 
                     legendShrink = 0.94,
                     ## The following arguments are now in inches
                     legendSpace = 0.5,   
                     legendWidth = 0.13,
                     legendGap = 0.09,
                     frame = TRUE,
                     frameTicks = FALSE, tickLen = 0.09,
                     ...)
{
 
  if (length(naColor)==0) naColor = 0;  ### Means transparent (as opposed to white) color.
  data = as.matrix(data); nCols = ncol(data); nRows = nrow(data);
  if (is.null(zlim)) 
  {
    zlim = range(data, na.rm = TRUE);
    if (signed) zlim = c(-max(abs(zlim)), max(abs(zlim)));
  }

  barplot(1, col = "white", border = "white", axisnames = FALSE,
                  axes = FALSE, ...);

  pin = par("pin");
  box = par("usr");
  xminAll = box[1]; 
  xmaxAll = box[2]; 
  yminAll = box[3]; 
  ymaxAll = box[4]; 

  legendSpace.usr = legendSpace/pin[1] * (xmaxAll-xminAll);
  legendWidth.usr = legendWidth/pin[1] * (xmaxAll-xminAll);
  legendGap.usr = legendGap/pin[1] * (xmaxAll-xminAll);
  tickLen.usr = tickLen/pin[1] * (xmaxAll-xminAll);

  if (!keepLegendSpace && !plotLegend)
  {
     legendSpace.usr = 0;
     legendWidth.usr = 0;
     legendGap.usr = 0;
  }

  ymin = yminAll; 
  ymax = ymaxAll; 
  xmin = xminAll; 
  xmax = xmaxAll - legendSpace.usr;
  if (xmax < xmin) stop("'legendSpace is too large, not enough space for the heatmap."); 

  xStep = (xmax - xmin)/nCols; 
  xLeft = xmin + c(0:(nCols-1)) * xStep;
  xRight = xLeft + xStep; 
  xMid = (xLeft + xRight)/2;

  yStep = (ymax - ymin)/nRows; yBot  = ymin + c(0:(nRows-1)) * yStep;
  yTop  = yBot + yStep; yMid = c(yTop+ yBot)/2;

  
  if (is.null(colorMatrix))
    colorMatrix = numbers2colors(data, signed, colors = colors, lim = zlim, naColor = naColor)
  dim(colorMatrix) = dim(data);
  if (reverseRows)
    colorMatrix = .reverseRows(colorMatrix);
  for (c in 1:nCols)
  {
    rect(xleft = rep(xLeft[c], nRows), xright = rep(xRight[c], nRows),
         ybottom = yBot, ytop = yTop, col = ifelse(colorMatrix[, c]==0, 0, colorMatrix[, c]), 
                border = ifelse(colorMatrix[, c]==0, 0, colorMatrix[, c]));
    ## Note: the ifelse seems superfluous here but it essentially converts a potentially character "0" to the number 0
    ## which the plotting system should understand as transparent color.
  }

  if (frame) lines( c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin) );

  if (plotLegend)
  {
      # Now plot the legend.
      .plotColorLegend(xmin = xmaxAll - (legendSpace.usr - legendGap.usr),
                       xmax = xmaxAll - (legendSpace.usr - legendGap.usr - legendWidth.usr),
                       ymin = yminAll + (1-legendShrink) * (ymaxAll - yminAll),
                       ymax =  ymaxAll - (1-legendShrink) * (ymaxAll - yminAll),
                       lim = zlim,
                       colors = colors,
                       tickLen.usr = tickLen.usr,
                       cex.legend = cex.legend);
    
  }

  list(xMid = xMid, yMid = if (reverseRows) rev(yMid) else yMid, 
       box = c(xmin, xmax, ymin, ymax), xLeft = xLeft, xRight = xRight,
       yTop = yTop, yBot = yBot);
  
}

