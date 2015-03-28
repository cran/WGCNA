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


.heatmapWithLegend = function(data, signed, colors, naColor = "grey", zlim = NULL, 
                     reverseRows = TRUE,
                     plotLegend = TRUE,
                     keepLegendSpace = plotLegend,
                     cex.legend = 1, 
                     legendShrink = 0.94,
                     legendSpace = 0.10,
                     legendWidth = 0.02,
                     legendGap = 0.02,
                     frame = TRUE,
                     frameTicks = FALSE, tickLen = 0.02,
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
  xmax = xmaxAll - legendSpace * (xmaxAll-xminAll);

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
      yminL = yminAll + (1-legendShrink) * (ymaxAll - yminAll);
      ymaxL = ymaxAll - (1-legendShrink) * (ymaxAll - yminAll);
      xminL = xmaxAll - (xmaxAll - xminAll) * (legendSpace - legendGap)
      xmaxL = xmaxAll - (xmaxAll - xminAll) * (legendSpace - legendGap - legendWidth)
    
      tickVal = .autoTicks(zlim[1], zlim[2]);
      tickY = (tickVal - zlim[1]) / (zlim[2] - zlim[1]) * (ymaxL - yminL) + yminL;
      nTicks = length(tickVal);
    
      # Ticks:
      for (t in 1:nTicks)
        lines(c(xmaxL, xmaxL + (xmaxAll - xminAll) * 0.8 * tickLen), c(tickY[t], tickY[t]));
      text(rep(xmaxL + (xmaxAll - xminAll) * tickLen), tickY, tickVal, adj = c(0, 0.5), cex = cex.legend,
           xpd = TRUE);
    
      # Fill with color:
      nColors = length(colors);
      ybl = (ymaxL-yminL)/nColors * (0:(nColors-1)) + yminL;
      ytl = (ymaxL-yminL)/nColors * (1:nColors) + yminL;
      rect(xleft = rep(xminL, nColors), xright = rep(xmaxL, nColors),
           ybottom = ybl, ytop = ytl, col = colors, border = colors);
    
      lines(c(xminL, xmaxL, xmaxL, xminL, xminL), c(yminL, yminL, ymaxL, ymaxL, yminL) );
  }

  list(xMid = xMid, yMid = if (reverseRows) .reverseVector(yMid) else yMid, 
       box = c(xmin, xmax, ymin, ymax), xLeft = xLeft, xRight = xRight,
       yTop = yTop, yBot = yBot);
  
}

