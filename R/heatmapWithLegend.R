# Replacement for the function image.plot

.autoTicks = function(min, max, maxTicks = 6, tickPos = c(1,2,5))
{
  if (max < min) { x = max; max = min; min = x }
  range = max - min;
  if (range==0) return(max);
  tick0 = range/(maxTicks+1-1e-6)
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
                            mar = c(0.5, 0.2, 0.5, 0.1),
                            lab = "",
                            horizontal = FALSE,
                            ...)
{
  par(mar = mar);
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "");
  box = par("usr");
  if (horizontal) box.eff = box[c(3,4,1,2)] else box.eff = box;
  tickVal = .autoTicks(lim[1], lim[2]);
  pin = par("pin");
  pin.eff = if (horizontal) pin[c(2,1)] else pin;
  wrange = box.eff[2] - box.eff[1];
  tickLen.usr = tickLen/pin.eff[1] * wrange
  tickGap.usr = tickGap/pin.eff[1] * wrange
  minBarWidth.usr = minBarWidth/pin.eff[1] * wrange
  maxBarWidth.usr = maxBarWidth/pin.eff[1] * wrange
  sizeFnc = if (horizontal) strheight else strwidth;
  maxTickWidth = max(sizeFnc(tickVal));
  if (maxTickWidth + tickLen.usr + tickGap.usr > box.eff[2]-box.eff[1]-minBarWidth.usr) 
     warning("Some tick labels will be truncated.");
  haveLab = length(lab) > 0
  if (haveLab && is.character(lab)) haveLab = lab!="";
  width = max(box.eff[2]-box.eff[1]-maxTickWidth - tickLen.usr - tickGap.usr- haveLab * 3*sizeFnc("M"), minBarWidth.usr);
  if (width > maxBarWidth.usr) width = maxBarWidth.usr;
  .plotColorLegend(box[1], if (horizontal) box[2] else box[1] + width,
                   if (horizontal) box[4]-width else box[3], box[4], 
                   colors = colors,
                   lim = lim,
                   tickLen.usr = tickLen.usr, horizontal = horizontal,
                   tickGap.usr = tickGap.usr, lab = lab, ...);
}

if (FALSE)
{
   source("~/Work/RLibs/WGCNA/R/heatmapWithLegend.R")
   .plotStandaloneLegend(colors = blueWhiteRed(10), lim = c(-25, 25))
   d = matrix(rnorm(100), 10, 10);
   par(mar = c(2,2,2,0));
   
   .heatmapWithLegend(d,
                     signed = TRUE,
                     colors = blueWhiteRed(20), 
                     plotLegend = TRUE,
                     cex.legendAxis = 1,
                     legendShrink = 0.94,
                     legendLabel = "",
                     cex.legendLabel = 1)
                     ## The following arguments are now in inches
                     #legendSpace = 0.5 + (legendLabel!="") * 1.5*strheight("M",units = "inch", cex = cex.legendLabel),
                     #legendWidth = 0.13,
                     #legendGap = 0.09,
                     #frame = TRUE,
                     #frameTicks = FALSE, tickLen = 0.09);

}

.plotColorLegend = function(xmin, xmax, ymin, ymax,
                            # colors can be a vector or a matrix (in which case a matrix of colors will be plotted)
                            colors,
                            horizontal = FALSE,
### FIXME: it would be good if these could respect settings in par("mgp")
                            tickLen.usr = 0.5* (if (horizontal) strheight("M") else strwidth("M")),
                            tickGap.usr = 0.5 * (if (horizontal) strheight("M") else strwidth("M")),
                            lim, cex.axis = 1, tickLabelAngle = if (horizontal) 0 else -90,
                            lab = "", cex.lab = 1, labAngle = 0, 
                            labGap = 0.6 * (if (horizontal) strheight("M") else strwidth("M"))
                            )
{
  tickVal = .autoTicks(lim[1], lim[2]);
  nTicks = length(tickVal);

  if (horizontal) {
    lmin = xmin; lmax = xmax; 
    tmin = ymin; tmax = ymax;
  } else {
    tmin = xmin; tmax = xmax; 
    lmin = ymin; lmax = ymax;
  }
  tickPos = (tickVal - lim[1]) / (lim[2] - lim[1]) * (lmax - lmin) + lmin;
  pin = par("pin");
  box = par("usr");
  asp = pin[2]/pin[1] * ( box[2]-box[1])/(box[4] - box[3]);
  # Ticks:
  
  if (horizontal) {
    angle0 = 0;
    angle = angle0 + tickLabelAngle;
    if (angle==0) adj = c(0.5, 1) else adj = c(1, 0.5);
    for (t in 1:nTicks) 
      lines(c(tickPos[t], tickPos[t]), c(ymin, ymin - tickLen.usr), xpd = TRUE);
    text(tickPos, rep(ymin - tickLen.usr - tickGap.usr), tickVal, adj = adj, cex = cex.axis,
           xpd = TRUE, srt = angle);
    tickLabelWidth = if (angle==0) max(strheight(tickVal)) else max(strwidth(tickVal))/asp;
  } else {
    angle0 = 90;
    angle = angle0 + tickLabelAngle;
    if (angle==0) adj = c(0, 0.5) else adj = c(0.5, 1);
    for (t in 1:nTicks) 
      lines(c(xmax, xmax + tickLen.usr), c(tickPos[t], tickPos[t]), xpd = TRUE);
    text(rep(xmax + tickLen.usr + tickGap.usr), tickPos, tickVal, adj = adj, cex = cex.axis,
         xpd = TRUE, srt = angle);
    tickLabelWidth = if (angle==0) max(strwidth(tickVal)) else max(strheight(tickVal)) * asp;
  }
  # Fill with color:
  colors = as.matrix(colors);
  nColumns = ncol(colors);
  nColors = nrow(colors);
  bl = (lmax-lmin)/nColors * (0:(nColors-1)) + lmin;
  tl = (lmax-lmin)/nColors * (1:nColors) + lmin;
  wi.all = tmax - tmin;
  wi1 = wi.all/nColumns
  if (horizontal) {
    for (col in 1:nColumns)
      rect(xleft = bl, xright = tl,
         ybottom = rep(tmin + (col-1) * wi1, nColors), ytop = rep(tmin + wi1*col, nColors), 
            col = colors[, col], border = colors[, col], xpd = TRUE);
  } else {
    for (col in 1:nColumns)
       rect(xleft = rep(tmin + (col-1) * wi1, nColors), xright = rep(tmin + wi1*col, nColors),
          ybottom = bl, ytop = tl, col = colors[, col], border = colors[, col], xpd = TRUE);
  }
  # frame for the legend
  lines(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin), xpd = TRUE );

  if (nColumns > 1) for (col in 2:nColumns) 
    if (horizontal) lines(c(xmin, xmax), c(tmin + (col-1) * wi1, tmin + (col-1) * wi1)) else 
                    lines(c(tmin + (col-1) * wi1, tmin + (col-1) * wi1), c(ymin, ymax));
  # Axis label
  if (length(lab)>0 && as.character(lab) != "")
  {
    if (horizontal)
    {
      y = ymin - tickLen.usr - tickGap.usr - tickLabelWidth - labGap;
      x = (xmin + xmax)/2;
      adj = if (labAngle==0) c(0.5, 1) else c(1, 0.5)
      angle = labAngle;
      text(x, y, lab, cex = cex.lab, srt = labAngle, xpd = TRUE, adj = adj);
    } else {
      y = (ymin + ymax)/2;
      x = xmax + tickLen.usr + tickGap.usr + tickLabelWidth + labGap;
      adj = if (labAngle==0) c(0.5, 1) else c(0, 0.5);
      angle = labAngle+90;
      text(x, y, lab, cex = cex.lab, srt = labAngle+90, xpd = TRUE, adj = adj);
    }
    height = strheight(lab);
    if (!horizontal) height = height * asp;
    labelInfo = list(x = x, y = y, angle = angle, adj = adj,
                     space.usr = height, gap.usr = labGap);
  } else labelInfo = list(space.usr = 0, gap.usr = 0);
  #### FIXME: also include a component named box that gives the outer coordinates of the area used by the legend, to the
  ###best approximation. Maybe include the padding around the color bar.
  invisible(list(bar = list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                            space.usr = tmax - tmin),
       ticks = list(length.usr = tickLen.usr, gap.usr = tickGap.usr, labelSpace.usr = tickLabelWidth),
       label = labelInfo));
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
                     cex.legendAxis = 1, 
                     legendShrink = 0.94,
                     legendPosition = 0.5, ## center; 1 means at the top, 0 means at the bottom
                     legendLabel = "",
                     cex.legendLabel = 1,
                     ## The following arguments are now in inches
                     legendSpace = 0.5 + (as.character(legendLabel)!="") * 1.5*
                            strheight("M",units = "inch", cex = cex.legendLabel),   
                     legendWidth = 0.13,
                     legendGap = 0.09,
                     maxLegendSize = 4,
                     legendLengthGap = 0.15,
                     frame = TRUE,
                     frameTicks = FALSE, tickLen = 0.09,
                     tickLabelAngle = 0,
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
  maxLegendSize.usr = maxLegendSize/pin[2] * (ymaxAll-yminAll);
  legendLengthGap.usr = legendLengthGap/pin[2] * (ymaxAll-yminAll)

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
      legendSize.usr = legendShrink * (ymaxAll - yminAll);
      if (legendSize.usr > maxLegendSize.usr) legendSize.usr = maxLegendSize.usr
      if (legendLengthGap.usr > 0.5*(ymaxAll - yminAll)*(1-legendShrink)) 
          legendLengthGap.usr = 0.5*(ymaxAll - yminAll)*(1-legendShrink);
      y0 = yminAll + legendLengthGap.usr;
      y1 = ymaxAll - legendLengthGap.usr;
      movementRange = (y1-y0 - legendSize.usr);
      if (movementRange < -1e-10) {browser(".heatmapWithLegend: movementRange is negative."); movementRange = 0;}
      ymin.leg = y0 + legendPosition * movementRange;
      ymax.leg = y0 + legendPosition * movementRange + legendSize.usr
      legendPosition = .plotColorLegend(xmin = xmaxAll - (legendSpace.usr - legendGap.usr),
                       xmax = xmaxAll - (legendSpace.usr - legendGap.usr - legendWidth.usr),
                       ymin = ymin.leg,
                       ymax =  ymax.leg,
                       lim = zlim,
                       colors = colors,
                       tickLen.usr = tickLen.usr,
                       cex.axis = cex.legendAxis,
                       lab = legendLabel,
                       cex.lab = cex.legendLabel,
                       tickLabelAngle = tickLabelAngle
                       );
    
  } else legendPosition = NULL

  invisible(list(xMid = xMid, yMid = if (reverseRows) rev(yMid) else yMid, 
       box = c(xmin, xmax, ymin, ymax), xLeft = xLeft, xRight = xRight,
       yTop = yTop, yBot = yBot,
       legendPosition = legendPosition));
  
}

