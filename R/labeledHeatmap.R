#---------------------------------------------------------------------------------------------------------
# labeledHeatmap.R
#---------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#
# .reverseRows = function(Matrix)
#
#--------------------------------------------------------------------------
#


.reverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,];
  #Matrix
}

.extend = function(x, n)
{
  nRep = ceiling(n/length(x));
  rep(x, nRep)[1:n];
}

# Adapt a numeric index to a subset
# Aim: if 'index' is a numeric index of special entries of a vector,
#    create a new index that references 'subset' elements of the vector  
.restrictIndex = function(index, subset)
{
  out = match(index, subset);
  out[!is.na(out)];
}

  
#--------------------------------------------------------------------------
#
# labeledHeatmap
#
#--------------------------------------------------------------------------
# This function plots a heatmap of the specified matrix 
# and labels the x and y axes wit the given labels.
# It is assumed that the number of entries in xLabels and yLabels is consistent 
# with the dimensions in.
# If colorLabels==TRUE, the labels are not printed and instead interpreted as colors --
#  -- a simple symbol with the appropriate color is printed instead of the label.
# The x,yLabels are expected to have the form "..color" as in "MEgrey" or "PCturquoise".
# xSymbol, ySymbols are additional markers that can be placed next to color labels

labeledHeatmap = function (
  Matrix, 
  xLabels, yLabels = NULL, 
  xSymbols = NULL, ySymbols = NULL, 
  colorLabels = NULL, 
  xColorLabels = FALSE, yColorLabels = FALSE,
  checkColorsValid = TRUE,
  invertColors = FALSE, 
  setStdMargins = TRUE,
  xLabelsPosition = "bottom",
  xLabelsAngle = 45,
  xLabelsAdj = 1,
  yLabelsPosition = "left",
  xColorWidth = 2*strheight("M"),
  yColorWidth = 2*strwidth("M"),
  xColorOffset = strheight("M")/3, 
  yColorOffset = strwidth("M")/3,
  # Content of heatmap
  colorMatrix = NULL,
  colors = NULL, 
  naColor = "grey",
  textMatrix = NULL, cex.text = NULL, 
  textAdj = c(0.5, 0.5),
  # labeling of rows and columns
  cex.lab = NULL, 
  cex.lab.x = cex.lab,
  cex.lab.y = cex.lab,
  colors.lab.x = 1,
  colors.lab.y = 1,
  font.lab.x = 1,
  font.lab.y = 1,
  bg.lab.x = NULL,
  bg.lab.y = NULL,
  x.adj.lab.y = 1,
  plotLegend = TRUE, 
  keepLegendSpace = plotLegend,
  # Separator line specification                   
  verticalSeparator.x = NULL,
  verticalSeparator.col = 1,  
  verticalSeparator.lty = 1,
  verticalSeparator.lwd = 1,
  verticalSeparator.ext = 0,
  verticalSeparator.interval = 0,

  horizontalSeparator.y = NULL,
  horizontalSeparator.col = 1,  
  horizontalSeparator.lty = 1,
  horizontalSeparator.lwd = 1,
  horizontalSeparator.ext = 0,
  horizontalSeparator.interval = 0,
  # optional restrictions on which rows and columns to actually show
  showRows = NULL,
  showCols = NULL,
  # Other arguments...
  ... ) 
{
  textFnc = match.fun("text");
  if (!is.null(colorLabels)) {xColorLabels = colorLabels; yColorLabels = colorLabels; }
  
  if (is.null(yLabels) & (!is.null(xLabels)) & (dim(Matrix)[1]==dim(Matrix)[2])) 
    yLabels = xLabels; 

  nCols = ncol(Matrix);
  nRows = nrow(Matrix);

  if (length(xLabels)!=nCols) 
    stop("Length of 'xLabels' must equal the number of columns in 'Matrix.'");

  if (length(yLabels)!=nRows)
    stop("Length of 'yLabels' must equal the number of rows in 'Matrix.'");

  if (is.null(showRows)) showRows = c(1:nRows);
  if (is.null(showCols)) showCols = c(1:nCols);

  nShowCols = length(showCols);
  nShowRows = length(showRows);

  if (nShowCols==0) stop("'showCols' is empty.");
  if (nShowRows==0) stop("'showRows' is empty.");

  if (checkColorsValid)
  {
    xValidColors = !is.na(match(substring(xLabels, 3), colors()));
    yValidColors = !is.na(match(substring(yLabels, 3), colors()));
  } else {
    xValidColors = rep(TRUE, length(xLabels));
    yValidColors = rep(TRUE, length(yLabels));
  }

  if (sum(xValidColors)>0) xColorLabInd = xValidColors[showCols]
  if (sum(!xValidColors)>0) xTextLabInd = !xValidColors[showCols]

  if (sum(yValidColors)>0) yColorLabInd = yValidColors[showRows]
  if (sum(!yValidColors)>0) yTextLabInd = !yValidColors[showRows]

  if (setStdMargins)
  {
    if (xColorLabels & yColorLabels)
    {
      par(mar=c(2,2,3,5)+0.2);
    } else {
      par(mar = c(7,7,3,5)+0.2);
    }
  }

  xLabels.show = xLabels[showCols];
  yLabels.show = yLabels[showRows];

  if (!is.null(xSymbols))
  {
     if (length(xSymbols)!=nCols)
       stop("When 'xSymbols' are given, their length must equal the number of columns in 'Matrix.'");
     xSymbols.show = xSymbols[showCols];
  } else 
     xSymbols.show = NULL;

  if (!is.null(ySymbols))
  {
     if (length(ySymbols)!=nRows)
       stop("When 'ySymbols' are given, their length must equal the number of rows in 'Matrix.'");
     ySymbols.show = ySymbols[showRows];
  } else 
     ySymbols.show = NULL;

  xLabPos = charmatch(xLabelsPosition, c("bottom", "top"));
  if (is.na(xLabPos))
    stop("Argument 'xLabelsPosition' must be (a unique abbreviation of) 'bottom', 'top'");

  yLabPos = charmatch(yLabelsPosition, c("left", "right"));
  if (is.na(yLabPos))
    stop("Argument 'yLabelsPosition' must be (a unique abbreviation of) 'left', 'right'");

  if (is.null(colors)) colors = heat.colors(30);
  if (invertColors) colors = rev(colors);

  labPos = .heatmapWithLegend(Matrix[showRows, showCols, drop = FALSE], 
              signed = FALSE, colorMatrix = colorMatrix, colors = colors, naColor = naColor, 
              cex.legend = cex.lab, plotLegend = plotLegend,  keepLegendSpace = keepLegendSpace, ...)
  plotbox = labPos$box;
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;
  # The positions below are for showCols/showRows-restriceted data
  xLeft = labPos$xLeft;
  xRight = labPos$xRight;
  yTop = labPos$yTop;
  yBot = labPos$yBot;

  xspacing = labPos$xMid[2] - labPos$xMid[1];
  yspacing = abs(labPos$yMid[2] - labPos$yMid[1]);

  offsetx = .extend(xColorOffset, nCols)[showCols]
  offsety = .extend(yColorOffset, nRows)[showRows]
  xColW = xColorWidth;
  yColW = yColorWidth;

  # Additional angle-dependent offsets for x axis labels
  textOffsetY = strheight("M") * cos(xLabelsAngle/180 * pi);

  if (any(xValidColors)) offsetx = offsetx + xColW;
  if (any(yValidColors)) offsety = offsety + yColW;

  # Create the background for column and row labels.

  extension.left = par("mai")[2] * # left margin width in inches
                   par("cxy")[1] / par("cin")[1]   # character size in user corrdinates/character size in inches

  extension.right = par("mai")[4] * # right margin width in inches
                   par("cxy")[1] / par("cin")[1]   # character size in user corrdinates/character size in inches

  extension.bottom = par("mai")[1] * 
                   par("cxy")[2] / par("cin")[2]- # character size in user corrdinates/character size in inches
                      offsetx   
                     
  extension.top = par("mai")[3] * 
                   par("cxy")[2] / par("cin")[2]-   # character size in user corrdinates/character size in inches
                     offsetx

  figureBox = par("usr");
  figXrange = figureBox[2] - figureBox[1];
  figYrange = figureBox[4] - figureBox[3];
  if (!is.null(bg.lab.x))
  {
    bg.lab.x = .extend(bg.lab.x, nCols)[showCols];
    if (xLabPos==1)
    {
      y0 = ymin;
      ext = extension.bottom;
      sign = 1;
    } else {
      y0 = ymax;
      ext = extension.top;
      sign = -1;
    }
    figureDims = par("pin");
    angle = xLabelsAngle/180*pi;
    ratio = figureDims[1]/figureDims[2] * figYrange/figXrange;
    ext.x = -sign * ext * 1/tan(angle)/ratio;
    ext.y = sign * ext * sign(sin(angle))

    #offset = (sum(xValidColors)>0) * xColW + offsetx + textOffsetY;
    offset = offsetx + textOffsetY;

    for (cc in 1:nShowCols)
       polygon(x = c(xLeft[cc], xLeft[cc], xLeft[cc] + ext.x, xRight[cc] + ext.x, xRight[cc], xRight[cc]),
               y = c(y0, y0-sign*offset[cc], y0-sign*offset[cc] - ext.y, y0-sign*offset[cc] - ext.y, 
                     y0-sign*offset[cc], y0), 
               border = bg.lab.x[cc], col = bg.lab.x[cc], xpd = TRUE);
  }

  if (!is.null(bg.lab.y))
  {
    bg.lab.y = .extend(bg.lab.y, nRows)
    reverseRows = TRUE;
    if (reverseRows) bg.lab.y = rev(bg.lab.y);
    bg.lab.y = bg.lab.y[showRows];

    if (yLabPos==1)
    {
      xl = xmin-extension.left;
      xr = xmin;
    } else {
      xl = xmax;
      xr = xmax + extension.right;
    }
    for (r in 1:nShowRows)
      rect(xl, yBot[r], xr, yTop[r],
           col = bg.lab.y[r], border = bg.lab.y[r], xpd = TRUE);
  }

  colors.lab.x = .extend(colors.lab.x, nCols)[showCols];
  font.lab.x = .extend(font.lab.x, nCols)[showCols];
  # Write out labels
  if (sum(!xValidColors)>0)
  {
    xLabYPos = if(xLabPos==1) ymin - offsetx- textOffsetY else ymax + offsetx + textOffsetY;
    if (is.null(cex.lab)) cex.lab = 1;
    mapply(textFnc, x = labPos$xMid[xTextLabInd], 
           y = xLabYPos, labels = xLabels.show[xTextLabInd],
           col = colors.lab.x[xTextLabInd],
           font = font.lab.x[xTextLabInd],
           MoreArgs = list(srt = xLabelsAngle, 
          adj = xLabelsAdj, xpd = TRUE, cex = cex.lab.x));
  }
  if (sum(xValidColors)>0)
  {
    baseY = if (xLabPos==1) ymin-offsetx else  ymax + offsetx;
    deltaY = if (xLabPos==1) xColW else -xColW;
    rect(xleft = labPos$xMid[xColorLabInd] - xspacing/2, ybottom = baseY[xColorLabInd],
         xright = labPos$xMid[xColorLabInd] + xspacing/2, ytop = baseY[xColorLabInd] + deltaY,
         density = -1,  col = substring(xLabels.show[xColorLabInd], 3), 
         border = substring(xLabels.show[xColorLabInd], 3), xpd = TRUE)
    if (!is.null(xSymbols))
      mapply(textFnc, x = labPos$xMid[xColorLabInd], 
             y = baseY[xColorLabInd] -textOffsetY - sign(deltaY)* strwidth("M")/3, 
             labels = xSymbols.show[xColorLabInd],
             col = colors.lab.x[xColorLabInd],
             font = font.lab.x[xColorLabInd],
              MoreArgs = list( adj = xLabelsAdj, 
             xpd = TRUE, srt = xLabelsAngle, cex = cex.lab.x));
  }
  x.adj.lab.y = .extend(x.adj.lab.y, nRows)[showRows]
  if (yLabPos==1)
  {
    marginWidth = par("mai")[2] / par("pin")[1] * xrange
  } else {
    marginWidth = par("mai")[4] / par("pin")[1] * xrange
  }
  xSpaceForYLabels = marginWidth-2*strwidth("M")/3 - ifelse(yValidColors[showRows], yColW, 0);
  xPosOfYLabels.relative = xSpaceForYLabels * (1-x.adj.lab.y) + offsety

  colors.lab.y = .extend(colors.lab.y, nRows)[showRows];
  font.lab.y = .extend(font.lab.y, nRows)[showRows];

  if (sum(!yValidColors)>0)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    if (yLabPos==1)
    {
      x = xmin - strwidth("M")/3 - xPosOfYLabels.relative[yTextLabInd]
      adj = x.adj.lab.y[yTextLabInd]
    } else {
      x = xmax + strwidth("M")/3 + xPosOfYLabels.relative[yTextLabInd];
      adj = 1-x.adj.lab.y[yTextLabInd];
    }
    mapply(textFnc, y = labPos$yMid[yTextLabInd], labels = yLabels.show[yTextLabInd],
           adj = lapply(adj, c, 0.5),
           x = x,
           col = colors.lab.y[yTextLabInd],
           font = font.lab.y[yTextLabInd],
           MoreArgs = list(srt = 0, xpd = TRUE, cex = cex.lab.y));
  } 
  if (sum(yValidColors)>0)
  {
    if (yLabPos==1)
    {
      xl = xmin-offsety;
      xr = xmin-offsety + yColW;
      xtext = xmin - strwidth("M")/3 - xPosOfYLabels.relative[yColorLabInd];
      adj = x.adj.lab.y[yColorLabInd]
    } else {
      xl = xmax + offsety - yColW;
      xr = xmax + offsety;
      xtext = xmin + strwidth("M")/3 + xPosOfYLabels.relative[yColorLabInd]
      adj = 1-x.adj.lab.y[yColorLabInd];
    }

    rect(xleft = xl[yColorLabInd], ybottom = rev(labPos$yMid[yColorLabInd]) - yspacing/2,
         xright = xr[yColorLabInd], ytop = rev(labPos$yMid[yColorLabInd]) + yspacing/2, 
         density = -1,  col = substring(rev(yLabels.show[yColorLabInd]), 3), 
         border = substring(rev(yLabels.show[yColorLabInd]), 3), xpd = TRUE)
    #for (i in yColorLabInd)
    #{
    #  lines(c(xmin- offsetx, xmin- offsetx+yColW), y = rep(labPos$yMid[i] - yspacing/2, 2), col = i, xpd = TRUE)
    #  lines(c(xmin- offsetx, xmin- offsetx+yColW), y = rep(labPos$yMid[i] + yspacing/2, 2), col = i, xpd = TRUE)
    #}
    if (!is.null(ySymbols))
      mapply(textFnc, y = labPos$yMid[yColorLabInd], labels = ySymbols.show[yColorLabInd],
             adj = lapply(adj, c, 0.5),
             x = xtext, col = colors.lab.y[yColorLabInd], 
             font = font.lab.y[yColorLabInd],
          MoreArgs = list(srt = 0, xpd = TRUE, cex = cex.lab.y));
  }

  # Draw separator lines, if requested

  showCols.ext = c(if (1 %in% showCols) 0 else NULL, showCols);
  showCols.shift = if (0 %in% showCols.ext) 1 else 0;

  if (length(verticalSeparator.x) > 0)
  {
    if (any(verticalSeparator.x < 0 | verticalSeparator.x > nCols))
      stop("If given. 'verticalSeparator.x' must all be between 0 and the number of columns.");
    colSepShowIndex = which(verticalSeparator.x %in% showCols.ext);
    verticalSeparator.x.show = .restrictIndex(verticalSeparator.x, showCols.ext)-showCols.shift;
  } else if (verticalSeparator.interval > 0)
  {
    verticalSeparator.x.show = verticalSeparator.x = 
           seq(from = verticalSeparator.interval, by = verticalSeparator.interval,
                                    length.out = floor(length(showCols)/verticalSeparator.interval));
    colSepShowIndex = 1:length(verticalSeparator.x);
  } else 
    verticalSeparator.x.show = NULL;

  if (length(verticalSeparator.x.show) > 0)
  {
    nLines = length(verticalSeparator.x);
    vs.col = .extend(verticalSeparator.col, nLines)[colSepShowIndex];
    vs.lty = .extend(verticalSeparator.lty, nLines)[colSepShowIndex];
    vs.lwd = .extend(verticalSeparator.lwd, nLines)[colSepShowIndex];
    vs.ext = .extend(verticalSeparator.ext, nLines)[colSepShowIndex];

    x.lines = ifelse(verticalSeparator.x.show>0, labPos$xRight[verticalSeparator.x.show], labPos$xLeft[1]);
    nLines.show = length(verticalSeparator.x.show);
    for (l in 1:nLines.show)
      lines(rep(x.lines[l], 2), c(ymin, ymax), col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l]);

    angle = xLabelsAngle/180*pi;
    if (angle==0) angle = pi/2;
    if (xLabelsPosition =="bottom") 
    {
      sign = 1;
      y0 = ymin;
      ext = extension.bottom;
    } else {
      sign = -1;
      y0 = ymax;
      ext = extension.top;
    }
    figureDims = par("pin");
    ratio = figureDims[1]/figureDims[2] * figYrange/figXrange;
    ext.x = -sign * ext * 1/tan(angle)/ratio;
    ext.y = sign * ext * sign(sin(angle))
    #offset = (sum(xValidColors)>0) * xColW + offsetx + textOffsetY;
    offset = offsetx + textOffsetY;
    for (l in 1:nLines.show)
         lines(c(x.lines[l], x.lines[l], x.lines[l] + vs.ext[l] * ext.x[l]), 
               c(y0, y0-sign*offset[l], y0-sign*offset[l] - vs.ext[l] * ext.y[l]),  
                 col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE);
  }

  showRows.ext = c(if (1 %in% showRows) 0 else NULL, showRows);
  showRows.shift = if (0 %in% showRows.ext) 1 else 0;

  if (length(horizontalSeparator.y) >0)
  {
    if (any(horizontalSeparator.y < 0 | horizontalSeparator.y > nRows))
      stop("If given. 'horizontalSeparator.y' must all be between 0 and the number of rows.");
    rowSepShowIndex = which( horizontalSeparator.y %in% showRows.ext);
    horizontalSeparator.y.show = .restrictIndex(horizontalSeparator.y, showRows.ext)-showRows.shift;
  } else if (horizontalSeparator.interval > 0)
  {
    horizontalSeparator.y.show = horizontalSeparator.y = 
            seq(from = horizontalSeparator.interval, by = horizontalSeparator.interval,
                                    length.out = floor(length(showRows)/horizontalSeparator.interval));
    rowSepShowIndex = 1:length(horizontalSeparator.y);
  } else 
    horizontalSeparator.y.show = NULL;
  
  if (length(horizontalSeparator.y.show) > 0)
  {
    reverseRows = TRUE;
    if (reverseRows) 
    {
      horizontalSeparator.y.show = nShowRows - horizontalSeparator.y.show+1;
      y.lines = ifelse( horizontalSeparator.y.show <=nShowRows, 
                               labPos$yBot[horizontalSeparator.y.show], labPos$yTop[nShowRows]);
    } else {
      y.lines = ifelse( horizontalSeparator.y.show > 0, labPos$yBot[horizontalSeparator.y.show], labPos$yTop[1]);
    }
    nLines = length(horizontalSeparator.y);
    vs.col = .extend(horizontalSeparator.col, nLines)[rowSepShowIndex];
    vs.lty = .extend(horizontalSeparator.lty, nLines)[rowSepShowIndex];
    vs.lwd = .extend(horizontalSeparator.lwd, nLines)[rowSepShowIndex];
    vs.ext = .extend(horizontalSeparator.ext, nLines)[rowSepShowIndex];
    nLines.show = length(horizontalSeparator.y.show);
    for (l in 1:nLines.show)
    {
      if (yLabPos==1)
      {
         xl = xmin-vs.ext[l]*extension.left;
         xr = xmax;
      } else {
         xl = xmin;
         xr = xmax + vs.ext[l]*extension.right;
      }
  
      lines(c(xl, xr), rep(y.lines[l], 2), 
            col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE);
    }
  }

  if (!is.null(textMatrix))
  {
    if (is.null(cex.text)) cex.text = par("cex");
    if (is.null(dim(textMatrix)))
      if (length(textMatrix)==prod(dim(Matrix))) dim(textMatrix)=dim(Matrix);
    if (!isTRUE(all.equal(dim(textMatrix), dim(Matrix))))
      stop("labeledHeatmap: textMatrix was given, but has dimensions incompatible with Matrix.");
    for (rw in 1:nShowRows)
      for (cl in 1:nShowCols)
      {
        text(labPos$xMid[cl], labPos$yMid[rw],
             as.character(textMatrix[showRows[rw],showCols[cl]]), xpd = TRUE, cex = cex.text, adj = textAdj);
      }
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
  invisible(labPos)
}

#===================================================================================================
#
# multi-page labeled heatmap
#
#===================================================================================================

labeledHeatmap.multiPage = function(
   # Input data and ornament[s
   Matrix,
   xLabels, yLabels = NULL,
   xSymbols = NULL, ySymbols = NULL,
   textMatrix = NULL,

   # Paging options
   rowsPerPage = NULL, maxRowsPerPage = 20,
   colsPerPage = NULL, maxColsPerPage = 10,
   addPageNumberToMain = TRUE,

   # Further arguments to labeledHeatmap
   zlim = NULL,
   signed = TRUE,
   main = "",

  verticalSeparator.x = NULL,
  verticalSeparator.col = 1,
  verticalSeparator.lty = 1,
  verticalSeparator.lwd = 1,
  verticalSeparator.ext = 0,

  horizontalSeparator.y = NULL,
  horizontalSeparator.col = 1,
  horizontalSeparator.lty = 1,
  horizontalSeparator.lwd = 1,
  horizontalSeparator.ext = 0,

   ...)
{

  nr = nrow(Matrix);
  nc = ncol(Matrix);

  if (is.null(rowsPerPage))
  {
    nPages.rows = ceiling(nr/maxRowsPerPage);
    rowsPerPage = allocateJobs(nr, nPages.rows);
  } else 
    nPages.rows = length(rowsPerPage);

  if (is.null(colsPerPage))
  {
    nPages.cols = ceiling(nc/maxColsPerPage);
    colsPerPage = allocateJobs(nc, nPages.cols);
  } else 
    nPages.cols = length(colsPerPage);

  if (is.null(zlim)) 
  {
    zlim = range(Matrix, na.rm = TRUE)
    if (signed) zlim = c(-max(abs(zlim)), max(abs(zlim)));
  }

  if (!is.null(verticalSeparator.x))
  {
    nvs = length(verticalSeparator.x);
    verticalSeparator.col= .extend(verticalSeparator.col, nvs);
    verticalSeparator.lty= .extend(verticalSeparator.lty, nvs);
    verticalSeparator.lwd= .extend(verticalSeparator.lwd, nvs);
    verticalSeparator.ext= .extend(verticalSeparator.ext, nvs);
  }
  
  if (!is.null(horizontalSeparator.y))
  {
    nhs = length(horizontalSeparator.y);
    horizontalSeparator.col= .extend(horizontalSeparator.col, nhs);
    horizontalSeparator.lty= .extend(horizontalSeparator.lty, nhs);
    horizontalSeparator.lwd= .extend(horizontalSeparator.lwd, nhs);
    horizontalSeparator.ext= .extend(horizontalSeparator.ext, nhs);
  }
  

  page = 1;
  multiPage = (nPages.cols > 1 | nPages.rows > 1)

  for (page.col in 1:nPages.cols) for (page.row in 1:nPages.rows)
  {
    rows = rowsPerPage[[page.row]];
    cols = colsPerPage[[page.col]];
    if (!is.null(verticalSeparator.x))
    {
      keep.vs = verticalSeparator.x %in% cols;
    } else 
      keep.vs = numeric(0);
    if (!is.null(horizontalSeparator.y))
    {
      keep.hs = horizontalSeparator.y %in% rows;
    } else 
      keep.hs = numeric(0);

    main.1 = main;
    if (addPageNumberToMain & multiPage) main.1 = spaste(main, "(page ", page, ")");
    labeledHeatmap(Matrix = Matrix[rows, cols, drop = FALSE],
                   xLabels = xLabels[cols], xSymbols = xSymbols[cols],
                   yLabels = yLabels[rows], ySymbols = ySymbols[rows],
                   textMatrix = textMatrix[rows, cols, drop = FALSE],
                   zlim = zlim, main = main.1, 
                   verticalSeparator.x = verticalSeparator.x[keep.vs] - min(cols) + 1,
                   verticalSeparator.col= verticalSeparator.col[keep.vs],
                   verticalSeparator.lty= verticalSeparator.lty[keep.vs],
                   verticalSeparator.lwd= verticalSeparator.lwd[keep.vs],
                   verticalSeparator.ext= verticalSeparator.ext[keep.vs],
 
                   horizontalSeparator.y = horizontalSeparator.y[keep.hs] - min(rows) + 1,
                   horizontalSeparator.col= horizontalSeparator.col[keep.hs],
                   horizontalSeparator.lty= horizontalSeparator.lty[keep.hs],
                   horizontalSeparator.lwd= horizontalSeparator.lwd[keep.hs],
                   horizontalSeparator.ext= horizontalSeparator.ext[keep.hs],
                   ...);
    page = page + 1;
  }
}
                   


