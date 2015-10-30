#=================================================================================================
#
# labelPoints: label points in a scatterplot while trying to avoid labels overlapping with one another
# and with points.
#
#=================================================================================================

labelPoints = function(x, y, labels, cex = 0.7, offs = 0.01, xpd = TRUE, jiggle = 0, 
                       protectEdges = TRUE, 
                       doPlot = TRUE, ...)
{
  nPts = length(labels);
  box = par("usr");
  dims = par("pin");
  scaleX = dims[1]/(box[2] - box[1]);
  scaleY = dims[2]/(box[4] - box[3]);

  #ish = charmatch(shape, .shapes);
  #if (is.na(ish))
  #  stop(paste("Unrecognized 'shape'. Recognized values are", paste(.shapes, collapse = ", ")));

  if (par("xlog"))
  {
     xx = log10(x);
  } else
     xx = x;

  if (par("ylog"))
  {
     yy = log10(y);
  } else
     yy = y;

  xx = xx * scaleX;
  yy = yy * scaleY;

  if (jiggle > 0)
  {
    rangeX = max(xx, na.rm = TRUE) - min(xx, na.rm = TRUE)
    jx = xx + jiggle * rangeX * (runif(nPts) - 0.5);  
    rangeY = max(yy, na.rm = TRUE) - min(yy, na.rm = TRUE)
    jy = yy + jiggle * rangeY * (runif(nPts) - 0.5);  
  } else {
    jx = xx;
    jy = yy;
  }
  dx = offs;
  dy = offs;
  labWidth = strwidth(labels, cex=cex) * scaleX; 
  labHeight = strheight(labels, cex=cex) * scaleY; 
  if (nPts==0) return(0);

  if (nPts==1) 
  {
    if (protectEdges)
    {
       shift = ifelse(x - labWidth/2/scaleX < box[1], box[1] - x + labWidth/2/scaleX, 
                       ifelse(x + labWidth/2/scaleX > box[2], box[2] - x - labWidth/2/scaleX, 0));
       x = x + shift;
       # Also check the top and bottom edges
       yShift = if (y + labHeight/scaleY  + offs/scaleY > box[4])  -(labHeight + 2*offs)/scaleY else 0;
       y = y + yShift
    } 
    text(x, y + labHeight/2/scaleY + offs/scaleY, labels, cex = cex, xpd = xpd, adj = c(0.5, 0.5), ...)
    return (0);
  }

  xMat = cbind(xx,yy);
  jxMat = cbind(jx, jy);
  distX = as.matrix(dist(jx));
  distY = as.matrix(dist(jy));

  dir = matrix(0, nPts, 2);

  d0SqX = (labWidth+2*offs)^2
  d0SqY = (labHeight + 2*offs)^2; 
  for (p in 1:nPts)
  {
    difs = matrix(jxMat[p, ], nPts, 2, byrow = TRUE) - jxMat;
    difSc = difs / sqrt(matrix(apply(difs^2, 1, sum, na.rm = TRUE), nPts, 2));
    difSx = rbind(difSc, c(0,1));
    difSx[p, ] = 0;
    w = c(exp(-distX[,p]^4 / d0SqX[p]^2 - distY[,p]^4/d0SqY^2));
    w[distX[, p]==0 & distY[,p]==0] = 0;
    w = c(w, 0.01);
    dir[p, ] = apply(difSx * matrix(w, (nPts+1), 2), 2, sum, na.rm = TRUE) / sum(w, na.rm = TRUE)
    
    if (sum(abs(dir[p, ]))==0) dir[p, ] = runif(2);
  }

  scDir = dir / sqrt(matrix(apply(dir^2, 1, sum, na.rm = TRUE), nPts, 2));
  offsMat = cbind(labWidth/2 + offs, labHeight/2 + offs)
  Rmat = abs(scDir / offsMat); 
  ind = Rmat[, 1] > Rmat[, 2]; # This is an indicator of whether the labels touch the vertical (TRUE ) or
                               # horizontal (FALSE) edge of the square around the point

  # These are preliminary text coordinates relative to their points.
  dx = offsMat[, 1] * sign(scDir[, 1])
  dx[!ind] = scDir[!ind, 1] * offsMat[!ind, 2]/abs(scDir[!ind,2]);
  dy = offsMat[, 2] * sign(scDir[, 2]);
  dy[ind] = scDir[ind, 2] * offsMat[ind, 1]/abs(scDir[ind,1]);

  # Absolute coordinates
  xt = (xx + dx)/scaleX;
  yt = (yy + dy)/scaleY;


  # Check if any of the points overlap with a label (of a different point)

  pointMaxx = matrix(xx + offs, nPts, nPts);
  pointMinx = matrix(xx - offs, nPts, nPts);
  pointMiny = matrix(yy - offs, nPts, nPts);
  pointMaxy = matrix(yy + offs, nPts, nPts);

  labelMinx = matrix(xt - labWidth/2, nPts, nPts, byrow = TRUE);
  labelMaxx = matrix(xt + labWidth/2, nPts, nPts, byrow = TRUE);
  labelMiny = matrix(yt - labHeight/2, nPts, nPts, byrow = TRUE);
  labelMaxy = matrix(yt + labHeight/2, nPts, nPts, byrow = TRUE);

  overlapF = function(x1min, x1max, x2min, x2max)
  {
     overlap = matrix(0, nPts, nPts);
     overlap[ x1max > x2min & x1max < x2max & x1min < x2min ] = 1;
     overlap[ x1max > x2min & x1max < x2max & x1min > x2min ] = 2;
     overlap[ x1max > x2max & x1min > x2min ] = 3;
     overlap;
  }

  overlapX = overlapF(pointMinx, pointMaxx, labelMinx, labelMaxx); 
  overlapY = overlapF(pointMiny, pointMaxy, labelMiny, labelMaxy); 

  indOvr = overlapX > 0 & overlapY >0;
  overlap = matrix(0, nPts, nPts);
  overlap[indOvr] = (overlapY[indOvr] - 1) * 3 + overlapX[indOvr];

  # For now try to fix cases of a single overlap.

  nOvrPerLabel = apply(overlap>0, 1, sum);

  #for (p in 1:nPts) if (nOverPerLabel[p]==1) 
  #{
     
  # Check if any of the labels extend past the left or right edge of the plot
  if (protectEdges)
  {
     shift = ifelse(xt - labWidth/2/scaleX < box[1], box[1] - xt + labWidth/2/scaleX, 
                     ifelse(xt + labWidth/2/scaleX > box[2], box[2] - xt - labWidth/2/scaleX, 0));
     xt = xt + shift;

     # Also check the top and bottom edges
     # Do labels overlap with points along the x coordinate?
     xOverlap = abs(xt-x) < (labWidth/2 + offs)/scaleX;

     yShift = ifelse(yt - labHeight/2/scaleY < box[3],  
                      ifelse(xOverlap, (labHeight + 2*offs)/scaleY, box[3] - yt + labHeight/2/scaleY),
                      ifelse(yt + labHeight/2/scaleY > box[4], -(labHeight + 2*offs)/scaleY, 0));
     yt = yt + yShift
  } 

  if (par("xlog")) xt = 10^xt;
  if (par("ylog")) yt = 10^yt;

  if (doPlot)
    text(xt, yt, labels, cex = cex, xpd = xpd, adj = c(0.5, 0.5), ...) 

  invisible(data.frame(x = xt, y= yt, label = labels));
}


