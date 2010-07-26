/*
 * Compilation:
 *  gcc --std=c99 -fPIC -O3 -o pivot.so -shared pivot.c
 *
 */


#include <stdlib.h>
#include <R.h>


void RprintV(double * v, int l)
{
  for (int i=0; i<l; i++) Rprintf("%5.3f ", v[i]);
  Rprintf("\n");
}

double vMax(double * v, int len)
{
  double mx = v[0];
  for (int i=1; i<len; i++)
    if (v[i] > mx) mx = v[i];
  return mx;
}

double vMin(double * v, int len)
{
  double mn = v[0];
  for (int i=1; i<len; i++)
    if (v[i] < mn) mn = v[i];
  return mn;
}


double pivot(double * v, int len, double target)
{
  // Rprintf("Entering pivot with len=%d and target=%f\n   ", len, target);
  // RprintV(v, len);

  if (len > 2)
  {
    // pick the pivot, say as the median of the first, middle and last
    int i1 = 0, i2 = len-1, i3 = (len-1)/2, ip;
    if (v[i1] <= v[i2])
    {
      if (v[i2] <= v[i3])
        ip = i2;
      else if (v[i3] >= v[i1])
         ip = i3;
      else 
         ip = i1;
    } else {
      if (v[i1] <= v[i3])
        ip = i1;
      else if (v[i2] <= v[i3])
        ip = i3;
      else
        ip = i2;
    }

    // put ip at the end
    double vp = v[ip];
    v[ip] = v[len-1];
    v[len-1] = vp;

    // Rprintf("   pivot value: %5.3f, index: %d\n", vp, ip);

    // pivot everything else
    int bound = 0;
    for (int i=0; i<len; i++) if (v[i] < vp)
    {
      double x = v[bound];
      v[bound] = v[i];
      v[i] = x;
      bound++;
    }

    v[len-1] = v[bound];
    v[bound] = vp;

    // Rprintf("   After pivoting: bound:%d and vector: ", bound); // RprintV(v, len);

    // Did we find the target?
    
    double crit = target - bound;
    // Rprintf("   crit: %5.3f\n", crit);
    if (fabs(crit) > 1.0)
    {
      if (crit < 0)
        return pivot(v, bound, target);
      else
        return pivot(v+bound+1, len-bound-1, target-bound-1);
    }
    // Rprintf("vMax(v, bound): %5.3f, vMin(v+bound+1, len-bound-1): %5.3f, vp: %5.3f\n", vMax(v, bound),
                // vMin(v+bound+1, len-bound-1), vp);
    if (crit < 0)
    {
       double v1 = vMax(v, bound);
       return (v1 *(-crit) + vp * (1+crit));
    } // else
    double v2 = vMin(v+bound+1, len-bound-1);
    return (vp * (1-crit) + v2 * crit);
  } 
  else if (len==2)
  {
      // Rprintf("  Short v, returning a direct value.\n"); 
      double v1 = vMin(v, 2);
      double v2 = vMax(v, 2);
      if (target < 0) return v1;
      else if (target > 1) return v2;
      else return (target * v2 + (1-target) * v1);
  }
  else 
  {
     // Rprintf("  length 1 v, returning a direct value.\n"); 
     return v[0];
  }
}


/*
 *
 * This isn't needed for now.
 *
 *
void testPivot(double * v, int * len, double * target, double * result)
{
   * result = pivot(v, *len, *target);
}
*/

