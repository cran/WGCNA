// The name of this file ends in .h so R CMD install does not attempt to compile it on its own.

class CLASS_NAME
{
  protected:

    TYPE 	* data_;
    int 	size_;
    int 	allocated;
    vector <int> dims;
    string	name_;

  public:

  #ifdef CheckDimensions
    TYPE value(int i)
       { if (i<dims[1]) return data_[i]; 
         else throw (Exception(string("Index out of range in variable" + name_).c_str())); 
       }
    TYPE value(int i, int j) 
       { if (dims.size()==2)
           if ( (i<dims[0]) && (j<dims[1]) )
              return data_[j*dims[0] + i]; 
           else
               throw (Exception(string("Index out of range in variable" + name_)));
         else
            throw (Exception(string("incorrect number of dimensions accessing variable" + name_)));
       }
    TYPE value(int i, int j, int k) 
       { if (dims.size()==3)
           if ((k < dims[2]) && (j<dims[1]) && (i<dims[0])) 
              return data_[(k*dims[1]+j)*dims[0] + i]; 
           else
               throw (Exception(string("Index out of range in variable" + name_)));
         else
            throw (Exception(string("Incorrect number of dimensions accessing variable" + name_)));
       }

    TYPE linValue(int i)
       {
         int max = 1;
         for (unsigned di=0; di < dims.size(); di++) max *= dims[di];
         if (i<max)
            return data_[i];
         else 
            throw (Exception(string("Linear index out of range in variable" + name_)));
       }
         
    void setValue(int i, TYPE r)
       { if (i<dims[0]) data_[i] = r;
         else throw (Exception(string("Index out of range in variable" + name_)));
       }
    void setValue(int i, int j, TYPE r)
       { if (dims.size()==2)
           if ( (i<dims[0]) && (j<dims[1]) )
              data_[j*dims[0] + i] = r;
           else
               throw (Exception(string("Index out of range in variable" + name_)));
         else
            throw (Exception(string("incorrect number of dimensions accessing variable" + name_)));
       }
    void setValue(int i, int j, int k, TYPE r)
       { if (dims.size()==3)
           if ((k < dims[2]) && (j<dims[1]) && (i<dims[0])) 
              data_[(k*dims[1]+j)*dims[0] + i] = r; 
           else
               throw (Exception(string("Index out of range in variable" + name_)));
         else
            throw (Exception(string("Incorrect number of dimensions accessing variable" + name_)));
       }

    void linValue(int i, TYPE r)
       {
         int max = 1;
         for (unsigned di=0; di < dims.size(); di++) max *= dims[di];
         if (i<max)
            data_[i] = r;
         else
            throw (Exception(string("Linear index out of range in variable" + name_)));
       }

    // void copyData(CLASS_NAME arr, int start = 0, int length = -1)

  #else

    TYPE linValue(int i)
       { return data_[i]; }
    TYPE value(int i)
       { return data_[i]; }
    TYPE value(int i, int j) 
       { return (data_[j*dims[0] + i]); }
    TYPE value(int i, int j, int k) 
       { return (data_[(k*dims[1]+j)*dims[0] + i]); }

    void linValue(int i, TYPE r)
       { data_[i] = r; }
    void value(int i, TYPE r)
       { data_[i]; }
    void value(int i, int j, TYPE r) 
       { data_[j*dims[0] + i] = r; }
    void value(int i, int j, int k, TYPE r) 
       { data_[(k*dims[1]+j)*dims[0] + i] = r; }
    

  #endif

    void name(string n) {name_ = n; }
    string name() {return name_; }

    void initData(int size);
    void initData(int size, TYPE val);

    void setDim(int length);
    void setDim(int nrow, int ncol);
    void setDim(int nrow, int ncol, int k);

    void setDim(vector <int> dims, int start=0);

    int nDim() { return dims.size(); }

    vector <int> dim() { return dims; }

    int size() { return size_; }
    int length() 
       { 
          if (dims.size()==0) return 0; 
          int prod = 1; for (unsigned i=0; i<dims.size(); i++) prod *= dims[i];
          return prod;
       }

    void wrap(TYPE * data, int len)	// Points the data_ pointer to given data. Will make sure the
                                        // data will not be deallocated in the destructor.
       {
          if (allocated) delete data_;
          allocated = 0;
          data_ = data;
          size_ = len;
          setDim(len);
       }

    void wrap(TYPE * data, int nrow, int ncol)	
       {
         wrap(data, nrow*ncol);
         setDim(nrow, ncol);
       }

    void wrap(TYPE * data, int nrow, int ncol, int k)	
       {
         wrap(data, nrow*ncol*k);
         setDim(nrow, ncol, k);
       }


    TYPE * data() { return data_; }

    TYPE max();
    TYPE min();
    vector <int> table();	// returns frequencies but no values
    vector <int> table(vector <TYPE> & values); // returns frequencies and values

    void copy2vector(int start, int length, vector <TYPE> & result);
    void colMWM(CLASS_NAME & minVal, INT_CLASS & which);
    void colQuantile(double q, CLASS_NAME & quantile);

    void sample(int size, CLASS_NAME & values, int replace = 0);

    // void sort();
    // vector <int> order();
    // vector <int> rank();

    CLASS_NAME() { allocated = 0; data_ = (TYPE *) NULL; dims.clear(); }

    CLASS_NAME(int size) { initData(size); setDim(size); }

    CLASS_NAME(int size, TYPE value) { initData(size, value); setDim(size); }

    // CLASS_NAME(CLASS_NAME arr);	// This constructor will copy the data from arr into *this

    ~CLASS_NAME() { if (allocated) { delete data_; allocated = 0; } }

};

void CLASS_NAME::initData(int size)
{
  size_ = size;
  data_ = new TYPE[size];
  allocated = 1;
  dims.clear();
  dims.push_back(size_);
}

void CLASS_NAME::initData(int size, TYPE val)
{
  initData(size);
  for (int i=0; i<size; i++) data_[i] = val;
}

void CLASS_NAME::setDim(int length)
{
  if (length > size_)
    throw (Exception("attempt to set linear dimension " + NumberToString(length)  + " higher than size "
                     + NumberToString(size()) + " in variable " + name()));
  else
  {
    dims.clear();
    dims.push_back(length);
  }
}

void CLASS_NAME::setDim(int nrow, int ncol)
{
  if (nrow*ncol > size())
    throw (Exception("attempt to set matrix dimensions " + NumberToString(nrow)  + ", " +
                     NumberToString(ncol) + " higher than size " +
                     NumberToString(size()) + " in variable " + name()));
  else
  {
    dims.clear();
    dims.push_back(nrow);
    dims.push_back(ncol);
  }
}

void CLASS_NAME::setDim(int nrow, int ncol, int k)
{
  if (nrow*ncol*k > size_)
    throw (Exception("attempt to set 3-dim CLASS_NAME dimensions " + NumberToString(nrow)  + ", " +
                     NumberToString(ncol) + ", " + NumberToString(k) + " higher than size " +
                     NumberToString(size()) + " in variable " + name()));
  else
  {
    dims.clear();
    dims.push_back(nrow);
    dims.push_back(ncol);
    dims.push_back(k);
  }
}

/*

void CLASS_NAME::copyData(CLASS_NAME arr, int start, int length)
{
  if (start >= arr.length())
    throw(Exception("attempt to copy non-existent data from variable" + arr.name()));
  if (length==-1) length = arr.length() - start;
  if (length > size())
    throw(Exception("attempt to copy data larger than target CLASS_NAME size."));
}

*/

TYPE CLASS_NAME::max()
{
  if (length()==0)
    throw(Exception(string("attempt to calculate max of an empty array.")));
  TYPE max = linValue(0);
  for (int i=1; i<length(); i++)
    if (!ISNAN(linValue(i)) && (linValue(i) > max)) max = linValue(i);
  return max;
}

TYPE CLASS_NAME::min()
{
  if (length()==0)
    throw(Exception(string("attempt to calculate min of an empty array.")));
  TYPE min = linValue(0);
  for (int i=1; i<length(); i++)
    if (!ISNAN(linValue(i)) && (linValue(i) < min)) min = linValue(i);
  return min;
}

vector <int> CLASS_NAME::table(vector <TYPE> & values)
{
  vector <int> counts;
  counts.clear();
  values.clear();
  for (int i=0; i<length(); i++)
  {
    TYPE v = linValue(i);
    unsigned j;
    for (j=0; (j<values.size()) && (values[j]!=v); j++)
    if (j==values.size())
    {
      values.push_back(v);
      counts.push_back(0);
    } else
      counts[j]++;
  }
  return counts;
}

vector <int> CLASS_NAME::table()
{
  vector <TYPE> values;    
  return table(values);
}

void CLASS_NAME::setDim(vector <int> dims, int start)
{
  int len = 1;
  for (unsigned i=start; i<dims.size(); i++) len *= dims[i];
  if (len > size())
    throw(Exception(string("setDim: not enough space to accomodate given dimensions.")));
  this->dims.clear();
  this->dims.reserve(dims.size()-start);
  for (unsigned i=start; i<dims.size(); i++) this->dims.push_back(dims[i]);
}
  

void CLASS_NAME::copy2vector(int start, int length, vector <TYPE> & result)
{
  if (start + length > this->length())
    throw(Exception(string("copy2vector: start+length exceed the actual length of the array.")));
  result.clear();
  // result.reserve(length);
  for (int i=start; i<start + length; i++)
    result.push_back(*(data_+i));
}

void CLASS_NAME::colMWM(CLASS_NAME & minVal, INT_CLASS & which)
{
  if (dim().size()==0)
    throw(Exception(string(
       "Attempt to calculate columnwise minimum of array that has no dimensions set.")));
  if (dim().size()==1)
  {
    minVal.setDim(1);
    which.setDim(1);
  } else {
    minVal.setDim(dim(), 1);
    which.setDim(dim(), 1);
  }

  int colLen = dim()[0], totLen = length();

  if (colLen==0)
    throw(Exception(string("colMWM: Column length is zero in variable") + name()));
  int col = 0; 
  for (int i=0; i<totLen; i+=colLen, col++)
  {
    TYPE cmin = linValue(i);
    int wmin = 0;
    for (int j=i+1; j<i+colLen; j++)
      if (linValue(j) < cmin)
      {
         cmin = linValue(j);
         wmin = j-i;
      }
    minVal.linValue(col, cmin);
    which.linValue(col, wmin);
  }
}

void CLASS_NAME::colQuantile(double q, CLASS_NAME & quantile)
{
  if (dim().size()==0)
    throw(Exception(string(
       "Attempt to calculate columnwise quantile of array that has no dimensions set.")));
  if (dim().size()==1)
    quantile.setDim(1);
  else
    quantile.setDim(dim(), 1);

  int colLen = dim()[0], totLen = length();

  if (colLen==0)
    throw(Exception(string("colQuantile: Column length is zero in variable") + name()));

  vector <TYPE> column;
  column.reserve(colLen);

  double index = q * (colLen-1);
  double i1 = floor(index), i2 = ceil(index);
  int ii1 = (int) i1, ii2 = (int) i2;
  // cout << "Quantile: q = " << q << ", ii1: " << ii1 << ", ii2: " << ii2 << ", index: "<< index << endl;
  for (int i=0, col=0; i<totLen; i+=colLen, col++)
  {
    // column.clear();
    // cout << "Sorting column starting at position" << i << endl;
    copy2vector(i, colLen, column);
    // cout << "original vector:" << endl;
    // for (int kk=0; kk<column.size(); kk++) cout << column[kk] << ", "; cout << endl;
    sort(column.begin(), column.end());
    // cout << "sorted vector:" << endl;
    // for (int kk=0; kk<column.size(); kk++) cout << column[kk] << ", "; cout << endl;
    double val;
    if (i1!=i2)
      val = (index - i1) * column[ii2] + (i2-index) * column[ii1];
    else
      val = column[ii1];
    // cout << "index: " << index << ", i1:" << i1 << ", i2: " << i2 << 
            // ", Quantile: " <<  val << endl;
    quantile.linValue(col, (TYPE) val);
  }
}

/*
   This function is designed to generate relatively small samples from large arrays.
   replace is honored exactly. Using this function on samples that are relatively large will be quite
   slow.
*/
void CLASS_NAME::sample(int size, CLASS_NAME & values, int replace)
{
  int len = length();
  if (replace)
  {
    if (size > length())
      throw(Exception(string("Attempt to sample too many samples without replacement.")));
    values.setDim(size);
    for (int i=0; i<size; i++)
    {
      int s = (int) (floor(unif_rand() * len));
      values.linValue(i, linValue(s));
    }
  } else {
    indArray taken(length(), false);
    values.setDim(size);
    for (int nSampled=0; nSampled < size; )
    {
      int s = (int) (floor(unif_rand() * len));
      if (!taken.value(s))
      {
        values.linValue(nSampled, linValue(s));
        taken.value(s, true);
        nSampled++;
      }
    }
  }
}

