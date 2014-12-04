// The name of this file ends in .h so R CMD install does not attempt to compile it on its own.

class CLASS_NAME
{
  protected:

    TYPE 	* data_;
    size_t 	size_;
    int 	allocated;
    vector <size_t> dims;
    string	name_;

  public:

  #ifdef CheckDimensions
    TYPE value(size_t i)
       { if (i<dims[1]) return data_[i]; 
         else throw (Exception(string("Index out of range in variable" + name_).c_str())); 
       }
    TYPE value(size_t i, size_t j) 
       { if (dims.size()==2)
           if ( (i<dims[0]) && (j<dims[1]) )
              return data_[j*dims[0] + i]; 
           else
               throw (Exception(string("Index out of range in variable" + name_)));
         else
            throw (Exception(string("incorrect number of dimensions accessing variable" + name_)));
       }
    TYPE value(size_t i, size_t j, size_t k) 
       { if (dims.size()==3)
           if ((k < dims[2]) && (j<dims[1]) && (i<dims[0])) 
              return data_[(k*dims[1]+j)*dims[0] + i]; 
           else
               throw (Exception(string("Index out of range in variable" + name_)));
         else
            throw (Exception(string("Incorrect number of dimensions accessing variable" + name_)));
       }

    TYPE linValue(size_t i)
       {
         size_t max = 1;
         for (size_t di=0; di < dims.size(); di++) max *= dims[di];
         if (i<max)
            return data_[i];
         else 
            throw (Exception(string("Linear index out of range in variable" + name_)));
       }
         
    void setValue(size_t i, TYPE r)
       { if (i<dims[0]) data_[i] = r;
         else throw (Exception(string("Index out of range in variable" + name_)));
       }
    void setValue(size_t i, size_t j, TYPE r)
       { if (dims.size()==2)
           if ( (i<dims[0]) && (j<dims[1]) )
              data_[j*dims[0] + i] = r;
           else
               throw (Exception(string("Index out of range in variable" + name_)));
         else
            throw (Exception(string("incorrect number of dimensions accessing variable" + name_)));
       }
    void setValue(size_t i, size_t j, size_t k, TYPE r)
       { if (dims.size()==3)
           if ((k < dims[2]) && (j<dims[1]) && (i<dims[0])) 
              data_[(k*dims[1]+j)*dims[0] + i] = r; 
           else
               throw (Exception(string("Index out of range in variable" + name_)));
         else
            throw (Exception(string("Incorrect number of dimensions accessing variable" + name_)));
       }

    void linValue(size_t i, TYPE r)
       {
         size_t max = 1;
         for (size_t di=0; di < dims.size(); di++) max *= dims[di];
         if (i<max)
            data_[i] = r;
         else
            throw (Exception(string("Linear index out of range in variable" + name_)));
       }

    // void copyData(CLASS_NAME arr, size_t start = 0, size_t length = -1)

  #else

    TYPE linValue(size_t i)
       { return data_[i]; }
    TYPE value(size_t i)
       { return data_[i]; }
    TYPE value(size_t i, size_t j) 
       { return (data_[j*dims[0] + i]); }
    TYPE value(size_t i, size_t j, size_t k) 
       { return (data_[(k*dims[1]+j)*dims[0] + i]); }

    void linValue(size_t i, TYPE r)
       { data_[i] = r; }
    void value(size_t i, TYPE r)
       { data_[i]; }
    void value(size_t i, size_t j, TYPE r) 
       { data_[j*dims[0] + i] = r; }
    void value(size_t i, size_t j, size_t k, TYPE r) 
       { data_[(k*dims[1]+j)*dims[0] + i] = r; }
    

  #endif

    void name(string n) {name_ = n; }
    string name() {return name_; }

    void initData(size_t size);
    void initData(size_t size, TYPE val);

    void setDim(size_t length);
    void setDim(size_t nrow, size_t ncol);
    void setDim(size_t nrow, size_t ncol, size_t k);

    void setDim(vector <size_t> dims, size_t start=0);

    size_t nDim() { return dims.size(); }

    vector <size_t> dim() { return dims; }

    size_t size() { return size_; }
    size_t length() 
       { 
          if (dims.size()==0) return 0; 
          size_t prod = 1; for (size_t i=0; i<dims.size(); i++) prod *= dims[i];
          return prod;
       }

    void wrap(TYPE * data, size_t len)	// Points the data_ pointer to given data. Will make sure the
                                        // data will not be deallocated in the destructor.
       {
          if (allocated) delete data_;
          allocated = 0;
          data_ = data;
          size_ = len;
          setDim(len);
       }

    void wrap(TYPE * data, size_t nrow, size_t ncol)	
       {
         wrap(data, nrow*ncol);
         setDim(nrow, ncol);
       }

    void wrap(TYPE * data, size_t nrow, size_t ncol, size_t k)	
       {
         wrap(data, nrow*ncol*k);
         setDim(nrow, ncol, k);
       }


    TYPE * data() { return data_; }

    TYPE max();
    TYPE min();
    vector <size_t> table();	// returns frequencies but no values
    vector <size_t> table(vector <TYPE> & values); // returns frequencies and values

    void copy2vector(size_t start, size_t length, vector <int> & result);
    void copy2vector(size_t start, size_t length, vector <double> & result);
    void colMWM(CLASS_NAME & minVal, INT_CLASS & which);
    void colQuantile(double q, dArray & quantile);
    void rowQuantile(double q, dArray & quantile);

    void sample(size_t size, CLASS_NAME & values, int replace = 0);

    // void sort();
    // vector <size_t> order();
    // vector <size_t> rank();

    CLASS_NAME() { allocated = 0; data_ = (TYPE *) NULL; dims.clear(); }

    CLASS_NAME(size_t size) { initData(size); setDim(size); }

    CLASS_NAME(size_t size, TYPE value) { initData(size, value); setDim(size); }

    // CLASS_NAME(CLASS_NAME arr);	// This constructor will copy the data from arr into *this

    ~CLASS_NAME() { if (allocated) { delete data_; allocated = 0; } }

};

void CLASS_NAME::initData(size_t size)
{
  size_ = size;
  data_ = new TYPE[size];
  allocated = 1;
  dims.clear();
  dims.push_back(size_);
}

void CLASS_NAME::initData(size_t size, TYPE val)
{
  initData(size);
  for (size_t i=0; i<size; i++) data_[i] = val;
}

void CLASS_NAME::setDim(size_t length)
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

void CLASS_NAME::setDim(size_t nrow, size_t ncol)
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

void CLASS_NAME::setDim(size_t nrow, size_t ncol, size_t k)
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

void CLASS_NAME::copyData(CLASS_NAME arr, size_t start, size_t length)
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
  for (size_t i=1; i<length(); i++)
    if (!ISNAN(linValue(i)) && (linValue(i) > max)) max = linValue(i);
  return max;
}

TYPE CLASS_NAME::min()
{
  if (length()==0)
    throw(Exception(string("attempt to calculate min of an empty array.")));
  TYPE min = linValue(0);
  for (size_t i=1; i<length(); i++)
    if (!ISNAN(linValue(i)) && (linValue(i) < min)) min = linValue(i);
  return min;
}

vector <size_t> CLASS_NAME::table(vector <TYPE> & values)
{
  vector <size_t> counts;
  counts.clear();
  values.clear();
  for (size_t i=0; i<length(); i++)
  {
    TYPE v = linValue(i);
    size_t j;
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

vector <size_t> CLASS_NAME::table()
{
  vector <TYPE> values;    
  return table(values);
}

void CLASS_NAME::setDim(vector <size_t> dims, size_t start)
{
  size_t len = 1;
  for (size_t i=start; i<dims.size(); i++) len *= dims[i];
  if (len > size())
    throw(Exception(string("setDim: not enough space to accomodate given dimensions.")));
  this->dims.clear();
  this->dims.reserve(dims.size()-start);
  for (size_t i=start; i<dims.size(); i++) this->dims.push_back(dims[i]);
}
  

void CLASS_NAME::copy2vector(size_t start, size_t length, vector <int> & result)
{
  if (start + length > this->length())
    throw(Exception(string("copy2vector: start+length exceed the actual length of the array.")));
  result.clear();
  // result.reserve(length);
  for (size_t i=start; i<start + length; i++)
    result.push_back((int) *(data_+i));
}

void CLASS_NAME::copy2vector(size_t start, size_t length, vector <double> & result)
{
  if (start + length > this->length())
    throw(Exception(string("copy2vector: start+length exceed the actual length of the array.")));
  result.clear();
  // result.reserve(length);
  for (size_t i=start; i<start + length; i++)
    result.push_back((double) *(data_+i));
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

  size_t colLen = dim()[0], totLen = length();

  if (colLen==0)
    throw(Exception(string("colMWM: Column length is zero in variable") + name()));
  size_t col = 0; 
  for (size_t i=0; i<totLen; i+=colLen, col++)
  {
    TYPE cmin = linValue(i);
    size_t wmin = 0;
    for (size_t j=i+1; j<i+colLen; j++)
      if (linValue(j) < cmin)
      {
         cmin = linValue(j);
         wmin = j-i;
      }
    minVal.linValue(col, cmin);
    which.linValue(col, wmin);
  }
}

/*
   This function is designed to generate relatively small samples from large arrays.
   replace is honored exactly. Using this function on samples that are relatively large will be quite
   slow.
*/
void CLASS_NAME::sample(size_t size, CLASS_NAME & values, int replace)
{
  size_t len = length();
  if (replace)
  {
    if (size > length())
      throw(Exception(string("Attempt to sample too many samples without replacement.")));
    values.setDim(size);
    for (size_t i=0; i<size; i++)
    {
      size_t s = (size_t) (floor(unif_rand() * len));
      values.linValue(i, linValue(s));
    }
  } else {
    indArray taken(length(), false);
    values.setDim(size);
    for (size_t nSampled=0; nSampled < size; )
    {
      size_t s = (size_t) (floor(unif_rand() * len));
      if (!taken.value(s))
      {
        values.linValue(nSampled, linValue(s));
        taken.value(s, true);
        nSampled++;
      }
    }
  }
}

