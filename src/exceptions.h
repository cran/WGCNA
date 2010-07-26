/*

   Exception handling. Just adding a bit more information to the standard exception.

*/

#ifndef __exception_h__

#define __exception_h__

#include <iostream>

using namespace std;

class Exception
{
  protected:

    string	_what;

  public:

    virtual string what() const throw() { return _what; }

    Exception(string wht) throw()
    {
      _what = wht;
    }

    ~Exception() throw() {}
};

#endif
