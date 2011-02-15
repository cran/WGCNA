/*

  Compiling:

  gcc --std=c99 -fPIC -O3 -o functionsThreaded.so -shared -lpthread -lm \
  -I/usr/local/lib/R-2.7.1-goto/include \
  -DWITH_THREADS \
  corFunctions-common.c corFunctions.c corFunctions-parallel.c networkFunctions.c pivot.c

  Without threads: 

  gcc --std=c99 -fPIC -O3 -o functions.so -shared -lpthread -lm \
  -I/usr/local/lib/R-2.7.1-goto/include \
  corFunctions-common.c corFunctions.c corFunctions-parallel.c networkFunctions.c pivot.c



  Home:

  gcc --std=c99 -fPIC -O3 -o functionsThreaded.so -shared -lpthread -lm \
  -I/usr/local/lib/R-2.8.0-patched-2008-12-06-Goto/include \
  -DWITH_THREADS \
  corFunctions-common.c corFunctions.c corFunctions-parallel.c networkFunctions.c pivot.c

  Without threads: 

  gcc --std=c99 -fPIC -O3 -o functions.so -shared -lpthread -lm \
  -I/usr/local/lib/R-2.8.0-patched-2008-12-06-Goto/include \
  corFunctions-common.c corFunctions.c corFunctions-parallel.c networkFunctions.c pivot.c



 
 gcc --std=gnu99 --shared -Ic:/PROGRA~1/R/R-27~0PA/include -o functions.dll functions.c
-Lc:/PROGRA~1/R/R-27~0PA/bin -lR -lRblas

"C:\Program Files\R\R-2.7.0pat\bin\R.exe" CMD SHLIB functions.c -lRblas

*/

