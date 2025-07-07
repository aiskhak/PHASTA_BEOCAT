//
// contains functions for dealing with reading and writing of
// Fortran unformatted file
//

class FortranIO {
public:
  FortranIO(const char *);
  ~FortranIO();

  void setCReadEndianFlag(int flag);

  // read into the array ar(n1,n2) skipping nth records
  // ar is read with n1 varying fastest
  FortranIO& read(int    *ar, int n1, int n2=1, int nth=0);
  FortranIO& read(double *ar, int n1, int n2=1, int nth=0);

  void readCInt(int    *ar, int n1, int n2, int nth);
  void readCDouble(double *ar, int n1, int n2, int nth);
  void readCData(void *ar, int n1, int n2, int nth, size_t type_size);

  // num  = global # of shape functions (nshg or numnp)
  // step = time step in restart file
  FortranIO& readRestartHeader(int &num);
  FortranIO& readRestartHeader(int &num, int &step);

  // write into the file
  // ar is written with n1 varying fastest
  FortranIO& write(int    *ar, int n1, int n2=1);
  FortranIO& write(double *ar, int n1, int n2=1);

  // num  = global # of shape functions (nshg or numnp)
  // step = time step in restart file
  FortranIO& writeRestartHeader(int num, int step=1);

  FortranIO& open();
  void close();

private:
  int length;
  char *fname;
  int unit;
  int CReadFlag;
  int EndianFlag;
  FILE *inFile;
  static int global_unit;
};
