#ifdef IBM
void flush(int* junk ){ return; }
#endif

#ifdef LINUX
void flush_(int* junk ){ return; }
#endif

#ifdef intel
void FLUSH(int* junk ){ return; }
#endif
