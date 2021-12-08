

#ifdef WIN32
//_declspec(dllexport) extern int number_of_thread;
_declspec(dllexport) void setNumber_of_thread(int num);
_declspec(dllexport) int getNumber_of_thread();
#endif

extern int number_of_thread;
void setNumber_of_thread(int num);
int getNumber_of_thread();