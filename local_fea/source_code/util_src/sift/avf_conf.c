#include "avf_conf.h"

//int number_of_thread =16;
int number_of_thread = 16;

void setNumber_of_thread(int num){
	number_of_thread = num;
	//std::cout<<number_of_thread<<std::endl;
}
int getNumber_of_thread(){
	return number_of_thread;
}
