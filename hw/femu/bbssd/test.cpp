// #include<stdio.h>
#include <iostream>
#include <unordered_map>
#include "lutl.hpp"

// using namespace std;


int test_cpp(ADD *obj) {
    ADD qwe;
    obj->add.push_back(1);
    obj->add.push_back(2);
    cout << "[test cpp]:hello world\n"  << endl;
    for (int i = 0; i < 2; i++) cout << obj->add[i] << " ";
    cout << endl;
    return 1;
}