#include <iostream>
#include <boost/version.hpp>
#include <LEDA/graph/graph.h>

int main(void){
    std::cout << "Using Boost " << BOOST_VERSION / 10000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
    return 0;
}

// : ^)