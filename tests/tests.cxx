#include <iostream>
#include <cassert>

// Include project headers
#include "../constants.hpp"
#include "../parameters.hpp"
#include "../array2d.hpp"
#include "../barycentric-fn.hpp"
#include "../sortindex.hpp"
#include "../utils.hpp"

// Define a simple test framework
#define TEST(name) void test_##name()
#define RUN_TEST(name) std::cout << "Running " << #name << "... "; test_##name(); std::cout << "PASSED\n"
#define ASSERT(condition) if (!(condition)) { std::cerr << "FAILED: " << __FILE__ << ":" << __LINE__ << "\n"; exit(1); }

// Test array2d functionality
TEST(array2d) {
    Array2D<int,3> a(10);
    a[0][0] = 1;
    a[0][1] = 2;
    a[0][2] = 3;
    
    ASSERT(a[0][0] == 1);
    ASSERT(a[0][1] == 2);
    ASSERT(a[0][2] == 3);
    
    // Test copy constructor
    Array2D<int,3> b(a);
    ASSERT(b[0][0] == 1);
    ASSERT(b[0][1] == 2);
    ASSERT(b[0][2] == 3);
    
    // Test assignment
    Array2D<int,3> c(5);
    c = a;
    ASSERT(c[0][0] == 1);
    ASSERT(c[0][1] == 2);
    ASSERT(c[0][2] == 3);
}

// Test sorting utility
TEST(sortindex) {
    double a[] = {3.0, 1.0, 4.0, 2.0};
    int indices[4];
    
    sortindex(a, a+4, indices);
    
    ASSERT(indices[0] == 1); // 1.0
    ASSERT(indices[1] == 3); // 2.0
    ASSERT(indices[2] == 0); // 3.0
    ASSERT(indices[3] == 2); // 4.0
}

// Add more tests here...

int main() {
    std::cout << "Running unit tests...\n";
    
    RUN_TEST(array2d);
    RUN_TEST(sortindex);
    
    std::cout << "All tests passed!\n";
    return 0;
}