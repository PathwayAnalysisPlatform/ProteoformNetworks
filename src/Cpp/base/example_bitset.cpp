#include "example_bitset.hpp"

int run_example() {
    base::dynamic_bitset<> b(50);

    std::cout << b.size() << " " << b.count() << "\n";
    for (int i = 0; i < 10; ++i) {
        b[std::rand() % 50] = true;
    }

    std::cout << b.size() << " " << b.count() << "\n";
    for (int i = 0; i < 50; ++i) {
        std::cout << b[i] << " ";
    }

    base::dynamic_bitset<> x1;
    base::dynamic_bitset<> x2;
    base::dynamic_bitset<> y1 = x1 & x2;
    base::dynamic_bitset<> y2 = x1 | x2;

    return 0;
}
