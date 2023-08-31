#include <classes.h>
#include <iostream>
#include <unordered_map>

int main() {
    Chemical c1,c2;

    ChemicalMap chems_in = {{c1, 2}};
    ChemicalMap chems_out = {{c2, 3}};

    Reaction r(chems_in, chems_out);

    std::cout << r.id << std::endl;

    return 0;
}