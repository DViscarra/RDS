#include <helpers.h>
#include <classes.h>
#include <sstream>
#include <iomanip>
#include <algorithm>

std::string double_to_string(double num, int precision) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << num;
    return stream.str();
}

bool is_subset(std::vector<Chemical> subset, std::vector<Chemical> superset) {
    for (Chemical chem : subset) if (std::find(superset.begin(), superset.end(), chem) == superset.end()) return false;

    return true;
}