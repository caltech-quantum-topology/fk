#include <string>

int string_to_int (std::string s) {
    int acc = 0;
    bool negate = false;
    int start = 0;
    if (s[0] == '-') {
        negate = true;
        start = 1;
    }
    for (int i = start; i < s.size(); i++) {
        acc += (s[i] - '0') * std::pow(10, s.size() - i - 1);
    }
    if (negate) {
        return -1 * acc;
    }
    return acc;
}

#include <iostream>
double string_to_double (std::string s) {
    double acc = 0;
    bool negate = false;
    int start = 0;
    if (s[0] == '-') {
        negate = true;
        start = 1;
    }
    int decindex = s.find('.');
    if (decindex == -1) {
        decindex = s.size();
    }
    for (int i = start; i < decindex; i++) {
        acc += (s[i] - '0') * std::pow(10, decindex - i - 1);
    }
    for (int i = decindex + 1; i < s.size(); i++) {
        acc += (s[i] - '0') * std::pow(10, -(i - decindex));
    }
    if (negate) {
        return -1 * acc;
    }
    return acc;
}