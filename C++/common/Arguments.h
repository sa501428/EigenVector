#ifndef POSSUM_ARGUMENTS_H
#define POSSUM_ARGUMENTS_H
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace possum {
inline int positiveInt(const char *text) {
    char *end = NULL;
    errno = 0;
    long value = std::strtol(text,&end,10);
    if (errno || end == text || *end || value <= 0 || value > INT_MAX)
        throw std::runtime_error(std::string("expected a positive integer, got: ")+text);
    return static_cast<int>(value);
}
}
#endif
