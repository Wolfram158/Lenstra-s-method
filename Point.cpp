#include <Point.hpp>

bool Point::get_is_o() {
    return is_o;
}

mpz_class Point::get_x() {
    return x;
}

mpz_class Point::get_y() {
    return y;
}

Point::Point(mpz_class x, mpz_class y, bool is_o) {
    this->x = x;
    this->y = y;
    this->is_o = is_o;
}