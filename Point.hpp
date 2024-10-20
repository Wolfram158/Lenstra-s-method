#include <Utils.cpp>

class Point {
public:
    Point(mpz_class x, mpz_class y, bool is_o = false);
    bool get_is_o();
    mpz_class get_x();
    mpz_class get_y();
private:
    mpz_class x;
    mpz_class y;
    bool is_o;
};