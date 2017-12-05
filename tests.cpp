#include <catch/catch.hpp>
#include "fuchsia.cpp"

TEST_CASE("eigenvectors") {
    symbol x("x");
    matrix m = {
        {x, 0, 0},
        {1, x+1, 0},
        {x, 0, x+1},
    };
    for (const ex evalue : vector<ex>({x, x+1})) {
        SECTION("eigenvectors right") {
            for (const matrix &ev : eigenvectors_right(m, evalue)) {
                ex mm = normal(m.mul(ev) - evalue*ev);
                REQUIRE(mm.is_zero_matrix());
                mm = normal(m.mul(ev) - 2*evalue*ev);
                REQUIRE(!mm.is_zero_matrix());
            }
        }
        SECTION("eigenvectors left") {
            for (const matrix &ev : eigenvectors_left(m, evalue)) {
                ex mm = normal(ev.mul(m) - evalue*ev);
                REQUIRE(mm.is_zero_matrix());
                mm = normal(ev.mul(m) - 2*evalue*ev);
                REQUIRE(!mm.is_zero_matrix());
            }
        }
    }
}

TEST_CASE("matrix save/load") {
    symbol x("x"), y("y");
    matrix m1 = {
        {x, 0, 2*y/3},
        {1, x+y, 0},
        {x, 0, (x-1)/(y*y*y+1)},
    };
    save_matrix("/tmp/fuchsia.test.m", m1);
    symtab s = {{"x", x}, {"y", y}};
    matrix m2;
    tie(m2, s) = load_matrix("/tmp/fuchsia.test.m", s);
    REQUIRE(m1 == m2);
}

TEST_CASE("poly_decompose") {
    symbol x("x"), y("y");
    ex e = power(1+2*x+3*y, 3);
    ex p = y*x + y + 1;
    lst c = poly_decompose(e, p, x);
    REQUIRE(c.nops() > 1);
    ex ee = 0;
    for (size_t i = 0; i < c.nops(); i++) {
        ee += c.op(i) * power(p, i);
    }
    REQUIRE(expand(e) == expand(ee));
}

TEST_CASE("partial_fraction, 1 variable") {
    auto check_pf = [](const ex &e, const symbol &x) {
        ex p = partial_fraction(e, x);
        REQUIRE(ratcan(e) == ratcan(p));
    };
    symbol x("x");
    SECTION("0") {
        check_pf(0, x);
    }
    SECTION("9") {
        check_pf(9, x);
    }
    SECTION("2*x+1") {
        check_pf(2*x+1, x);
    }
    SECTION("x*x+2*x+3") {
        check_pf(x*x+2*x+3, x);
    }
    SECTION("1/x") {
        check_pf(1/x, x);
    }
    SECTION("1/x/x") {
        check_pf(1/x/x, x);
    }
    SECTION("1/(x*x-1)") {
        check_pf(1/(x*x-1), x);
    }
    SECTION("x/(x*x+2)/(x-1)") {
        check_pf(x/(x*x+2)/(x-1), x);
    }
    SECTION("(x*x*x+2*x+3)/(3*x*x+x+1)") {
        check_pf((x*x*x+2*x+3)/(3*x*x+x+1), x);
    }
    SECTION("x*x+2*(x-1)+(2*x-1)/(2*x*x-1)^3") {
        check_pf(x*x+2*(x-1)+(2*x-1)/pow(2*x*x-1, 3), x);
    }

}

TEST_CASE("partial_fraction, 2 variables") {
    auto check_pf = [](const ex &e, const symbol &x) {
        ex p = partial_fraction(e, x);
        REQUIRE(ratcan(e) == ratcan(p));
    };
    symbol x("x"), y("y");
    SECTION("1/(x-1)/(y-1)") {
        check_pf(1/(x-1)/(y-1), x);
    }
    SECTION("1/(x*x-y)") {
        check_pf(1/(x*x-y), x);
    }
    SECTION("1/(x*x-y*y)") {
        check_pf(1/(x*x-y*y), x);
    }
    SECTION("(x+1)/(y*x*x-1)/(x-y)") {
        check_pf((x+1)/(y*x*x-1)/(x-y), x);
    }
    SECTION("(x*x*x+x*x+x+1)/(y*x-1)") {
        check_pf((x*x*x+x*x+x+1)/(y*x-1), x);
    }
    SECTION("(x*x*x+x*x+x+1)/(y*x*x-1)") {
        check_pf((x*x*x+x*x+x+1)/(y*x*x-1), x);
    }
}

TEST_CASE("pfmatrix") {
    symbol x("x"), y("y");
    matrix m = {
        {(y-1)/(x-1)/(y*x-2), (x*x-y*y)/(y-1)},
        {(x*x+2)/(y*x+1/y)+pow(y-1, 2)*x, (x-y)*(y-1/x)}
    };
    pfmatrix pfm(m, x);
    matrix mm = pfm.to_matrix();
    REQUIRE(ratcan(m) == ratcan(mm));
}

TEST_CASE("eigenvalues & eigenvectors") {
    symbol a("a"), b("b"), c("c"), x("x"), y("y"), z("z");
    matrix m = {
        {1, 0, 8, 0, 0, c},
        {0, 2, a, 0, 0, 0},
        {1, 0, 3, 0, b, 0},
        {0, 0, y, x, 0, 0},
        {0, 0, 0, 0, 5, 0},
        {0, 0, 0, 0, 0, z},
    };
    for (auto &eval : eigenvalues(m)) {
        for (auto &evec : eigenvectors_right(m, eval.first)) {
            matrix z = ex_to<matrix>((m*evec - eval.first*evec).evalm());
            for (unsigned i = 0; i < z.rows(); i++) {
                for (unsigned j = 0; j < z.cols(); j++) {
                    REQUIRE(normal(z(i, j)) == 0);
                }
            }
        }
    }
}

TEST_CASE("pfmatrix add_*") {
    symbol x("x"), y("y");
    matrix m = {
        {y+1, y/x, x+1/pow(x-1, 2), 0},
        {0, x/(x-1), 1/x-1, (y-1)/x/x},
        {2/x/y, x/(x-1), x*(x-1)-y/x, y/x/(x-1)},
        {x*x-1/x, 2/(x-1), x+y/x, 1-2/(x-1)}
    };
    matrix c = {
        {2-y, 2, 3*y, 1},
        {y*y, 0, -y, 3},
        {y/2, 2, y+1, 0},
        {0, 2*y, 1/y, 1}
    };
    auto test_add_mul = [&](auto c, auto pi, auto ki, auto x0) {
        matrix mm1 = ex_to<matrix>((m + c*pow(x - pi, ki)*(x - x0)).evalm());
        pfmatrix pfm(m, x);
        pfm.add_mul(c, pi, ki, x0);
        matrix mm2 = pfm.to_matrix();
        REQUIRE(ratcan(mm1) == ratcan(mm2));
    };
    SECTION("add_mul, 1") { test_add_mul(c, 1, -2, -y*y); }
    SECTION("add_mul, 2") { test_add_mul(c, 1, -1, y/2); }
    SECTION("add_mul, 3") { test_add_mul(c, 0, -1, 1); }
    SECTION("add_mul, 4") { test_add_mul(c, 0, 0, 3*y); }
    SECTION("add_mul, 5") { test_add_mul(c, 0, 1, 0); }
    SECTION("add_mul, 6") { test_add_mul(c, 0, 2, 2*y-1); }
    auto test_add_div = [&](auto c, auto pi, auto ki, auto x0) {
        matrix mm1 = ex_to<matrix>((m + c*pow(x - pi, ki)/(x - x0)).evalm());
        pfmatrix pfm(m, x);
        pfm.add_div(c, pi, ki, x0);
        matrix mm2 = pfm.to_matrix();
        REQUIRE(ratcan(mm1) == ratcan(mm2));
    };
    SECTION("add_div, 1") { test_add_div(c, 1, -2, -y*y); }
    SECTION("add_div, 2") { test_add_div(c, 1, -1, y/2); }
    SECTION("add_div, 3") { test_add_div(c, 0, -1, 1); }
    SECTION("add_div, 4") { test_add_div(c, 0, 0, 3*y); }
    SECTION("add_div, 5") { test_add_div(c, 0, 1, 0); }
    SECTION("add_div, 6") { test_add_div(c, 0, 2, 2*y-1); }
    auto test_add_pow = [&](auto c, auto p1, auto k1, auto p2, auto k2) {
        matrix mm1 = ex_to<matrix>((m + c*pow(x - p1, k1)*pow(x - p2, k2)).evalm());
        pfmatrix pfm(m, x);
        pfm.add_pow(c, p1, k1, p2, k2);
        matrix mm2 = pfm.to_matrix();
        REQUIRE(ratcan(mm1) == ratcan(mm2));
    };
    SECTION("add_pow, 1") { test_add_pow(c, 1, -2, 0, 2); }
    SECTION("add_pow, 2") { test_add_pow(c, 1, -1, 1, -1); }
    SECTION("add_pow, 3") { test_add_pow(c, 0, -1, 0, 0); }
    SECTION("add_pow, 4") { test_add_pow(c, 0, 0, 1, -1); }
    SECTION("add_pow, 5") { test_add_pow(c, 0, 1, 0, 2); }
    SECTION("add_pow, 6") { test_add_pow(c, 0, 2, 0, -1); }
    SECTION("add_pow, 7") { test_add_pow(c, 0, 2, y, -2); }
    SECTION("add_pow, 8") { test_add_pow(c, 1, -2, 1/y, -3); };
    SECTION("add_pow, 9") { test_add_pow(c, y, -1, 0, -2); };
}

TEST_CASE("pfmatrix with_*") {
    symbol x("x"), y("y");
    matrix m = {
        {x+2+3/(x-y)+4/x/x, 5+6/x, 7*x*x, 8/x/x},
        {9*x+10/x/(x-y), 11*x*x+12/x, 13+14/(x-y)+15/x/x, 0},
        {16/x+17/x/x, 18*x*x+19, 20*x+21/x, 22+24*x+25/x},
        {26, 27*x*x + 28/x, 29/x/(x-y), 30/x+31*x},
    };
    pfmatrix pfm(m, x);
    REQUIRE(ratcan(m) == ratcan(pfm.to_matrix()));
    SECTION("with_constant_t") {
        matrix t = {
            {1, 2, 3, 4},
            {5, 4, 3, 2},
            {0, 2, 4, 8},
            {9, 5, 2, 1}
        };
        matrix mt = t.inverse().mul(m).mul(t);
        matrix m2 = pfm.with_constant_t(t).to_matrix();
        matrix m3 = pfm.with_constant_t(t.inverse(), t).to_matrix();
        REQUIRE(ratcan(mt) == ratcan(m2));
        REQUIRE(ratcan(mt) == ratcan(m3));
    }
    SECTION("with_balance_t") {
        matrix u = {
            {0, ex(3)/5, ex(4)/5, 0},
            {ex(5)/13, 0, 0, ex(12)/13},
        };
        matrix p = u.transpose().mul(u);
        REQUIRE(p.mul(p) == p);
        SECTION("0 --> y") {
            matrix mt = with_balance_t(m, p, 0, y, x);
            matrix m2 = pfm.with_balance_t(p, 0, y).to_matrix();
            REQUIRE(ratcan(mt) == ratcan(m2));
        }
        SECTION("y --> infinity") {
            matrix mt = with_balance_t(m, p, y, infinity, x);
            matrix m2 = pfm.with_balance_t(p, y, infinity).to_matrix();
            REQUIRE(ratcan(mt) == ratcan(m2));
        }
        SECTION("infinity --> 0") {
            matrix mt = with_balance_t(m, p, infinity, 0, x);
            matrix m2 = pfm.with_balance_t(p, infinity, 0).to_matrix();
            REQUIRE(ratcan(mt) == ratcan(m2));
        }
    }
}

TEST_CASE("jordan") {
    matrix j0 = {
        {2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3},
    };
    matrix p = {
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 6, 7, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 8, 9, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 3, 4, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 5, 6, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 7, 8, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 9, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
    };
    matrix m = p.mul(j0).mul(p.inverse());
    matrix q = jordan(m).first;
    REQUIRE(q.rank() == j0.rows());
    matrix j = q.inverse().mul(m).mul(q);
    for (unsigned r = 0; r < j.rows(); r++) {
        for (unsigned c = 0; c < j.cols(); c++) {
            const ex &x = j(r, c);
            if (r == c) {
                REQUIRE((x - 2)*(x - 3) == 0);
            }
            else if (r + 1 == c) {
                REQUIRE((x - 0)*(x - 1) == 0);
            }
            else {
                REQUIRE(x == 0);
            }
        }
    }
}

TEST_CASE("jordan 2") {
    symbol x("x");
    matrix m = {
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,-18*x*x*x+15*x*x-3*x,0,0,0,0},
        {0,-12*x*x*x+10*x*x-2*x,0,0,0,0},
    };
    auto qcs = jordan(m);
    auto q = qcs.first;
    auto cs = qcs.second;
    for (auto n : cs) {
        REQUIRE(n > 0);
    }
    REQUIRE(q.rank() == m.rows());
}
