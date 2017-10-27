#include "fuchsia.cpp"
#include <assert.h>
#include <chrono>

int
main(int argc, char *argv[])
{
    symtab s;
    matrix m;
    tie(m, s) = load_matrix((argc > 1) ? argv[1] : "git_409.m", s);
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    symbol x = ex_to<symbol>(s["x"]);
    pfmatrix pfm(m, x);
    matrix pfm_mm = pfm.to_matrix(x);
    for (unsigned i = 0; i < m.rows(); i++) {
        for (unsigned j = 0; j < m.cols(); j++) {
            assert(normal(m(i,j) - pfm_mm(i,j)).is_zero());
        }
    }
    cout << "mx size = " << m.rows() << " x " << m.cols() << endl;
    for (const auto kv : pfm.residues) {
        const ex &x0 = kv.first.first;
        const int &k = kv.first.second;
        const matrix &c = kv.second;
        cout << "residue at x=" << x0 << ", k=" << k << endl;
        for (auto ev : eigenvalues(c)) {
            cout << "  ev: " << ev << endl;
        }
    }
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    std::cout << "time taken: " << time_span.count() << " seconds" << endl;
    return 0;
}
