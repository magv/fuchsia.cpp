#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <fstream>
#include <tuple>
#include <assert.h>

using namespace GiNaC;
using namespace std;

/* Load matrix from a file in Mathematica format.
 */
pair<matrix, symtab>
load_matrix(const char *filename, const symtab &table)
{
    parser reader(table);
    ifstream i(filename);
    ex x = reader(i);
    i.close();
    // TODO: check that the input is indeed a matrix of rational
    // expressions, with no imaginary or irrational numbers; signal
    // errors if this is not the case.
    matrix m = ex_to<matrix>(lst_to_matrix(ex_to<lst>(x)));
    return make_pair(m, reader.get_syms());
}

/* Save matrix to a file in Mathematica format.
 */
void
save_matrix(const char *filename, const matrix &m)
{
    ofstream f(filename);
    f << "{";
    for (unsigned i = 0; i < m.rows(); i++) {
        if (i != 0) f << ",\n";
        f << "{";
        for (unsigned j = 0; j < m.cols(); j++) {
            if (j != 0) f << ",";
            f << m(i, j);
        }
        f << "}";
    }
    f << "}";
    f.close();
}

/* Find the set of right eigenvectors for a given eigenvalue.
 */
vector<matrix>
eigenvectors_right(matrix m, ex eigenvalue)
{
    unsigned n = m.cols();
    // Construct eigenvector of symbol temporaries.
    std::vector<symbol> tmp;
    matrix ev(n, 1);
    for (unsigned i = 0; i < n; i++) {
        symbol t;
        tmp.push_back(t);
        ev(i, 0) = t;
    }
    // Construct mm = m - eigenvalue*I
    matrix mm(n, n);
    for (unsigned i = 0; i < n; i++) {
        mm(i, i) = eigenvalue;
    }
    mm = m.sub(mm);
    // Solve mm*v = 0
    matrix rhs(n, 1);
    matrix s = mm.solve(ev, rhs);
    // Find eigenvectors
    vector<matrix> evectors;
    for (unsigned k = 0; k < n; k++) {
        matrix ev(n, 1);
        for (unsigned i = 0; i < n; i++) {
            ev(i, 0) = normal(s(i, 0).coeff(tmp[k]));
        }
        // is_zero_matrix() is only sufficient here because
        // of the normal() call above.
        if (not ev.is_zero_matrix()) {
            evectors.push_back(ev);
        }
    }
    return evectors;
}

/* Find the set of left eigenvectors for a given eigenvalue.
 */
vector<matrix>
eigenvectors_left(matrix m, ex eigenvalue)
{
    auto v = eigenvectors_right(m.transpose(), eigenvalue);
    for (auto &e : v) e = e.transpose();
    return v;
}

/* Divide one polynomial in x by another, return the quotient
 * and the remainder.
 */
pair<ex, ex>
poly_divmod(const ex &a, const ex &b, const symbol &x)
{
    ex aa = a.expand().collect(x);
    ex bb = b.expand().collect(x);
    int deg_b = bb.degree(x);
    if (deg_b == 0)
        return make_pair(a/b, ex(0));
    ex b0 = normal(bb.coeff(x, deg_b));
    ex q = 0;
    for (;;) {
        int deg_a = aa.degree(x);
        if (deg_a < deg_b) {
            return make_pair(q, aa);
        }
        ex a0 = aa.coeff(x, deg_a);
        ex a0_b0 = normal(a0/b0);
        if (a0_b0.is_zero()) {
            // a0 is an un-normal zero; let's drop it from aa and try again.
            aa -= a0*pow(x, deg_a);
            continue;
        }
        ex k = a0_b0*pow(x, deg_a - deg_b);
        aa = (aa - bb*k).expand().collect(x);
        q += k;
    }
}

/* Decompose polynomial over rationals e(x) into {c_i(x)},
 * such that e(x) = \sum_i c_i(x) p(x)^i.
 */
lst
poly_decompose(const ex &e, const ex &p, const symbol &x)
{
    lst c;
    ex q = e;
    while (!q.is_zero()) {
        auto qr = poly_divmod(q, p, x);
        c.append(qr.second);
        // XXX: is normal() call needed here?
        q = qr.first;
    }
    return c;
}

/* Iterate through factors of e, call yield(f, k) for each
 * factor of the form f^k.
 *
 * Note that this function doesn't factor e itself, only iterates
 * through the factors already explicitly present.
 */
template <typename F> void
factor_iter(const ex &e, F yield)
{
    if (is_a<mul>(e)) {
        for (const auto &f : e) {
            if (is_a<power>(f)) {
                yield(f.op(0), ex_to<numeric>(f.op(1)).to_int());
            } else {
                yield(f, 1);
            }
        }
    } else {
        if (is_a<power>(e)) {
            yield(e.op(0), ex_to<numeric>(e.op(1)).to_int());
        } else {
            yield(e, 1);
        }
    }
}

/* Iterate through terms of a partial fraction decomposition
 * of a rational expression in x. Call yield(p, q, k) for
 * each p*q^k term.
 */
template<typename F> void
partial_fraction_iter(const ex &e, const symbol &x, F yield)
{
    ex nd = numer_denom(e);
    ex numer = nd.op(0);
    ex denom = nd.op(1);
    auto qr = poly_divmod(numer, denom, x);
    ex q = qr.first;
    ex r = qr.second;
    denom = factor(denom);
    struct Term { ex f; int n; };
    vector<Term> factors;
    factor_iter(denom,
        [&](const ex &f, int i) {
            factors.push_back(Term{f, i});
        });
    lst clist;
    ex eq = -r;
    // For each factor...
    for (const auto &f : factors) {
        int deg = f.f.degree(x);
        // For each power of factor...
        for (int k = 1; k <= f.n; k++) {
            ex pfnum = 0;
            // For each power of the parfrac term numerator...
            for (int l = 0; l < deg; l++) {
                symbol c;
                clist.append(c);
                pfnum += c*pow(x, l);
            }
            eq += expand(pfnum*(denom/pow(f.f, k)));
        }
    }
    lst eqlist;
    int deg_lo = eq.ldegree(x);
    int deg_hi = eq.degree(x);
    for (int k = deg_lo; k <= deg_hi; k++) {
        eqlist.append(eq.coeff(x, k) == 0);
    }
    ex sol = lsolve(eqlist, clist);
    q = q.expand().collect(x);
    for (int k = q.ldegree(x); k <= q.degree(x); k++) {
        ex c = q.coeff(x, k);
        if (!c.is_zero()) yield(c, x, k);
    }
    // For each factor...
    int idx = 0;
    for (const auto &f : factors) {
        int deg = f.f.degree(x);
        // For each power of factor...
        for (int k = 1; k <= f.n; k++) {
            ex pfnum = 0;
            // For each power of the parfrac term numerator...
            for (int l = 0; l < deg; l++) {
                pfnum += sol.op(idx).rhs()*pow(x, l);
                idx++;
            }
            if (!pfnum.is_zero()) yield(pfnum, f.f, -k);
        }
    }
}

/* Transform a rational expression into partial fraction form in x.
 */
ex
partial_fraction(const ex &e, const symbol &x)
{
    ex pf = 0;
    partial_fraction_iter(e, x, [&](const auto &p, const auto &q, int n) {
        pf += p*pow(q, n);
    });
    return pf;
}

/* Matrix in a partial fraction form:
 *     M(x) = \sum_i R_i (x-p_i)^{k_i},
 *     where p_i = 0, if k_i >= 0.
 *
 * Note that this class expects its p_i arguments to be in some
 * sort of canonical form. Honor this expectation. Or else.
 */
struct pfmatrix {
    typedef pair<ex, int> key;
    struct key_is_less {
        bool operator ()(const key &k1, const key &k2) const;
    };

    map<key, matrix, key_is_less> residues;
    unsigned nrows;
    unsigned ncols;

    pfmatrix(const matrix &m, const symbol &x);
    matrix &operator ()(const ex &p, int k);
    matrix to_matrix(const ex &x);
    // M += C*(x-pi)^ki
    void add(const matrix &C, const ex &pi, int ki);
    // M += C*(x-pi)^ki/(x-x0)
    void add_div(const matrix &C, const ex &pi, int ki, const ex &x0);
    // M += C*(x-pi)^ki*(x-x0)
    void add_mul(const matrix &C, const ex &pi, int ki, const ex &x0);
    // M += C*(x-p1)^k1*(x-p2)^k2
    void add_pow(const matrix &C, const ex &p1, int k1, const ex &p2, int k2);
};

bool
pfmatrix::key_is_less::operator ()(const key &k1, const key &k2) const
{
    return (k1.second < k2.second) || ex_is_less()(k1.first, k2.first);
}

pfmatrix::pfmatrix(const matrix &m, const symbol &x)
    : nrows(m.rows()), ncols(m.cols())
{
    for (unsigned i = 0; i < nrows; i++) {
        for (unsigned j = 0; j < ncols; j++) {
            partial_fraction_iter(m(i, j), x,
                [&](const auto &p, const auto &q, int k) {
                    int deg = q.degree(x);
                    if (k >= 0) {
                        // q == x
                        (*this)(0, k)(i, j) = p;
                    }
                    else if (deg == 1) {
                        // q == c0 + c1*x
                        ex c0 = q.coeff(x, 0);
                        ex c1 = q.coeff(x, 1);
                        (*this)(normal(-c0/c1), k)(i, j) = p*pow(c1, k);
                    }
                    else {
                        throw runtime_error("pfmatrix(): can't solve equations of 2nd degree or higher");
                    }
                });
        }
    }
}

matrix &
pfmatrix::operator ()(const ex &p, int k)
{
    auto key = make_pair(p, k);
    decltype(residues)::iterator it;
    tie(it, ignore) = residues.emplace(
            piecewise_construct,
            make_tuple(key),
            std::make_tuple(nrows, ncols));
    return (*it).second;
}

matrix
pfmatrix::to_matrix(const ex &x)
{
    ex m = matrix(nrows, ncols);
    for (const auto &k : residues) {
        m += k.second*pow(x - k.first.first, k.first.second);
    }
    return ex_to<matrix>(m.evalm());
}

void
pfmatrix::add(const matrix &C, const ex &pi, int ki)
{
    matrix &R = (*this)(pi, ki);
    R = R.add(C);
}

void
pfmatrix::add_div(const matrix &C, const ex &pi, int ki, const ex &x0)
{
    if (pi == x0) {
        add(C, pi, ki - 1);
    }
    else if (ki >= 0) {
        assert(pi == 0);
        for (int k = 0; k < ki; k++) {
            add(C.mul_scalar(pow(x0, ki - 1 - k)), pi, k);
        }
        add(C.mul_scalar(pow(x0, ki)), x0, -1);
    }
    else {
        ex dx = x0 - pi;
        for (int k = ki; k < 0; k++) {
            add(C.mul_scalar(-pow(dx, ki - 1 - k)), pi, k);
        }
        add(C.mul_scalar(pow(dx, ki)), x0, -1);
    }
}

void
pfmatrix::add_mul(const matrix &C, const ex &pi, int ki, const ex &x0)
{
    add(C.mul_scalar(pi - x0), pi, ki);
    if (ki == -1) {
        add(C, 0, ki + 1);
    }
    else {
        add(C, pi, ki + 1);
    }
}

void
pfmatrix::add_pow(const matrix &C, const ex &p1, int k1, const ex &p2, int k2)
{
    if (k1 >= 0) assert(p1 == 0);
    if (k2 >= 0) assert(p2 == 0);
    if (p1 == p2) {
        add(C, p1, k1 + k2);
    }
    else if (k1 == 0) {
        add(C, p2, k2);
    }
    else if (k2 == 0) {
        add(C, p1, k1);
    }
    else if ((k1 >= 0) && (k2 >= 0)) {
        add(C, p1, k1 + k2);
    }
    else if ((k1 < 0) && (k2 < 0)) {
        ex d = p1 - p2;
        add(C.mul_scalar(pow(d, k2)), p1, k1);
        for (int i = k2; i < 0; i++) {
            add_pow(C.mul_scalar(-pow(d, i)),
                    (k1 + 1 == 0) ? 0 : p1, k1 + 1,
                    p2, k2 - i - 1);
        }
    }
    else if ((k1 < 0) && (k2 > 0)) {
        long b = 1;
        for (int i = 0; i < min(-k1, k2 + 1); i++) {
            add(C.mul_scalar(b*pow(p1, k2 - i)), p1, k1 + i);
            b = b*(k2 - i)/(i + 1);
        }
        for (int i = min(-k1, k2 + 1); i < k2 + 1; i++) {
            long b2 = b;
            for (int j = 0; j < k1 + i + 1; j++) {
                add(C.mul_scalar(((k1 + i - j) % 2 == 0 ? b2 : -b2)*pow(p1, k1 + k2 - j)), 0, j);
                b2 = b2*(k1 + i - j)/(j + 1);
            }
            b = b*(k2 - i)/(i + 1);
        }
    }
    else if ((k1 > 0) && (k2 < 0)) {
        add_pow(C, p2, k2, p1, k1);
    }
    else {
        assert(false);
    }
}

/* Block-triangular permutation of a matrix.
 * This class computes it, and keeps the results.
 *
 * This one is a class instead of a function, because
 * recursive closures are awkward (impossible?) in C++.
 */
class block_triangular_permutation {
    public:
    block_triangular_permutation(const matrix &m);
    const matrix & t();
    const vector<int> &block_size();

    private:
    struct node {
        int index;
        int link;
        bool onstack;
        node() : index(-1) {}
    };
    const matrix &m;
    vector<node> nodes;
    stack<int> st;
    matrix result_t;
    int index_counter;
    int output_counter;
    vector<int> result_block_size;

    void visit(int i);
};

block_triangular_permutation::
block_triangular_permutation(const matrix &m)
    : m(m),
        nodes(m.rows()),
        result_t(m.rows(), m.cols()),
        index_counter(0),
        output_counter(0)
{
    for (int i = 0; i < (int)m.rows(); i++) {
        if (nodes[i].index == -1)
            visit(i);
    }
}

const matrix &
block_triangular_permutation::t()
{ return result_t; }

const vector<int> &
block_triangular_permutation::block_size()
{ return result_block_size; }

void
block_triangular_permutation::visit(int i)
{
    nodes[i].index = nodes[i].link = index_counter++;
    st.push(i);
    nodes[i].onstack = true;
    for (int j = 0; j < (int)m.rows(); j++) {
        if (i == j) continue;
        if (m(i, j).is_zero()) continue;
        if (nodes[j].index == -1) {
            visit(j);
            nodes[i].link = min(nodes[i].link, nodes[j].link);
        }
        else if (nodes[j].onstack) {
            nodes[i].link = min(nodes[i].link, nodes[j].index);
        }
    }
    if (nodes[i].link == nodes[i].index) {
        int size = 0;
        for (;;) {
            int j = st.top();
            st.pop();
            nodes[j].onstack = false;
            result_t(j, output_counter++) = 1;
            size++;
            if (j == i) break;
        }
        result_block_size.push_back(size);
    }
}

/* Compute and return characteristic polynomial of a matrix by
 * first permuting it into block-triangular shape.
 *
 * This is significantly faster than the naive approach for
 * larger (e.g. 10x10), sparse matrices.
 */
ex
charpoly_by_blocks(const matrix &m, const ex &lambda)
{
    block_triangular_permutation btp(m);
    matrix mm = ex_to<matrix>((btp.t().transpose()*m*btp.t()).evalm());
    ex cp = 1;
    int o = 0;
    for (int size : btp.block_size()) {
        matrix b = ex_to<matrix>(sub_matrix(mm, o, size, o, size));
        cp *= b.charpoly(lambda);
        o += size;
    }
    return cp;
}

template<typename F> void
factor_twice_and_iter(const ex &e, F yield)
{
    factor_iter(factor(e),
        [&](const ex &f1, int k1) {
            factor_iter(factor(f1),
                [&](const ex &f2, int k2) {
                    yield(f2, k1*k2);
                });
        });
}

/* Find and return eigenvalues of a matrix.
 * Throw an error, if unable.
 */
vector<ex>
eigenvalues(const matrix &m)
{
    vector<ex> eigenvalues;
    symbol lambda("λ");
    ex charpoly = charpoly_by_blocks(m, lambda);
    factor_iter(factor(charpoly),
        [&](const ex &f, int k) {
            int deg = f.degree(lambda);
            if (deg == 0) {
                // No roots in λ here.
            }
            else if (deg == 1) {
                // f == c0 + c1*λ
                ex c0 = f.coeff(lambda, 0);
                ex c1 = f.coeff(lambda, 1);
                for (int i = 0; i < k; i++)
                    eigenvalues.push_back(normal(-c0/c1));
            }
            else {
                throw runtime_error("eigenvalues(): can't solve equations of 2nd degree or higher");
            }
        });
    return eigenvalues;
}
