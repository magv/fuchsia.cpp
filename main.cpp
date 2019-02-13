#include "fuchsia.cpp"

static const char usagetext[] = R"(
Ss{NAME}
    Nm{fuchsia} -- transform linear differential equations into an epsilon
    form.

Ss{SYNOPSYS}
    Nm{fuchsia} [options] Cm{command} Ar{args} ...

Ss{DESCRIPTION}
    Nm{fuchsia} transforms systems of linear differential equations,
        Dl{∂/∂x I(x,ε) = M(x,ε) I(x,ε),}

    into an epsilon form,
        Dl{∂/∂x J(x,ε) = ε S(x) J(x,ε),}

    where Ql{I} and Ql{J} are column vectors of functions in the original
    and epsilon basis, Ql{M} is the original matrix, Ql{ε*S} is the
    matrix in an epsilon form, and Ql{I} is related to Ql{J} via the
    transformation matrix Ql{T(x,ε)} such that
        Dl{I = T J.}

    In all cases Ql{M} can depend on additional symbolic variables, which
    are treated as independent constants.

Ss{EXAMPLES}
    To reduce a single-variable differential system of equations to an
    epsilon form, use this:

    Nm{fuchsia} Cm{reduce} Fl{-x} Ar{x} Fl{-e} Ar{eps} Ar{matrix.orig} Fl{-m} Ar{matrix.ep} Fl{-t} Ar{matrix.ep.t} \
        Fl{-C} 2>&1 | tee Ar{matrix.ep.log}

    For differential equations in multiple variables this is the usage:

    Nm{fuchsia} Cm{reduce} Ar{matrix.x} Ar{matrix.y} Fl{-x} Ar{x} Fl{-x} Ar{y} Fl{-e} Ar{eps} \
        Fl{-m} Ar{matrix.ep.x} Fl{-m} Ar{matrix.ep.y} Fl{-t} Ar{matrix.ep.t} \
        Fl{-C} 2>&1 | tee Ar{matrix.ep.log}

Ss{COMMANDS}
    Cm{show} [Fl{-x} Ar{name}] Ar{matrix}
        Show a description of a given matrix.

    Cm{reduce} [Fl{-x} Ar{name}] [Fl{-e} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix}
        Find an epsilon form of the given matrix. Internally
        this is a combination of Cm{reduce-diagonal-blocks},
        Cm{fuchsify-off-diagonal-blocks} and Cm{factorize}.

    Cm{reduce} [Fl{-x} Ar{name}] ... [Fl{-e} Ar{name}] [Fl{-m} Ar{path}] ... [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix} ...
        Find an epsilon form of a given multivariate differential equation
        system. A matching number of Ar{matrix} arguments, Fl{-x}, and
        Fl{-m} flags is required.

        The matrices are reduced one by one, and a single transformation
        is computed that simultaneously transforms all of them into an
        epsilon form. It may be best to list the simplest matrix first.

        NOTE: this command is under development, and may fail when it
        shouldn't.

    Cm{reduce-diagonal-blocks} [Fl{-x} Ar{name}] [Fl{-e} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix}
        Transform the matrix into a block-triangular form and reduce the
        diagonal blocks into an epsilon form.

    Cm{fuchsify-off-diagonal-blocks} [Fl{-x} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix}
        Transform the off-diagonal blocks of a block-triangular matrix
        into a Fuchsian form, assuming the diagonal blocks are already in
        an epsilon form, thus making the whole matrix normalized Fuchsian.

    Cm{factorize} [Fl{-x} Ar{name}] [Fl{-e} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix}
        Find a transformation that will make a given normalized Fuchsian
        matrix proportional to the infinitesimal parameter.

    Cm{fuchsify} [Fl{-x} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix}
        Find a transformation that will transform a given matrix into a
        Fuchsian form. This is less efficient than block-based commands,
        because it effectively treats the whole matrix as one big block.

    Cm{normalize} [Fl{-x} Ar{name}] [Fl{-e} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix}
        Find a transformation that will transform a given Fuchsian matrix
        into a normalized form. This is less efficient than block-based
        commands, because it effectively treats the whole matrix as one
        big block.

    Cm{sort} [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix}
        Find a block-triangular form of the given matrix by shuffling.

    Cm{transform} [Fl{-x} Ar{name}] [Fl{-m} Ar{path}] Ar{matrix} Ar{transform} ...
        Transform a given matrix using a given transformation.

    Cm{changevar} [Fl{-x} Ar{name}] [Fl{-y} Ar{name}] [Fl{-m} Ar{path}] Ar{matrix} Ar{expr}
        Perform a change of variable from x to y, such that x=expr(y).

    Cm{suggest-changevar} [Fl{-x} Ar{name}] [Fl{-y} Ar{name}] Ar{matrix}
        Suggest a rational change of variable that will transform residue
        eigenvalues of the form n/2+k*eps into n+k*eps, thus making it
        possible to find an epsilon form of the matrix.

    Cm{simplify} [Fl{-x} Ar{name}] [Fl{-m} Ar{path}] [Fl{-t} Ar{path}] [Fl{-i} Ar{path}] Ar{matrix}
        Try to find a transformation that makes a given matrix
        simpler, for some definition of "simple".

Ss{OPTIONS}
    Fl{-x} Ar{name}    Use this name for the free variable (default: x).
    Fl{-0} Ar{expr}    Set this value for x during multivariate reduction (default: 0).
    Fl{-y} Ar{name}    Use this name for the new free variable (default: y).
    Fl{-e} Ar{name}    Use this name for the infinitesimal parameter (default: eps).
    Fl{-m} Ar{path}    Save the resulting matrix into this file.
    Fl{-t} Ar{path}    Save the resulting transformation into this file.
    Fl{-i} Ar{path}    Save the inverse transformation into this file.
    Fl{-C}         Force colored output even if stdout is not a tty.
    Fl{-P}         Paranoid mode: spend more time checking internal invariants.
    Fl{-q}         Print a more quiet log.
    Fl{-h}         Show this help message.
    Fl{-V}         Print version information.

Ss{ARGUMENTS}
    Ar{matrix}     Read the input matrix from this file.
    Ar{transform}  Read the transformation matrix from this file.
    Ar{expr}       Arbitrary expression.

Ss{AUTHORS}
    Vitaly Magerya <vitaly.magerya@tx97.net>
)";

/* Return a mapping of the form (A*x+B)/(C*x+D) that maps
 * x=0 into a, x=1 into b, x=infinity into c.
 */
ex
invmoebius(const ex &x, const ex &a, const ex &b, const ex &c)
{
    if (a == infinity) return c + (b - c)/x;
    if (b == infinity) return (c*x - a)/(x - 1);
    if (c == infinity) return a + (b - a)*x;
    return ((c - b)*a + (b - a)*c*x)/((c - b) + (b - a)*x);
}

void
print_matrix_shape(ostream &f, const matrix &m, const char *ident)
{
    for (unsigned i = 0; i < m.rows(); i++) {
        f << ident;
        for (unsigned j = 0; j < m.cols(); j++) {
            if (COLORS && (i == j)) f << "\033[1;34m";
            f << (m(i, j).is_zero() ? '.' : '#');
            if (COLORS && (i == j)) f << "\033[0m";
        }
        f << endl;
    }
}

void
usage()
{
    const char *p = strchr(usagetext, '\n') + 1;
    for (;;) {
        const char *l = strchr(p + 2, '{');
        if (l == NULL) break;
        const char *r = strchr(l, '}');
        if (r == NULL) break;
        const char *a = "", *b = "\033[0m";
        if (l[-2] == 'S' && l[-1] == 's') { a = "\033[1m"; }
        if (l[-2] == 'N' && l[-1] == 'm') { a = "\033[1;35m"; }
        if (l[-2] == 'F' && l[-1] == 'l') { a = "\033[33m"; }
        if (l[-2] == 'C' && l[-1] == 'm') { a = "\033[1m"; }
        if (l[-2] == 'A' && l[-1] == 'r') { a = "\033[32m"; }
        if (l[-2] == 'D' && l[-1] == 'l') { a = "\033[36m"; }
        if (l[-2] == 'Q' && l[-1] == 'l') { a = "\033[36m"; }
        cout.write(p, l - p - 2);
        if (COLORS) cout << a;
        cout.write(l + 1, r - l - 1);
        if (COLORS) cout << b;
        p = r + 1;
    }
    cout << p;
}

#define IFCMD(name, condition) \
    if ((argc >= 1) && !strcmp(argv[0], name)) \
        if (!(condition)) { \
            cerr << "fuchsia: malformed '" << argv[0] \
                 << "' invocation (use -h to see usage)" << endl; \
            return 1; \
        } else

int
main(int argc, char *argv[])
{
    vector<const char *> var_x_names;
    vector<const char *> var_x_values;
    const char *var_y_name = "y";
    const char *var_eps_name = "eps";
    vector <const char *>matrix_m_paths;
    const char *matrix_t_path = NULL;
    const char *matrix_i_path = NULL;
    for (int opt; (opt = getopt(argc, argv, "hq0:x:e:y:m:t:i:s:CPV")) != -1;) {
        switch (opt) {
        case 'h': usage(); return 0;
        case 'V': cout << VERSION; return 0;
        case 'q': VERBOSE = false; break;
        case 'x': var_x_names.push_back(optarg); break;
        case '0': var_x_values.push_back(optarg); break;
        case 'y': var_y_name = optarg; break;
        case 'e': var_eps_name = optarg; break;
        case 'm': matrix_m_paths.push_back(optarg); break;
        case 't': matrix_t_path = optarg; break;
        case 'i': matrix_i_path = optarg; break;
        case 'P': PARANOID = true; break;
        case 'C': COLORS = true; break;
        default: return 1;
        }
    }
    argc -= optind;
    argv += optind;
    vector<matrix> matrices_m;
    matrix matrix_m(0, 0);
    matrix matrix_t(0, 0);
    matrix matrix_i(0, 0);
    symbol y(var_y_name);
    symbol eps(var_eps_name);
    parser reader(symtab({{var_y_name, y}, {var_eps_name, eps}}));
    vector<symbol> vars_x;
    if (var_x_names.size() == 0) {
        var_x_names.push_back("x");
    };
    for (auto &&name : var_x_names) {
        symbol x(name);
        reader.get_syms()[name] = x;
        vars_x.push_back(x);
    }
    const symbol &x = vars_x[vars_x.size() - 1];
    IFCMD("help", argc == 1) {
        usage();
        return 0;
    }
    else IFCMD("show", argc == 2) {
        matrix m = load_matrix(argv[1], reader);
        cout << "Matrix size: " << m.rows() << "x" << m.cols() << endl;
        pfmatrix pfm(m, x);
        cout << "Matrix complexity: " << complexity(pfm) << endl;
        cout << "Matrix expansion:" << endl;
        for (const auto kv : pfm.residues) {
            const auto &xi = kv.first.first;
            const auto &ki = kv.first.second;
            const auto &c = kv.second;
            if (c.is_zero_matrix()) continue;
            cout << "  pole of power " << ki << " at " << x << "=" << (ki >= 0 ? infinity : xi) << endl;
            cout << "    complexity: " << complexity(c) << endl;
            for (auto ev : eigenvalues(c, true)) {
                cout << "    e-value^" << ev.second << ": " << ev.first << endl;
            }
            cout << "    shape: " << endl;
            print_matrix_shape(cout, c, "      ");
        }
        matrix c = normal(c0_infinity(pfm));
        if (!c.is_zero_matrix()) {
            cout << "  effective pole of power " << -1 << " at " << x << "=" << infinity << endl;
            cout << "    complexity: " << complexity(c) << endl;
            for (auto ev : eigenvalues(c, true)) {
                cout << "    e-value^" << ev.second << ": " << ev.first << endl;
            }
            cout << "    shape: " << endl;
            print_matrix_shape(cout, c, "      ");
        }
    }
    else IFCMD("sort", argc == 2) {
        matrix m = load_matrix(argv[1], reader);
        block_triangular_permutation btp(m);
        logi("Block sizes: {}", btp.block_size());
        matrix_m = btp.t().transpose().mul(m).mul(btp.t());
        matrix_t = btp.t();
        matrix_i = btp.t().transpose();
    }
    else IFCMD("fuchsify", argc == 2) {
        pfmatrix pfm = load_pfmatrix(argv[1], x, reader);
        auto r = fuchsify(pfm);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else IFCMD("normalize", argc == 2) {
        pfmatrix pfm = load_pfmatrix(argv[1], x, reader);
        auto r = normalize(pfm, eps);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else IFCMD("factorize", argc == 2) {
        pfmatrix pfm = load_pfmatrix(argv[1], x, reader);
        auto r = factorize(pfm, eps);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else IFCMD("reduce", argc >= 2) {
        if ((size_t)(argc - 1) != var_x_names.size()) {
            cerr << "fuchsia: got " << argc - 1
                 << " matrices and " << var_x_names.size()
                 << " variable names -- these should match" << endl;
            return 1;
        }
        vector<matrix> mlist;
        for (int i = 1; i < argc; i++) {
            matrix m = load_matrix(argv[i], reader);
            mlist.push_back(m);
        }
        vector<ex> vals_x;
        for (auto &&value : var_x_values) {
            vals_x.push_back(reader(value));
        }
        auto r = reduce_multivar(mlist, vars_x, vals_x, eps);
        for (auto &&pfm : r.first) {
            matrices_m.push_back(pfm.to_matrix());
        }
        matrix_t = r.second.first;
        matrix_i = r.second.second;
    }
    else IFCMD("reduce-diagonal-blocks", argc == 2) {
        pfmatrix pfm = load_pfmatrix(argv[1], x, reader);
        auto r = reduce_diagonal_blocks(pfm, eps);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else IFCMD("fuchsify-off-diagonal-blocks", argc == 2) {
        pfmatrix pfm = load_pfmatrix(argv[1], x, reader);
        auto r = fuchsify_off_diagonal_blocks(pfm);
        matrix_m = r.first.to_matrix();
        matrix_t = r.second.to_matrix();
        matrix_i = r.second.to_inverse_matrix();
    }
    else IFCMD("transform", argc >= 3) {
        matrix m = load_matrix(argv[1], reader);
        for (int i = 2; i < argc; i++) {
            matrix t = load_matrix(argv[i], reader);
            m = t.inverse().mul(m.mul(t).sub(ex_to_matrix(t.diff(x))));
        }
        pfmatrix pfm(m, x);
        matrix_m = pfm.to_matrix();
    }
    else IFCMD("changevar", argc == 3) {
        matrix m = load_matrix(argv[1], reader);
        auto xsubs = reader(argv[2]);
        matrix m2 = ex_to_matrix(m.subs(exmap{{x, xsubs}})).mul_scalar(xsubs.diff(y));
        matrix_m = pfmatrix(m2, y).to_matrix();
    }
    else IFCMD("suggest-changevar", argc == 2) {
        pfmatrix pfm = load_pfmatrix(argv[1], x, reader);
        exset halfpoints;
        for (auto &&kv : pfm.residues) {
            const auto &pi = kv.first.first;
            const auto &ci = kv.second;
            for (const auto &ev : eigenvalues(ci)) {
                const auto &eval = ev.first;
                ex ev0 = eval.subs(exmap{{eps, 0}}, subs_options::no_pattern);
                assert(is_a<numeric>(ev0));
                numeric ev0n = ex_to<numeric>(ev0);
                numeric den = ev0n.denom();
                if (den == 1) { continue; }
                else if (den == 2) {
                    logi("Found a half-integer point at {}={} with eigenvalue {}", x, pi, eval);
                    halfpoints.insert(pi);
                }
                else assert(false);
            }
        }
        matrix c0inf = c0_infinity(pfm);
        int infinity_too = 0;
        for (const auto &ev : eigenvalues(c0inf)) {
            const auto &eval = ev.first;
            ex ev0 = eval.subs(exmap{{eps, 0}}, subs_options::no_pattern);
            assert(is_a<numeric>(ev0));
            numeric ev0n = ex_to<numeric>(ev0);
            numeric den = ev0n.denom();
            if (den == 1) { }
            else if (den == 2) {
                logi("Found a half-integer point at {}={} with eigenvalue {}", x, infinity, eval);
                infinity_too = 1;
            }
            else assert(false);
        }
        logi("The matrix has {} half-integer points in total", halfpoints.size() + infinity_too);
        if (halfpoints.size() + infinity_too == 0) {
            logi("No variable change is necessary");
        }
        else if (halfpoints.size() + infinity_too == 1) {
            assert(!"Is only one half-integer point even possible?");
        }
        else if (halfpoints.size() == 1 && infinity_too) {
            auto &&it = halfpoints.begin();
            ex a = *it++;
            logi("This variable change from {} to {} should help:\n{} = {}",
                    x, y,
                    x, a + y*y); // invmoebius(y*y, a, a + 1, infinity)
        }
        else if (halfpoints.size() == 2 && !infinity_too) {
            auto &&it = halfpoints.begin();
            ex a = *it++, b = *it++;
            logi("Any of these variable changes from {} to {} should help:\n{} = {}\n{} = {}",
                    x, y,
                    x, (a + b*y*y)/(y*y + 1), // invmoebius(y*y, a, (a + b)/2, b)
                    x, (b + a*y*y)/(y*y + 1)); // invmoebius(y*y, b, (a + b)/2, a)
        }
        else if (halfpoints.size() == 2 && infinity_too) {
            auto &&it = halfpoints.begin();
            ex a = *it++, b = *it++;
            logi("Any of these variable changes from {} to {} should help:\n{} = {}\n{} = {}",
                    x, y,
                    x, invmoebius(pow((y*y + 1)/(2*y), 2), a, b, infinity),
                    x, invmoebius(pow((y*y + 1)/(2*y), 2), b, a, infinity));
        }
        else if (halfpoints.size() == 3 && !infinity_too) {
            auto &&it = halfpoints.begin();
            ex a = *it++, b = *it++, c = *it++;
            logi("Any of these variable changes from {} to {} should help:\n"
                    "{} = {}\n{} = {}\n{} = {}\n{} = {}\n{} = {}\n{} = {}",
                    x, y,
                    x, invmoebius(pow((y*y+1)/(2*y), 2), a, b, c),
                    x, invmoebius(pow((y*y+1)/(2*y), 2), c, a, b),
                    x, invmoebius(pow((y*y+1)/(2*y), 2), b, c, a),
                    x, invmoebius(pow((y*y+1)/(2*y), 2), a, c, b),
                    x, invmoebius(pow((y*y+1)/(2*y), 2), b, a, c),
                    x, invmoebius(pow((y*y+1)/(2*y), 2), c, b, a));
        }
        else {
            loge("The matrix has {} half-integer points in total, "
                    "no rational variable change can cure that",
                    halfpoints.size() + infinity_too);
        }
    }
    else IFCMD("simplify", argc == 2) {
        pfmatrix pfm = load_pfmatrix(argv[1], x, reader);
        auto r1 = simplify_off_diagonal_blocks(pfm);
        auto r2 = simplify_by_rescaling(r1.first);
        r1.second.add(r2.second);
        matrix_m = r2.first.to_matrix();
        matrix_t = r1.second.to_matrix();
        matrix_i = r1.second.to_inverse_matrix();
    }
    else if (argc == 0) {
        cerr << "fuchsia: no command provided (use -h to see usage)" << endl;
        return 1;
    }
    else {
        cerr << "fuchsia: unrecognized command '"
             << argv[0]
             << "' (use -h to see usage)" << endl;
        return 1;
    }
    if (matrix_m.nops() != 0) matrices_m.push_back(matrix_m);
    // Saving the resulting diff. eq. matrices.
    size_t nsave = min(matrix_m_paths.size(), matrices_m.size());
    for (size_t i = 0; i < nsave; i++) {
        logi("Saving a diff. eq. matrix to {}", matrix_m_paths[i]);
        save_matrix(matrix_m_paths[i], matrices_m[i]);
    }
    for (size_t i = nsave; i < matrices_m.size(); i++) {
        logi("Printing a diff. eq. matrix to the standard output (no -m argument)");
        save_matrix(cout, matrices_m[i]);
        cout << endl;
    }
    for (size_t i = nsave; i < matrix_m_paths.size(); i++) {
        loge("No matrix to save into {}", matrix_m_paths[i]);
    }
    // Saving the transformation matrix.
    if (matrix_t.nops() != 0 && matrix_t_path != NULL) {
        logi("Saving the (unsimplified) transformation to {}", matrix_t_path);
        save_matrix(matrix_t_path, matrix_t);
    }
    if (matrix_t.nops() != 0 && matrix_t_path == NULL) {
        logw("Not saving the transformation matrix (no -t argument)");
    }
    if (matrix_t.nops() == 0 && matrix_t_path != NULL) {
        loge("No transformation matrix to save into {}", matrix_t_path);
    }
    // Saving the inverse transformation matrix.
    if (matrix_i.nops() != 0 && matrix_i_path != NULL) {
        logi("Saving the (unsimplified) inverse transformation to {}", matrix_i_path);
        save_matrix(matrix_i_path, matrix_i);
    }
    if (matrix_i.nops() != 0 && matrix_i_path == NULL) {
        logd("Not saving the inverse transformation matrix (no -i argument)");
    }
    if (matrix_i.nops() == 0 && matrix_i_path != NULL) {
        loge("No inverse transformation matrix to save into {}", matrix_i_path);
    }
    // Re-saving the simplified transformation matrices.
    //
    // This final simplification can take a while, and the idea
    // is that the user can Ctrl-C out of Fuchsia without waiting
    // for it, while still having the unsimplified versions
    // saved.
    if (matrix_t.nops() != 0 && matrix_t_path != NULL) {
        logd("Simplifying the transformation...");
        matrix t = normal(matrix_t);
        logi("Saving the simplified transformation to {}", matrix_t_path);
        save_matrix(matrix_t_path, t);
    }
    if (matrix_i.nops() != 0 && matrix_i_path != NULL) {
        logd("Simplifying the inverse transformation...");
        matrix i = normal(matrix_i);
        logi("Saving the simplified inverse transformation to {}", matrix_i_path);
        save_matrix(matrix_i_path, i);
    }
    return 0;
}
