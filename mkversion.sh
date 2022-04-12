gitid() {
    cidlong=$(git --git-dir "$1/.git" describe --always --abbrev=0 --match '' 2>/dev/null)
    if [ $? -ne 0 ]; then return 1; fi
    cidshort=$(git --git-dir "$1/.git" describe --always --abbrev=12 --match '' --dirty=+)
    date=$(git --git-dir "$1/.git" log -1 $cidlong --date=format:'%Y-%m-%d' --format=format:'%cd')
    echo "commit $cidshort from $date"
}

pkgid_rpm() {
    info=$(rpm --queryformat="%{VERSION}-%{RELEASE}\n" -q "$1" 2>/dev/null)
    if [ $? -ne 0 ]; then return 1; fi
    printf "%s" "$info" | head -1
}

pkgid_dpkg() {
    info=$(dpkg -s "$1" 2>/dev/null)
    if [ $? -ne 0 ]; then return 1; fi
    printf "%s" "$info" | sed -n '/Version: /{s/^.* //;p}'
}

pkgid_pkgconf() {
    info=$(pkg-config --modversion "$1" 2>/dev/null)
    if [ $? -ne 0 ]; then return 1; fi
    printf "%s" "$info" | head -1
}

pkgid() {
    local pkg
    for pkg in "$@"; do
        pkgid_rpm "$pkg" && return
        pkgid_dpkg "$pkg" && return
        pkgid_pkgconf "$pkg" && return
    done
    echo "unknown"
}

cat <<EOF
static const char VERSION[] = R"(\
Fuchsia, $(gitid .)
Libraries:
  GiNaC $(pkgid ginac libginac libginac-dev)
  CLN $(pkgid cln libcln libcln-dev)
  GMP $(pkgid gmp libgmp-dev)
  glibc $(pkgid glibc libc6)
  libstdc++ $(pkgid libstdc++ libstdc++6)
Compiler:
  $(${CXX:-c++} --version | head -1)
Compressor:
  UPX $(pkgid upx upx-ucl)
)";
EOF
