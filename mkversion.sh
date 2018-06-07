gitver() {
    git -C "$1" describe --tags --abbrev=0 | sed '{s/^[^0-9.-]*//;s/-/./g}'
}

gitid() {
    git -C "$1" log -1 --date=format:'%Y-%m-%d' --format=format:'commit %h from %cd'
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

pkgid() {
    local pkg
    for pkg in "$@"; do
        pkgid_rpm "$pkg" && return
        pkgid_dpkg "$pkg" && return
    done
    echo "unknown"
}

cat <<EOF
Fuchsia, $(gitid .)
Libraries:
  GiNaC $(gitver "$GINAC"), $(gitid "$GINAC")
  CLN $(gitver "$CLN"), $(gitid "$CLN")
  GMP $(pkgid gmp libgmp-dev)
  glibc $(pkgid glibc libc6)
  libstdc++ $(pkgid libstdc++ libstdc++6)
Compiler:
  GCC $(pkgid gcc-c++ g++)
Compressor:
  UPX $(pkgid upx upx-ucl)
EOF
