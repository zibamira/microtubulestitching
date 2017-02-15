#!/bin/bash
# vim: sw=4
set -o errexit -o nounset -o pipefail -o noglob

githubUrl='https://github.com/zibamira/microtubulestitching.git'
githubLocal='/tmp/github-zibamira-microtubulestitching'
pkgs='dai mpir hxalignmicrotubules'
logPaths="${pkgs} microtubulestitching"

toplevelFiles='
AUTHORS.md
LICENSE.md
README-technical.md
README.md
export-to-github.sh
'

main() {
    warn=

    echo "~~~ FETCH ${githubLocal}"
    if [ -e ${githubLocal} ]; then
        [ "$(git -C ${githubLocal} config remote.origin.url)" = "${githubUrl}" ] ||
            die "Existing ${githubLocal} with wrong remote.origin.url; expected ${githubUrl}."
        gitgh fetch
        gitgh checkout master
        gitgh reset --hard origin/master
    else
        git clone ${githubUrl} ${githubLocal}
    fi
    gitgh show -s
    echo

    echo "~~~ CONFIG"
    prev=$(gitgh show -s --pretty=%s | cut -d '@' -f 2)
    echo "Previous version at GitHub: ${prev}"
    git rev-parse -q --verify refs/tags/${prev} >/dev/null 2>&1 ||
        die "Missing tag ${prev}"

    this=$(git describe | cut -d - -f 1,2 | egrep 'zibedition-\d{4}.\d{2}')
    echo "Latest local version: ${this}"
    git rev-parse -q --verify refs/tags/${this} >/dev/null 2>&1 ||
        die "Missing tag ${this}"

    echo
    echo "Relevant commits: previous..latest -- ${logPaths}"
    git -C .. --no-pager log --oneline --reverse ${prev}..${this} -- ${logPaths}
    echo
    echo "shortlog:"
    git -C .. log ${prev}..${this} -- ${logPaths} | git --no-pager shortlog
    echo

    echo "~~~ COPY ${githubLocal}"
    (
        cd ${githubLocal} &&
        git ls-files -z | xargs -0 rm -f
    )
    for pkg in ${pkgs}; do
        dst="${githubLocal}/${pkg}"
        echo "${dst}..."
        mkdir -p "${dst}"
        git -C .. archive ${this}:${pkg} | gtar -C "${dst}" -xf -
    done
    for top in ${toplevelFiles}; do
        dst="${githubLocal}/${top}"
        echo "mv ${dst}"
        mv "${githubLocal}/hxalignmicrotubules/${top}" "${dst}" || {
            echo "WARNING: ignoring mv error."
            warn=t
        }
    done
    echo

    echo "~~~ COMMIT"
    (
        cd ${githubLocal} &&
        echo '*.am filter=silo -text' >'.gitattributes' &&
        git silo init
    )
    gitgh add -- .
    gitgh commit -s -S -m \
"Update to zib-amira@${this}

Relevant commits since the previous export:

$(git -C .. log --oneline --reverse ${prev}..${this} -- ${logPaths} | sed -e 's/^/    /')
"

    if test ${warn}; then
        echo
        echo 'There have been warnings!'
        echo
    fi
}

gitgh() {
    git -C ${githubLocal} "$@"
}

die() {
    printf >&2 'Error: %s\n' "$1"
    exit 1
}

main "$@"
