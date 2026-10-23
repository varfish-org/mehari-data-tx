#!/usr/bin/env python3
"""Check for a cdot data release newer than the pinned one.

If there is one, bumps cdot, mehari to the newest version on bioconda and the
HGNC complete set to the newest quarterly set, in config/config.yaml,
workflow/envs/mehari.yaml and MEHARI_VERSION of .github/workflows/build.yml.
A new mehari or HGNC version alone changes nothing. Writes a summary for the
pull request body. Run it from the repository root. With --dry-run, it only
prints the summary and the diff. Uses GITHUB_TOKEN for the GitHub API if set.
"""

import argparse
import difflib
import json
import os
import re
import urllib.request
from pathlib import Path

CONFIG = Path("config/config.yaml")
MEHARI_ENV = Path("workflow/envs/mehari.yaml")
BUILD_WORKFLOW = Path(".github/workflows/build.yml")

BIOCONDA_MEHARI = "https://api.anaconda.org/package/bioconda/mehari"
HGNC_QUARTERLY_LISTING = (
    "https://storage.googleapis.com/storage/v1/b/public-download-files/o"
    "?prefix=hgnc/archive/archive/quarterly/json/hgnc_complete_set_"
)
# a source entry of config.yaml, e.g. GRCh38-ensembl with its cdot release and file name
SOURCE_RE = re.compile(
    r"^  (?P<assembly>\w+)-(?P<kind>ensembl|refseq):\n"
    r"    cdot:\n"
    r"      release: (?P<release>\S+)\n"
    r"      custom: (?P<custom>\S+)$",
    re.MULTILINE,
)


def get_json(url: str):
    request = urllib.request.Request(url)
    if url.startswith("https://api.github.com/"):
        request.add_header("Accept", "application/vnd.github+json")
        if token := os.environ.get("GITHUB_TOKEN"):
            request.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def github(path: str):
    return get_json(f"https://api.github.com/{path}")


def version_key(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def configured_releases(config: str) -> dict[tuple[str, str], list[str]]:
    """The releases lists of the reference section, by assembly and source."""
    releases = {}
    in_reference = False
    assembly = kind = None
    for line in config.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        if indent == 0:
            in_reference = line.startswith("reference:")
        elif in_reference and indent == 2:
            assembly = stripped.rstrip(":")
        elif in_reference and indent == 4:
            kind = stripped.rstrip(":")
        elif in_reference and indent == 6 and stripped.startswith("releases:"):
            items = stripped.split(":", 1)[1].strip().strip("[]").split(",")
            releases[(assembly, kind)] = [
                item.strip().strip("\"'") for item in items if item.strip()
            ]
    return releases


def check_cdot(files: dict[Path, str], summary: list[str]) -> bool:
    """Bump cdot if there is a newer data release that offers all configured releases."""
    config = files[CONFIG]
    sources = list(SOURCE_RE.finditer(config))
    if not sources:
        raise ValueError(f"no cdot sources found in {CONFIG}")
    pinned = max((source["release"] for source in sources), key=version_key)
    data_releases = {
        match[1]: release
        for release in github("repos/SACGF/cdot/releases?per_page=100")
        if (match := re.fullmatch(r"data_v(\d+(?:\.\d+)*)", release["tag_name"]))
        and not release["draft"]
        and not release["prerelease"]
    }
    version = max(data_releases, key=version_key)
    assets = {asset["name"] for asset in data_releases[version]["assets"]}
    listed = configured_releases(config)
    notes = []
    missing = []
    for source in sources:
        assembly, kind = source["assembly"], source["kind"]
        pattern = re.escape(f"cdot-{version}.{source['custom']}.json.gz").replace(
            re.escape(f"{{{kind}_release}}"), r"([A-Za-z0-9_.]+)"
        )
        offered = sorted(
            {match[1] for name in assets if (match := re.fullmatch(pattern, name))}
        )
        configured = listed.get((assembly, kind), [])
        missing += [f"{assembly}-{kind}-{r}" for r in configured if r not in offered]
        if unlisted := [r for r in offered if r not in configured]:
            notes.append(
                f"  - {assembly} {kind}: offered but not listed in config: {', '.join(unlisted)}"
            )
    bumped = False
    if version_key(version) <= version_key(pinned):
        summary.append(f"- cdot: {pinned} is the newest data release.")
    elif missing:
        summary.append(
            f"- cdot: {version} is out, but lacks the configured {', '.join(missing)}; not bumped."
        )
    else:
        files[CONFIG] = re.sub(
            r"^(    cdot:\n      release: )\S+$",
            rf"\g<1>{version}",
            config,
            flags=re.MULTILINE,
        )
        summary.append(f"- cdot: {pinned} -> {version}")
        bumped = True
    summary += notes
    return bumped


def check_mehari(files: dict[Path, str], summary: list[str]):
    match = re.search(r"^mehari:\n  version: (\S+)$", files[CONFIG], re.MULTILINE)
    if not match:
        raise ValueError(f"no mehari version found in {CONFIG}")
    pinned = match[1]
    versions = [
        version
        for version in get_json(BIOCONDA_MEHARI)["versions"]
        if re.fullmatch(r"\d+(?:\.\d+)*", version)
    ]
    version = max(versions, key=version_key)
    if version_key(version) <= version_key(pinned):
        summary.append(f"- mehari: {pinned} is the newest version on bioconda.")
        return
    summary.append(f"- mehari: {pinned} -> {version} (newest on bioconda)")
    files[CONFIG] = re.sub(
        r"^(mehari:\n  version: )\S+(\n  docker: \"docker://ghcr.io/varfish-org/mehari:)[^\"]+",
        rf"\g<1>{version}\g<2>{version}",
        files[CONFIG],
        flags=re.MULTILINE,
    )
    files[MEHARI_ENV] = re.sub(r"(mehari ==)\S+", rf"\g<1>{version}", files[MEHARI_ENV])
    files[BUILD_WORKFLOW] = re.sub(
        r'(MEHARI_VERSION: ")[^"]+', rf"\g<1>{version}", files[BUILD_WORKFLOW]
    )


def check_hgnc(files: dict[Path, str], summary: list[str]):
    pattern = r'^(hgnc:\n(?:  .*\n)*?  version: ")(\d{4}-\d\d-\d\d)(")'
    match = re.search(pattern, files[CONFIG], re.MULTILINE)
    if not match:
        raise ValueError(f"no HGNC version found in {CONFIG}")
    pinned = match[2]
    items = get_json(HGNC_QUARTERLY_LISTING).get("items", [])
    dates = [
        found[1]
        for item in items
        if (
            found := re.search(
                r"hgnc_complete_set_(\d{4}-\d\d-\d\d)\.json$", item["name"]
            )
        )
    ]
    newest = max(dates)
    if newest <= pinned:
        summary.append(f"- HGNC complete set: {pinned} is the newest quarterly set.")
        return
    config = re.sub(pattern, rf"\g<1>{newest}\g<3>", files[CONFIG], flags=re.MULTILINE)
    files[CONFIG] = re.sub(
        r"^(hgnc:\n(?:  .*\n)*?  archive: )\w+",
        r"\g<1>quarterly",
        config,
        flags=re.MULTILINE,
    )
    summary.append(f"- HGNC complete set: {pinned} -> {newest} (quarterly)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the updates and the diff, change no files",
    )
    parser.add_argument("--summary", type=Path, help="write the summary to this file")
    args = parser.parse_args()

    original = {path: path.read_text() for path in (CONFIG, MEHARI_ENV, BUILD_WORKFLOW)}
    files = dict(original)
    summary = ["Result of the weekly update check:", ""]
    # only a new cdot data release triggers an update
    if check_cdot(files, summary):
        check_mehari(files, summary)
        check_hgnc(files, summary)
    text = "\n".join(summary) + "\n"

    if args.dry_run:
        print(text)
        for path, content in files.items():
            print(
                "".join(
                    difflib.unified_diff(
                        original[path].splitlines(keepends=True),
                        content.splitlines(keepends=True),
                        fromfile=f"a/{path}",
                        tofile=f"b/{path}",
                    )
                ),
                end="",
            )
    else:
        for path, content in files.items():
            if content != original[path]:
                path.write_text(content)
    if args.summary:
        args.summary.write_text(text)
    if step_summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(step_summary, "a") as file:
            file.write(text)


if __name__ == "__main__":
    main()
