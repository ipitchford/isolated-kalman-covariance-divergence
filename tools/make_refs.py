#!/usr/bin/env python3
"""Build refs.bib from Crossref and arXiv records saved in tools/records/.

Every title, author list, year, venue, volume, issue, pages, DOI and arXiv
identifier is copied from a saved record.  If a record is missing it is fetched
once (api.crossref.org/works/DOI or export.arxiv.org/api/query?id_list=ID) and
the raw response is saved before use; existing records are never overwritten.
Only mechanical transformations are applied:
  * HTML entities are unescaped, inline HTML tags removed, and non-ASCII
    letters converted to LaTeX accents;
  * family (given) names entirely in capitals are converted to title case;
  * an author record with no given name is set as one braced name;
  * titles are wrapped in braces so that their capitalisation is kept;
  * page ranges use an en dash;
  * year = year of print publication if present, otherwise of first issue
    (online, posting, or as a last resort record creation, in that order).
The only hand-written entry is the software reference ikf_lib (a @misc with the
pinned commit URL), as permitted for software.  The script stops with an error
if a record is missing and cannot be fetched, or a character cannot be converted.

Usage: python3 -I tools/make_refs.py refs.bib tools/records
"""
import html
import json
import os
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET

CROSSREF = {
    "10.1007/s10514-025-10191-x": "JungLuftWeiss2025",
    "10.15607/RSS.2016.XII.016": "Luft2016",
    "10.1177/0278364918760698": "Luft2018",
    "10.1098/rsos.172348": "Georgescu2018",
    "10.1109/ACC.1997.609105": "JulierUhlmann1997",
    "10.1109/TRA.2002.803461": "RoumeliotisBekey2002",
    "10.1016/B978-1-4831-6716-9.50011-4": "Schmidt1966",
    "10.1109/IVS.2012.6232155": "LiNashashibi2012",
    "10.1109/TAC.1984.1103465": "ChanGoodwinSin1984",
    "10.1162/089976601750541769": "WeissFreeman2001",
    "10.1109/LCSYS.2025.3576671": "DeCarliDimarogonas2025",
}
ARXIV = {
    "2408.12360": "Jung2023thesis",
    "1111.4580": "KhanJadbabaie2011",
    "2603.01404": "DGVIO2026",
}
SOFTWARE = r"""@misc{ikflib,
  author = {{Control of Networked Systems, University of Klagenfurt}},
  title = {{ikf\_lib}},
  year = {2025},
  howpublished = {GitHub repository, pinned commit 4c0d0ec4600b663cb116698d514e055fb371919d},
  url = {https://github.com/aau-cns/ikf_lib/tree/4c0d0ec4600b663cb116698d514e055fb371919d},
  note = {Source inspected 2 November 2026; not compiled or run}
}
"""

ACCENTS = {
    "á": r"{\'a}", "Á": r"{\'A}", "é": r"{\'e}", "É": r"{\'E}", "í": r"{\'i}",
    "ó": r"{\'o}", "ú": r"{\'u}", "ñ": r"{\~n}", "Ñ": r"{\~N}", "ï": r'{\"i}',
    "ü": r'{\"u}', "ö": r'{\"o}', "ä": r'{\"a}', "è": r"{\`e}", "ç": r"{\c c}",
    "Ö": r'{\"O}', "Ü": r'{\"U}', "ß": r"{\ss}",
    "‐": "-", "‑": "-", "–": "--", "—": "---",
    "“": "``", "”": "''", "‘": "`", "’": "'", " ": "~",
}
UA = {"User-Agent": "make_refs (mailto:anonymous@example.org)"}


def tex(s):
    s = re.sub(r"<[^>]+>", "", html.unescape(s))
    s = " ".join(s.split())
    out = []
    for ch in s:
        if ch in ACCENTS:
            out.append(ACCENTS[ch])
        elif ord(ch) < 128:
            out.append({"&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_"}.get(ch, ch))
        else:
            raise SystemExit(f"unconvertible character {ch!r} in {s!r}")
    return "".join(out)


def family(name):
    name = html.unescape(name)
    if name.isupper():
        name = name.capitalize()
    return tex(name)


def given(name):
    name = html.unescape(name)
    if name.isupper():
        name = name.title()
    return tex(name)


def fname(ident):
    return "".join(c if c.isalnum() or c in ".-" else "_" for c in ident)


def fetch(url, path):
    if os.path.exists(path):
        return
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read()
    with open(path, "xb") as fh:
        fh.write(data)


def year(m):
    for k in ("published-print", "issued", "published-online", "posted", "created"):
        dp = m.get(k, {}).get("date-parts")
        if dp and dp[0] and dp[0][0]:
            return dp[0][0]
    raise SystemExit("no year")


def fmt(kind, key, f):
    body = ",\n".join(f"  {k} = {v}" if v.startswith("{{") else f"  {k} = {{{v}}}"
                      for k, v in f.items())
    return f"@{kind}{{{key},\n{body}\n}}\n"


def crossref_entry(doi, key, m):
    authors = " and ".join(
        f"{family(a['family'])}, {given(a['given'])}" if a.get("given")
        else "{" + family(a["family"]) + "}" for a in m["author"])
    f = {"author": authors, "title": "{{" + tex(m["title"][0]) + "}}",
         "year": str(year(m))}
    typ = m["type"]
    container = tex(m["container-title"][0]) if m.get("container-title") else ""
    if typ == "journal-article":
        kind = "article"
        f["journal"] = container
        if m.get("volume"):
            f["volume"] = m["volume"]
        if m.get("issue"):
            f["number"] = m["issue"]
        if m.get("page"):
            f["pages"] = m["page"].replace("-", "--")
        elif m.get("article-number"):
            f["pages"] = m["article-number"]
    elif typ in ("proceedings-article", "book-chapter"):
        kind = "inproceedings" if typ == "proceedings-article" else "incollection"
        f["booktitle"] = container
        if m.get("volume"):
            f["volume"] = m["volume"]
        if m.get("page"):
            f["pages"] = m["page"].replace("-", "--")
        if m.get("publisher"):
            f["publisher"] = tex(m["publisher"])
    else:
        raise SystemExit(f"unhandled type {typ}")
    f["doi"] = m["DOI"]
    return fmt(kind, key, f)


def arxiv_entry(aid, key, raw):
    ns = {"a": "http://www.w3.org/2005/Atom", "x": "http://arxiv.org/schemas/atom"}
    e = ET.fromstring(raw).find("a:entry", ns)
    if e is None or aid not in e.find("a:id", ns).text:
        raise SystemExit(f"arXiv record mismatch for {aid}")
    names = [a.find("a:name", ns).text.strip() for a in e.findall("a:author", ns)]

    def bibname(n):
        parts = n.split()
        return f"{tex(parts[-1])}, {tex(' '.join(parts[:-1]))}"
    f = {"author": " and ".join(bibname(n) for n in names),
         "title": "{{" + tex(e.find("a:title", ns).text) + "}}",
         "year": e.find("a:published", ns).text[:4],
         "eprint": aid, "archivePrefix": "arXiv",
         "howpublished": f"arXiv:{aid}"}
    c = e.find("x:comment", ns)
    if c is not None and c.text:
        f["note"] = "arXiv comment: " + tex(c.text)
    return fmt("misc", key, f)


def main():
    out, rec = sys.argv[1], sys.argv[2]
    os.makedirs(rec, exist_ok=True)
    parts = ["% Generated by tools/make_refs.py from saved Crossref/arXiv records. Do not edit by hand.\n"]
    for doi, key in CROSSREF.items():
        path = os.path.join(rec, "crossref_" + fname(doi) + ".json")
        fetch("https://api.crossref.org/works/" + doi, path)
        with open(path) as fh:
            m = json.load(fh)["message"]
        if m["DOI"].lower() != doi.lower():
            raise SystemExit(f"DOI mismatch {doi} vs {m['DOI']}")
        parts.append(crossref_entry(doi, key, m))
    for aid, key in ARXIV.items():
        path = os.path.join(rec, "arxiv_" + fname(aid) + ".xml")
        fetch("https://export.arxiv.org/api/query?id_list=" + aid, path)
        with open(path, "rb") as fh:
            parts.append(arxiv_entry(aid, key, fh.read()))
    parts.append(SOFTWARE)
    with open(out, "w") as fh:
        fh.write("\n".join(parts))
    print(f"wrote {len(CROSSREF) + len(ARXIV) + 1} entries to {out}")


if __name__ == "__main__":
    main()
