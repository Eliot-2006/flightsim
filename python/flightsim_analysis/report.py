"""Write out/metrics.json, out/index.html (static, Plotly JS inlined) and docs/results.md.

    python -m flightsim_analysis.report [--out out] [--png-dir docs/img]
"""

import argparse
import html
import json
import re
from pathlib import Path

from plotly.offline import get_plotlyjs

from . import metrics as metrics_mod
from . import plots

ROOT = Path(__file__).resolve().parents[2]
PLACEHOLDER = re.compile(r"\{\{\s*([\w.\-]+)\s*(?::([^{}]+?))?\s*\}\}")


def render(template: str, m: dict) -> str:
    """Fill {{key}} or {{key:format}} from the flat metrics dict; a missing key raises KeyError."""

    def fill(match):
        key, fmt = match.group(1), match.group(2)
        value = m[key]
        return format(value, fmt) if fmt else str(value)

    return PLACEHOLDER.sub(fill, template)


def write_results(m: dict, template: Path, dest: Path) -> None:
    dest.write_text(render(template.read_text(), m))


# ---------------------------------------------------------------------------------------------- HTML
CSS = """
:root{--bg:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--rule:#dcdad4;--grid:rgba(0,0,0,.09);--accent:#2a78d6;
--serif:Charter,"Iowan Old Style","Palatino Linotype",Georgia,serif;--mono:ui-monospace,"SF Mono",Menlo,monospace}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])){--bg:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;
--rule:#3a3935;--grid:rgba(255,255,255,.12);--accent:#3987e5}}
:root[data-theme="dark"]{--bg:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--rule:#3a3935;--grid:rgba(255,255,255,.12);--accent:#3987e5}
*{box-sizing:border-box}
html{color-scheme:light dark}
body{margin:0;background:var(--bg);color:var(--ink);font:1.0625rem/1.6 var(--serif)}
::selection{background:var(--accent);color:var(--bg)}
main{max-width:64rem;margin:0 auto;padding:3rem 1rem 5rem}
.text{max-width:68ch}
h1{font-size:2rem;line-height:1.2;margin:0 0 .5rem;text-wrap:balance}
h2{font-size:1.3rem;line-height:1.3;margin:3.5rem 0 .75rem;text-wrap:balance}
p{margin:0 0 1rem}
.meta,figcaption,.note{color:var(--ink2);font-size:.9375rem}
.meta{font-family:var(--mono);font-size:.8125rem}
a{color:var(--accent);text-underline-offset:.2em}
a:focus-visible,summary:focus-visible{outline:2px solid var(--accent);outline-offset:3px}
nav ul{list-style:none;padding:0;margin:1.5rem 0 0;display:flex;flex-wrap:wrap;gap:.25rem 1.25rem;font-size:.9375rem}
figure{margin:1rem 0 0}
figcaption{max-width:68ch;margin-top:.5rem}
figcaption b{color:var(--ink);font-weight:600}
.plot{width:100%}
table{border-collapse:collapse;font-size:.9375rem;margin:1rem 0;width:100%}
caption{text-align:left;color:var(--ink2);padding-bottom:.4rem}
th,td{padding:.35rem .75rem .35rem 0;border-bottom:1px solid var(--rule);text-align:left;vertical-align:top}
th{font-weight:600}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums;font-family:var(--mono);font-size:.875rem;padding-right:1rem}
.scroll{overflow-x:auto}
details{margin:1rem 0}
summary{cursor:pointer}
@media (max-width:40rem){h1{font-size:1.6rem}}
"""

# Re-colour Plotly traces and axes for the active colour scheme. Traces tagged meta=<slot> take the scheme's
# categorical colour; meta="ref" traces (reference slopes, markers) keep their grey.
JS = """
const SCHEMES={light:{cat:%s,ink:'#52514e',grid:'rgba(0,0,0,.09)'},dark:{cat:%s,ink:'#c3c2b7',grid:'rgba(255,255,255,.12)'}};
function scheme(){const t=document.documentElement.dataset.theme;
 return (t||(matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light'))==='dark'?SCHEMES.dark:SCHEMES.light}
function recolor(){const s=scheme();
 document.querySelectorAll('.plotly-graph-div').forEach(gd=>{
  const lay={'font.color':s.ink};
  Object.keys(gd._fullLayout||{}).filter(k=>/^[xy]axis\\d*$/.test(k)).forEach(k=>{lay[k+'.gridcolor']=s.grid;lay[k+'.linecolor']=s.grid;lay[k+'.zerolinecolor']=s.grid});
  Plotly.relayout(gd,lay);
  (gd.data||[]).forEach((tr,i)=>{if(typeof tr.meta==='number'){Plotly.restyle(gd,{'line.color':s.cat[tr.meta],'marker.color':s.cat[tr.meta]},[i])}})})}
addEventListener('load',recolor);
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',recolor);
"""
DARK = ["#3987e5", "#d95926", "#199e70", "#c98500", "#9085e9"]


def _table(caption, header, rows):
    th = "".join(f'<th scope="col"{" class=n" if i else ""}>{html.escape(h)}</th>' for i, h in enumerate(header))
    body = "".join("<tr>" + "".join(
        f'<{"th scope=row" if i == 0 else "td"}{" class=n" if i else ""}>{c}</{"th" if i == 0 else "td"}>'
        for i, c in enumerate(r)) + "</tr>" for r in rows)
    return f'<div class="scroll"><table><caption>{caption}</caption><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>'


def build_html(out: Path, m: dict) -> str:
    sci = lambda v: f"{v:.2e}"
    drift_rows = [[plots.LABEL[n], sci(m[f"r3.{n}.dE"]), sci(m[f"r3.{n}.dH"]), sci(m[f"r3.{n}.dq"])]
                  for n in metrics_mod.INTEGRATORS]
    bench_rows = [[plots.LABEL[n], f"{m[f'bench.step.{n}.ns']:.0f}", f"{m[f'bench.r3run.{n}.us']:.0f}"]
                  for n in ("euler", "rk4", "rkmk4")]
    figs = []
    for f in plots.FIGURES:
        div = f.build(out).to_html(full_html=False, include_plotlyjs=False, div_id=f"fig-{f.id}",
                                   config={"responsive": True, "displaylogo": False})
        figs.append(f'<section aria-labelledby="h-{f.id}"><h2 id="h-{f.id}">{html.escape(f.title)}</h2>'
                    f'<figure><div class="plot">{div}</div><figcaption><b>What it shows.</b> {html.escape(f.caption)}'
                    f'<br><b>Asserted by.</b> {html.escape(f.test_ref)}</figcaption></figure>{_extra(f.id, drift_rows)}'
                    "</section>")
    toc = "".join(f'<li><a href="#h-{f.id}">{html.escape(f.title.split(":")[0])}</a></li>' for f in plots.FIGURES)
    plotly_js = get_plotlyjs()
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>flightsim validation report</title><meta name="color-scheme" content="light dark">
<style>{CSS}</style><script>{plotly_js}</script></head>
<body><main>
<header class="text"><h1>flightsim validation report</h1>
<p>Rigid-body core: {m['orders.in_bounds']} of {m['orders.asserted']} asserted convergence orders are inside their bounds
(first-order methods {m['orders.first.min']:.3f}–{m['orders.first.max']:.3f}, fourth-order {m['orders.fourth.min']:.3f}–{m['orders.fourth.max']:.3f}).
Every figure below is drawn from the output of the <code>flightsim</code> command line; no physics is computed in Python.</p>
<p class="meta">flightsim {m['env.version']} · {html.escape(m['env.git_sha'])} · {html.escape(m['env.compiler'])} · {m['env.build_type']}</p>
<nav aria-label="Figures"><ul>{toc}<li><a href="#h-bench">Benchmarks</a></li><li><a href="#h-aircraft">Aircraft</a></li></ul></nav></header>
{"".join(figs)}
<section aria-labelledby="h-bench"><h2 id="h-bench">Benchmarks</h2>
{_table(f"Release build, single thread, {m['bench.cpus']} logical CPUs. Targets are reported, not gated.", ["Integrator", "Rigid-body step [ns]", "R3 10 s run [µs]"], bench_rows)}</section>
<section aria-labelledby="h-aircraft"><h2 id="h-aircraft">Aircraft</h2>
<p class="text note">F-16 results (trim, linearization, JSBSim cross-validation, autopilot) are coming in M4–M7.</p></section>
</main><script>{JS % (json.dumps(plots.LIGHT), json.dumps(DARK))}</script></body></html>"""


def _extra(fig_id, drift_rows):
    if fig_id != "drift":
        return ""
    return _table("Maximum over the run, R3, h = 0.02 s (values behind the lines above).",
                  ["Integrator", "Energy drift [-]", "|H_N| drift [-]", "|‖q‖ − 1| [-]"], drift_rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=ROOT / "out")
    ap.add_argument("--results", default=ROOT / "docs" / "results.md")
    ap.add_argument("--png-dir", help="also export headline figures as PNG via kaleido")
    args = ap.parse_args(argv)
    out = Path(args.out)
    m = metrics_mod.build(out)
    (out / "metrics.json").write_text(json.dumps(m, indent=2, sort_keys=True) + "\n")
    (out / "index.html").write_text(build_html(out, m))
    write_results(m, ROOT / "docs" / "results.template.md", Path(args.results))
    if args.png_dir:
        Path(args.png_dir).mkdir(parents=True, exist_ok=True)
        for f in plots.FIGURES:
            if f.headline:
                f.build(out).write_image(Path(args.png_dir) / f"{f.id}.png", width=1100, scale=2)
    print(f"report: {out / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
