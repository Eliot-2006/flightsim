"""Write out/metrics.json, out/index.html (static, Plotly JS inlined), docs/results.md and README.md.

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
html{color-scheme:light dark;scrollbar-color:var(--rule) var(--bg);scroll-behavior:smooth}
@media (prefers-reduced-motion:reduce){html{scroll-behavior:auto}}
body{margin:0;background:var(--bg);color:var(--ink);font:1.0625rem/1.6 var(--serif)}
::selection{background:var(--accent);color:var(--bg)}
main{max-width:64rem;margin:0 auto;padding:3.5rem 1rem 5rem}
.text{max-width:68ch}
h1{font-size:clamp(2.25rem,5vw,3.25rem);line-height:1.05;letter-spacing:-.02em;margin:0 0 1rem}
.lede{font-size:1.3rem;line-height:1.45;margin:0 0 1rem;text-wrap:pretty}
h2{font-size:1.75rem;line-height:1.2;letter-spacing:-.01em;margin:5rem 0 .75rem;padding-top:1.25rem;border-top:1px solid var(--rule);text-wrap:balance}
h3{font-size:1.2rem;line-height:1.3;margin:3rem 0 .5rem;text-wrap:balance}
p{margin:0 0 1rem}
.meta,figcaption,.note{color:var(--ink2);font-size:.9375rem}
.meta{font-family:var(--mono);font-size:.8125rem}
a{color:var(--accent);text-underline-offset:.2em;text-decoration-thickness:1px}
a:hover{text-decoration-thickness:2px}
a:focus-visible,summary:focus-visible{outline:2px solid var(--accent);outline-offset:3px;border-radius:2px}
.links{display:flex;flex-wrap:wrap;gap:.25rem 1.5rem;margin:1.25rem 0 .75rem;font-size:1rem}
.skip{position:absolute;left:-999px}.skip:focus{left:1rem;top:1rem;background:var(--bg);padding:.25rem .5rem}
.toc{display:grid;grid-template-columns:repeat(auto-fit,minmax(12.5rem,1fr));gap:1.5rem 2rem;margin:1rem 0 0;font-size:.9375rem}
.toc ul{list-style:none;padding:0;margin:.35rem 0 0}
.toc li{margin:.15rem 0}
.toc>div>a{font-weight:600;font-size:1rem}
.chapter-intro{max-width:68ch}
.toc-h{margin-top:3rem}
.top{font-size:.875rem;margin-top:2rem}
figure{margin:1rem 0 0}
figcaption{max-width:68ch;margin-top:.5rem;overflow-wrap:anywhere}
p,caption,td,th{overflow-wrap:break-word}
figcaption b{color:var(--ink);font-weight:600}
code{font-family:var(--mono);font-size:.85em}
.plot{width:100%;max-width:100%;overflow:hidden}
.plot .plotly-graph-div{max-width:100%}
figure,section{min-width:0}
table{border-collapse:collapse;font-size:.9375rem;margin:1rem 0;width:100%}
caption{text-align:left;color:var(--ink2);padding-bottom:.4rem}
th,td{padding:.45rem .75rem .45rem 0;border-bottom:1px solid var(--rule);text-align:left;vertical-align:top}
thead th{font-weight:600;border-bottom-color:var(--ink2)}
tbody th{font-weight:400}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums;font-family:var(--mono);font-size:.875rem;padding-right:1rem}
.glance td.v{font-variant-numeric:tabular-nums;white-space:nowrap}
.scroll{overflow-x:auto}
details{margin:1rem 0}
summary{cursor:pointer}
@media (max-width:40rem){.lede{font-size:1.15rem}h2{margin-top:3.5rem}.glance td.v{white-space:normal}}
"""

# Re-colour Plotly traces and axes for the active colour scheme. Traces tagged meta=<slot> take the scheme's
# categorical colour; meta="ref" traces (reference slopes, markers) keep their grey. 3D (WebGL) scenes ignore alpha in
# grid colours, so they get solid equivalents of the rule colour.
JS = """
const SCHEMES={light:{cat:%s,ink:'#52514e',grid:'rgba(0,0,0,.09)',grid3:'#e4e2dc'},dark:{cat:%s,ink:'#c3c2b7',grid:'rgba(255,255,255,.12)',grid3:'#3a3935'}};
function scheme(){const t=document.documentElement.dataset.theme;
 return (t||(matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light'))==='dark'?SCHEMES.dark:SCHEMES.light}
function recolor(){const s=scheme();
 document.querySelectorAll('.plotly-graph-div').forEach(gd=>{
  const lay={'font.color':s.ink};
  Object.keys(gd._fullLayout||{}).filter(k=>/^[xy]axis\\d*$/.test(k)).forEach(k=>{lay[k+'.gridcolor']=s.grid;lay[k+'.linecolor']=s.grid;lay[k+'.zerolinecolor']=s.grid});
  Object.keys(gd._fullLayout||{}).filter(k=>/^scene\\d*$/.test(k)).forEach(k=>{['xaxis','yaxis','zaxis'].forEach(a=>{lay[k+'.'+a+'.gridcolor']=s.grid3;lay[k+'.'+a+'.linecolor']=s.grid3;lay[k+'.'+a+'.color']=s.ink})});
  Plotly.relayout(gd,lay);
  (gd.data||[]).forEach((tr,i)=>{if(typeof tr.meta==='number'){Plotly.restyle(gd,{'line.color':s.cat[tr.meta],'marker.color':s.cat[tr.meta]},[i])}})})}
addEventListener('load',recolor);
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',recolor);
"""
DARK = ["#3987e5", "#d95926", "#199e70", "#c98500", "#9085e9"]

REPO = "https://github.com/Eliot-2006/flightsim"
SITE = "https://eliot-2006.github.io/flightsim/"

# Page chapters: (anchor, title, intro template, [(figure id, nav label)]). Every figure in plots.FIGURES must appear
# exactly once (asserted in build_html), so a new figure cannot silently go missing from the page.
CHAPTERS = [
    ("f16-jsbsim", "F-16 against JSBSim",
     "JSBSim, an independent open-source flight dynamics engine, flies an aircraft generated from the same data file flightsim "
     "reads, on a matched non-rotating planet with the same gravity and atmosphere. Both simulators trim and fly the same three "
     "maneuvers. The worst difference over every signal is {{jsbsim.j3.worst:.2f}} % of that signal's peak (criterion 2 %), and "
     "the discrepancy budget accounts for what remains, source by source. This checks the implementation (equations of motion, "
     "frames, interpolation, integration), not how faithfully the tables describe a real F-16.",
     [("j3-elevator", "Elevator doublet"), ("j3-aileron", "Aileron pulse"), ("j3-rudder", "Rudder doublet"),
      ("flight-3d", "3D flight path"), ("j3-budget", "Discrepancy budget")]),
    ("f16-linear", "Trim and linear model",
     "flightsim solves for the angle of attack, throttle and elevator that hold steady level flight, then linearizes the nonlinear "
     "model about that point to get the aircraft's modes (short period, phugoid, Dutch roll, roll, spiral). Against published "
     "values, the trim agrees within {{trim.max_rel_pct:.3f}} % and the eigenvalues within {{linear.modes.max_rel_pct:.2f}} %; those "
     "references come from a single secondary source, as the full results explain. The figures check the linearization itself: "
     "a ten-times smaller input shrinks the linear–nonlinear gap about a hundred times ({{linear.ln2.ratio:.0f}}×), and the "
     "finite-difference step sits in its flat region.",
     [("ln2", "Linear vs nonlinear response"), ("ln3", "Finite-difference step sweep")]),
    ("autopilot", "Autopilot",
     "A cascaded altitude and airspeed autopilot (pitch-rate, pitch and altitude loops on the elevator; airspeed on the "
     "throttle), sampled at 50 Hz through actuator lag, rate and position limits. The gains were designed on the linear model, "
     "with loop margins of at least {{control.c5.worst_gm:.0f}} dB and {{control.c5.worst_pm:.0f}}°, and were not tuned on the "
     "nonlinear aircraft; these figures verify them there.",
     [("c1-overlay", "+100 ft altitude step"), ("c2-overlay", "+1,000 ft altitude step"), ("c2-windup", "Anti-windup"),
      ("c3-overlay", "+50 ft/s airspeed step"), ("c4-disturbance", "Gust and headwind"), ("control-actuators", "Actuator usage")]),
    ("numerics", "Numerical core",
     "The rigid-body engine and its five integrators are checked against closed-form solutions. On these log–log plots, slopes "
     "of 1 and 4 are the measured convergence orders, and {{orders.in_bounds}} of {{orders.asserted}} of them are asserted in CI. "
     "Energy and angular momentum are checked on torque-free rotation, including the intermediate-axis flip.",
     [("convergence", "Convergence orders"), ("work-precision", "Work–precision"), ("drift", "Energy and momentum drift"),
      ("qnorm", "Quaternion norm"), ("r4", "Intermediate-axis flip"), ("t4", "Harmonic oscillator")]),
]


def _table(caption, header, rows, cls=""):
    th = "".join(f'<th scope="col"{" class=n" if i else ""}>{html.escape(h)}</th>' for i, h in enumerate(header))
    body = "".join("<tr>" + "".join(
        f'<{"th scope=row" if i == 0 else "td"}{" class=n" if i else ""}>{c}</{"th" if i == 0 else "td"}>'
        for i, c in enumerate(r)) + "</tr>" for r in rows)
    return (f'<div class="scroll"><table{f" class={cls}" if cls else ""}><caption>{caption}</caption><thead><tr>{th}</tr>'
            f"</thead><tbody>{body}</tbody></table></div>")


def _code(text: str) -> str:
    """Escape caption text and render `backticked` spans as <code>."""
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", html.escape(text))


def _glance(m: dict) -> str:
    """Results at a glance: what was checked, the measured value, the requirement, and where the evidence is."""
    rows = [
        ("Agrees with JSBSim flying the same F-16 data", f"{m['jsbsim.j3.worst']:.2f} % of peak, worst signal", "≤ 2 %", "#f16-jsbsim"),
        ("Reproduces the published F-16 state-derivative check case", f"{m['checks.a1.max_rel']:.1e} relative", "≤ 1e-3", REPO + "/blob/main/docs/results.md"),
        ("Trim matches published values", f"{m['trim.max_rel_pct']:.3f} %", "≤ 0.5 %", "#f16-linear"),
        ("Linear-model eigenvalues match published values", f"{m['linear.modes.max_rel_pct']:.2f} %", "≤ 1 %", "#f16-linear"),
        ("Autopilot, +100 ft altitude step", f"{m['control.c1.overshoot']:.1f} % overshoot, settles in {m['control.c1.settling']:.0f} s", "≤ 10 %, ≤ 25 s", "#h-c1-overlay"),
        ("Anti-windup, +1,000 ft step", f"settles in {m['control.c2.settling']:.0f} s (without: {m['control.c2n.settling']:.0f} s)", "≤ 120 s", "#h-c2-windup"),
        ("Autopilot loop margins (linear design)", f"≥ {m['control.c5.worst_gm']:.0f} dB, ≥ {m['control.c5.worst_pm']:.0f}°", "≥ 6 dB, ≥ 45°", "#autopilot"),
        ("Integrator convergence orders", f"{m['orders.in_bounds']} of {m['orders.asserted']} within theory", "order ±0.05 / ±0.1", "#numerics"),
        ("Speed, 60 s closed-loop F-16 flight", f"{m['bench.f16.realtime_factor']:,.0f}× real time", "≥ 1,000×", "#h-bench"),
    ]
    body = "".join(f'<tr><th scope="row">{html.escape(a)}</th><td class="v">{html.escape(b)}</td><td class="v">{html.escape(c)}</td>'
                   f'<td><a href="{d}">{"Evidence" if d.startswith("#") else "Results table"}</a></td></tr>' for a, b, c, d in rows)
    return ('<div class="scroll"><table class="glance"><caption>Results at a glance. Every value is computed by this build; '
            'requirements are enforced by automated tests.</caption><thead><tr><th scope="col">Check</th><th scope="col">Measured</th>'
            f'<th scope="col">Requirement</th><th scope="col"><span class="sr">Link</span></th></tr></thead><tbody>{body}</tbody></table></div>')


def build_html(out: Path, m: dict) -> str:
    sci = lambda v: f"{v:.2e}"
    drift_rows = [[plots.LABEL[n], sci(m[f"r3.{n}.dE"]), sci(m[f"r3.{n}.dH"]), sci(m[f"r3.{n}.dq"])]
                  for n in metrics_mod.INTEGRATORS]
    bench_rows = [[plots.LABEL[n], f"{m[f'bench.step.{n}.ns']:.0f}", f"{m[f'bench.r3run.{n}.us']:.0f}"]
                  for n in ("euler", "rk4", "rkmk4")]
    f16_rows = [[label, target, f"{m[key] / div:.4g} {unit}", "met" if m[key] <= limit else "missed"]
                for key, label, limit, target, div, unit in metrics_mod.BENCH_TARGETS]
    by_id = {f.id: f for f in plots.FIGURES}
    placed = [fid for _, _, _, figs in CHAPTERS for fid, _ in figs]
    assert sorted(placed) == sorted(by_id), f"figures and chapters disagree: {set(placed) ^ set(by_id)}"

    def figure(f):
        div = f.build(out).to_html(full_html=False, include_plotlyjs=False, div_id=f"fig-{f.id}",
                                   config={"responsive": True, "displaylogo": False})
        return (f'<section aria-labelledby="h-{f.id}"><h3 id="h-{f.id}">{html.escape(f.title)}</h3>'
                f'<figure><div class="plot">{div}</div><figcaption><b>What it shows.</b> {_code(f.caption)}'
                f'<br><b>Checked by.</b> {_code(f.test_ref)}</figcaption></figure>{_extra(f.id, drift_rows)}</section>')

    chapters = "".join(
        f'<section aria-labelledby="{cid}"><h2 id="{cid}">{html.escape(title)}</h2>'
        f'<p class="chapter-intro">{html.escape(render(intro, m))}</p>'
        + "".join(figure(by_id[fid]) for fid, _ in figs)
        + '<p class="top"><a href="#contents">Back to contents</a></p></section>'
        for cid, title, intro, figs in CHAPTERS)
    toc = "".join(
        f'<div><a href="#{cid}">{html.escape(title)}</a><ul>'
        + "".join(f'<li><a href="#h-{fid}">{html.escape(label)}</a></li>' for fid, label in figs) + "</ul></div>"
        for cid, title, _, figs in CHAPTERS) + '<div><a href="#h-bench">Benchmarks</a></div>'
    desc = ("A C++20 six-degree-of-freedom flight dynamics simulator with an F-16 model from NASA wind-tunnel data, "
            "validated against closed-form physics, published reference values and JSBSim.")
    plotly_js = get_plotlyjs()
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>flightsim: a validated 6-DOF flight dynamics simulator</title><meta name="color-scheme" content="light dark">
<meta name="description" content="{desc}">
<meta property="og:title" content="flightsim: a validated 6-DOF flight dynamics simulator">
<meta property="og:description" content="{desc}"><meta property="og:type" content="website">
<meta property="og:url" content="{SITE}"><meta property="og:image" content="{SITE}og.png">
<meta name="twitter:card" content="summary_large_image">
<style>{CSS}.sr{{position:absolute;width:1px;height:1px;overflow:hidden;clip-path:inset(50%)}}</style><script>{plotly_js}</script></head>
<body><a class="skip" href="#contents">Skip to contents</a><main>
<header class="text"><h1>flightsim</h1>
<p class="lede">A six-degree-of-freedom flight dynamics simulator in C++20, with a nonlinear F-16 model built from NASA
wind-tunnel data, a trim solver, numerical linearization and an altitude and airspeed autopilot.</p>
<p>Everything on this page is evidence. Each figure is drawn from the simulator's own output and named with the automated
test that checks it. The aircraft is cross-validated against JSBSim, an independent flight dynamics engine flying the same data.</p>
<p class="links"><a href="{REPO}">Source code on GitHub</a><a href="{REPO}#readme">README</a>
<a href="{REPO}/blob/main/docs/design.md">Design and decision log</a><a href="{REPO}/blob/main/docs/results.md">All result tables</a></p>
<p class="meta">Built from commit {html.escape(m['env.git_sha'])} · flightsim {m['env.version']} · {html.escape(m['env.compiler'])} · {m['env.build_type']}</p></header>
{_glance(m)}
<p class="note text">The published reference values (check case, trim, eigenvalues) come from secondary reproductions of the
textbook tables, a single lineage, because the book itself was not available. The physics-law checks, the analytic cases and
the JSBSim comparison do not depend on them.</p>
<nav id="contents" aria-label="Contents"><h2 class="toc-h">Contents</h2><div class="toc">{toc}</div></nav>
{chapters}
<section aria-labelledby="h-bench"><h2 id="h-bench">Benchmarks</h2>
<p class="text note">{html.escape(m['host.cpu'])}, {html.escape(m['host.os'])}, {html.escape(m['host.compiler'])}, {m['host.build_type']}; median of {m['host.bench_repetitions']} repetitions. Full methodology in <a href="{REPO}/blob/main/docs/results.md">docs/results.md</a>.</p>
{_table("F-16 against the PRD §11 targets. Targets are reported, not gated.", ["Benchmark", "Target", "Measured", "Status"], f16_rows)}
{_table("Integrator cost.", ["Integrator", "Rigid-body step [ns]", "R3 10 s run [µs]"], bench_rows)}
<p class="top"><a href="#contents">Back to contents</a></p></section>
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
    write_results(m, ROOT / "docs" / "README.template.md", ROOT / "README.md")
    if args.png_dir:
        Path(args.png_dir).mkdir(parents=True, exist_ok=True)
        for f in plots.FIGURES:
            if f.headline:
                f.build(out).write_image(Path(args.png_dir) / f"{f.id}.png", width=1100, scale=2)
    print(f"report: {out / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
