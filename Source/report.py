"""Build the Overleaf project and compile its PDF using pdfLaTeX and BibTeX."""
from pathlib import Path
import shutil
import subprocess
import tempfile
from zipfile import ZipFile, ZIP_DEFLATED
from model import ROOT

template = Path(__file__).resolve().parents[1] / "Report" / "LaTeX"
report = ROOT / "Report"
source = report / "LaTeX"
images = source / "figures"
images.mkdir(parents=True, exist_ok=True)
for name in ("main.tex", "references.bib"):
    if (template/name).resolve() != (source/name).resolve():
        shutil.copy2(template/name, source/name)
shutil.copy2(report/"figures"/"cad_cutaway.png", images/"cad_cutaway.png")
assets = [source/"main.tex", source/"references.bib", images/"cad_cutaway.png"]
assets.extend(images/(name+".pdf") for name in ("mesh","mach_field","four_cases","shock_detail","residuals"))
package = report/"Air_Nozzle_Overleaf.zip"
with ZipFile(package, "w", ZIP_DEFLATED) as archive:
    for asset in assets:
        archive.write(asset, asset.relative_to(source))
# Compile from the actual archive, independently of all files in the source tree.
with tempfile.TemporaryDirectory() as tmp:
    with ZipFile(package) as archive:
        archive.extractall(tmp)
    for command in (["pdflatex","-halt-on-error","-interaction=nonstopmode","main.tex"],
                    ["bibtex","main"],
                    ["pdflatex","-halt-on-error","-interaction=nonstopmode","main.tex"],
                    ["pdflatex","-halt-on-error","-interaction=nonstopmode","main.tex"]):
        subprocess.run(command, cwd=tmp, check=True)
    log=(Path(tmp)/"main.log").read_text(errors="replace")
    if "undefined references" in log or "undefined citations" in log:
        raise RuntimeError("Unresolved cross-references in the final LaTeX pass.")
    shutil.copyfile(Path(tmp)/"main.pdf", report/"Air_Nozzle_Report.pdf")
