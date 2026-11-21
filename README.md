# XPPAUT 8.4.2

XPPAUT is a tool for simulating, animating, and analyzing dynamical systems:
differential equations, difference equations, delay equations, functional
equations, boundary value problems, stochastic equations and differential-algebraic
equations.  It combines a graphical interface, many numerical methods,
two- and three-dimensional graphics, and AUTO for bifurcation analysis.

Home page (downloads, tutorial, and more): <https://sites.pitt.edu/~phase/>

## Building

XPPAUT needs a C compiler and the X11 development files.

**Linux (Debian/Ubuntu)**

    sudo apt install build-essential libx11-dev
    make
    ./xppaut ode/lecar.ode

`make install` (as root) installs the program, the manual, the man page and the
examples.  `make uninstall` removes them.

**macOS**: install [XQuartz](https://www.xquartz.org) and the command line tools,
then run `make`.  The Makefile already looks for the X11 headers in `/opt/X11/include`.

**Windows**: see the home page for binaries.

If the zoom rectangle is not drawn, start the program with `-xorfix`:

    ./xppaut -xorfix ode/lecar.ode

## Documentation

| File | What it is |
|------|------------|
| `xpp_doc.pdf` (`xpp_doc.tex`) | the manual |
| `tree.pdf` (`tree.tex`) | every command, menu by menu |
| `xpp_sum.pdf` | a short summary |
| `install.pdf` | installation notes |
| `xppaut.1` | the man page (`man ./xppaut.1`) |
| `help/` | HTML help |

Build the PDFs with `make docs` (needs LaTeX and the figure `DutchBoy.png`).

## Examples

* `ode/` -- several hundred ODE files, with animations (`.ani`), sets, and scripts (`.xps`)
* `tstauto/` -- bifurcation examples for AUTO
* `canonical/` -- canonical models

Try `./xppaut ode/lorenz.ode`.  Examples that go with the manual's AUTO and
scripting chapters: `ode/hom.ode` and `hom.xps` (homoclinic orbit), `ode/het.ode` and
`het.xps` (heteroclinic orbit), `ode/lorenz_lyap.xps` (Lyapunov exponents).

## Scripts

A script is a text file of commands that XPPAUT carries out for you, and XPPAUT
can record one while you work (File, recOrd script):

    ./xppaut ode/hom.ode -script ode/hom.xps

See the scripting chapter of the manual.  With `-nopause -scriptexit` a script
runs unattended; under `xvfb-run` this works in a batch job with no display.

## What is new

See `HISTORY`.  Version 8.4.2 adds scripting and script recording, a new Lyapunov
exponent code (the full spectrum, by the QR method), `Initialconds Periodic`,
homoclinic and heteroclinic continuation in AUTO, phase-reduction tools, and many
bug fixes; it also builds with current compilers.

## License

GNU General Public License, version 2.  See `LICENSE`.
