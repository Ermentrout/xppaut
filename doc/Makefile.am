LATEX = latex
DVIPS = dvips
PS2PDF = ps2pdf

pkgxbmdir = $(docdir)/xbm

dist_man_MANS = \
	xppaut.1

dist_doc_DATA = \
	changes.txt \
	install.pdf \
	tree.pdf \
	xpp_doc.pdf \
	xpp_sum.pdf

EXTRA_DIST = \
	install.tex \
	tree.tex \
	xpp_doc.tex \
	xpp_sum.tex

dist_pkgxbm_DATA = \
	bitmap/feather.xbm \
	bitmap/galahs.xbm \
	bitmap/heart.xbm \
	bitmap/kitty.xbm \
	bitmap/mf.xbm \
	bitmap/pg.xbm \
	bitmap/santa.xbm \
	bitmap/ste.xbm \
	bitmap/triskele.xbm \
	bitmap/woodgrain.xbm

CLEANFILES = \
	install.aux \
	install.log \
	install.pdf \
	xpp_doc.aux \
	xpp_doc.log \
	xpp_doc.out \
	xpp_doc.pdf \
	xpp_doc.toc \
	xpp_sum.aux \
	xpp_sum.log \
	xpp_sum.pdf


.tex.dvi:
	$(LATEX) $<
	$(LATEX) $<

.dvi.ps:
	$(DVIPS) -o $@ $<

.ps.pdf:
	$(PS2PDF) $<
