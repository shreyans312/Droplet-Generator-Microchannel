# Figure 6b data provenance

`fig6_digitized.csv` contains all 18 markers in Figure 6b of the retained
paper. The values were extracted from the PDF's vector paths rather than read
by eye.

The page-space calibration was:

$$
K_p=\frac{x-407.785159}{46.98103},
$$

$$
\frac{D_d}{w}=1.76+\frac{215.230468-y}{191.186534}.
$$

Marker centers give `kp` and `d_over_w`. The upper and lower cap coordinates
of each vector error bar give `uncertainty`. The paper states that these error
bars are experimental standard deviations.

The PDF artwork and extracted data are not author-supplied raw measurements.
Small coordinate-rounding uncertainty remains, but it is much smaller than
the plotted experimental error bars.

Figure 6's embedded raster is 563 pixels wide and is placed at a scale of
0.239753 PDF points per pixel. Its vector scale bar is 9.875 points for
200 um, or about 41.2 raster pixels. The widened droplet-channel gap is about
38-41 pixels, so the implemented region-B width is rounded to 200 um. This is
an image-derived assumption, not an author-supplied dimension.
