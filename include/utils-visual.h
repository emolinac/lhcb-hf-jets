#ifndef UTILS_VISUAL_H
#define UTILS_VISUAL_H

#include "TLatex.h"
#include "TH1.h"

void set_histogram_style(TH1* h, int color, int line_width, int marker, int marker_size);

void draw_lhcb_tag(TLatex* latex);

void set_lhcb_watermark_properties(TLatex* latex);

#endif