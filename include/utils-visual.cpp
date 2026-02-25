#include "utils-visual.h"
#include "TLatex.h"
#include "TH1.h"

void set_histogram_style(TH1* h, int color, int line_width, int marker, double marker_size)
{
        h->SetLineColor(color);
        h->SetLineWidth(line_width);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(marker);
        h->SetMarkerSize(marker_size);
}

void draw_lhcb_tag(TLatex* latex)
{
        latex->SetLineWidth(1);
        latex->DrawLatexNDC(0.7, 0.90, "#font[132]{LHCb}");
        latex->DrawLatexNDC(0.7, 0.85, "#font[132]{pp collisions}");
        latex->DrawLatexNDC(0.7, 0.80, "#font[132]{#sqrt{s} = 13 TeV}");
        latex->DrawLatexNDC(0.7, 0.75, "#font[132]{AK5 Z-Tagged Jets}");
}

void set_lhcb_watermark_properties(TLatex* latex)
{
        latex->SetTextColorAlpha(16, 0.3);
        latex->SetTextSize(0.1991525);
        latex->SetTextAngle(26.15998);
        latex->SetLineWidth(2);
}
