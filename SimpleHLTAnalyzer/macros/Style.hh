#ifndef STYLE_HH
#define STYLE_HH

#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>

TCanvas* MakeCanvas   (const char* name, const char *title, int dX = 500, int dY = 500);
void     InitSubPad   (TPad* pad, int i);
void     InitHist     (TH1 *hist, const char *xtit, const char *ytit  = "Number of Entries",
                       int color = kBlack, int style = 0);
void     InitGraph    (TGraph *gr, const char *title, const char *xtit, const char *ytit, double xmin, double xmax, double ymin, double ymax, int color = kBlack);
void     SetStyle     ();

TCanvas* MakeCanvas(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);
  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.18);
  canvas->SetRightMargin    (0.05);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.15);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void InitSubPad(TPad* pad, int i)
{
  pad->cd(i);
  TPad *tmpPad = (TPad*) pad->GetPad(i);
  tmpPad->SetLeftMargin  (0.18);
  tmpPad->SetTopMargin   (0.05);
  tmpPad->SetRightMargin (0.07);
  tmpPad->SetBottomMargin(0.15);
  return;
}

void InitHist(TH1 *hist, const char *xtit, const char *ytit, int color, int style)
{
  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(    3.);
  hist->SetFillColor(color );
  hist->SetFillStyle(style );
  hist->SetTitleSize  (0.055,"Y");
  hist->SetTitleOffset(1.600,"Y");
  hist->SetLabelOffset(0.014,"Y");
  hist->SetLabelSize  (0.050,"Y");
  hist->SetLabelFont  (42   ,"Y");
  hist->SetTitleSize  (0.055,"X");
  hist->SetTitleOffset(1.300,"X");
  hist->SetLabelOffset(0.014,"X");
  hist->SetLabelSize  (0.050,"X");
  hist->SetLabelFont  (42   ,"X");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize (0.6);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->SetTitle("");  
  return;
}

void InitGraph(TGraph *gr, const char *title, const char *xtit, const char *ytit, double xmin, double xmax, double ymin, double ymax, int color)
{
  gr->SetTitle(title);
  gr->GetXaxis()->SetTitle(xtit);
  gr->GetYaxis()->SetTitle(ytit);
  gr->SetLineColor(color);
  gr->SetLineWidth(    3.);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(color);
  gr->SetMarkerSize (0.6);
  gr->GetYaxis()->SetTitleSize  (0.055);
  gr->GetYaxis()->SetTitleOffset(1.600);
  gr->GetYaxis()->SetLabelOffset(0.014);
  gr->GetYaxis()->SetLabelSize  (0.050);
  gr->GetYaxis()->SetLabelFont  (42   );
  gr->GetXaxis()->SetTitleSize  (0.055);
  gr->GetXaxis()->SetTitleOffset(1.300);
  gr->GetXaxis()->SetLabelOffset(0.014);
  gr->GetXaxis()->SetLabelSize  (0.050);
  gr->GetXaxis()->SetLabelFont  (42   );
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetLimits(xmin, xmax);
  gr->GetYaxis()->SetLimits(ymin, ymax);
  return;
}

void SetStyle()
{
  TStyle *MyStyle = new TStyle("New Style","");
  gStyle = MyStyle;

  // Canvas
  MyStyle->SetCanvasColor     (0);
  MyStyle->SetCanvasBorderSize(10);
  MyStyle->SetCanvasBorderMode(0);
  MyStyle->SetCanvasDefH      (700);
  MyStyle->SetCanvasDefW      (700);
  MyStyle->SetCanvasDefX      (100);
  MyStyle->SetCanvasDefY      (100);
  //
  // color palette for 2D temperature plots
  MyStyle->SetPalette(1,0);
  //
  // Pads
  MyStyle->SetPadColor       (0);
  MyStyle->SetPadBorderSize  (10);
  MyStyle->SetPadBorderMode  (0);
  MyStyle->SetPadBottomMargin(0.13);
  MyStyle->SetPadTopMargin   (0.08);
  MyStyle->SetPadLeftMargin  (0.15);
  MyStyle->SetPadRightMargin (0.05);
  MyStyle->SetPadGridX       (0);
  MyStyle->SetPadGridY       (0);
  MyStyle->SetPadTickX       (1);
  MyStyle->SetPadTickY       (1);
  //
  // Frames
  MyStyle->SetLineWidth(3);
  MyStyle->SetFrameFillStyle ( 0);
  MyStyle->SetFrameFillColor ( 0);
  MyStyle->SetFrameLineColor ( 1);
  MyStyle->SetFrameLineStyle ( 0);
  MyStyle->SetFrameLineWidth ( 2);
  MyStyle->SetFrameBorderSize(10);
  MyStyle->SetFrameBorderMode( 0);
  //
  // Histograms
  MyStyle->SetHistFillColor(2);
  MyStyle->SetHistFillStyle(0);
  MyStyle->SetHistLineColor(1);
  MyStyle->SetHistLineStyle(0);
  MyStyle->SetHistLineWidth(3);
  //MyStyle->SetNdivisions(101);
  //
  // Functions
  MyStyle->SetFuncColor(1);
  MyStyle->SetFuncStyle(0);
  MyStyle->SetFuncWidth(2);
  //
  // Various
  MyStyle->SetMarkerStyle(20);
  MyStyle->SetMarkerColor(kBlack);
  MyStyle->SetMarkerSize (1.4);
  //
  MyStyle->SetTitleBorderSize(0);
  MyStyle->SetTitleFillColor (0);
  MyStyle->SetTitleX         (0.2);
  //
  MyStyle->SetTitleSize  (0.055,"X");
  MyStyle->SetTitleOffset(1.200,"X");
  MyStyle->SetLabelOffset(0.005,"X");
  MyStyle->SetLabelSize  (0.050,"X");
  MyStyle->SetLabelFont  (42   ,"X");
  //
  MyStyle->SetStripDecimals(kFALSE);
  MyStyle->SetLineStyleString(11,"20 10");
  //
  MyStyle->SetTitleSize  (0.055,"Y");
  MyStyle->SetTitleOffset(1.600,"Y");
  MyStyle->SetLabelOffset(0.010,"Y");
  MyStyle->SetLabelSize  (0.050,"Y");
  MyStyle->SetLabelFont  (42   ,"Y");
  //
  MyStyle->SetTextSize   (0.055);
  MyStyle->SetTextFont   (42);
  //
  MyStyle->SetStatFont   (42);
  MyStyle->SetTitleFont  (42);
  MyStyle->SetTitleFont  (42,"X");
  MyStyle->SetTitleFont  (42,"Y");
  //
  MyStyle->SetOptStat    (0);

  return;
}

#endif
