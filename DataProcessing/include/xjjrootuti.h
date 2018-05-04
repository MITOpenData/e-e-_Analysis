#ifndef _XJJROOT_UTI_H_
#define _XJJROOT_UTI_H_

/* xjjrootuti.h */

#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TH1.h>

namespace xjjroot
{
  const float margin_pad_left = 0.18;
  const float margin_pad_right = 0.043;
  const float margin_pad_bottom = 0.145;
  const float margin_pad_top = 0.1;

  void setgstyle(Int_t padtick=0);
  template <class T> void sethempty(T* hempty, Float_t xoffset=0, Float_t yoffset=0);
  template <class T> void setthgr(T* hempty, Float_t xoffset=0, Float_t yoffset=0);
  template <class T> void setthgrstyle(T* h, Color_t mcolor=-1, Style_t mstyle=-1, Size_t msize=-1, Color_t lcolor=-1, Style_t lstyle=-1, Width_t lwidth=-1, Color_t fcolor=-1, Float_t falpha=-1, Style_t fstyle=-1);
  void drawCMS(TString collision="", TString snn="5.02", Float_t xpos=0, Float_t ypos=0, Bool_t drawenergy=true);
  void settex(TLatex* tex, Float_t tsize=0.04, Short_t align=12, Style_t font=42);
  void drawtex(Double_t x, Double_t y, const char *text, Float_t tsize=0.04, Short_t align=12, Style_t font=42);
  void setleg(TLegend* leg, Float_t tsize=0.04);
  void setlegndraw(TLegend* leg, Float_t tsize=0.04);
  void setline(TLine* l, Color_t lcolor=kBlack, Style_t lstyle=1, Width_t lwidth=2);
  void drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor=kBlack, Style_t lstyle=1, Width_t lwidth=2);
  void drawbox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t fcolor=kGray, Float_t falpha=0.4, Style_t fstyle=1001, Color_t lcolor=0, Style_t lstyle=1, Width_t lwidth=0);

  void dividebinwid(TH1* h);

  void setbranchaddress(TTree* nt, const char* bname, void* addr);
  template <class T> T* copyobject(const T* obj, TString objname);
}

/* ---------- */

void xjjroot::setgstyle(Int_t padtick/*=0*/)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(xjjroot::margin_pad_right);
  gStyle->SetPadLeftMargin(xjjroot::margin_pad_left);
  gStyle->SetPadTopMargin(xjjroot::margin_pad_top);
  gStyle->SetPadBottomMargin(xjjroot::margin_pad_bottom);
  gStyle->SetTitleX(.0f);
  if(padtick)
    {
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
    }
}

template <class T>
void xjjroot::sethempty(T* hempty, Float_t xoffset/*=0*/, Float_t yoffset/*=0*/)
{
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.10+xoffset);
  hempty->GetYaxis()->SetTitleOffset(1.30+yoffset);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.045);
  hempty->GetYaxis()->SetLabelSize(0.045);
  hempty->SetStats(0);
}

template <class T>
void xjjroot::setthgr(T* hempty, Float_t xoffset/*=0*/, Float_t yoffset/*=0*/)
{
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.10+xoffset);
  hempty->GetYaxis()->SetTitleOffset(1.30+yoffset);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.045);
  hempty->GetYaxis()->SetLabelSize(0.045);
  hempty->SetStats(0);
  hempty->SetMarkerStyle(20);
  hempty->SetMarkerSize(1.1);
  hempty->SetMarkerColor(kBlack);
  hempty->SetLineColor(kBlack);
}

template <class T>
void xjjroot::setthgrstyle(T* h, Color_t mcolor/*=-1*/, Style_t mstyle/*=-1*/, Size_t msize/*=-1*/, Color_t lcolor/*=-1*/, Style_t lstyle/*=-1*/, Width_t lwidth/*=-1*/, Color_t fcolor/*=-1*/, Float_t falpha/*=-1*/, Style_t fstyle/*=-1*/)
{
  if(mcolor>=0) h->SetMarkerColor(mcolor);
  if(mstyle>=0) h->SetMarkerStyle(mstyle);
  if(msize>=0)  h->SetMarkerSize(msize);
  if(lcolor>=0) h->SetLineColor(lcolor);
  if(lstyle>=0) h->SetLineStyle(lstyle);
  if(lwidth>=0) h->SetLineWidth(lwidth);
  if(fcolor>=0) h->SetFillColor(fcolor);
  if(falpha>=0) h->SetFillColorAlpha(fcolor, falpha);
  if(fstyle>=0) h->SetFillStyle(fstyle);
}

void xjjroot::drawCMS(TString collision/*=""*/, TString snn/*="5.02"*/, Float_t xpos/*=0*/, Float_t ypos/*=0*/, Bool_t drawenergy/*=true*/)
{
  TLatex* texCms = new TLatex(0.18+xpos,0.93+ypos, "#scale[1.25]{CMS} Preliminary");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();
  if(!drawenergy) return;
  TLatex* texCol = new TLatex(0.96+xpos,0.93+ypos, Form("%s #sqrt{s_{NN}} = %s TeV", collision.Data(), snn.Data()));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);
  texCol->Draw();
}

void xjjroot::settex(TLatex* tex, Float_t tsize/*=0.04*/, Short_t align/*=12*/, Style_t font/*=42*/)
{
  tex->SetNDC();
  tex->SetTextFont(font);
  tex->SetTextAlign(align);
  tex->SetTextSize(tsize);
}

void xjjroot::drawtex(Double_t x, Double_t y, const char* text, Float_t tsize/*=0.04*/, Short_t align/*=12*/, Style_t font/*=42*/)
{
  TLatex* tex = new TLatex(x, y, text);
  xjjroot::settex(tex, tsize, align, font);
  tex->Draw();
}

void xjjroot::setleg(TLegend* leg, Float_t tsize/*=0.04*/)
{
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(tsize);
}

void xjjroot::setlegndraw(TLegend* leg, Float_t tsize/*=0.04*/)
{
  xjjroot::setleg(leg, tsize);
  leg->Draw();
}

void xjjroot::setline(TLine* l, Color_t lcolor/*=kBlack*/, Style_t lstyle/*=1*/, Width_t lwidth/*=2*/)
{
  l->SetLineColor(lcolor);
  l->SetLineStyle(lstyle);
  l->SetLineWidth(lwidth);
}

void xjjroot::drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor/*=kBlack*/, Style_t lstyle/*=1*/, Width_t lwidth/*=2*/)
{
  TLine* l = new TLine(x1, y1, x2, y2);
  xjjroot::setline(l, lcolor, lstyle, lwidth);
  l->Draw();
}

void xjjroot::drawbox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t fcolor/*=kGray*/, Float_t falpha/*=0.4*/, Style_t fstyle/*=1001*/, Color_t lcolor/*=0*/, Style_t lstyle/*=1*/, Width_t lwidth/*=0*/)
{
  TBox* b = new TBox(x1, y1, x2, y2);
  b->SetFillColor(fcolor);
  b->SetFillColorAlpha(fcolor, falpha);
  b->SetFillStyle(fstyle);
  b->SetLineColor(lcolor);
  b->SetLineStyle(lstyle);
  b->SetLineWidth(lwidth);
  b->Draw();
}

/* ----- */

void xjjroot::setbranchaddress(TTree* nt, const char* bname, void* addr)
{
  nt->SetBranchStatus(bname, 1);
  nt->SetBranchAddress(bname, addr);
}

template <class T>
T* xjjroot::copyobject(const T* obj, TString objname)
{
  T* newobj = new T(*obj);
  newobj->SetName(objname);
  return newobj;
}

void xjjroot::dividebinwid(TH1* h)
{
  for(int i=0;i<h->GetNbinsX();i++)
    {
      Float_t val = h->GetBinContent(i+1)/h->GetBinWidth(i+1);
      Float_t valErr = h->GetBinError(i+1)/h->GetBinWidth(i+1);
      h->SetBinContent(i+1,val);
      h->SetBinError(i+1,valErr);
    }
}

#endif
