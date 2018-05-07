//
//  formatMultipanelViewer.h
//  
//
//  Created by Anthony Badea on 5/5/18.
//
//

#ifndef formatMultipanelViewer_h
#define formatMultipanelViewer_h

// C headers
#include <stdio.h>

// root headers
#include <TPad.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLatex.h>

// local headers
#include "styleUtil.h"

void formatMultipanelViewer(TCanvas *multiPanel, TPad *pads[], int rows, int columns)
{
    // Parameters for the canvas division
    static const float leftMargin = 0.35;
    static const float rightMargin = 0.18;
    static const float bottomMargin = 0.35;
    static const float topMargin = 0.22;
    static const float xMargin = 0.0; // for multipanel put xMargin/yMargin = 0.0
    static const float yMargin = 0.0;
    static const float frameWidth = 0.8;
    static const float frameHeight = 0.8;
    static const float yMinOffset = 0.05;
    
    divideCanvas(multiPanel,pads,rows,columns,leftMargin,rightMargin,bottomMargin,topMargin,xMargin,yMargin,frameWidth,frameHeight,yMinOffset);
}
#endif /* formatMultipanelViewer_h */
