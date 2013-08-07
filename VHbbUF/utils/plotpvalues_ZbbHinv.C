{
///
/// macro adpated from Michele de Gruttola
///

    gROOT->SetStyle("Plain");
    gROOT->LoadMacro("tdrstyle.C");
    setTDRStyle();
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetLineStyleString(11,"16 16");

    bool unblind = false;
    int n_points = 6;
    double x_vals[n_points] = {105, 115, 125, 135, 145, 150};


    // -------------------------------------------------------------------------

    // ZbbHinv (2013-07-31)
    
    double y_observed[n_points]=    {0.141247,0.161316,0.179391,0.199058,0.214818,0.228364,};
    // using --toysFreq
    double y_vals[n_points]=        {0.141247,0.161316,0.179391,0.199058,0.214818,0.228364,};
    // not using --toyFreq
    //double y_vals[n_points]=        {0.141247,0.161316,0.179391,0.199058,0.214818,0.228364,};
    
    // -------------------------------------------------------------------------
    
    // Expected
    m_y_line_graph = new TGraph(n_points, x_vals, y_vals);
    m_y_line_graph->SetLineWidth(2);
    m_y_line_graph->SetLineStyle(11);
    m_y_line_graph->SetLineColor(kBlue);
    //m_y_line_graph->SetFillColor(kWhite);

    // Observed
    m_y_lineObs_graph = new TGraph(n_points, x_vals, y_observed);
    m_y_lineObs_graph->SetLineWidth(1);
    //m_y_lineObs_graph->SetLineColor(kBlack);
    //m_y_lineObs_graph->SetFillColor(kWhite);
    m_y_lineObs_graph->SetMarkerStyle(20);
    m_y_lineObs_graph->SetMarkerSize(1.1);

    // evil axes
    m_y_line_graph->GetXaxis()->SetNdivisions(505);
    m_y_line_graph->GetXaxis()->SetTitleSize(0.05);
    m_y_line_graph->GetXaxis()->SetTitleOffset(1.0);
    m_y_line_graph->GetXaxis()->SetTitle("m_{H} [GeV]");
    //m_y_line_graph->GetXaxis()->SetLimits(103, 152);
    m_y_line_graph->GetXaxis()->SetLimits(x_vals[0]-2, x_vals[n_points-1]+2);
    
    //m_y_line_graph->GetYaxis()->SetNdivisions(510);
    m_y_line_graph->GetYaxis()->SetTitleSize(0.048);
    m_y_line_graph->GetYaxis()->SetTitleOffset(1.25);
    m_y_line_graph->GetYaxis()->SetTitle("Local p-value");
    m_y_line_graph->GetYaxis()->SetLabelSize(0.045);
    
    m_y_line_graph->SetMaximum(1);
    m_y_line_graph->SetMinimum(8e-4 * 0.90);

    // y = n sigma lines
    double y_sigmas[5]={0.158655254,0.022750132,0.001349898,0.000031671,0.00000028665};
    m_one_line   = new TLine(x_vals[0],y_sigmas[0],x_vals[n_points-1],y_sigmas[0]);
    m_two_line   = new TLine(x_vals[0],y_sigmas[1],x_vals[n_points-1],y_sigmas[1]);
    m_three_line = new TLine(x_vals[0],y_sigmas[2],x_vals[n_points-1],y_sigmas[2]);
    m_four_line  = new TLine(x_vals[0],y_sigmas[3],x_vals[n_points-1],y_sigmas[3]);
    m_five_line  = new TLine(x_vals[0],y_sigmas[4],x_vals[n_points-1],y_sigmas[4]);
    m_one_line->SetLineColor(2);
    m_one_line->SetLineWidth(2);
    m_two_line->SetLineColor(2);
    m_two_line->SetLineWidth(2);
    m_three_line->SetLineColor(2);
    m_three_line->SetLineWidth(2);
    m_four_line->SetLineColor(2);
    m_four_line->SetLineWidth(2);
    m_five_line->SetLineColor(2);
    m_five_line->SetLineWidth(2);

    // Legend
    TLegend* m_legend = 0;
    if (unblind) {
        m_legend = new TLegend(0.64,0.24,0.96,0.36);
        m_legend->SetTextFont(42);
        m_legend->SetFillStyle(0);
        //m_legend->SetFillColor(0);
        //m_legend->SetLineColor(kBlack);
        m_legend->SetBorderSize(0);
        //m_legend->SetBorderSize(1);
        m_legend->AddEntry(m_y_lineObs_graph, "Observed", "lp");
        m_legend->AddEntry(m_y_line_graph,"Expected", "l");
        
    } else {
        m_legend = new TLegend(0.64,0.24,0.96,0.30);
        m_legend->SetTextFont(42);
        m_legend->SetFillStyle(0);
        //m_legend->SetFillColor(0);
        //m_legend->SetLineColor(kBlack);
        m_legend->SetBorderSize(0);
        //m_legend->SetBorderSize(1);
        m_legend->AddEntry(m_y_line_graph,"Expected", "l");
    }


    // -------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas();
    //TCanvas *c1 = new TCanvas("c1","c1",300,300,700,500);
    //c1->SetGridx();
    //c1->SetGridy();
    c1->SetLogy();
    
    // Draw
    m_y_line_graph->Draw("AL");
    if (unblind) {
        m_y_lineObs_graph->Draw("LP");
    }
     
    m_one_line->Draw();
    m_two_line->Draw();
    m_three_line->Draw();
    m_four_line->Draw();
    m_five_line->Draw();
    
    m_legend->Draw("0");

    // Stamp
    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(12);
    latex->SetTextFont(42);
    //latex->SetTextSize(0.052);
    //latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
    //latex->SetTextSize(0.04);
    //latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 19.0 fb^{-1}");
    latex->SetTextSize(0.032);
    latex->DrawLatex(0.16, 0.97, "CMS Preliminary  #sqrt{s} = 8 TeV, L = 19.0 fb^{-1}");
    latex->SetTextSize(0.045);
    latex->DrawLatex(0.19, 0.90, "Z(#rightarrow b#bar{b}) + H(#rightarrow E_{T}^{miss})");

    // righthand y label
    TLatex * latex_y = new TLatex();
    //latex_y->SetNDC();
    latex_y->SetTextAlign(12);
    latex_y->SetTextFont(42);
    latex_y->SetTextSize(0.04);
    latex_y->SetTextColor(2);
    latex_y->DrawLatex(x_vals[n_points-1]+2.5,m_one_line->GetY1()  ,"1#sigma");
    latex_y->DrawLatex(x_vals[n_points-1]+2.5,m_two_line->GetY1()  ,"2#sigma");
    latex_y->DrawLatex(x_vals[n_points-1]+2.5,m_three_line->GetY1(),"3#sigma");
    latex_y->DrawLatex(x_vals[n_points-1]+2.5,m_four_line->GetY1() ,"4#sigma");
    latex_y->DrawLatex(x_vals[n_points-1]+2.5,m_five_line->GetY1() ,"5#sigma");

    gPad->RedrawAxis();
    gPad->Print("pvalue_ZbbHinv_BDT_20130731.pdf");
    gPad->Print("pvalue_ZbbHinv_BDT_20130731.png");
}
