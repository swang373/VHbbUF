{
///
/// macro adpated from Michele de Gruttola
///

    gROOT->SetStyle("Plain");
    gROOT->LoadMacro("tdrstyle.C");
    setTDRStyle();
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetLineStyleString(11,"16 16");
    gStyle->SetLineStyleString(12,"24 24");
    gStyle->SetLineStyleString(13,"32 32");

    bool unblind = false;
    int n_points = 5;
    double x_vals[n_points] = {105, 115, 125, 135, 145};
    //int n_points = 6;
    //double x_vals[n_points] = {105, 115, 125, 135, 145, 150};
    //int n_points = 9;
    //double x_vals[n_points] = {110, 115, 120, 125, 130, 135, 140, 145, 150};
    //int n_points = 6;
    //double x_vals[n_points] = {110, 115, 120, 125, 130, 135};


    // -------------------------------------------------------------------------

    // ZHinv (2013-08-09)
    double y_injected[n_points]=    {  0, 0, 0, 0, 0 };
    double y_observed[n_points]=    {  0.7209, 0.7683, 0.8483, 0.8974, 0.9449 };
    double y_down_points2[n_points]={  0.3716, 0.3943, 0.4354, 0.463, 0.4836 };
    double y_down_points1[n_points]={  0.5039, 0.5357, 0.5916, 0.6282, 0.657 };
    double y_vals[n_points]=        {  0.7207, 0.7676, 0.8477, 0.8945, 0.9414 };
    double y_up_points1[n_points]=  {  1.0511, 1.1194, 1.2295, 1.2974, 1.3729 };
    double y_up_points2[n_points]=  {  1.4743, 1.5702, 1.7297, 1.8143, 1.9258 };
    double y_vals_old[n_points]=    {  0.7441, 0.793, 0.9258, 0.9883, 1.0586 };
    double y_vals_older[n_points]=  {  1.7266, 1.875, 2.0391, 2.2266, 2.3828 };

    // -------------------------------------------------------------------------

    // Prepare error bars
    double* y_down_bars2 = new double[n_points];
    double* y_down_bars1 = new double[n_points];
    double* y_up_bars1 = new double[n_points];
    double* y_up_bars2 = new double[n_points];

    for (int i=0;i<n_points;++i){
        y_down_bars2[i]=y_vals[i]-y_down_points2[i];
        y_down_bars1[i]=y_vals[i]-y_down_points1[i];   
        y_up_bars2[i]=y_up_points2[i]-y_vals[i];
        y_up_bars1[i]=y_up_points1[i]-y_vals[i];
    }

    // Expected
    m_y_line_graph = new TGraph(n_points, x_vals, y_vals);
    m_y_line_graph->SetLineWidth(2);
    m_y_line_graph->SetLineStyle(11);
    //m_y_line_graph->SetFillColor(kWhite);
    
    // Expected, old
    m_y_lineOld_graph = new TGraph(n_points, x_vals, y_vals_old);
    m_y_lineOld_graph->SetLineWidth(2);
    m_y_lineOld_graph->SetLineStyle(12);
    //m_y_lineOld_graph->SetFillColor(kWhite);
    m_y_lineOld_graph->SetLineColor(kBlue-1);
    
    // Expected, older
    m_y_lineOlder_graph = new TGraph(n_points, x_vals, y_vals_older);
    m_y_lineOlder_graph->SetLineWidth(2);
    m_y_lineOlder_graph->SetLineStyle(13);
    //m_y_lineOlder_graph->SetFillColor(kWhite);
    m_y_lineOlder_graph->SetLineColor(kBlue-2);

    // Observed
    m_y_lineObs_graph = new TGraph(n_points, x_vals, y_observed);
    m_y_lineObs_graph->SetLineWidth(1);
    //m_y_lineObs_graph->SetLineColor(kBlack);
    //m_y_lineObs_graph->SetFillColor(kWhite);
    m_y_lineObs_graph->SetMarkerStyle(20);
    m_y_lineObs_graph->SetMarkerSize(1.1);

    // Signal injected
    m_y_lineSI_graph = new TGraph(n_points, x_vals, y_injected);
    m_y_lineSI_graph->SetLineWidth(2);
    m_y_lineSI_graph->SetLineStyle(11);
    m_y_lineSI_graph->SetLineColor(kRed);
    //m_y_lineSI_graph->SetFillColor(kWhite);

    // y band 1 sigma
    m_y_band_graph_1sigma = new TGraphAsymmErrors(n_points, x_vals, y_vals, 0, 0, y_down_bars1, y_up_bars1);
    m_y_band_graph_1sigma->SetFillColor(kGreen);
    //m_y_band_graph_1sigma->SetLineColor(kGreen);
    //m_y_band_graph_1sigma->SetMarkerColor(kGreen);

    // y band 2 sigma
    m_y_band_graph_2sigma = new TGraphAsymmErrors(n_points, x_vals, y_vals, 0, 0, y_down_bars2, y_up_bars2);
    m_y_band_graph_2sigma->SetFillColor(kYellow);
    //m_y_band_graph_2sigma->SetLineColor(kYellow);
    //m_y_band_graph_2sigma->SetMarkerColor(kYellow);

    // evil axes
    m_y_band_graph_2sigma->GetXaxis()->SetNdivisions(505);
    m_y_band_graph_2sigma->GetXaxis()->SetTitleSize(0.05);
    m_y_band_graph_2sigma->GetXaxis()->SetTitleOffset(1.0);
    m_y_band_graph_2sigma->GetXaxis()->SetTitle("m_{H} [GeV]");
    //m_y_band_graph_2sigma->GetXaxis()->SetLimits(103, 152);
    m_y_band_graph_2sigma->GetXaxis()->SetLimits(x_vals[0]-2, x_vals[n_points-1]+2);
    
    m_y_band_graph_2sigma->GetYaxis()->SetNdivisions(505);
    m_y_band_graph_2sigma->GetYaxis()->SetTitleSize(0.048);
    m_y_band_graph_2sigma->GetYaxis()->SetTitleOffset(1.25);
    m_y_band_graph_2sigma->GetYaxis()->SetTitle("95% CL limit on #sigma_{ZH} x BR_{inv}/#sigma_{ZH,SM}");
    //m_y_band_graph_2sigma->GetYaxis()->SetTitle("95% CL limit on #sigma_{ZH} x BR_{inv} [pb]");
    //m_y_band_graph_2sigma->GetYaxis()->SetTitle("Asymptotic 95% CL Limit on #sigma/#sigma_{SM}");

    m_y_band_graph_2sigma->SetMaximum(3.5);
    m_y_band_graph_2sigma->SetMinimum(0.0);

    // y = 1 line
    m_one_line = new TLine(x_vals[0], 1, x_vals[n_points-1], 1);
    //m_one_line = new TLine(110, 1, 150, 1);
    //m_one_line = new TLine(110, 1, 135, 1);
    m_one_line->SetLineColor(2);
    m_one_line->SetLineWidth(2);

    // Legend
    TLegend* m_legend = 0;
    if (unblind) {
        m_legend = new TLegend(0.64,0.72,0.96,0.92);
        m_legend->SetTextFont(42);
        m_legend->SetFillStyle(0);
        //m_legend->SetFillColor(0);
        //m_legend->SetLineColor(kBlack);
        m_legend->SetBorderSize(0);
        //m_legend->SetBorderSize(1);
        m_legend->AddEntry(m_y_lineObs_graph, "Observed", "lp");
        m_legend->AddEntry(m_y_lineSI_graph, "Signal injected", "l");
        m_legend->AddEntry(m_y_line_graph, "Expected", "l");
        m_legend->AddEntry(m_y_band_graph_1sigma, "Expected #pm 1 #sigma", "f");
        m_legend->AddEntry(m_y_band_graph_2sigma, "Expected #pm 2 #sigma", "f");
        
    } else {
        m_legend = new TLegend(0.64,0.72,0.925,0.92);
        m_legend->SetTextFont(42);
        m_legend->SetFillStyle(0);
        //m_legend->SetFillColor(0);
        //m_legend->SetLineColor(kBlack);
        m_legend->SetBorderSize(0);
        //m_legend->SetBorderSize(1);
        m_legend->AddEntry(m_y_lineOlder_graph, "Expected from Z(bb)", "l");
        m_legend->AddEntry(m_y_lineOld_graph, "Expected from Z(ll)", "l");
        m_legend->AddEntry(m_y_line_graph, "Expected", "l");
        m_legend->AddEntry(m_y_band_graph_1sigma, "Expected #pm 1 #sigma", "f");
        m_legend->AddEntry(m_y_band_graph_2sigma, "Expected #pm 2 #sigma", "f");
    }


    delete[] y_down_bars2;
    delete[] y_down_bars1;
    delete[] y_up_bars2;
    delete[] y_up_bars1;


    // -------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas();
    //TCanvas *c1 = new TCanvas("c1","c1",300,300,700,500);
    //c1->SetGridx();
    //c1->SetGridy();

    // Draw
    m_y_band_graph_2sigma->Draw("A3");
    m_y_band_graph_1sigma->Draw("3");
    m_y_lineOlder_graph->Draw("L");
    m_y_lineOld_graph->Draw("L");
    m_y_line_graph->Draw("L");
    if (unblind) {
        m_y_lineSI_graph->Draw("L");
        m_y_lineObs_graph->Draw("LP");
    }
    
    m_one_line->Draw();
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
    latex->DrawLatex(0.16, 0.97, "CMS Preliminary  #sqrt{s} = 8 TeV, L #leq 19.6 fb^{-1}");
    latex->SetTextSize(0.045);
    latex->DrawLatex(0.19, 0.90, "Z(#rightarrow l^{+} l^{-}, b#bar{b}) H(#rightarrow E_{T}^{miss})");

    gPad->RedrawAxis();
    gPad->Print("Limit_ZHinv_BDT_20130809.pdf");
    gPad->Print("Limit_ZHinv_BDT_20130809.png");
}
