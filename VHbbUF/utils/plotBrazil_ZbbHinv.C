{
///
/// macro adpated from Michele de Gruttola
///

    gROOT->SetStyle("Plain");
    gROOT->LoadMacro("tdrstyle.C");
    setTDRStyle();
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetLineStyleString(11,"16 16");

    bool unblind = true;
    bool isMJJ = false;
    int n_points = 6;
    double x_vals[n_points] = {105, 115, 125, 135, 145, 150};
    //int n_points = 9;
    //double x_vals[n_points] = {110, 115, 120, 125, 130, 135, 140, 145, 150};
    //int n_points = 6;
    //double x_vals[n_points] = {110, 115, 120, 125, 130, 135};


    // -------------------------------------------------------------------------

    // ZbbHinv (2013-07-31)
    //double y_injected[n_points]=    {  1.8672, 2.0391, 2.2109, 2.3984, 2.5703, 2.7422 };
    //double y_observed[n_points]=    {  1.8698, 2.0352, 2.2061, 2.4033, 2.5697, 2.7441 };
    //double y_down_points2[n_points]={  0.9883, 1.0793, 1.1659, 1.2695, 1.3605, 1.4461 };
    //double y_down_points1[n_points]={  1.3282, 1.4429, 1.5619, 1.7061, 1.8283, 1.9473 };
    //double y_vals[n_points]=        {  1.8672, 2.0391, 2.2109, 2.3984, 2.5703, 2.7422 };
    //double y_up_points1[n_points]=  {  2.6636, 2.8925, 3.1363, 3.4214, 3.6666, 3.9117 };
    //double y_up_points2[n_points]=  {  3.6408, 3.965, 4.2992, 4.707, 5.0443, 5.3816 };

    // ZbbHinv (2013-08-07)
    //double y_injected[n_points]=    {  1.7266, 1.8828, 2.0547, 2.2422, 2.3984, 2.5703 };
    //double y_observed[n_points]=    {  1.7241, 1.8727, 2.0377, 2.2258, 2.3873, 2.5483 };
    //double y_down_points2[n_points]={  0.9037, 0.9851, 1.0713, 1.1698, 1.2519, 1.3322 };
    //double y_down_points1[n_points]={  1.2155, 1.3258, 1.4418, 1.5702, 1.6892, 1.8002 };
    //double y_vals[n_points]=        {  1.7266, 1.875, 2.0391, 2.2266, 2.3828, 2.5547 };
    //double y_up_points1[n_points]=  {  2.463, 2.6896, 2.925, 3.194, 3.4181, 3.6646 };
    //double y_up_points2[n_points]=  {  3.4102, 3.7132, 4.0381, 4.4095, 4.7488, 5.0753 };

    // ZbbHinv (2013-09-06)
    //double y_injected[n_points]=    {  1.6641, 1.8203, 1.9609, 2.1328, 2.2578, 2.4141 };
    //double y_observed[n_points]=    {  1.4109, 1.5401, 1.728, 1.828, 2.0541, 2.2494 };
    //double y_down_points2[n_points]={  0.8497, 0.9317, 0.9974, 1.0795, 1.1452, 1.2191 };
    //double y_down_points1[n_points]={  1.1465, 1.2507, 1.3459, 1.4566, 1.5452, 1.6406 };
    //double y_vals[n_points]=        {  1.6172, 1.7734, 1.8984, 2.0547, 2.1797, 2.3203 };
    //double y_up_points1[n_points]=  {  2.3069, 2.5298, 2.7081, 2.9474, 3.1267, 3.3284 };
    //double y_up_points2[n_points]=  {  3.1738, 3.458, 3.7257, 4.0433, 4.3166, 4.5951 };

    // ZbbHinv (2013-09-06 mT)
    //double y_injected[n_points]=    {  2.1953, 2.3516, 2.5078, 2.7422, 2.9453, 3.1406 };
    //double y_observed[n_points]=    {  1.5053, 1.5684, 1.7188, 1.8894, 2.0629, 2.2055 };
    //double y_down_points2[n_points]={  1.1616, 1.2437, 1.334, 1.4489, 1.5604, 1.6665 };
    //double y_down_points1[n_points]={  1.5633, 1.6782, 1.8, 1.9551, 2.1085, 2.2486 };
    //double y_vals[n_points]=        {  2.2109, 2.3672, 2.5391, 2.7578, 2.9922, 3.1719 };
    //double y_up_points1[n_points]=  {  3.1539, 3.3957, 3.6625, 3.956, 4.2922, 4.5753 };
    //double y_up_points2[n_points]=  {  4.3669, 4.688, 5.0734, 5.4616, 5.9257, 6.3379 };

    // ZbbHinv (2013-09-24)
    //double y_injected[n_points]=    {  1.7266, 1.8984, 2.0703, 2.2578, 2.3984, 2.5703 };
    //double y_observed[n_points]=    {  1.4711, 1.615, 1.823, 1.9387, 2.1964, 2.4127 };
    //double y_down_points2[n_points]={  0.8825, 0.9728, 1.0506, 1.1452, 1.2227, 1.3012 };
    //double y_down_points1[n_points]={  1.1908, 1.3126, 1.4147, 1.5452, 1.6445, 1.7557 };
    //double y_vals[n_points]=        {  1.6797, 1.8516, 1.9922, 2.1797, 2.3359, 2.4766 };
    //double y_up_points1[n_points]=  {  2.4095, 2.6413, 2.8578, 3.1267, 3.3322, 3.5526 };
    //double y_up_points2[n_points]=  {  3.3054, 3.6337, 3.9203, 4.2893, 4.5843, 4.9046 };

    // ZbbHinv (2013-09-24 mT)
    //double y_injected[n_points]=    {  2.2891, 2.4609, 2.6328, 2.9141, 3.1406, 3.3594 };
    //double y_observed[n_points]=    {  1.5679, 1.6444, 1.8137, 2.003, 2.203, 2.3635 };
    //double y_down_points2[n_points]={  1.2109, 1.3012, 1.4079, 1.5392, 1.6665, 1.7814 };
    //double y_down_points1[n_points]={  1.6253, 1.7557, 1.8898, 2.0661, 2.2486, 2.4037 };
    //double y_vals[n_points]=        {  2.3047, 2.4766, 2.6797, 2.9297, 3.1719, 3.3906 };
    //double y_up_points1[n_points]=  {  3.2877, 3.5526, 3.844, 4.2026, 4.5753, 4.8908 };
    //double y_up_points2[n_points]=  {  4.5375, 4.9046, 5.3405, 5.8019, 6.3379, 6.775 };

    // ZbbHinv (2013-11-11) 
    double y_injected[n_points]=    {  1.6484, 1.8047, 1.9453, 2.1016, 2.2266, 2.3828 };
    double y_observed[n_points]=    {  1.4863, 1.6251, 1.8216, 1.9198, 2.1662, 2.3655 };
    double y_down_points2[n_points]={  0.8414, 0.9221, 0.9855, 1.0673, 1.137, 1.2027 };
    double y_down_points1[n_points]={  1.1295, 1.2393, 1.3255, 1.4355, 1.5261, 1.6228 };
    double y_vals[n_points]=        {  1.6016, 1.7422, 1.8828, 2.0391, 2.1641, 2.2891 };
    double y_up_points1[n_points]=  {  2.2846, 2.4991, 2.6858, 2.9087, 3.087, 3.2836 };
    double y_up_points2[n_points]=  {  3.1431, 3.4284, 3.6951, 4.0017, 4.2743, 4.5332 };


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
    
    m_y_band_graph_2sigma->GetYaxis()->SetNdivisions(510);
    m_y_band_graph_2sigma->GetYaxis()->SetTitleSize(0.048);
    m_y_band_graph_2sigma->GetYaxis()->SetTitleOffset(1.25);
    m_y_band_graph_2sigma->GetYaxis()->SetTitle("95% CL limit on #sigma_{ZH} x BR_{inv}/#sigma_{ZH,SM}");
    //m_y_band_graph_2sigma->GetYaxis()->SetTitle("95% CL limit on #sigma_{ZH} x BR_{inv} [pb]");
    //m_y_band_graph_2sigma->GetYaxis()->SetTitle("Asymptotic 95% CL Limit on #sigma/#sigma_{SM}");

    m_y_band_graph_2sigma->SetMaximum(7.0);
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
        //m_legend->AddEntry(m_y_lineSI_graph, "Signal injected", "l");
        m_legend->AddEntry(m_y_line_graph, "Expected", "l");
        m_legend->AddEntry(m_y_band_graph_1sigma, "Expected #pm 1 #sigma", "f");
        m_legend->AddEntry(m_y_band_graph_2sigma, "Expected #pm 2 #sigma", "f");
        
    } else {
        m_legend = new TLegend(0.64,0.80,0.96,0.92);
        m_legend->SetTextFont(42);
        m_legend->SetFillStyle(0);
        //m_legend->SetFillColor(0);
        //m_legend->SetLineColor(kBlack);
        m_legend->SetBorderSize(0);
        //m_legend->SetBorderSize(1);
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
    m_y_line_graph->Draw("L");
    if (unblind) {
        //m_y_lineSI_graph->Draw("L");
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
    latex->DrawLatex(0.16, 0.97, "CMS Preliminary  #sqrt{s} = 8 TeV, L = 18.9 fb^{-1}");
    latex->SetTextSize(0.045);
    latex->DrawLatex(0.19, 0.90, "Z(#rightarrow b#bar{b}) H(#rightarrow inv)");
    if (isMJJ) {
        latex->DrawLatex(0.19, 0.84, "#color[4]{m_{T} analysis}");
    }

    gPad->RedrawAxis();
    TString postfix = "_20131111";
    if (!isMJJ) {
        gPad->Print("Limit_ZbbHinv_BDT"+postfix+".pdf");
        gPad->Print("Limit_ZbbHinv_BDT"+postfix+".png");
    } else {
        gPad->Print("Limit_ZbbHinv_MJJ"+postfix+".pdf");
        gPad->Print("Limit_ZbbHinv_MJJ"+postfix+".png");
    }
}

