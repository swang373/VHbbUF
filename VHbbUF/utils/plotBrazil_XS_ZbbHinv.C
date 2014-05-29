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
    //double y_injected[n_points]=    {  0.1906, 0.1578, 0.1318, 0.1115, 0.0942, 0.0895 };
    //double y_observed[n_points]=    {  0.1908, 0.1575, 0.1315, 0.1117, 0.0942, 0.0896 };
    //double y_down_points2[n_points]={  0.1009, 0.0835, 0.0695, 0.059, 0.0499, 0.0472 };
    //double y_down_points1[n_points]={  0.1356, 0.1116, 0.0931, 0.0793, 0.067, 0.0636 };
    //double y_vals[n_points]=        {  0.1906, 0.1578, 0.1318, 0.1115, 0.0942, 0.0895 };
    //double y_up_points1[n_points]=  {  0.2718, 0.2238, 0.187, 0.159, 0.1344, 0.1277 };
    //double y_up_points2[n_points]=  {  0.3716, 0.3068, 0.2563, 0.2188, 0.1849, 0.1757 };

    // ZbbHinv (2013-08-07)
    //double y_injected[n_points]=    {  0.1762, 0.1457, 0.1225, 0.1042, 0.0879, 0.0839 };
    //double y_observed[n_points]=    {  0.176, 0.1449, 0.1215, 0.1035, 0.0875, 0.0832 };
    //double y_down_points2[n_points]={  0.0922, 0.0762, 0.0639, 0.0544, 0.0459, 0.0435 };
    //double y_down_points1[n_points]={  0.1241, 0.1026, 0.086, 0.073, 0.0619, 0.0588 };
    //double y_vals[n_points]=        {  0.1762, 0.1451, 0.1216, 0.1035, 0.0873, 0.0834 };
    //double y_up_points1[n_points]=  {  0.2514, 0.2081, 0.1744, 0.1485, 0.1253, 0.1196 };
    //double y_up_points2[n_points]=  {  0.348, 0.2873, 0.2407, 0.2049, 0.174, 0.1657 };

    // ZbbHinv (2013-09-06)
    //double y_injected[n_points]=    {  0.1767, 0.1475, 0.1231, 0.1051, 0.0882, 0.0842 };
    //double y_observed[n_points]=    {  0.1498, 0.1248, 0.1085, 0.0901, 0.0802, 0.0785 };
    //double y_down_points2[n_points]={  0.0902, 0.0755, 0.0626, 0.0532, 0.0447, 0.0425 };
    //double y_down_points1[n_points]={  0.1217, 0.1013, 0.0845, 0.0718, 0.0603, 0.0573 };
    //double y_vals[n_points]=        {  0.1717, 0.1437, 0.1192, 0.1012, 0.0851, 0.081 };
    //double y_up_points1[n_points]=  {  0.2449, 0.2049, 0.1701, 0.1452, 0.1221, 0.1162 };
    //double y_up_points2[n_points]=  {  0.337, 0.2801, 0.2339, 0.1992, 0.1686, 0.1604 };

    // ZbbHinv (2013-09-06 mT)
    //double y_injected[n_points]=    {  0.2331, 0.1905, 0.1575, 0.1351, 0.115, 0.1096 };
    //double y_observed[n_points]=    {  0.1598, 0.1271, 0.1079, 0.0931, 0.0806, 0.077 };
    //double y_down_points2[n_points]={  0.1233, 0.1008, 0.0838, 0.0714, 0.0609, 0.0582 };
    //double y_down_points1[n_points]={  0.166, 0.136, 0.113, 0.0963, 0.0823, 0.0785 };
    //double y_vals[n_points]=        {  0.2347, 0.1918, 0.1594, 0.1359, 0.1169, 0.1107 };
    //double y_up_points1[n_points]=  {  0.3349, 0.2751, 0.23, 0.1949, 0.1676, 0.1597 };
    //double y_up_points2[n_points]=  {  0.4636, 0.3798, 0.3186, 0.2691, 0.2314, 0.2212 };

    // ZbbHinv (2013-09-24)
    //double y_injected[n_points]=    {  0.1762, 0.1469, 0.1234, 0.1049, 0.0879, 0.0839 };
    //double y_observed[n_points]=    {  0.1501, 0.125, 0.1087, 0.0901, 0.0805, 0.0788 };
    //double y_down_points2[n_points]={  0.0901, 0.0753, 0.0626, 0.0532, 0.0448, 0.0425 };
    //double y_down_points1[n_points]={  0.1215, 0.1016, 0.0843, 0.0718, 0.0603, 0.0573 };
    //double y_vals[n_points]=        {  0.1714, 0.1433, 0.1188, 0.1013, 0.0856, 0.0808 };
    //double y_up_points1[n_points]=  {  0.2459, 0.2044, 0.1704, 0.1453, 0.1221, 0.116 };
    //double y_up_points2[n_points]=  {  0.3373, 0.2811, 0.2337, 0.1994, 0.168, 0.1601 };
    
    // ZbbHinv (2013-09-24 no XS uncert)
    //double y_injected[n_points]=    {  0.1698, 0.1408, 0.1178, 0.0999, 0.0845, 0.0798 };
    //double y_observed[n_points]=    {  0.1481, 0.1235, 0.1072, 0.0889, 0.0794, 0.0776 };
    //double y_down_points2[n_points]={  0.0892, 0.0745, 0.0619, 0.0529, 0.0444, 0.0423 };
    //double y_down_points1[n_points]={  0.1204, 0.1002, 0.0835, 0.071, 0.0599, 0.0568 };
    //double y_vals[n_points]=        {  0.1698, 0.1408, 0.1178, 0.0999, 0.0845, 0.0798 };
    //double y_up_points1[n_points]=  {  0.2409, 0.1998, 0.1672, 0.1424, 0.1198, 0.1139 };
    //double y_up_points2[n_points]=  {  0.3281, 0.2739, 0.2276, 0.1947, 0.1643, 0.1557 };

    // ZbbHinv (2013-09-24 mT)
    //double y_injected[n_points]=    {  0.2336, 0.1904, 0.157, 0.1354, 0.1151, 0.1097 };
    //double y_observed[n_points]=    {  0.16, 0.1272, 0.1081, 0.0931, 0.0807, 0.0772 };
    //double y_down_points2[n_points]={  0.1236, 0.1007, 0.0839, 0.0715, 0.0611, 0.0582 };
    //double y_down_points1[n_points]={  0.1659, 0.1358, 0.1127, 0.096, 0.0824, 0.0785 };
    //double y_vals[n_points]=        {  0.2352, 0.1916, 0.1598, 0.1362, 0.1163, 0.1107 };
    //double y_up_points1[n_points]=  {  0.3355, 0.2749, 0.2292, 0.1953, 0.1677, 0.1597 };
    //double y_up_points2[n_points]=  {  0.4631, 0.3795, 0.3184, 0.2697, 0.2323, 0.2212 };

    // ZbbHinv (2013-09-24 mT, no XS uncert)
    //double y_injected[n_points]=    {  0.2304, 0.1892, 0.1579, 0.134, 0.1151, 0.1097 };
    //double y_observed[n_points]=    {  0.1582, 0.1258, 0.1069, 0.0921, 0.0798, 0.0763 };
    //double y_down_points2[n_points]={  0.122, 0.1001, 0.083, 0.0709, 0.0605, 0.0576 };
    //double y_down_points1[n_points]={  0.1639, 0.1342, 0.1113, 0.0953, 0.0816, 0.0773 };
    //double y_vals[n_points]=        {  0.2304, 0.1892, 0.1579, 0.134, 0.1151, 0.1097 };
    //double y_up_points1[n_points]=  {  0.3287, 0.2699, 0.2252, 0.1911, 0.1642, 0.1564 };
    //double y_up_points2[n_points]=  {  0.4522, 0.3713, 0.3099, 0.263, 0.2259, 0.2152 };

    // ZbbHinv (2013-11-11)
    //double y_injected[n_points]=    {  0.1667, 0.1399, 0.1163, 0.0989, 0.0833, 0.0791 };
    //double y_observed[n_points]=    {  0.1558, 0.1298, 0.1129, 0.0934, 0.0835, 0.0814 };
    //double y_down_points2[n_points]={  0.0883, 0.074, 0.0615, 0.0522, 0.0439, 0.0415 };
    //double y_down_points1[n_points]={  0.1186, 0.0995, 0.0827, 0.0703, 0.0592, 0.0561 };
    //double y_vals[n_points]=        {  0.1667, 0.1399, 0.1163, 0.0989, 0.0833, 0.0791 };
    //double y_up_points1[n_points]=  {  0.2378, 0.1984, 0.1649, 0.1411, 0.1188, 0.1128 };
    //double y_up_points2[n_points]=  {  0.3251, 0.2702, 0.2261, 0.1929, 0.1624, 0.1542 };

    // ZbbHinv (2013-11-15)
    //double y_injected[n_points]=    {  0.1674, 0.1396, 0.116, 0.0991, 0.0833, 0.0793 };
    //double y_observed[n_points]=    {  0.1557, 0.1301, 0.1132, 0.0936, 0.0837, 0.0817 };
    //double y_down_points2[n_points]={  0.0886, 0.0739, 0.0614, 0.0521, 0.0439, 0.0417 };
    //double y_down_points1[n_points]={  0.1191, 0.0993, 0.0825, 0.0703, 0.0592, 0.0559 };
    //double y_vals[n_points]=        {  0.1674, 0.1396, 0.116, 0.0991, 0.0833, 0.0793 };
    //double y_up_points1[n_points]=  {  0.2375, 0.1981, 0.1654, 0.1406, 0.1189, 0.1125 };
    //double y_up_points2[n_points]=  {  0.3256, 0.2715, 0.2261, 0.1928, 0.1625, 0.1542 };

    // ZbbHinv (2013-11-18)
    //double y_injected[n_points]=    {  0.1682, 0.1396, 0.1169, 0.0991, 0.0833, 0.0793 };
    //double y_observed[n_points]=    {  0.156, 0.1299, 0.113, 0.0936, 0.0836, 0.0815 };
    //double y_down_points2[n_points]={  0.0884, 0.0739, 0.0614, 0.0525, 0.0441, 0.0417 };
    //double y_down_points1[n_points]={  0.1193, 0.0996, 0.0829, 0.0705, 0.0593, 0.0562 };
    //double y_vals[n_points]=        {  0.1682, 0.1396, 0.1169, 0.0991, 0.0833, 0.0793 };
    //double y_up_points1[n_points]=  {  0.2387, 0.1992, 0.1658, 0.1414, 0.1189, 0.1131 };
    //double y_up_points2[n_points]=  {  0.3271, 0.2723, 0.2273, 0.1933, 0.1635, 0.1547 };

    // ZbbHinv (2013-11-18, fixed -- not including 0.1512)
    double y_injected[n_points]=    {  1.1127, 0.9235, 0.7732, 0.6556, 0.5511, 0.5246 };
    double y_observed[n_points]=    {  1.0321, 0.8593, 0.7473, 0.6193, 0.5528, 0.5391 };
    double y_down_points2[n_points]={  0.5846, 0.4888, 0.4062, 0.347, 0.2917, 0.2756 };
    double y_down_points1[n_points]={  0.7888, 0.6586, 0.5482, 0.4664, 0.392, 0.3719 };
    double y_vals[n_points]=        {  1.1127, 0.9235, 0.7732, 0.6556, 0.5511, 0.5246 };
    double y_up_points1[n_points]=  {  1.5784, 1.3173, 1.0968, 0.9353, 0.7861, 0.7483 };
    double y_up_points2[n_points]=  {  2.1636, 1.8007, 1.5035, 1.2784, 1.0815, 1.0229 };

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
    
    m_y_band_graph_2sigma->GetYaxis()->SetNdivisions(506);
    m_y_band_graph_2sigma->GetYaxis()->SetTitleSize(0.048);
    m_y_band_graph_2sigma->GetYaxis()->SetTitleOffset(1.25);
    //m_y_band_graph_2sigma->GetYaxis()->SetTitle("95% CL limit on #sigma_{ZH} x BR_{inv}/#sigma_{ZH,SM}");
    //m_y_band_graph_2sigma->GetYaxis()->SetTitle("95% CL limit on #sigma_{ZH} x BR_{inv} [pb]");
    m_y_band_graph_2sigma->GetYaxis()->SetTitle("95% CL limit on #sigma_{ZH} x BR(H #rightarrow inv) [pb]");
    //m_y_band_graph_2sigma->GetYaxis()->SetTitle("Asymptotic 95% CL Limit on #sigma/#sigma_{SM}");

    m_y_band_graph_2sigma->SetMaximum(3.0);
    m_y_band_graph_2sigma->SetMinimum(0.0);

    // y = 1 line
    m_one_line = new TLine(x_vals[0], 1, x_vals[n_points-1], 1);
    //m_one_line = new TLine(110, 1, 150, 1);
    //m_one_line = new TLine(110, 1, 135, 1);
    m_one_line->SetLineColor(4);
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

    //TCanvas *c1 = new TCanvas();
    TCanvas *c1 = new TCanvas("c1","c1",600,400);
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
    
    //m_one_line->Draw();
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
    //latex->SetTextSize(0.032);
    //latex->DrawLatex(0.16, 0.97, "CMS Preliminary  #sqrt{s} = 8 TeV, L = 18.9 fb^{-1}");
    //latex->SetTextSize(0.045);
    //latex->DrawLatex(0.19, 0.90, "Z(#rightarrow b#bar{b}) H(#rightarrow inv)");
    latex->SetTextSize(0.045);
    latex->DrawLatex(0.19, 0.90, "CMS Preliminary, Z(#rightarrow b#bar{b}) H(#rightarrow inv)");
    latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 18.9 fb^{-1}");
    if (isMJJ) {
        latex->DrawLatex(0.19, 0.84, "#color[4]{m_{T} analysis}");
    }

    gPad->RedrawAxis();
    TString postfix = "_20131118";
    if (!isMJJ) {
        gPad->Print("Limit_XS_ZbbHinv_BDT"+postfix+".pdf");
        gPad->Print("Limit_XS_ZbbHinv_BDT"+postfix+".png");
    } else {
        gPad->Print("Limit_XS_ZbbHinv_MJJ"+postfix+".pdf");
        gPad->Print("Limit_XS_ZbbHinv_MJJ"+postfix+".png");
    }
}

