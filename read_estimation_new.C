#include <Rtypes.h>
#if !defined(__CINT__) && !defined(__CLING__)

// ROOT includes
#include "TString.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TChain.h"
#include <TMath.h>
#include "TObject.h"
#include "TParticle.h"
#include "TVector3.h"

#ifndef ROOT_TParticlePDG
 #include "TParticlePDG.h"
#endif
#ifndef ROOT_TDatabasePDG
 #include "TDatabasePDG.h"
#endif

// Fair includes
#include "FairRunAna.h"
#include "FairFileSource.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTask.h"
#include "FairField.h"
#include "FairTimeStamp.h"              // for FairTimeStamp
#include "FairHit.h"            

// MPD includes
#include "MpdTpcHitProducer.h"
#include "MpdTpcClusterFinderTask.h"
#include "MpdTpcDigitizerAZ.h"
#include "MpdTpcClusterFinderAZ.h"
#include "MpdTpcClusterFinderMlem.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanTrack.h"
#include "MpdVertexZfinder.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdKfPrimaryVertexFinder.h"
#include "MpdTofHitProducer.h"
#include "MpdTofPoint.h"
#include "MpdEtofHitProducer.h"
#include "MpdEctTrackFinderTpc.h"
#include "MpdEctTrackFinderTof.h"
#include "MpdEctTrackFinderCpc.h"
#include "MpdTofMatching.h"
#include "MpdZdcDigiProducer.h"
#include "MpdEtofMatching.h"
#include "MpdFillDstTask.h"
#include "MpdGetNumEvents.h"
#include "MpdEmcHitCreation.h"
#include "MpdTofMatching.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdTofHit.h"
#include "MpdKalmanTrack.h"

#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

using namespace std;
#endif

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

/*
 * 
 */

int read_estimation_new(TString inFile = "mpddst_tof_method_t_ev_10000_plug.root") 
{

    //.........................................................................................................................................

    // load libries for mpd
    mpdloadlibs();

    // timer
    TStopwatch timer;
    timer.Start();
    
    ////////////////////////////////////////////////////////////////////////////
    // Init Input file
    Printf("Open file %s\n", inFile.Data());

    // Read file
    TFile * RootFileIn = new TFile(inFile.Data(), "READ");

    // check file is exist
    if (!RootFileIn->IsOpen()) {
        printf("File %s does not find\n", inFile.Data());
        return -1;
    }

    RootFileIn->ls();

    TH2F* h1 = (TH2F*)RootFileIn->Get("h1");//для estimation of efficiency(Xi2)
    TH2F* h1_1 = (TH2F*)RootFileIn->Get("h1_1");//для estimation of efficiency(elements of interval)
    TH2F* h1_2 = (TH2F*)RootFileIn->Get("h1_2");//для estimation of efficiency(elements of event)

    TH2F* h1_mod = (TH2F*)RootFileIn->Get("h1_mod");//для estimation of efficiency(Xi2) for |Pdg_code|
    TH2F* h1_1_mod = (TH2F*)RootFileIn->Get("h1_1_mod");//для estimation of efficiency(elements of interval) for |Pdg_code|
    TH2F* h1_2_mod = (TH2F*)RootFileIn->Get("h1_2_mod");//для estimation of efficiency(elements of event) for |Pdg_code|

    TH2F* h2 = (TH2F*)RootFileIn->Get("h2"); //для Realtime, s (number of elements of tofmatching)
    TH2F* h3 = (TH2F*)RootFileIn->Get("h3"); //для CPUtime, s (number of elements of tofmatching)

    TH2F* h4 = (TH2F*)RootFileIn->Get("h4"); //для time_best_ev, ns вычесленного нашим алгоритмом от (number of elements of tofmatching: N_tofmatching)
    TH2F* h5 = (TH2F*)RootFileIn->Get("h5"); //для time_ev (детектор T0), ns (number of elements of tofmatching: N_tofmatching)

    TH2F* h6 = (TH2F*)RootFileIn->Get("h6");//для time_best_ev_sigma, ps (number of elements of tofmatching: N_tofmatching)
    TH2F* h7 = (TH2F*)RootFileIn->Get("h7");//для time_ev_sigma_interval, ps (1 / sqrt(number of elements of interval) )
    TH2F* h_sigma_ev = (TH2F*)RootFileIn->Get("h_sigma_ev");//для time_ev_sigma_interval, ns (number of elements of interval)

    TH2F* h8 = (TH2F*)RootFileIn->Get("h8");//для nsigma: eff(number of elements of intetval)
    TH2F* h9 = (TH2F*)RootFileIn->Get("h9");//для nsigma: eff(number of elements of tofmatching: N_tofmatching)

    TH2F* h_pi_nsigma = (TH2F*)RootFileIn->Get("h_pi_nsigma");//если nsigma<3 то отмечаем пион на плоте m^2[GeV^2](momentum, GeV/c)
    TH2F* h_k_nsigma = (TH2F*)RootFileIn->Get("h_k_nsigma");//если nsigma<3 то отмечаем каон на плоте m^2[GeV^2](momentum GeV/c)
    TH2F* h_p_nsigma = (TH2F*)RootFileIn->Get("h_p_nsigma");//если nsigma<3 то отмечаем протон на плоте m^2[GeV^2](momentum GeV/c)


//...................................   estimation of efficiency of particle determination
    auto Canvas = new TCanvas("estimation of efficiency of particle determination", "estimation of efficiency of particle determination");
    Canvas->Divide(1,3);
    
    Canvas->cd(1);
    gStyle->SetGridColor(kGray);
    Canvas->cd(1)->SetGrid();
    h1->Draw("col");
    h1->GetYaxis()->SetTitle("efficiency");
    h1->GetXaxis()->SetTitle("Xi2");
    h1->GetYaxis()->SetAxisColor(12);
    h1->GetXaxis()->SetAxisColor(12);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetXaxis()->SetLabelSize(0.05);
    h1->SetTitle("efficiency(Xi2)");
    Canvas->cd(1)->Update();

    Canvas->cd(2);
    gStyle->SetGridColor(kGray);
    Canvas->cd(2)->SetGrid();
    h1_1->SetMarkerColor(4); //Blue
    h1_1->SetMarkerStyle(8); //круг
    h1_1->SetMarkerSize(0.3);
    h1_1->Draw();
    h1_1->GetYaxis()->SetTitle("efficiency");
    h1_1->GetXaxis()->SetTitle("elements of interval");
    h1_1->GetYaxis()->SetAxisColor(12);
    h1_1->GetXaxis()->SetAxisColor(12);
    h1_1->GetYaxis()->SetLabelSize(0.05);
    h1_1->GetXaxis()->SetLabelSize(0.05);
    h1_1->SetTitle("efficiency(elements of interval)");
    Canvas->cd(2)->Update();

    Canvas->cd(3);
    gStyle->SetGridColor(kGray);
    Canvas->cd(3)->SetGrid();
    h1_2->Draw("col");
    h1_2->GetYaxis()->SetTitle("efficiency");
    h1_2->GetXaxis()->SetTitle("elements of event");
    h1_2->GetYaxis()->SetAxisColor(12);
    h1_2->GetXaxis()->SetAxisColor(12);
    h1_2->GetYaxis()->SetLabelSize(0.05);
    h1_2->GetXaxis()->SetLabelSize(0.05);
    h1_2->SetTitle("efficiency(elements of event)");
    Canvas->cd(3)->Update();

//...................................   estimation of efficiency for condition |Pdg_code|
    auto Canvas_mod = new TCanvas("estimation of efficiency of particle determination for condition |Pdg_code|", "estimation of efficiency of particle determination for condition |Pdg_code|");
    Canvas_mod->Divide(1,3);

    Canvas_mod->cd(1);
    gStyle->SetGridColor(kGray);
    Canvas_mod->cd(1)->SetGrid();
    h1_mod->Draw("col");
    h1_mod->GetYaxis()->SetTitle("efficiency");
    h1_mod->GetXaxis()->SetTitle("Xi2");
    h1_mod->GetYaxis()->SetAxisColor(12);
    h1_mod->GetXaxis()->SetAxisColor(12);
    h1_mod->GetYaxis()->SetLabelSize(0.05);
    h1_mod->GetXaxis()->SetLabelSize(0.05);
    h1_mod->SetTitle("efficiency(Xi2)");
    Canvas_mod->cd(1)->Update();

    Canvas_mod->cd(2);
    gStyle->SetGridColor(kGray);
    Canvas_mod->cd(2)->SetGrid();
    h1_1_mod->SetMarkerColor(4); //Blue
    h1_1_mod->SetMarkerStyle(8); //круг
    h1_1_mod->SetMarkerSize(0.3);
    h1_1_mod->Draw();
    h1_1_mod->GetYaxis()->SetTitle("efficiency");
    h1_1_mod->GetXaxis()->SetTitle("elements of interval");
    h1_1_mod->GetYaxis()->SetAxisColor(12);
    h1_1_mod->GetXaxis()->SetAxisColor(12);
    h1_1_mod->GetYaxis()->SetLabelSize(0.05);
    h1_1_mod->GetXaxis()->SetLabelSize(0.05);
    h1_1_mod->SetTitle("efficiency(elements of interval)");
    Canvas_mod->cd(2)->Update();

    Canvas_mod->cd(3);
    gStyle->SetGridColor(kGray);
    Canvas_mod->cd(3)->SetGrid();
    h1_2_mod->Draw("col");
    h1_2_mod->GetYaxis()->SetTitle("efficiency");
    h1_2_mod->GetXaxis()->SetTitle("elements of event");
    h1_2_mod->GetYaxis()->SetAxisColor(12);
    h1_2_mod->GetXaxis()->SetAxisColor(12);
    h1_2_mod->GetYaxis()->SetLabelSize(0.05);
    h1_2_mod->GetXaxis()->SetLabelSize(0.05);
    h1_2_mod->SetTitle("efficiency(elements of event)");
    Canvas_mod->cd(3)->Update();

//...................................   resource
    auto Canvas_resource = new TCanvas("resource", "resource");
    Canvas_resource->Divide(1,2);

    h2->SetMarkerColor(4); //Blue
    h2->SetMarkerStyle(8); //круг
    h2->SetMarkerSize(0.5);

    h3->SetMarkerColor(4); //Blue
    h3->SetMarkerStyle(8); //круг
    h3->SetMarkerSize(0.5);

    Canvas_resource->cd(1);
    gStyle->SetGridColor(kGray);
    Canvas_resource->cd(1)->SetGrid();
    h2->Draw();
    h2->GetYaxis()->SetTitle("Realtime, s");
    h2->GetXaxis()->SetTitle("number of elements of tofmatching");
    h2->GetYaxis()->SetAxisColor(12);
    h2->GetXaxis()->SetAxisColor(12);
    h2->GetYaxis()->SetLabelSize(0.05);
    h2->GetXaxis()->SetLabelSize(0.05);
    h2->SetTitle("Realtime(number of elements of tofmatching)");
    Canvas_resource->cd(1)->Update();
  
    Canvas_resource->cd(2);
    gStyle->SetGridColor(kGray);
    Canvas_resource->cd(2)->SetGrid();
    h3->Draw();
    h3->GetYaxis()->SetTitle("CPUtime, s");
    h3->GetXaxis()->SetTitle("number of elements of tofmatching");
    h3->GetYaxis()->SetAxisColor(12);
    h3->GetXaxis()->SetAxisColor(12);
    h3->GetYaxis()->SetLabelSize(0.05);
    h3->GetXaxis()->SetLabelSize(0.05);
    h3->SetTitle("CPUtime(number of elements of tofmatching)");
    Canvas_resource->cd(2)->Update();

//...................................   time_ev
    auto Canvas_time = new TCanvas("time_ev", "time_ev");
    Canvas_time->Divide(1,2);

    Canvas_time->cd(1);
    gStyle->SetGridColor(kGray);
    Canvas_time->cd(1)->SetGrid();
    h4->Draw("col");
    h4->GetYaxis()->SetTitle("time_TOF_ev, ns");
    h4->GetXaxis()->SetTitle("number of elements of tofmatching");
    h4->GetYaxis()->SetAxisColor(12);
    h4->GetXaxis()->SetAxisColor(12);
    h4->GetYaxis()->SetLabelSize(0.05);
    h4->GetXaxis()->SetLabelSize(0.05);
    h4->SetTitle("time_TOF_ev(number of elements of tofmatching)");
    Canvas_time->cd(1)->Update();
 
    Canvas_time->cd(2);
    gStyle->SetGridColor(kGray);
    Canvas_time->cd(2)->SetGrid();
    h5->Draw("col");
    h5->GetYaxis()->SetTitle("time_T0_ev, ns");
    h5->GetXaxis()->SetTitle("number of elements of tofmatching");
    h5->GetYaxis()->SetAxisColor(12);
    h5->GetXaxis()->SetAxisColor(12);
    h5->GetYaxis()->SetLabelSize(0.05);
    h5->GetXaxis()->SetLabelSize(0.05);
    h5->SetTitle("time_T0_ev(number of elements of tofmatching)");
    Canvas_time->cd(2)->Update();

//...................................   time_TOF_ev_sigma
    auto Canvas_sigma = new TCanvas("time_TOF_ev_sigma", "time_TOF_ev_sigma");
    Canvas_sigma->Divide(1,3);

    Canvas_sigma->cd(1);
    gStyle->SetGridColor(kGray);
    Canvas_sigma->cd(1)->SetGrid();
    h6->Draw("col");
    h6->GetYaxis()->SetTitle("time_TOF_ev_sigma, ps");
    h6->GetXaxis()->SetTitle("number of elements of tofmatching");
    h6->GetYaxis()->SetAxisColor(12);
    h6->GetXaxis()->SetAxisColor(12);
    h6->GetYaxis()->SetLabelSize(0.05);
    h6->GetXaxis()->SetLabelSize(0.05);
    h6->SetTitle("time_TOF_ev_sigma(number of elements of tofmatching)");
    Canvas_sigma->cd(1)->Update();

    Canvas_sigma->cd(2);
    gStyle->SetGridColor(kGray);
    Canvas_sigma->cd(2)->SetGrid();
    h_sigma_ev->SetMarkerColor(4); //Blue
    h_sigma_ev->SetMarkerStyle(8); //круг
    h_sigma_ev->SetMarkerSize(0.3);
    h_sigma_ev->Draw();
    h_sigma_ev->GetYaxis()->SetTitle("time_TOF_ev_sigma, ps");
    h_sigma_ev->GetXaxis()->SetTitle("number of elements of interval");
    h_sigma_ev->GetYaxis()->SetAxisColor(12);
    h_sigma_ev->GetXaxis()->SetAxisColor(12);
    h_sigma_ev->GetYaxis()->SetLabelSize(0.05);
    h_sigma_ev->GetXaxis()->SetLabelSize(0.05);
    h_sigma_ev->SetTitle("time_TOF_ev_sigma(number of elements of interval)");
    Canvas_sigma->cd(2)->Update();

    Canvas_sigma->cd(3);
    gStyle->SetGridColor(kGray);
    Canvas_sigma->cd(3)->SetGrid();
    h7->SetMarkerColor(4); //Blue
    h7->SetMarkerStyle(8); //круг
    h7->SetMarkerSize(0.3);
    h7->Draw();
    h7->GetYaxis()->SetTitle("time_TOF_ev_sigma, ps");
    h7->GetXaxis()->SetTitle("1 / sqrt(number of elements of interval)");
    h7->GetYaxis()->SetAxisColor(12);
    h7->GetXaxis()->SetAxisColor(12);
    h7->GetYaxis()->SetLabelSize(0.05);
    h7->GetXaxis()->SetLabelSize(0.05);
    h7->SetTitle("time_TOF_ev_sigma(1 / sqrt(number of elements of interval) )");
    Canvas_sigma->cd(3)->Update();

//...................................   estimation of efficiency for nsigma
    auto Canvas_nsigma = new TCanvas("estimation of efficiency of particle determination for nsigma", "estimation of efficiency for nsigma");
    Canvas_nsigma->Divide(1,2);

    Canvas_nsigma->cd(1);
    gStyle->SetGridColor(kGray);
    Canvas_nsigma->cd(1)->SetGrid();
    h8->SetMarkerColor(4); //Blue
    h8->SetMarkerStyle(8); //круг
    h8->SetMarkerSize(0.3);
    h8->Draw();
    h8->GetYaxis()->SetTitle("efficiency");
    h8->GetXaxis()->SetTitle("number of elements of intetval");
    h8->GetYaxis()->SetAxisColor(12);
    h8->GetXaxis()->SetAxisColor(12);
    h8->GetYaxis()->SetLabelSize(0.05);
    h8->GetXaxis()->SetLabelSize(0.05);
    h8->SetTitle("estimation of efficiency for nsigma(number of elements of intetval)");
    Canvas_nsigma->cd(1)->Update();
    
    Canvas_nsigma->cd(2); 
    gStyle->SetGridColor(kGray);
    Canvas_nsigma->cd(2)->SetGrid();
    h9->Draw("col");
    h9->GetYaxis()->SetTitle("efficiency");
    h9->GetXaxis()->SetTitle("number of elements of tofmatching");
    h9->GetYaxis()->SetAxisColor(12);
    h9->GetXaxis()->SetAxisColor(12);
    h9->GetYaxis()->SetLabelSize(0.05);
    h9->GetXaxis()->SetLabelSize(0.05);
    h9->SetTitle("estimation of efficiency for nsigma(number of elements of tofmatching)");
    Canvas_nsigma->cd(2)->Update();

//...................................   nsigma < 3: m^2[GeV^2](momentum, GeV/c) for Pi, K, proton
    auto Canvas_m2_p = new TCanvas("m^2(momentum) for Pi, K, p", "m^2(momentum) for Pi, K, p");
    Canvas_m2_p->Divide(1,1);

    h_pi_nsigma->SetMarkerColor(2); //Red
    h_pi_nsigma->SetMarkerStyle(7); //точка -пион
    h_pi_nsigma->SetMarkerSize(0.2);

    h_k_nsigma->SetMarkerColorAlpha(3, 0.35); //Green
    h_k_nsigma->SetMarkerStyle(5); //крест - каон
    h_k_nsigma->SetMarkerSize(0.2);

    h_p_nsigma->SetMarkerColorAlpha(4, 0.35); //Blue
    h_p_nsigma->SetMarkerStyle(2); //точка - протон
    h_p_nsigma->SetMarkerSize(0.2);
    
    Canvas_m2_p->cd(1); 
    gStyle->SetGridColor(kGray);
    Canvas_m2_p->SetGrid();
    h_pi_nsigma->Draw(); 
    h_pi_nsigma->GetYaxis()->SetTitle("m, GeV");
    h_pi_nsigma->GetXaxis()->SetTitle("momentum, GeV/c");
    h_pi_nsigma->GetYaxis()->SetAxisColor(12);
    h_pi_nsigma->GetXaxis()->SetAxisColor(12);
    h_pi_nsigma->GetYaxis()->SetLabelSize(0.05);
    h_pi_nsigma->GetXaxis()->SetLabelSize(0.05);
    h_pi_nsigma->SetTitle("Pi- red, K- green, proton- blue");
    h_k_nsigma->Draw("same"); // "same"
    h_p_nsigma->Draw("same");
    Canvas_m2_p->cd(1)->Update();

/*
    Canvas_m2_p->cd(1); 
    gStyle->SetGridColor(kGray);
    Canvas_m2_p->SetGrid();
    h_pi_nsigma->Draw("col");
    h_pi_nsigma->GetYaxis()->SetTitle("m^2, GeV^2");
    h_pi_nsigma->GetXaxis()->SetTitle("momentum, GeV/c");
    h_pi_nsigma->SetTitle("Pi");
    Canvas_m2_p->cd(1)->Update();
     
    Canvas_m2_p->cd(2); 
    gStyle->SetGridColor(kGray);
    Canvas_m2_p->SetGrid();
    h_k_nsigma->Draw("col");
    h_k_nsigma->GetYaxis()->SetTitle("m^2, GeV^2");
    h_k_nsigma->GetXaxis()->SetTitle("momentum, GeV/c");
    h_k_nsigma->SetTitle("K"); 
    Canvas_m2_p->cd(2)->Update();

    Canvas_m2_p->cd(3);
    gStyle->SetGridColor(kGray);
    Canvas_m2_p->SetGrid();
    h_p_nsigma->Draw("col");
    h_p_nsigma->GetYaxis()->SetTitle("m^2, GeV^2");
    h_p_nsigma->GetXaxis()->SetTitle("momentum, GeV/c");
    h_p_nsigma->SetTitle("proton");  
    Canvas_m2_p->cd(3)->Update();
*/
    // -----   Finish   -------------------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    cout << endl << endl;
    cout << "Macro finished successfully." << endl; // marker of successful execution
    cout << "Input file is " << inFile << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
    cout << endl;
    return 1;

}
