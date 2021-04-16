/*
 * 
 * Author: Viktor Baryshnikov
 *
 * 
 */

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
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TRandom2.h>
#include <TClonesArray.h>


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


#include <assert.h>
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <vector>
#include <fstream> // работа с файлами
#include <iomanip> // манипуляторы ввода/вывода

using namespace std;
#endif

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"
#include <macro/mpd/library_estimation.h> //Здесь функции
/*
 * 
 */
//mpddst_1000.root     mpddst_probparticles.root<-вот здесь нет листьев timeEv
//.........................................................................................................................................
int read_mpddst_tof_method_t_ev_new(TString inFile = "mpddst_10000.root", TString outFile = "mpddst_tof_method_t_ev_10000_plug.root", Int_t nStartEvent = 0, Int_t nEvents4Read = 5831)// 5831
{

    //  для гистограммы........................................................................................................................
    auto h1 = new TH2F("h1", "", 300, 0, 3000, 20, 0, 1); //для estimation of efficiency(Xi2)
    auto h1_1 = new TH2F("h1_1", "", 21, 0, 21, 20, 0, 1); //для estimation of efficiency(elements of interval)
    auto h1_2 = new TH2F("h1_2", "", 300, 0, 300, 20, 0, 1); //для estimation of efficiency(elements of event)

    auto h1_mod = new TH2F("h1_mod", "", 300, 0, 3000, 20, 0, 1); //для estimation of efficiency(Xi2) for |Pdg_code|
    auto h1_1_mod = new TH2F("h1_1_mod", "", 21, 0, 21, 10, 0, 1); //для estimation of efficiency(elements of interval) for |Pdg_code|
    auto h1_2_mod = new TH2F("h1_2_mod", "", 300, 0, 300, 20, 0, 1); //для estimation of efficiency(elements of event) for |Pdg_code|

    auto h2 = new TH2F("h2", "", 300, 0, 300, 400, 0, 20); //для Realtime, s (number of elements of tofmatching)
    auto h3 = new TH2F("h3", "", 300, 0, 300, 400, 0, 20); //для CPUtime, s (number of elements of tofmatching)

    auto h4 = new TH2F("h4", "", 300, 0, 300, 200, -0.3, 0.3); //для time_best_ev, ns вычесленного нашим алгоритмом от (number of elements of tofmatching: N_tofmatching)
    auto h5 = new TH2F("h5", "", 300, 0, 300, 200, -1, 1); //для time_ev (детектор T0), ns (number of elements of tofmatching: N_tofmatching)

    auto h6 = new TH2F("h6", "", 300, 0, 300, 100, 0, 100); //для time_best_ev_sigma, ps (number of elements of tofmatching: N_tofmatching)
    auto h7 = new TH2F("h7", "", 100, 0, 1, 100, 0, 100); //для time_ev_sigma_interval, ps (1 / sqrt(number of elements of interval) )
    auto h_sigma_ev = new TH2F("h_sigma_ev", "", 21, 0, 21, 100, 0, 100); //для time_ev_sigma_interval, ps (number of elements of interval)

    auto h8 = new TH2F("h8", "", 21, 0, 21, 20, 0, 1); //для nsigma: eff(number of elements of intetval)
    auto h9 = new TH2F("h9", "", 300, 0, 300, 20, 0, 1); //для nsigma: eff(number of elements of tofmatching: N_tofmatching)

    auto h_pi_nsigma = new TH2F("h_pi_nsigma", "", 300, 0, 3, 300, 0, 3); //если nsigma>3 то отмечаем пион на плоте m^2[GeV^2](momentum, GeV/c)
    auto h_k_nsigma = new TH2F("h_k_nsigma", "", 300, 0, 3, 300, 0, 3); //если nsigma>3 то отмечаем каон на плоте m^2[GeV^2](momentum GeV/c)
    auto h_p_nsigma = new TH2F("h_p_nsigma", "", 300, 0, 3, 300, 0, 3); //если nsigma>3 то отмечаем протон на плоте m^2[GeV^2](momentum GeV/c)
    //.........................................................................................................................................

    // load libries for mpd
    mpdloadlibs();

    // timer
    TStopwatch timer;
    timer.Start();

    // массив таймеров для оценки времязатратности в зависимости от N_tofmatching
    TStopwatch* timers = new TStopwatch[12]; // max(N_tofmatching) / 50 = 600 / 50 = 12

    ////////////////////////////////////////////////////////////////////////////
    // Init Input file
    Printf("\nOpen file %s for read.\n", inFile.Data());

    // Read file
    TFile * RootFileIn = new TFile(inFile.Data(), "READ");

    // check file is exist
    if (!RootFileIn->IsOpen()) {
        printf("File %s does not find\n", inFile.Data());
        return -1;
    }
    
    // date file for t_ev and t_ev_sigma
    ofstream fout("new_mpddst_tof_method_t_ev_1000.txt", ios_base::out | ios_base::trunc);

    // check file is exist
    if (!fout.is_open()) 
    {
       cout << "File new_mpddst_tof_method_t_ev_1000.txt does not find\n" << endl; 
    }

    // Get TTree object from file
    TTree *InTree = (TTree *) RootFileIn->Get("mpdsim");
    //    InTree->Print();

    // Create array for TpcKalmanTrack branch
    TClonesArray *tpckalmantrack = new TClonesArray("MpdTpcKalmanTrack");
    // connect pointer with branch of the tree
    InTree->SetBranchAddress("TpcKalmanTrack", &tpckalmantrack);

    // The same for TofHit branch
    TClonesArray *tofhit = new TClonesArray("MpdTofHit");
    InTree->SetBranchAddress("TOFHit", &tofhit);

    // The same for TofMatching branch
    TClonesArray *tofmatching = new TClonesArray("MpdTofMatchingData");
    InTree->SetBranchAddress("TOFMatching", &tofmatching);
 
    // The same for MCTrack branch 
    TClonesArray *mctrack = new TClonesArray("MpdMCTrack");
    InTree->SetBranchAddress("MCTrack", &mctrack);

    //Get number of events in the tree
    Long64_t nEvents = InTree->GetEntries();

    //
    if (nEvents < (nStartEvent + nEvents4Read)) nEvents4Read = nEvents - nStartEvent;
    if (nEvents4Read <= 0) {
        nEvents4Read = nEvents;
        nStartEvent = 0;
    }
    printf("Total number of events = %lld; read %d events from %d to %d\n", nEvents, nEvents4Read - nStartEvent, nStartEvent, nEvents4Read);

    //
    Double_t timeEvSigma = 0.050; // Uncertainties of timeEv, gaus sigma [ns]

    //
    Double_t timeSigma = 0.060; // Uncertainties of time, gaus sigma [ns]
    
    //
    Double_t speed_light = TMath::C() * pow(10, -7); //29.9792458;  // [cm/nc]


    // гипотезы масс
    /*
    Double_t M_e_ = 0.511/1000; // [GeV] -
    Double_t M_u = 105.6583745/1000; // [GeV] -
    Double_t M_pi = 139.57/1000; // [GeV]
    Double_t M_K = 493.7/1000; // [GeV]
    Double_t M_proton = 938.27/1000; // [GeV]
    Double_t M_deuteron = 1875.612/1000; // [GeV]
    */

    // гипотезы в виде pdgCode
    Int_t pdg_photon = 22;       // photon
    Int_t pdg_anti_neutron = -2112;    // anti-neutron
    Int_t pdg_e = -11;      // e+
    Int_t pdg_anti_Lambda = -3122;    // anti-Lambda
    Int_t pdg_e_ = 11;       // e-
    Int_t pdg_Sigma_ = -3222;    // Sigma- 
    Int_t pdg_e_neutrino = 12;       // e-neutrino (NB: flavour undefined by Geant)
    Int_t pdg_anti_Sigma0 = -3212;    // anti Sigma0 
    Int_t pdg_mu = -13;      // mu+
    Int_t pdg_Sigma_PB = -3112;    // Sigma+ (PB)*/
    Int_t pdg_mu_ = 13;       // mu-
    Int_t pdg_anti_Xi0 = -3322;    // anti Xi0
    Int_t pdg_pi0 = 111;      // pi0
    Int_t pdg_Xi = -3312;    // Xi+
    Int_t pdg_pi = 211;      // pi+
    Int_t pdg_Omega = -3334;    // Omega+ (PB)
    Int_t pdg_pi_ = -211;     // pi-
    Int_t pdg_tau = -15;      // tau+
    Int_t pdg_K_long = 130;      // K long
    Int_t pdg_tau_ = 15;       // tau-
    Int_t pdg_K = 321;      // K+
    Int_t pdg_D = 411;      // D+
    Int_t pdg_K_ = -321;     // K-
    Int_t pdg_D_ = -411;     // D-
    Int_t pdg_neutron = 2112;     // n
    Int_t pdg_D0 = 421;      // D0
    Int_t pdg_proton = 2212;     // p
    Int_t pdg_anti_D0 = -421;     // anti D0
    Int_t pdg_anti_proton = -2212;    // anti-proton
    Int_t pdg_Ds = 431;      // Ds+
    Int_t pdg_K_short = 310;      // K short
    Int_t pdg_anti_Ds_ = -431;     // anti Ds-
    Int_t pdg_eta = 221;      // eta
    Int_t pdg_Lamba_c = 4122;     // Lamba_c+
    Int_t pdg_Lambda = 3122;     // Lambda
    Int_t pdg_W = 24;       // W+
    Int_t pdg_Sigma = 3222;     // Sigma+
    Int_t pdg_W_ = -24;      // W-
    Int_t pdg_Sigma0 = 3212;     // Sigma0
    Int_t pdg_Z = 23;       // Z
    Int_t pdg_Sigma_PB_ = 3112;     // Sigma- (PB)/
    Int_t pdg_deuteron = 0;        // deuteron ????
    Int_t pdg_Xi0 = 3322;     // Xi0
    Int_t pdg_triton = 0;        // triton
    Int_t pdg_Xi_ = 3312;     // Xi-
    Int_t pdg_alpha = 0;        // alpha
    Int_t pdg_Omega_ = 3334;     // Omega- (PB)
    Int_t pdg_K0 = 311; // K0
    Int_t pdg_anti_K0 = -311; // anti-K0
    // Int_t pdg_undefined = 0;        // G nu ? PDG ID 0 is undefined

    // гипотезы масс
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg_e_);
    Double_t M_e_ = particle->Mass();  // particle mass in GeV -

    particle = TDatabasePDG::Instance()->GetParticle(pdg_mu_);
    Double_t M_mu_ = particle->Mass();  // particle mass in GeV - 

    particle = TDatabasePDG::Instance()->GetParticle(pdg_pi);
    Double_t M_pi = particle->Mass();  // particle mass in GeV 

    particle = TDatabasePDG::Instance()->GetParticle(pdg_K); 
    Double_t M_K = particle->Mass();  // particle mass in GeV 

    particle = TDatabasePDG::Instance()->GetParticle(pdg_proton);
    Double_t M_proton = particle->Mass();  // particle mass in GeV 

    //почему то выдает массу дейтрона как 0
    particle = TDatabasePDG::Instance()->GetParticle(pdg_deuteron);
    Double_t M_deuteron = particle->Mass();  // particle mass in GeV 
    M_deuteron  = 1875.612/1000; // [GeV]

    cout << "\n" << endl;

    Double_t M_i[] = { M_pi, M_K, M_proton}; // i = pi, K, proton|     deuteron; , M_deuteron
    Int_t Pdg_i[] = { pdg_pi, pdg_K, pdg_proton}; // i = pi, K, proton|        , deuteron; , pdg_deuteron
    string Name_i[] = { "pi", "K", "proton"}; // , "deuteron"

    Int_t sizeM_i = sizeof(M_i) / sizeof(Double_t);
    Int_t sizePdg_i = sizeof(Pdg_i) / sizeof(Int_t);
    printf("size M_i: %d\n", sizeM_i);
    printf("size Pdg_i: %d\n", sizePdg_i);
    printf("\n");

    for (Int_t ii = 0; ii < sizeM_i; ii++)
    {
    cout << "M_["<<Name_i[ii]<<", "<<ii<<"]: "<< M_i[ii]*1000 <<", MeV; " << "Pdg_["<<Name_i[ii]<<", "<<ii<<"]: "<< Pdg_i[ii] <<";"<< endl;
    }
/*
    //заранее создаем массив сметченных хитов
    MpdTofMatchingData *Matching = new MpdTofMatchingData(); //пустой констурктор

    //заранее создаем массив треков
    MpdTpcKalmanTrack *track = new MpdTpcKalmanTrack(); //пустой констурктор

    //заранее создае массив хитов
    MpdTofHit *Hit = new MpdTofHit(); //пустой констурктор

    //заранее создае массив mctrack
    MpdMCTrack *MCtrack = new MpdMCTrack(); //пустой констурктор
*/ 
    int count_pi = 0;
    int count_k = 0;
    int count_p = 0;
    int count_d = 0;
    int count_anti_pi = 0;
    int count_anti_k = 0;
    int count_anti_p = 0;
    int count_anti_d = 0;

    int c_pi = 0;
    int c_k = 0;
    int c_p = 0;
    int c_d = 0;
    int c_anti_pi = 0;
    int c_anti_k = 0;
    int c_anti_p = 0;
    int c_anti_d = 0;

    int neutral_particles = 0;
    int total_particles = 0;
    int photons = 0;
    int ions = 0;

    vector<TSIGMA> TSIGMA_n;
    //TSIGMA_n.reserve(nEvents4Read - nStartEvent);
    Int_t involved = 0;

    //main loop over events in the tree
    for (Long64_t iEv = nStartEvent; iEv < nEvents4Read; iEv++) {

        //Clear TClonesArray
        tpckalmantrack->Clear();
        tofhit->Clear();
        tofmatching->Clear();
        mctrack->Clear();

        //TRandom2* generator = new TRandom2();
 
        //Get one event from the tree
        InTree->GetEntry(iEv);

        //print event number
        cout << "Event# " << iEv << endl;

        //Кол.во сметченных треков в кажом ивенте
        Int_t old_N_tofmatching = tofmatching->GetEntriesFast();

        // проверка N_tofmatching / 50 (целочисленное) ->индекс таймера.  Запуск таймера[индекс таймера]
        Int_t Index_timers = Int_t(old_N_tofmatching / 50);
        timers[Index_timers].Start();
 
        if (old_N_tofmatching > 0)
        {//НАЧАЛО old_N_tofmatching
        
        Int_t skipping = 0;
        // loop over tofmatching array для того что бы выкинуть неинтересующие/невозможные треки
        for (Int_t i = 0; i < old_N_tofmatching; i++) 
        {
            total_particles = total_particles + 1;
      
            //берем один элемент из массива сметченных хитов
            MpdTofMatchingData *Matching = (MpdTofMatchingData*) tofmatching->UncheckedAt(i);

            //из метченного хита вытаскиваем индех tpc трека (в массиве треков)
            //в хедере tof/MpdTofMatchingData.h находим метод 
            //int			GetKFTrackIndex(void)const{ return fKFTrackIndex;};
            // эта строчка говорит нам что у объекта класса "MpdTofMatchingData" есть метод "GetKFTrackIndex" который возвращает нам целочисленное значение типа int
            // это число есть номер в масиве треков. трек под этим номером использовался для метчинга
            Int_t TrackIndex = Matching->GetKFTrackIndex();

            //из массива треков вытаскиваем трек, индекс которого нашли в сметченном хите
            MpdTpcKalmanTrack *track = (MpdTpcKalmanTrack*) tpckalmantrack->UncheckedAt(TrackIndex);

            //из метченного хита вытаскиваем индекс tof хита (в массиве треков)
            Int_t HitIndex = Matching->GetTofHitIndex();

            //из массива треков вытаскиваем хит, индекс которого нашли в сметченном хите
            MpdTofHit *Hit = (MpdTofHit*) tofhit->UncheckedAt(HitIndex);

            //из хита вытаскиваем RefIndex
            Int_t RefIndex = Hit->GetRefIndex();

            //из массива треков вытаскиваем mctrack, RefIndex которого нашли в хите
            MpdMCTrack *MCtrack = (MpdMCTrack*) mctrack->UncheckedAt(RefIndex);

            //из трека вытаскиваем импульс
            Double_t p = track->Momentum(); //Gev/c
            
            //из трека вытаскиваем информацию о частица/античастица?
            Int_t charge = track->Charge(); //1 - частица, -1 - античастица

            //из MCtrack вытаскиваем MotherId нашего трека 
            Int_t MotherId = MCtrack->GetMotherId();

            //из MCtrack вытаскиваем PdgCode нашего трека 
            Int_t PdgCode = MCtrack->GetPdgCode();

            //берем время начала события
            Double_t timeEv = Hit->GetTimeEv();

            //берем время пролета частицы
            Double_t time = Hit->GetTime(); //[nc]

            //берем длину трэка (всего) 
            Double_t length = Matching->GetTrackLength(); //[cm]

            //Double_t M_prop = p*sqrt(pow(((time - timeEv)*speed_light/length),2) - 1);
            Double_t M_prop_2 = p*p*abs(pow(((time - timeEv)*speed_light/length),2) - 1);

                                                    //                pi0                     K0                      
//D0                   Sigma0                    other                           ion                  photon           
            if( sqrt(M_prop_2) > (2 * M_proton) || abs(PdgCode) == 111 || abs(PdgCode) == 311 || abs(PdgCode) == 2112 || abs(PdgCode) == 421 || abs(PdgCode) == 3212 || abs(PdgCode) == 0 || MotherId != -1 || abs(PdgCode) > 10000 || PdgCode == pdg_photon)
            {

                  if(PdgCode == pdg_photon) photons = photons + 1;
                  if(abs(PdgCode) == 111 || abs(PdgCode) == 311 || abs(PdgCode) == 2112 || abs(PdgCode) == 421 || abs(PdgCode) == 3212)
                  { 
                     neutral_particles = neutral_particles + 1;         
                     //cout << "charge" << charge << endl;
                  }
                  if(abs(PdgCode) > 10000 || sqrt(M_prop_2) > (2 * M_proton) ) ions = ions + 1;

                  skipping = skipping + 1;
                  continue;
            }
           
        }//end loop over tofmatching array для того что бы выкинуть неинтересующие/невозможные треки

        Int_t N_tofmatching = old_N_tofmatching - skipping;

        cout <<"\nold_N_tofmatching = " << old_N_tofmatching << "; N_tofmatching = old_N_tofmatching -  skipping = " << N_tofmatching << "; skipping = " << skipping << ";" << endl;


 
        if (N_tofmatching > 3)
        {//НАЧАЛО N_tofmatching

        cout << "\nN_tofmatching > 3" << endl;

        skipping = 0;

        //массив timeExp_i и timeExpSigma_i для каждого сметченного хита при разных гипотезах масс
        Double_t** timeExp_i =  new Double_t *[N_tofmatching]; // [i]{0, 0, 0, 0, 0, 0};  //[nc]
        Double_t** timeExpSigma_i =  new Double_t *[N_tofmatching]; // [i]{0, 0, 0, 0, 0, 0}; //[nc]

        //массив SigmaPID2_i nSigmaTOF_i для каждого сметченного хита при разных гипотезах масс
        Double_t* SigmaPID2_i =  new Double_t [N_tofmatching];
        Double_t* nSigmaTOF_i =  new Double_t [N_tofmatching];

        //двумерный массив весов для сметченных трэков в выбранном ивенте
        Double_t** Weight_i = new Double_t *[N_tofmatching]; 

        //двумерный массив для сметченных трэков в выбранном ивенте [i][0] = индекс сметченного трэка, [i][1] = p - импульс сметченного трэка, [i]21] = pSigma - сигма импульса сметченного трэка 
        Double_t** Data_n = new Double_t *[N_tofmatching]; 

        //массив упрощенных гипотез (2 гипотезы) масс для каждого сметченного трека в ивенте
        Double_t** Mass_i = new Double_t *[N_tofmatching];

        //массив упрощенных гипотез (2 гипотезы) пдг кодов для упрощенных гипотез (2 гипотезы) масс для каждого сметченного трека в ивенте
        Int_t** Pdg_Code_i = new Int_t *[N_tofmatching];
 
        //время пролета частицы для каждого сметченного трэка
        Double_t* time_n = new Double_t [N_tofmatching];

        Int_t step_matching = 0;

        // loop over tofmatching array
        for (Int_t i = 0; i < old_N_tofmatching; i++) // cout i для проверки выхода за N_tofmatching 
        {

            //берем один элемент из массива сметченных хитов
            MpdTofMatchingData *Matching = (MpdTofMatchingData*) tofmatching->UncheckedAt(i);

            //из метченного хита вытаскиваем индех tpc трека (в массиве треков)
            //в хедере tof/MpdTofMatchingData.h находим метод 
            //int			GetKFTrackIndex(void)const{ return fKFTrackIndex;};
            // эта строчка говорит нам что у объекта класса "MpdTofMatchingData" есть метод "GetKFTrackIndex" который возвращает нам целочисленное значение типа int
            // это число есть номер в масиве треков. трек под этим номером использовался для метчинга
            Int_t TrackIndex = Matching->GetKFTrackIndex();

            //из массива треков вытаскиваем трек, индекс которого нашли в сметченном хите
            MpdTpcKalmanTrack *track = (MpdTpcKalmanTrack*) tpckalmantrack->UncheckedAt(TrackIndex);

            //из метченного хита вытаскиваем индекс tof хита (в массиве треков)
            Int_t HitIndex = Matching->GetTofHitIndex();

            //из массива треков вытаскиваем хит, индекс которого нашли в сметченном хите
            MpdTofHit *Hit = (MpdTofHit*) tofhit->UncheckedAt(HitIndex);

            //из хита вытаскиваем RefIndex
            Int_t RefIndex = Hit->GetRefIndex();

            //из массива треков вытаскиваем mctrack, RefIndex которого нашли в хите
            MpdMCTrack *MCtrack = (MpdMCTrack*) mctrack->UncheckedAt(RefIndex);

            //из трека вытаскиваем импульс
            Double_t p = track->Momentum(); //Gev/c
            
            //из трека вытаскиваем информацию о частица/античастица?
            Int_t charge = track->Charge(); //1 - частица, -1 - античастица

            //из MCtrack вытаскиваем MotherId нашего трека и сохраняем его в массив Data_n
            Int_t MotherId = MCtrack->GetMotherId();

            //из MCtrack вытаскиваем PdgCode нашего трека и сохраняем его в массив Data_n
            Int_t PdgCode = MCtrack->GetPdgCode();

            //берем время начала события
            Double_t timeEv = Hit->GetTimeEv();

            //для time_ev (детектор T0) (number of elements of tofmatching: N_tofmatching)
            h5->Fill(N_tofmatching, timeEv);

            //берем время пролета частицы
            Double_t time = Hit->GetTime(); //[nc]

            //берем длину трэка (всего) 
            Double_t length = Matching->GetTrackLength(); //[cm]
            Double_t lengthSigma = 1; //cm, можно считать что сигма вычисления длины стремится к 0, но нам нужно зазадать какой нибудь минимум что бы корректно вычислять t_exp_sigma

            Double_t M_prop_2 = p*p*abs(pow(((time - timeEv)*speed_light/length),2) - 1);

                                                    //                pi0                     Xi0                      
//D0                   Sigma0                    other                           ion                  photon           
            if( sqrt(M_prop_2) > (2* M_proton) || abs(PdgCode) == 111 || abs(PdgCode) == 311 || abs(PdgCode) == 2112 || abs(PdgCode) == 421 || abs(PdgCode) == 3212 || abs(PdgCode) == 0 || MotherId != -1 || abs(PdgCode) > 10000 || PdgCode == pdg_photon)
            {
                  skipping = skipping + 1;
                  continue;
            }


            time_n[step_matching] = time;          

            //двумерный массив весов для упрощенных гипотез масс сметченных трэков
            Weight_i[step_matching] = new Double_t [2];

            //массив timeExp_i и timeExpSigma_i для каждого сметченного хита для упрощенных гипотезах масс
            timeExp_i[step_matching] = new Double_t [2];
            timeExpSigma_i[step_matching] = new Double_t [2];

            //массив упрощенных гипотез (2 гипотезы) масс и пдг кодов для них для каждого сметченного трека в ивенте
            Mass_i[step_matching] = new Double_t [2];
            Pdg_Code_i[step_matching] = new Int_t [2];


            //досоздаем двумерный массив для сметченных трэков в выбранном ивенте [i][0] = индекс сметченного трэка, [i][1] = p - импульс сметченного трэка и заполняем его
            Data_n[step_matching] = new Double_t [6]; 
            Data_n[step_matching][0] = step_matching;
            Data_n[step_matching][1] = p;
            Data_n[step_matching][2] = PdgCode; // pdgcode
            Data_n[step_matching][3] = 0; // mass_pdg
            Data_n[step_matching][4] = MotherId; // MotherId
            Data_n[step_matching][5] = length;

            //printf("\n tofmatching# %d; PdgCode: %d; \n", step_matching, PdgCode);

            //если все сделали правильно то лезем в TDatabasePDG и пдгкод используем как индекс для того чтобы из этого класса вытащить массу
            if ( TDatabasePDG::Instance() ) 
            {
               particle = TDatabasePDG::Instance()->GetParticle(PdgCode);
 
               if (particle) 
               {             
               Double_t mass_pdg = particle->Mass();  // particle mass in GeV
               Data_n[step_matching][3] = mass_pdg;  // particle mass in GeV
                   
               //printf("\n tofmatching# %d; PdgCode: %d; mass_pdg = %.5f;\n", step_matching, PdgCode, mass_pdg);
               }
               
            }

            //cout << "'Proposed' mass^2 : " << M_prop_2 << ", GeV^2;" << endl;

            //заполянем упрощенный массив гипотез масс
            {
            Int_t k = 1; 
            Int_t st = 0;
            if ( sizePdg_i == sizeM_i )
            {//если sizePdg_i == sizeM_i

            while (k < sizeM_i && st == 0) 
            {
            
                if ( M_prop_2 < M_i[k]*M_i[k] )
                { 
                   Mass_i[step_matching][0] = M_i[k-1];
                   Mass_i[step_matching][1] = M_i[k];
                   if(charge == 1)
                   {
                      Pdg_Code_i[step_matching][0] = Pdg_i[k-1];
                      Pdg_Code_i[step_matching][1] = Pdg_i[k];
                   }
                   else
                   {
                      Pdg_Code_i[step_matching][0] = -Pdg_i[k-1];
                      Pdg_Code_i[step_matching][1] = -Pdg_i[k]; 
                   }
                   
                st = 1;
                }
                    
                k++;
            }
            if (M_prop_2 < M_i[0]*M_i[0])
            {
               Mass_i[step_matching][0] = M_i[0];
               Mass_i[step_matching][1] = M_i[1];

               if(charge == 1)
               {
                   Pdg_Code_i[step_matching][0] = Pdg_i[0];
                   Pdg_Code_i[step_matching][1] = Pdg_i[1];
               }
               else
               {
                   Pdg_Code_i[step_matching][0] = -Pdg_i[0];
                   Pdg_Code_i[step_matching][1] = -Pdg_i[1];
               }
            
            }
            if (M_prop_2 > M_i[sizeM_i - 1]*M_i[sizeM_i - 1])
            {
               Mass_i[step_matching][0] = M_i[sizeM_i - 1];
               Mass_i[step_matching][1] = M_i[sizeM_i - 2];

               if(charge == 1)
               {
                   Pdg_Code_i[step_matching][0] = Pdg_i[sizePdg_i - 1];
                   Pdg_Code_i[step_matching][1] = Pdg_i[sizePdg_i - 2];
               }
               else
               {
                   Pdg_Code_i[step_matching][0] = -Pdg_i[sizePdg_i - 1];
                   Pdg_Code_i[step_matching][1] = -Pdg_i[sizePdg_i - 2];
               }
            }
   
            }//конец если sizePdg_i == sizeM_i
            else
            {
               printf("Error: sizePdg_i != sizeM_i;\n");
               return -1;
            }

            }

            //уже есть заполяненный упрощенный массив гипотез масс и по нему теперь считаем и заполянем массивы для timeExp_i,timeExpSigma_i,Weight_i
            for (Int_t k = 0; k < 2; k++) 
            {

               timeExp_i[step_matching][k] = length*sqrt(p*p + Mass_i[step_matching][k] * Mass_i[step_matching][k])/p/speed_light; // timeExp_i = СУММ(dLk/c/Pk*(Pk^2+Mi^2)^1/2)
               timeExpSigma_i[step_matching][k] = sqrt( pow(Sigma_p(p), 2)*pow( (length/speed_light/sqrt(p*p + Mass_i[step_matching][k]*Mass_i[step_matching][k]) - timeExp_i[step_matching][k]/p) , 2) + pow(lengthSigma, 2)*pow(timeExp_i[step_matching][k]/length ,2) );

                //cout << "\n timeExp_i[" << step_matching << "][" << k << "]" << timeExp_i[step_matching][k] << "; p = " << p << endl;
                //cout << "\n timeExpSigma_i[" << step_matching << "][" << k << "]" << timeExpSigma_i[step_matching][k] << "; time = " << time << endl;


            //заполняем двумерный массив весов для сметченных трэков в выбранном ивенте для выбранного сметченного хита
            Weight_i[step_matching][k] = 1/(timeSigma*timeSigma + timeExpSigma_i[step_matching][k]*timeExpSigma_i[step_matching][k]);

            //cout << "\n Weight_i[" << step_matching << "][" << k << "]" << Weight_i[step_matching][k] << endl;

            //printf("tofmatching# %d; mass proposed: %.3f, MeV; mass^2 proposed: %.6f, GeV^2; time: %.3f; timeEv: %.3f; timeExp_i: %.3f; timeExpSigma_i = %.3f; Weight for prop_mass: %.3f;\n", step_matching, Mass_i[step_matching][k]*1000, Mass_i[step_matching][k]*Mass_i[step_matching][k], time, timeEv, timeExp_i[step_matching][k], timeExpSigma_i[step_matching][k], Weight_i[step_matching][k]);
            }

        step_matching = step_matching + 1;
            
        }// end loop over tofmatching

        //подсчет t_tof_ev: для начала надо взять двумерный массив для сметченных трэков в выбранном ивенте [i][0] = индекс сметченного трэка, [i][1] = p - импульс, ..., и отсортировать по p, поделить отсортированный массив по имупульсам на интервалы и t_tof_ev для 1 интервала = СУММ(t_tof_ev для 2: t_tof_ev для полследнего)/ (последний - 1) и так для остальных
        qsotr_Data_n(Data_n, 0, (N_tofmatching -1 ) ); //работает
/*
        for (Int_t gg = 0; gg < N_tofmatching; gg++)
        {
        cout <<"\n"<< endl;
        cout <<"Array[index][p]: "<<"[index: "<<Data_n[gg][0]<<"]"<<"[p: "<<Data_n[gg][1]<<"], Gev/c"<<";  "<<"[PdgCode: "<<Data_n[gg][2]<<"]"<<";  "<<"[Mass_pdg: " << Data_n[gg][3]<<"], GeV" << ";  " << "[MotherId: " << Data_n[gg][4]<<"]" << ";\n  " << endl;
        }
*/       


        //разбиение сметченых треков по группам, так что бы в каждой было >3 и <=30, интервалы сохраняем в виде индексов от 0 до N_tofmatching в массиве: например 10 треков -> массив[2] = {5, 10}.
        Int_t *intervals_p;
        Int_t N_intervals = 0;
        {
        Int_t h = 20;  
        Int_t status_arr = 0;
        Int_t subtrahend = -1;

        while ( subtrahend < 4 && status_arr == 0 )
        { //начало0

        subtrahend = subtrahend + 1;
        Int_t N_tofmatching_new = N_tofmatching - subtrahend;
        
        if (N_tofmatching_new % 4 == 0 || N_tofmatching_new % 5 == 0 || N_tofmatching_new % 6 == 0 || N_tofmatching_new % 7 == 0 || N_tofmatching_new % 8  == 0 || N_tofmatching_new % 9 == 0 || N_tofmatching_new % 10  == 0 || N_tofmatching_new % 11 == 0 || N_tofmatching_new % 12 == 0 || N_tofmatching_new % 13 == 0 || N_tofmatching_new % 14 == 0 || N_tofmatching_new % 15 == 0 || N_tofmatching_new % 16 == 0 || N_tofmatching_new % 17 == 0 || N_tofmatching_new % (18 - subtrahend) == 0 || N_tofmatching_new % (19 - subtrahend) == 0 || N_tofmatching_new % (20 - subtrahend) == 0/* || N_tofmatching_new % 21 == 0 || N_tofmatching_new % 22 == 0 || N_tofmatching_new % 23 == 0 || N_tofmatching_new % 24 == 0 || N_tofmatching_new % 25 == 0 || N_tofmatching_new % 26 == 0 || N_tofmatching_new % 27 == 0 || N_tofmatching_new % (28 - subtrahend) == 0 || N_tofmatching_new % (29  - subtrahend) == 0 || N_tofmatching_new % (30 - subtrahend) == 0*/)  
        { //начало1
         while (h > (3 + subtrahend) && status_arr == 0)
         { //начало while
          
           Int_t h_new = h - subtrahend;
           if (N_tofmatching_new % h_new == 0)
           {
              //cout << "\nN_intervals: " << N_intervals << ";" << endl;
              N_intervals = N_tofmatching_new/h_new;
              if (N_intervals == 1) 
              {
              //cout << "\nN_intervals = 1 =: " << N_intervals << ";" << endl;
                 if (h_new/2 >= 10)
                 {
                    intervals_p = new Int_t [3];
                    intervals_p[0] = 0;
                    intervals_p[1] = 10;
                    intervals_p[2] = h_new;
                    N_intervals = 2;
                    status_arr = 1;
                    //cout << "\nN_intervals: " << N_intervals << ";" << endl;
                    continue;
                 }
                 if (h_new/2 >= 8)
                 {
                    intervals_p = new Int_t [3];
                    intervals_p[0] = 0;
                    intervals_p[1] = 8;
                    intervals_p[2] = h_new;
                    N_intervals = 2;
                    status_arr = 1;
                    //cout << "\nN_intervals: " << N_intervals << ";" << endl;
                    continue;
                 }
                 if (h_new/2 >= 6)
                 {
                    intervals_p = new Int_t [3];
                    intervals_p[0] = 0;
                    intervals_p[1] = 6;
                    intervals_p[2] = h_new;
                    N_intervals = 2;
                    status_arr = 1;
                    //cout << "\nN_intervals: " << N_intervals << ";" << endl;
                    continue;
                 }
                 if (h_new/2 >= 4)
                 {
                    intervals_p = new Int_t [3];
                    intervals_p[0] = 0;
                    intervals_p[1] = 4;
                    intervals_p[2] = h_new;
                    N_intervals = 2;
                    status_arr = 1;
                    //cout << "\nN_intervals: " << N_intervals << ";" << endl;
                    continue;
                 }
                 if (h_new/2 < 4) 
                 {                    
                    intervals_p = new Int_t [2]; 
                    intervals_p[0] = 0;
                    intervals_p[1] = h_new;
                    N_intervals = 1;
                    status_arr = 1;
                    //cout << "\nN_intervals: " << N_intervals << ";" << endl;
                    continue;
                 }
              }
              else 
              {
                 //cout << "\nN_intervals: " << N_intervals << ";" << endl;
                 intervals_p = new Int_t [N_intervals + 1];
                 intervals_p[0] = 0;
                 for (Int_t j = 0; j < (N_intervals - 1); j++) 
                 {
                    intervals_p[j+1] = (j+1)*h_new;
                 }
                 intervals_p[N_intervals] = N_intervals * h_new + subtrahend;
                 status_arr = 1;
              }
           }
         h = h - 1;  
         } //конец while

        } //конец1
        } //конец0
        }
        
        //cout <<"N_tofmatching: "<< N_tofmatching <<";  "<< endl;
/*
        printf("\n");
        
        for (Int_t gg = 0; gg < N_intervals + 1; gg++)
        {
        cout <<"intervals_p["<< gg << "]: "<< intervals_p[gg]<<";  "<< endl;
        }
        cout <<"\n"<< endl;
*/
        
        cout << "N_intervals: " << N_intervals << ";" << endl;

        //класс в котором будем хранить t_best_Ev, t_best_Ev_sigma, Xi2 и гипотезу для выбраного интверала
        Result* result_n = new Result[N_intervals];
 
        Double_t efficiency_hypotez_event = 0;
        Double_t efficiency_hypotez_event_mod = 0;
        Double_t num_plug_event = 0;
        
        //определение TimeTofEv для каждого из ИМПУЛЬСНЫХ интервалов
        for (Int_t j = 0; j < N_intervals; j++)
        {
           result_n[j] = Result();

           Double_t Xi2_min = 100000000;
           Int_t number_elem_interval = intervals_p[j+1] - intervals_p[j];

           Int_t* hipotez = new Int_t [number_elem_interval]; //внутри цикла будем работать с массивом 0 и 1, а сохранять ифнормация о гипотезе будем в виде? "001010111" с помощью некого преобразования

           memset( hipotez, 0, (sizeof(Int_t)*number_elem_interval) );

           for (Int_t count = 0; count < pow( 2, number_elem_interval ); count++) //перебор всевозможных комбинация масс для выбранного интервала
           {  
              //printf_hipotez (hipotez, number_elem_interval);

              Double_t time_tof_ev = 0;
              Double_t Cymm_weight = 0;
              Double_t time_tof_ev_sigma = 0;
              Double_t Xi2 = 0;
              Int_t hh = 0;

              for (Int_t n_tracks = intervals_p[j]; n_tracks < intervals_p[j+1]; n_tracks++) // 0, 1, 2 , 3, 4;   5, 6, 7, 8, 9
              {                                                                                                 //0, 1, 2, 3, 4
                 //cout << "..." << n_tracks << endl;
                 Int_t id_track = Int_t(Data_n[n_tracks][0]); 

                 time_tof_ev = time_tof_ev + Weight_i[id_track][ hipotez[hh] ] * (time_n[id_track] - timeExp_i[id_track][ hipotez[hh] ]); 
                 Cymm_weight = Cymm_weight + Weight_i[id_track][ hipotez[hh] ];
 
                 hh = hh + 1;
              }

              time_tof_ev = time_tof_ev / Cymm_weight;
              time_tof_ev_sigma = sqrt(1 / Cymm_weight);
   
              hh = 0;

              for (Int_t n_tracks = intervals_p[j]; n_tracks < intervals_p[j+1]; n_tracks++)
              { 
                 Int_t id_track = Int_t(Data_n[n_tracks][0]);
                 
                 Xi2 = Xi2 + pow( (time_n[id_track] - timeExp_i[id_track][ hipotez[hh] ] - time_tof_ev) , 2) * Weight_i[id_track][ hipotez[hh] ]; 

                 hh = hh + 1;
              }
              
              if (Xi2 < Xi2_min) //сохраняем t_best_Ev, t_best_Ev_sigma, Xi2 и гипотезу и number_elem_interval
              {

                 Result Z (time_tof_ev, time_tof_ev_sigma, Xi2, int_hipotez(hipotez, number_elem_interval),  number_elem_interval);
                 result_n[j] = Z;

                 Xi2_min = Xi2; 
              }

              increase_elem_arr (hipotez, number_elem_interval); // каким либо образом приращать гипотезу масс в зависимости от этапа цикла

           } // конец перебора всевозможных комбинация масс для выбранного интервала

           //гистограмма sigma_t_ev по кол.ва треков в интервале
           h_sigma_ev->Fill(number_elem_interval, 1000*result_n[j].time_tof_ev_sigma);

           h7->Fill(1/sqrt(number_elem_interval), 1000*result_n[j].time_tof_ev_sigma);

           // 1 метод оценки эффективности определения сорта частиц
           Double_t efficiency_hypotez_Xi2min = 0;
           Double_t efficiency_hypotez_Xi2min_mod = 0;
           memset( hipotez, 0, (sizeof(Int_t)*number_elem_interval) );

           //функция которая int_hipotez(hipotez, number_elem_interval) переведет обратно в одномерный массив hipotez
           (result_n[j]).int_in_hipotez( hipotez ); 

           //printf_hipotez( hipotez, number_elem_interval); //выведем нашу гипотезу и посмотрим как мы ее преобразовали

           Int_t hh = 0;

           //система заглушек: для нейтарльных частиц[0], для ионов[1], для фотонов[2]
           int plug[3] = {1, 1, 1}; // 000- включаем все заглушки, если != 000 - выключаем все заглушки

           Int_t num_plug = 0; //кол.во сработавших систем заглушек

           for (Int_t n_tracks = intervals_p[j]; n_tracks < intervals_p[j+1]; n_tracks++)
           {
              Int_t Track_Index = Int_t(Data_n[n_tracks][0]);

              Int_t Pdg_Code = Int_t(Data_n[n_tracks][2]);
              Double_t Pdg_Mass = Data_n[n_tracks][3];

              Int_t Hypotez_Pdg = Pdg_Code_i[Track_Index][hipotez[hh]]; 
              Double_t Hypotez_Mass = Mass_i[Track_Index][hipotez[hh]];

              Int_t MotherId = Int_t(Data_n[n_tracks][4]);

              // printf("\nIndex interval %d; Pdg_Code = %d; Hypotez_Pdg = %d; Pdg_Mass = %f, GeV; Hypotez_Mass = %f, GeV; \n", j, Pdg_Code, Hypotez_Pdg, Pdg_Mass, Hypotez_Mass);
 
//                                       \/          \/
              if (plug_motherid(MotherId, 1) == 1) // 0 - включаем заглушку
              {
                 //printf("\nMotherId == -1");
                 if(Pdg_Code == 211) count_pi = count_pi + 1;
                 if(Pdg_Code == 321) count_k = count_k + 1;
                 if(Pdg_Code == 2212) count_p = count_p + 1;
                 if(Pdg_Code == -211) count_anti_pi = count_pi + 1;
                 if(Pdg_Code == -321) count_anti_k = count_k + 1;
                 if(Pdg_Code == -2212) count_anti_p = count_p + 1;

                 //система заглушек: для нейтарльных частиц[0], для ионов[1], для фотонов[2] 
                 if (plug_for_efficiency(Pdg_Code, plug) != 0)
                 {
                    if (Pdg_Code == Hypotez_Pdg  && Pdg_Mass == Hypotez_Mass) 
                    { 
                       if(Pdg_Code == 211) c_pi = c_pi + 1;
                       if(Pdg_Code == 321) c_k = c_k + 1;
                       if(Pdg_Code == 2212) c_p = c_p + 1;
                       if(Pdg_Code == -211) c_anti_pi = c_pi + 1;
                       if(Pdg_Code == -321) c_anti_k = c_k + 1;
                       if(Pdg_Code == -2212) c_anti_p = c_p + 1;
   
                       efficiency_hypotez_Xi2min = efficiency_hypotez_Xi2min + 1;
                       efficiency_hypotez_event = efficiency_hypotez_event + 1;

                       //printf("\nIndex interval %d; Pdg_Code = %d; Hypotez_Pdg = %d; Pdg_Mass = %f, GeV; Hypotez_Mass = %f, GeV; MotherId = %d; \n", j, Pdg_Code, Hypotez_Pdg, Pdg_Mass, Hypotez_Mass, MotherId);

                    }

                    if (abs(Pdg_Code) == abs(Hypotez_Pdg)  && Pdg_Mass == Hypotez_Mass) 
                    { 
                       efficiency_hypotez_Xi2min_mod = efficiency_hypotez_Xi2min_mod + 1;
                       efficiency_hypotez_event_mod = efficiency_hypotez_event_mod + 1;
                    }
                 }
                 else
                 {
                    num_plug = num_plug + 1;
                    num_plug_event = num_plug_event + 1;
                    //printf("num_plug %d; \n", num_plug );
                 } 
              }
              else
              {
                 num_plug = num_plug + 1;
                 num_plug_event = num_plug_event + 1;
                 //printf("num_plug %d; \n", num_plug );
              }

              hh = hh + 1;
           }

           efficiency_hypotez_Xi2min = efficiency_hypotez_Xi2min / (number_elem_interval - num_plug) ;
           efficiency_hypotez_Xi2min_mod = efficiency_hypotez_Xi2min_mod / (number_elem_interval - num_plug) ;

           if( isnan(efficiency_hypotez_Xi2min)== false )   
           {
              // записываем efficiency_hypotez_Xi2min и Xi2min в гистрограмму
              h1->Fill(result_n[j].Xi2, efficiency_hypotez_Xi2min);
              
              // записываем efficiency_hypotez_Xi2min и n_elements в гистрограмму
              h1_1->Fill(result_n[j].number_elem_interval, efficiency_hypotez_Xi2min);
               
              //result_n[j].printf_xi2(j, efficiency_hypotez_Xi2min);
           }

           if( isnan(efficiency_hypotez_Xi2min_mod)== false )   
           {
              // записываем efficiency_hypotez_Xi2min и Xi2min в гистрограмму
              h1_mod->Fill(result_n[j].Xi2, efficiency_hypotez_Xi2min_mod);
              
              // записываем efficiency_hypotez_Xi2min и n_elements в гистрограмму
              h1_1_mod->Fill(result_n[j].number_elem_interval, efficiency_hypotez_Xi2min_mod);
               
              //result_n[j].printf_xi2(j, efficiency_hypotez_Xi2min);
           }

           //удаление объектов 

           delete[] hipotez;

        } // конец определение TimeTofEv для каждого из интервалов

        efficiency_hypotez_event = efficiency_hypotez_event / (N_tofmatching - num_plug_event);
        efficiency_hypotez_event_mod = efficiency_hypotez_event_mod / (N_tofmatching - num_plug_event);

        if( isnan(efficiency_hypotez_event)== false )   
        {   
           // записываем efficiency_hypotez_event и N_tofmatching в гистрограмму
           h1_2->Fill(N_tofmatching, efficiency_hypotez_event);
        }

        if( isnan(efficiency_hypotez_event_mod)== false )   
        {   
           // записываем efficiency_hypotez_event и N_tofmatching в гистрограмму
           h1_2_mod->Fill(N_tofmatching, efficiency_hypotez_event_mod);
        }
 
        //вывод результатов для каждого интервала t_best_Ev, t_best_Ev_sigma, Xi2 и гипотезу, number_elem_interval   
        //printf_result (result_n, N_intervals);

        if (N_intervals > 1)
        {
        // пересчитываем time_tof_ev и time_tof_ev_sigma для каждого интервала | t_TOFev вычисляется для каждого импульсного интервала с использованием только треков, принадлежащих другим импульсным интервалам.
           Double_t* time_tof_ev = new Double_t[N_intervals]; 
           Double_t* time_tof_ev_sigma = new Double_t[N_intervals];
           Double_t* ves = new Double_t[N_intervals];

           for (Int_t j0 = 0; j0 < N_intervals; j0++)
           {
              Double_t cymm = 0;
              Double_t cymm_sigma = 0;
              Double_t weigh = 0;

              for (Int_t j1 = 0; j1 < N_intervals; j1++)
              {
                 Double_t p0 = (1 / result_n[j1].Xi2);
                 if ( j1 != j0)
                 {
                    cymm = cymm + (result_n[j1].time_tof_ev * p0);
                    cymm_sigma = cymm_sigma + (result_n[j1].time_tof_ev_sigma * p0);
                    weigh = weigh + p0;               
                 }
              }
           
              time_tof_ev[j0] = cymm / weigh;
              time_tof_ev_sigma[j0] = cymm_sigma / weigh;
              ves[j0] = sqrt( 1 / weigh);
           }

           // считаем уже конечные time_best_ev и time_best_ev_sigma для каждого ивента
           Double_t time_best_ev = 0; 
           Double_t time_best_ev_sigma = 0;
           Double_t weigh = 0;

           for (Int_t j = 0; j < N_intervals; j++)
           {
              time_best_ev = time_best_ev + time_tof_ev[j]*(1/ves[j]);
              time_best_ev_sigma = time_best_ev_sigma + time_tof_ev_sigma[j]*(1/ves[j]);
              weigh = weigh + (1/ves[j]);
           }

           time_best_ev = time_best_ev / weigh;
           time_best_ev_sigma = time_best_ev_sigma / weigh;

           //удаление объектов 

           delete[] time_tof_ev;
           delete[] time_tof_ev_sigma;
           delete[] ves;
 
           //cout << "\niEv: " << iEv << "; time_best_ev = " << time_best_ev << ", ns; time_best_ev_sigma = " << time_best_ev_sigma << ", ns;\n" << endl;
           fout << "\niEv: " << iEv << "; time_best_ev = " << time_best_ev << ", ns; time_best_ev_sigma = " << time_best_ev_sigma << ", ns;\n" << endl;

           //забиваем в гистрограмму time_best_ev(N_tofmatching)
           if( isnan(time_best_ev)== false)   
           {
              h4->Fill(N_tofmatching, time_best_ev);
           }

           //забиваем в гистрограмму time_best_ev_sigma(N_tofmatching) и  гистрограмму time_best_ev_sigma(1/sqrt(N_tofmatching))
           if( isnan(time_best_ev_sigma)== false)   
           {
              h6->Fill(N_tofmatching, time_best_ev_sigma*1000);
           }
           
           Double_t eff_nsigma = 0; 
           Double_t eff_nsigma_event = 0;

           //определение TimeTofEv для каждого из ИМПУЛЬСНЫХ интервалов
           for (Int_t j = 0; j < N_intervals; j++)
           {
              Int_t number_elem_interval = intervals_p[j+1] - intervals_p[j];
              Int_t* hipotez = new Int_t [number_elem_interval]; 
              memset( hipotez, 0, (sizeof(Int_t)*number_elem_interval) );

              //функция которая int_hipotez(hipotez, number_elem_interval) переведет обратно в одномерный массив hipotez
              (result_n[j]).int_in_hipotez( hipotez );

              Int_t hh = 0;

              for (Int_t n_tracks = intervals_p[j]; n_tracks < intervals_p[j+1]; n_tracks++)
              {
                 Int_t Track_Index = Int_t(Data_n[n_tracks][0]);
                 
                 //массивы SigmaPID2_i nSigmaTOF_i
                 SigmaPID2_i[Track_Index] = timeSigma*timeSigma + pow(time_best_ev_sigma, 2) + timeExpSigma_i[Track_Index][hipotez[hh]]*timeExpSigma_i[Track_Index][hipotez[hh]];
                 nSigmaTOF_i[Track_Index] = (time_n[Track_Index] - time_best_ev - timeExp_i[Track_Index][hipotez[hh]]) / sqrt(SigmaPID2_i[Track_Index]);

                 //cout << "N_intervals =" << j << "; SigmaPID2_i[" << Track_Index << "] = " << SigmaPID2_i[Track_Index] << "; nSigmaTOF_i["<< Track_Index << "] = " << nSigmaTOF_i[Track_Index] << endl;
                 if(abs(nSigmaTOF_i[Track_Index]) < 3)
                 {
                    //cout << "\n nSigmaTOF_i < 3" << endl;
                    eff_nsigma = eff_nsigma + 1;
                    eff_nsigma_event = eff_nsigma_event + 1;

                    Int_t Pdg_Code = Int_t(Data_n[n_tracks][2]);
                    Double_t Pdg_Mass = Data_n[n_tracks][3];

                    Double_t M_propagetion_2 = Data_n[n_tracks][1]*Data_n[n_tracks][1]*( (time_n[Track_Index] - time_best_ev)*(time_n[Track_Index] - time_best_ev)*speed_light*speed_light/Data_n[n_tracks][5]/Data_n[n_tracks][5] - 1);

                    if(M_propagetion_2 > 0)
                    {

                       if(abs(Pdg_Code) == pdg_pi)
                       {
                          h_pi_nsigma->Fill(Data_n[n_tracks][1], M_propagetion_2);
                          cout << "N_intervals =" << j << "; SigmaPID2_i[" << Track_Index << "] = " << SigmaPID2_i[Track_Index] << "; nSigmaTOF_i["<< Track_Index << "] = " << nSigmaTOF_i[Track_Index] << endl;
                          cout << "\n Pi, M_prop =" << sqrt(M_propagetion_2) <<"; M_pi =" << M_pi << endl;
                       }
                       if(abs(Pdg_Code) == pdg_K) 
                       {
                          h_k_nsigma->Fill(Data_n[n_tracks][1], M_propagetion_2);
                          cout << "N_intervals =" << j << "; SigmaPID2_i[" << Track_Index << "] = " << SigmaPID2_i[Track_Index] << "; nSigmaTOF_i["<< Track_Index << "] = " << nSigmaTOF_i[Track_Index] << endl;
                          cout << "\n K, M_prop =" << sqrt(M_propagetion_2) <<"; M_K =" << M_K << endl;
                       }
                       if(abs(Pdg_Code) == pdg_proton) 
                       {
                          h_p_nsigma->Fill(Data_n[n_tracks][1], M_propagetion_2);
                          cout << "N_intervals =" << j << "; SigmaPID2_i[" << Track_Index << "] = " << SigmaPID2_i[Track_Index] << "; nSigmaTOF_i["<< Track_Index << "] = " << nSigmaTOF_i[Track_Index] << endl;
                          cout << "\n proton, M_prop =" << sqrt(M_propagetion_2) <<"; M_proton =" << M_proton << endl;
                       }
                    }
                 }
                 hh = hh + 1;
              }
              
              eff_nsigma = eff_nsigma / number_elem_interval;

              if( isnan(eff_nsigma)== false)   
              {   
                 h8->Fill(number_elem_interval, eff_nsigma);
              }
           }

           eff_nsigma_event = eff_nsigma_event / N_tofmatching;
           
           if( isnan(eff_nsigma_event)== false)   
           {   
              h9->Fill(N_tofmatching, eff_nsigma_event);
           }

        }
              
        //удаление объектов 

        delete[] result_n;
        delete[] intervals_p;


        // останавливаем таймер[индекс таймера]
        timers[Index_timers].Stop();
/*        cout << endl;
        cout << "Real time[" << Index_timers << "]: "<< timers[Index_timers].RealTime() << " s, CPU time[" << Index_timers << "]: "<< timers[Index_timers].CpuTime() << " s" << endl;
*/     
        //забиваем во 2 и 3 гистрограммы timers[Index_timers].RealTime() и timers[Index_timers].CpuTime() в зависимости от N_tofmatching
        h2->Fill(N_tofmatching, timers[Index_timers].RealTime());
        h3->Fill(N_tofmatching, timers[Index_timers].CpuTime());

        //сбрасываем таймер[индекс таймера]
        timers[Index_timers].Reset();
        
    //удаление объектов 

    for (int k = 0; k < N_tofmatching; k++)
    {
       delete[] timeExp_i[k];
       delete[] timeExpSigma_i[k];

       delete[] Weight_i[k];

       delete[] Data_n[k];

       delete[] Mass_i[k];
       delete[] Pdg_Code_i[k];       

    }

    delete[] timeExp_i;
    delete[] timeExpSigma_i;

    delete[] SigmaPID2_i;
    delete[] nSigmaTOF_i;

    delete[] Weight_i;

    delete[] Data_n;

    delete[] Mass_i;
    delete[] Pdg_Code_i; 
 
    delete[] time_n;
            
    
    }//КОНЕЦ N_tofmatching  
    else cout << "number of New tofmatching: " << N_tofmatching << " <= 3; Skipping the entire algorithm" << endl;     
    }//КОНЕЦ old_N_tofmatching    
    else 
    {
    cout << "number of old tofmatching: " << old_N_tofmatching << "; Skipping the entire algorithm" << endl;

    // останавливаем таймер[индекс таймера]
    timers[Index_timers].Stop();

    /*cout << endl;
    cout << "Real time[" << Index_timers << "]: "<< timers[Index_timers].RealTime() << " s, CPU time[" << Index_timers << "]: "<< timers[Index_timers].CpuTime() << " s" << endl;*/
/*
    //забиваем во 2 и 3 гистрограммы timers[Index_timers].RealTime() и timers[Index_timers].CpuTime() в зависимости от N_tofmatching
    h2->Fill(N_tofmatching, timers[Index_timers].RealTime());
    h3->Fill(N_tofmatching, timers[Index_timers].CpuTime());
*/
    //сбрасываем таймер[индекс таймера]
    timers[Index_timers].Reset();

    }
    }// end main loop over events

    //закрываем файл для чтения и записи данных t_ev, t_ev_sigma
    RootFileIn->Close();
    fout.close();

    //сохраняем гистограммы

    // Init Output file
    Printf("\nOpen file %s for writing data.\n", outFile.Data());

    // Read file
    TFile * RootFileOut = new TFile(outFile.Data(), "recreate");

    // check file is exist
    if (!RootFileOut->IsOpen()) {
        printf("File %s does not find\n", outFile.Data());
        return -1;
    }
    h1->Write();
    h1_1->Write();
    h1_2->Write();

    h1_mod->Write();
    h1_1_mod->Write();
    h1_2_mod->Write();

    h2->Write();
    h3->Write();
 
    h4->Write();
    h5->Write();

    h6->Write();
    h7->Write();
    h_sigma_ev->Write();

    h8->Write();
    h9->Write();

    h_pi_nsigma->Write();
    h_k_nsigma->Write();
    h_p_nsigma->Write();

    //закрываем файл для записи
    RootFileOut->Close();

    printf("\ncount_pi = %d/%d; count_k = %d/%d; count_p = %d/%d; count_anti_pi = %d/%d; count_anti_k = %d/%d; count_anti_p = %d/%d;\n", c_pi, count_pi, c_k, count_k, c_p, count_p, c_anti_pi, count_anti_pi, c_anti_k, count_anti_k, c_anti_p, count_anti_p);

    printf("\nTOTAL for all reading events:   (neutral_particles %d/ total_particles %d) = %f procent;  \n", neutral_particles, total_particles, double(100 * neutral_particles / total_particles) );

    printf("\nTOTAL for all reading events:   (photons %d/ total_particles %d) = %f procent;  \n", photons, total_particles, double(100 * photons / total_particles) );

    printf("\nTOTAL for all reading events:   (ions %d/ total_particles %d) = %f procent;  \n", ions, total_particles, double(100 * ions / total_particles) );

    // -----   Finish   -------------------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    cout << endl << endl;
    cout << "Macro finished successfully." << endl; // marker of successful execution
    cout << "Input file is " << inFile << endl;
    cout << "Output file is " << outFile << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
    cout << endl;
    return 1;

}
