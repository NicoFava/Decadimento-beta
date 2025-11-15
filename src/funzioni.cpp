#include "funzioni.h"

bool IsInPhysicalRegion(double Ee, double Ev) {
    if ( Ee < M_ELECTRON ) 
        return false;

    double k = sqrt(Ee*Ee - M_ELECTRON*M_ELECTRON);
    double k_prime = Ev;

    if (k <= 0 || k_prime <= 0)
        return false;
    
    // Energia cinetica del protone
    double Tp = M_NEUTRON - Ee - Ev - M_PROTON;
    
    if (Tp < 0)
        return false;

    double Ep = Tp + M_PROTON;
    
    // Calcolo del coseno dell'angolo tra elettrone e neutrino
    double cos_theta_ev = (Ep*Ep - M_PROTON*M_PROTON - k*k - k_prime*k_prime) / (2*k*k_prime);

    return (cos_theta_ev >= -1.0 && cos_theta_ev <= 1.0);
}

// Genera una direzione casuale nello spazio
TVector3 GenerateRandomDirection(TRandom3* gRandom) {
    double theta = acos(gRandom->Uniform(-1.0, 1.0));  // Angolo polare
    double phi = gRandom->Uniform(0.0, 2*M_PI);         // Angolo azimutale
    
    TVector3 dir;
    dir.SetMagThetaPhi(1.0, theta, phi);
    return dir;
}

// Genera un evento di decadimento completo
Event GenerateDecayEvent(TRandom3* gRandom, double Ee, double Ev) {
    Event event;
    
    // Calcoli cinematici nel piano del decadimento
    double Te = Ee - M_ELECTRON;  // Energia cinetica dell'elettrone
    double ke = sqrt(Ee*Ee - M_ELECTRON*M_ELECTRON);  // Momento dell'elettrone
    double kv = Ev;  // Momento del neutrino (approssimato come energia)
    
    // Energia cinetica del protone
    double Tp = M_NEUTRON - Ee - Ev - M_PROTON;
    double Ep = Tp + M_PROTON;
    double kp = sqrt(Ep*Ep - M_PROTON*M_PROTON);  // Momento del protone
    
    // Calcolo dell'angolo tra elettrone e neutrino
    double cos_theta_ev = (Ep*Ep - M_PROTON*M_PROTON - ke*ke - kv*kv) / (2*ke*kv);
    cos_theta_ev = max(-1.0, min(1.0, cos_theta_ev));  // Clamp tra -1 e 1
    double theta_ev = acos(cos_theta_ev);
    
    // Nel piano del decadimento (asse z)
    // Elettrone lungo l'asse x
    event.p_electron.SetVect(TVector3(ke, 0, 0));
    event.p_electron.SetE(Ee);
    
    // Neutrino ad angolo theta_ev
    double kv_x = kv * cos(theta_ev);
    double kv_y = kv * sin(theta_ev);
    event.p_neutrino.SetVect(TVector3(kv_x, kv_y, 0));
    event.p_neutrino.SetE(Ev);
    
    // Protone per conservazione del momento
    TVector3 p_proton_vec = -(event.p_electron.Vect() + event.p_neutrino.Vect());
    event.p_proton.SetVect(p_proton_vec);
    event.p_proton.SetE(Ep);
    
    // Rotazione casuale nel piano xy
    double phi_plane = gRandom->Uniform(0.0, 2*M_PI);
    event.p_electron.RotateZ(phi_plane);
    event.p_neutrino.RotateZ(phi_plane);
    event.p_proton.RotateZ(phi_plane);
    
    // Genera una direzione casuale nello spazio
    TVector3 random_direction = GenerateRandomDirection(gRandom);
    
    // Ruota i vettori in modo che giacciano nel piano perpendicolare alla direzione
    // Creiamo una base ortonormale con random_direction come normale
    TVector3 axis1 = random_direction.Orthogonal();
    axis1 = axis1.Unit();
    TVector3 axis2 = random_direction.Cross(axis1);
    axis2 = axis2.Unit();
    
    // Applichiamo la rotazione a tutti i vettori
    auto rotate_to_plane = [&](TLorentzVector& vec) {
        TVector3 v = vec.Vect();
        double x = v.Dot(axis1);
        double y = v.Dot(axis2);
        double new_mag = sqrt(x*x + y*y);
        
        if (new_mag > 1e-10) {
            TVector3 new_v = (x/new_mag)*axis1 + (y/new_mag)*axis2;
            new_v = new_v * new_mag;
            vec.SetVect(new_v);
        }
    };
    
    rotate_to_plane(event.p_electron);
    rotate_to_plane(event.p_neutrino);
    rotate_to_plane(event.p_proton);
    
    return event;
}

void run_dalitz() {
    TRandom3 *gRandom = new TRandom3(0);
    const double E_max_box = M_NEUTRON - M_PROTON;
    
    TH2F *h_dalitz = new TH2F("h_dalitz", "Dalitz Plot (Ee vs Ev);E_e (MeV);E_nu (MeV)",
                                100, M_ELECTRON, E_max_box,
                                100, 0, E_max_box);

    // Istogrammi per i 4-momenti
    TH1F *h_pe = new TH1F("h_pe", "Momento elettrone;p (MeV);Counts", 50, 0, 1.5);
    TH1F *h_pv = new TH1F("h_pv", "Momento neutrino;p (MeV);Counts", 50, 0, 1.5);
    TH1F *h_pp = new TH1F("h_pp", "Momento protone;p (MeV);Counts", 50, 0, 1.5);
    
    TH1F *h_angle_ev = new TH1F("h_angle_ev", "Angolo e-nu;#theta (rad);Counts", 50, 0, M_PI);

    long N_TRIALS = 5000000;
    long accepted_events = 0;

    cout << "Generating Dalitz plot with " << N_TRIALS << " trials..." << endl;
    cout << "Generating kinematic 4-vectors..." << endl;

    for (long i = 0; i < N_TRIALS; ++i) {
        double Ee_trial = gRandom->Uniform(M_ELECTRON, E_max_box);
        double Ev_trial = gRandom->Uniform(0.0, E_max_box);

        if (IsInPhysicalRegion(Ee_trial, Ev_trial)) {
            h_dalitz->Fill(Ee_trial, Ev_trial);
            
            // Genera l'evento completo
            Event event = GenerateDecayEvent(gRandom, Ee_trial, Ev_trial);
            
            // Fill degli istogrammi dei momenti
            h_pe->Fill(event.p_electron.Vect().Mag());
            h_pv->Fill(event.p_neutrino.Vect().Mag());
            h_pp->Fill(event.p_proton.Vect().Mag());
            
            // Angolo tra elettrone e neutrino
            double cos_angle = (event.p_electron.Vect().Unit()).Dot(event.p_neutrino.Vect().Unit());
            cos_angle = max(-1.0, min(1.0, cos_angle));
            h_angle_ev->Fill(acos(cos_angle));
            
            accepted_events++;
        }
    }
    
    cout << "Dalitz plot generation completed." << endl;
    cout << "Accepted events: " << accepted_events << endl;
    cout << "Efficiency: " << (double)accepted_events / N_TRIALS * 100 << "%" << endl;

    TCanvas *c1 = new TCanvas("c1", "Dalitz Plot Decadimento Beta", 800, 600);
    h_dalitz->SetStats(0);
    h_dalitz->Draw("COLZ");
    c1->SaveAs("dalitz_plot.png");
    
    TCanvas *c2 = new TCanvas("c2", "Momenti e Angoli", 1200, 800);
    c2->Divide(2, 2);
    
    c2->cd(1);
    h_pe->Draw();
    h_pe->SetLineColor(kBlue);
    
    c2->cd(2);
    h_pv->Draw();
    h_pv->SetLineColor(kRed);
    
    c2->cd(3);
    h_pp->Draw();
    h_pp->SetLineColor(kGreen);
    
    c2->cd(4);
    h_angle_ev->Draw();
    h_angle_ev->SetLineColor(kMagenta);
    
    c2->SaveAs("momenti_e_angoli.png");
    
    delete gRandom;
    delete h_dalitz;
    delete h_pe;
    delete h_pv;
    delete h_pp;
    delete h_angle_ev;
}