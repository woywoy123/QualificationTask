#ifndef EVALUATE_H
#define EVALUATE_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/Statistics.h>
#include<numeric>

void CompareToTruth(TString dir);
void CompileCout(std::vector<TString> JetEnergy, std::vector<TString> Algo_Strings, std::vector<std::map<TString, std::map<TString, std::vector<float>>>> Collector); 
void MultiTrackTruthComparison(TString dir); 




static void CompileHeading(std::vector<TString> Algo_Strings, TString sep, int margin, TString* out)
{
  CoutText(&(*out), margin-1, " "); 
  for (int i(0); i < Algo_Strings.size(); i++)
  {
    *out += (sep); 
    *out += Algo_Strings[i]; 
    CoutText(&(*out), margin - Algo_Strings[i].Sizeof(), " "); 
  }
  *out += (sep); 
} 

static void InjectString(std::vector<TString> Algo_Strings, TString sep, int margin, TString* out, TString inject)
{
  CoutText(&(*out), margin-1, " "); 
  for (int i(0); i < Algo_Strings.size(); i++)
  {
    *out += (sep); 
    *out += inject; 
    CoutText(&(*out), margin - inject.Sizeof(), " "); 
  }
  *out += (sep); 
}

static void Spacer(TString* out, TString inj, int margin)
{
  *out += inj; 
  CoutText(&(*out), margin - inj.Sizeof(), " "); 
  *out += (" | ");
}

static void CompileStatusCodeTable(MMMVF Errors, VS* Output)
{
  std::set<TString> Algos;
  
  // S: Status, ntru: ntrk_ntru, s: Status Code, count: Count of status code 
  MMII S_ntrk_Alg_s_count; 

  // S: Status, LE: Layer energy, Alg: Algorithm, ntrkS: ntrk status
  MMMSI S_LE_Alg_ntrk_s; 

  for (MMMVFi x = Errors.begin(); x != Errors.end(); x++)
  {
    MMVF Er = Errors[x -> first]; 
    for (MMVFi p = Er.begin(); p != Er.end(); p++)
    {
      if ((p -> first).Contains("Truth")){ continue; }

      if (!(Algos.find(p -> first) != Algos.end())){ Algos.insert(p -> first); }
      
      MVF entry = p -> second; 
      for (MVFi m = entry.begin(); m != entry.end(); m++)
      {
        if (!(m -> first).Contains("status")){ continue; }
        TString fit_n = (m -> first)(0,6);

        if (S_ntrk_Alg_s_count[fit_n][p -> first][(m -> second)[0]] == 0){S_ntrk_Alg_s_count[fit_n][p -> first][(m -> second)[0]] = 1;}
        else{S_ntrk_Alg_s_count[fit_n][p -> first][(m -> second)[0]]++;}


        S_LE_Alg_ntrk_s[x -> first][p -> first][fit_n] = (m -> second)[0]; 
      }
    } 
  }

  int margin = 30; 
  TString out; 
  
  VS Algos_V(Algos.begin(), Algos.end()); 
  CompileHeading(Algos_V, " | ", margin, &out); 
  (*Output).push_back(out); 
  out.Clear(); 
  
  VS meaning(Algos_V.size(), "status (#status)"); 
  CompileHeading(meaning, " | ", margin, &out); 
  (*Output).push_back(out); 
  out.Clear(); 

  for (MMIIi x = S_ntrk_Alg_s_count.begin(); x != S_ntrk_Alg_s_count.end(); x++)
  {
    Spacer(&out, x -> first, margin); 
    MII Alg_s = x -> second;
    
    for (TString i : Algos_V)
    {
      TString temp; 
      for (IIi s = Alg_s[i].begin(); s != Alg_s[i].end(); s++)
      {
        temp += (s -> first); 
        temp += (" ("); 
        temp += (s -> second); 
        temp += (") "); 
      }
      Spacer(&out, temp, margin);
    }

    (*Output).push_back(out); 
    out.Clear();
  }
  (*Output).push_back(""); 
 
  VS meaning2(Algos_V.size(), "ntrk [status]"); 
  CompileHeading(meaning2, " | ", margin, &out); 
  (*Output).push_back(out); 
  out.Clear(); 

  for (MMMSIi x = S_LE_Alg_ntrk_s.begin(); x != S_LE_Alg_ntrk_s.end(); x++)
  {
    Spacer(&out, x -> first, margin); 
    MMSI LE_s = x -> second;
    for (TString i : Algos_V)
    {
      TString temp; 
      MI p = LE_s[i];
      for (MIi t  = p.begin(); t != p.end(); t++)
      {
        temp += (t -> first)(5,6); 
        temp += (" ["); 
        temp += (t -> second); 
        temp += ("] "); 
      }
      Spacer(&out, temp, margin);
    }
    
    (*Output).push_back(out); 
    out.Clear();
  }
}


static void CompileTruthTrackEvaluation(MMVT Results, VS* Output)
{
  auto Compare =[&] (VT ntruP, VT ntruT)
  {
    VF Errors(ntruT.size(), 0); 
    float igrl = 0; 
    for (TH1F* H : ntruT){ igrl += H -> Integral(); }
    for (int i(0); i < ntruP.size(); i++){Errors[i] = ((ntruP[i] -> Integral()) / igrl) * ErrorByIntegral(ntruP[i], ntruT[i]);}
    return Errors; 
  };
  
  std::set<TString> Algos; 
  MMMVF ErrorAlgorithmMap; 
  for (MMVi i = Results.begin(); i != Results.end(); i++)
  {
    TString JEnergy = i -> first; 
    MVT Alg_ntrkntru = i -> second; 
    
    // Collect the algorithms in the container 
    for (MVi a = Alg_ntrkntru.begin(); a != Alg_ntrkntru.end(); a++)
    { 
      TString alg = (a -> first)(0, (a -> first).Sizeof()-8); 
      if (!(Algos.find(alg) != Algos.end()) && !alg.Contains("ntrk_")){Algos.insert(alg);}
    }

    // Loop over the algorithms and evaluate the fits to truth
    VT ntrk_1_T = Results[JEnergy]["ntrk_1_Truth"]; 
    VT ntrk_2_T = Results[JEnergy]["ntrk_2_Truth"]; 
    VT ntrk_3_T = Results[JEnergy]["ntrk_3_Truth"]; 
    VT ntrk_4_T = Results[JEnergy]["ntrk_4_Truth"]; 
    VVT All_Truth = {ntrk_1_T, ntrk_2_T, ntrk_3_T, ntrk_4_T};     
    
    for (TString a : Algos)
    {
      VT ntrk_1_Algo = Results[JEnergy][a + "_ntrk_1"];
      VT ntrk_2_Algo = Results[JEnergy][a + "_ntrk_2"];
      VT ntrk_3_Algo = Results[JEnergy][a + "_ntrk_3"];
      VT ntrk_4_Algo = Results[JEnergy][a + "_ntrk_4"];
      VF trk1_E = Compare(ntrk_1_Algo, ntrk_1_T); 
      VF trk2_E = Compare(ntrk_2_Algo, ntrk_2_T); 
      VF trk3_E = Compare(ntrk_3_Algo, ntrk_3_T); 
      VF trk4_E = Compare(ntrk_4_Algo, ntrk_4_T);
      
      ErrorAlgorithmMap[JEnergy]["ntrk_1"][a] = trk1_E; 
      ErrorAlgorithmMap[JEnergy]["ntrk_2"][a] = trk2_E; 
      ErrorAlgorithmMap[JEnergy]["ntrk_3"][a] = trk3_E; 
      ErrorAlgorithmMap[JEnergy]["ntrk_4"][a] = trk4_E; 
    }
  }
  
  TString out; 
  int margin = 25; 
  VS Algos_V(Algos.begin(), Algos.end()); 
  CompileHeading(Algos_V, " | ", margin, &out);
  (*Output).push_back(out); 
  out.Clear(); 

  InjectString(Algos_V, " | ", margin, &out, "Bin by Bin Err (%)"); 
  (*Output).push_back(out); 
  out.Clear(); 

  MMSI BestMap; 
  for (TString i : Algos_V)
  {
    BestMap["ntrk_1"][i] = 0; 
    BestMap["ntrk_2"][i] = 0; 
    BestMap["ntrk_3"][i] = 0; 
    BestMap["ntrk_4"][i] = 0; 
  }

  for (MMMVFi i = ErrorAlgorithmMap.begin(); i != ErrorAlgorithmMap.end(); i++)
  {

    out.Clear(); 
    MMVF alg_trk = i -> second;   
    bool skip = false; 
    VS temp_V; 
    for (MMVFi nt = alg_trk.begin(); nt != alg_trk.end(); nt++)
    {
      
      TString temp; 
      float mini_er = 0; 
      TString mini_alg_err; 
      MF alg_err; 
      for (TString alg : Algos_V)
      {
        VF r = (nt->second)[alg];
        float f = 0; 
        for (float k : r){ f += k;} 
        if (f == 0 || std::isnan(f)){ alg_err[alg] = -1; continue; }
        if (mini_er == 0){ mini_er = f; mini_alg_err = alg; }
        if (f < mini_er) { mini_er = f; mini_alg_err = alg; }
        alg_err[alg] = f; 
      }
      
      Spacer(&temp, "--> " + nt-> first, margin); 
      int all_false = 0; 
      for (TString alg : Algos_V)
      {
        if (alg == mini_alg_err){ Spacer(&temp, "(+1) " + PrecisionString(alg_err[alg] * 100, 3, false), margin); all_false++;}
        else 
        {
          if (alg_err[alg] < 0){ Spacer(&temp, "Failed", margin); continue;}
          Spacer(&temp, PrecisionString(alg_err[alg] * 100, 3,false), margin);
        }
      }
      if (all_false == 0){ continue;}


      BestMap[nt -> first][mini_alg_err]++; 
      temp_V.push_back(temp); 
    }
    
    if (skip){ continue; }
    Spacer(&out, (i -> first), margin); 
    (*Output).push_back(out); 
    for (TString t : temp_V){(*Output).push_back(t);} 
  }
  (*Output).push_back("");   
  (*Output).push_back(""); 
  (*Output).push_back(""); 

  out.Clear(); 
  CompileHeading(Algos_V, " | ", margin, &out);
  (*Output).push_back(out); 
  out.Clear(); 

  for (MMSIi x = BestMap.begin(); x != BestMap.end(); x++)
  {
    out.Clear(); 
    Spacer(&out, (x -> first), margin); 
    MI r = x -> second; 
    for (TString i : Algos_V)
    {
      TString sc; sc += (r[i]); 
      Spacer(&out, sc, margin); 
    }

    (*Output).push_back(out); 
  }
}

static void CompileFLostTable(MMVT Results, MMMVF Errors, VS* Output, TString FL)
{
  
  MMVF Fl; 
  std::set<TString> Algos;
  for (MMVi i = Results.begin(); i != Results.end(); i++)
  {
    MMVF Er = Errors[i -> first];
    std::cout << "----" << FL << " " << i -> first << std::endl;
    for (MMVFi p = Er.begin(); p != Er.end(); p++)
    {
      TString Alg = p -> first;
      std::cout << Alg << std::endl;
      VT ntrk1 = (i -> second)[Alg + "_ntrk_1"]; 
      VT ntrk2 = (i -> second)[Alg + "_ntrk_2"]; 
      VT ntrk3 = (i -> second)[Alg + "_ntrk_3"]; 
      VT ntrk4 = (i -> second)[Alg + "_ntrk_4"]; 

      VVF Errors_Vec; 
      if (!(Algos.find(Alg) !=  Algos.end())){ Algos.insert( Alg ); }

      for (MVFi g = (p -> second).begin(); g != (p -> second).end(); g++)
      {if ((g -> first).Contains("Normalization_Error")){ Errors_Vec.push_back(g -> second); }}

      if (Alg.Contains("Truth"))
      {
        ntrk1 = (i -> second)["ntrk_1_" + Alg]; 
        ntrk2 = (i -> second)["ntrk_2_" + Alg]; 
        ntrk3 = (i -> second)["ntrk_3_" + Alg]; 
        ntrk4 = (i -> second)["ntrk_4_" + Alg]; 

        TH1F* FL2_H = SumHists(ntrk2, "FL2_H"); 
        TH1F* FL3_H = SumHists(ntrk3, "FL3_H"); 
        if (FL.Contains("FL2") && FL2_H -> Integral() < 20000 ){ continue; }
        if (FL.Contains("FL3") && FL3_H -> Integral() < 20000 ){ continue; }
      } 
      VVT ntrks = {ntrk1, ntrk2, ntrk3, ntrk4};  
     
      VF FLost_Vec; 
      if ( FL == "FL2" ){ FLost_Vec = Flost2(ntrks, Errors_Vec); }
      if ( FL == "FL3" ){ FLost_Vec = Flost3(ntrks, Errors_Vec); }
        

      Fl[i -> first][Alg] = FLost_Vec; 
    }
  }
  TString out; 
  int margin = 35; 
  
  VS Algos_V; 
  for (TString t : Algos){ if (t.Contains("Truth")){continue;} Algos_V.push_back(t); }
  
  CompileHeading(Algos_V, " | ", margin, &out); 
  Spacer(&out, "Truth", 6);
  (*Output).push_back(out);  
  out.Clear(); 
  
  VS meaning(Algos_V.size(), FL + "   ((p-t)/t)%" + "   [err/p]%"); 
  CompileHeading(meaning, " | ", margin, &out); 
  (*Output).push_back(out); 
  out.Clear(); 

  MI FLost_pred;
  for (TString i : Algos_V){ FLost_pred[i] = 0; }
  
  for (MMVFi i = Fl.begin(); i != Fl.end(); i++)
  {
    Spacer(&out, i -> first, margin); 
    VF t = Fl[i -> first]["Truth"];  
    if (t.size() == 0){out.Clear(); continue;} 

    float best = -1; 
    TString best_alg; 
    for (TString a : Algos_V)
    {
      VF x = Fl[i -> first][a]; 
      float dif = std::abs(((x[0] - t[0]) / t[0])*100); 
      if (best == -1 && !std::isnan(dif)){ best = dif; best_alg = a; }
      if (std::abs(dif) < best && !std::isnan(dif)){ best = dif; best_alg = a; }
    }
    FLost_pred[best_alg]++; 

    for (TString a : Algos_V)
    {
      VF x = Fl[i -> first][a]; 
      TString g = PrecisionString(x[0], 3, false); 
      float dif = ((x[0] - t[0]) / t[0])*100; 
      TString dif_s = PrecisionString(dif, 3, false); 
      
      float er_p = (x[1]/x[0]); 
      TString er_s = PrecisionString(er_p, 3, true); 

      if (std::isnan(x[0])){ dif_s = ""; g = "Failed";  er_s = ""; }
      if (std::isnan(x[1])){ er_s = "x"; }
     
      TString inj;
      if (a == best_alg){inj = g + "  (" + dif_s + ")  [" + er_s + "] (+1)";}
      else {inj = g + "  (" + dif_s + ")  [" + er_s + "]";}
      Spacer(&out, inj, margin);
    }
    Spacer(&out, PrecisionString(t[0], 3,false), 6); 
  
    (*Output).push_back(out); 
    out.Clear(); 
  }
  
  (*Output).push_back(""); 
  (*Output).push_back(""); 
  (*Output).push_back(""); 

  CompileHeading(Algos_V, " | ", margin, &out); 
  (*Output).push_back(out); 
  out.Clear();
  
  if (FL == "FL2"){Spacer(&out, "Flost 2 Scores", margin);}
  if (FL == "FL3"){Spacer(&out, "Flost 3 Scores", margin);}
  for (TString i : Algos_V)
  { 
    TString x; x += FLost_pred[i]; 
    Spacer(&out, x, margin);}
  (*Output).push_back(out);
}; 



#endif
