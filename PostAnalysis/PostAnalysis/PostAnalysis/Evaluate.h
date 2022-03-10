#ifndef EVALUATE_H
#define EVALUATE_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>


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
#endif
