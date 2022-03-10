#include<TFile.h>
#include<TKey.h>
#include<TString.h>
#include<iostream>
#include<fstream>

std::vector<TString> ReturnCurrentDirs(bool FolderOnly = true)
{
  std::vector<TString> Output;
  TDirectory* dir = gDirectory;
  for (TObject* key : *dir -> GetListOfKeys())
  {
    auto k = dynamic_cast<TKey*>(key);
    TString dir = (TString)k->GetName(); 
    Output.push_back(dir);
  }
  return Output;
}

int main()
{
  TFile* file = new TFile("Fit_Tracks.root", "READ"); 
  std::vector<TString> dir = ReturnCurrentDirs(); 
 
  bool Rerun = false;  
  if (dir.size() == 0){ Rerun = true; }
  for (TString L : dir)
  {
    gDirectory -> cd(L); 
    std::vector<TString> dir2 = ReturnCurrentDirs(); 
    if (dir2.size() == 0){ Rerun = true;} 
  }
  if (Rerun == true)
  {
    std::ofstream print;
    print.open("empty.txt"); 
    print.close(); 
  }

  return 0; 
}
