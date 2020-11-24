#ifndef IO_H
#define IO_H
#include<TFile.h>
#include<TKey.h>
#include<iostream>
#include<TString.h>
#include<TH1F.h>

std::map<TString, std::vector<TH1F*>> ReadEntries(TFile* F); 


#endif 
