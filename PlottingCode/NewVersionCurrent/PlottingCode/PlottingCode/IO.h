#ifndef IO_H
#define IO_H

#include "../PlottingCode/SharedFunctions.h"

std::vector<_TS> ReturnCurrentDirs();
std::map<_TS, std::vector<TH1F*>> ReadCTIDE(_TS dir);
std::map<_TS, std::map<_TS, std::map<_TS, std::map<_TS, TH1F*>>>> ReadDebugging(_TS dir);
void ReadPostAnalysis(_TS dir);




#endif
