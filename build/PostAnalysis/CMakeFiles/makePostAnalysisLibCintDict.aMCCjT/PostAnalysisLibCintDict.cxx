// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PostAnalysisLibCintDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "PostAnalysis/AlgorithmFunctions.h"
#include "PostAnalysis/BaseFunctions.h"
#include "PostAnalysis/Convolution.h"
#include "PostAnalysis/Debug.h"
#include "PostAnalysis/DistributionGenerator.h"
#include "PostAnalysis/Evaluate.h"
#include "PostAnalysis/ExperimentalBench.h"
#include "PostAnalysis/IO.h"
#include "PostAnalysis/Plotting.h"
#include "PostAnalysis/ReportFigures.h"
#include "PostAnalysis/RooFitBaseFunctions.h"
#include "PostAnalysis/SimpleExampleFits.h"
#include "PostAnalysis/Statistics.h"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_libPostAnalysisLib_Impl() {
    static const char* headers[] = {
"PostAnalysis/AlgorithmFunctions.h",
"PostAnalysis/BaseFunctions.h",
"PostAnalysis/Convolution.h",
"PostAnalysis/Debug.h",
"PostAnalysis/DistributionGenerator.h",
"PostAnalysis/Evaluate.h",
"PostAnalysis/ExperimentalBench.h",
"PostAnalysis/IO.h",
"PostAnalysis/Plotting.h",
"PostAnalysis/ReportFigures.h",
"PostAnalysis/RooFitBaseFunctions.h",
"PostAnalysis/SimpleExampleFits.h",
"PostAnalysis/Statistics.h",
0
    };
    static const char* includePaths[] = {
"/home/tnom6927/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/QualificationTask/PostAnalysis/PostAnalysis",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/include",
"/home/tnom6927/QualificationTask/build/PostAnalysis/CMakeFiles/makePostAnalysisLibCintDict.aMCCjT/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libPostAnalysisLib dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libPostAnalysisLib dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef HAVE_PRETTY_FUNCTION
  #define HAVE_PRETTY_FUNCTION 1
#endif
#ifndef HAVE_64_BITS
  #define HAVE_64_BITS 1
#endif
#ifndef __IDENTIFIER_64BIT__
  #define __IDENTIFIER_64BIT__ 1
#endif
#ifndef ATLAS
  #define ATLAS 1
#endif
#ifndef ROOTCORE
  #define ROOTCORE 1
#endif
#ifndef XAOD_STANDALONE
  #define XAOD_STANDALONE 1
#endif
#ifndef XAOD_ANALYSIS
  #define XAOD_ANALYSIS 1
#endif
#ifndef ROOTCORE_RELEASE_SERIES
  #define ROOTCORE_RELEASE_SERIES 25
#endif
#ifndef PACKAGE_VERSION
  #define PACKAGE_VERSION "PostAnalysis-00-00-00"
#endif
#ifndef PACKAGE_VERSION_UQ
  #define PACKAGE_VERSION_UQ PostAnalysis-00-00-00
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "PostAnalysis/AlgorithmFunctions.h"
#include "PostAnalysis/BaseFunctions.h"
#include "PostAnalysis/Convolution.h"
#include "PostAnalysis/Debug.h"
#include "PostAnalysis/DistributionGenerator.h"
#include "PostAnalysis/Evaluate.h"
#include "PostAnalysis/ExperimentalBench.h"
#include "PostAnalysis/IO.h"
#include "PostAnalysis/Plotting.h"
#include "PostAnalysis/ReportFigures.h"
#include "PostAnalysis/RooFitBaseFunctions.h"
#include "PostAnalysis/SimpleExampleFits.h"
#include "PostAnalysis/Statistics.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libPostAnalysisLib",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libPostAnalysisLib_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libPostAnalysisLib_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libPostAnalysisLib() {
  TriggerDictionaryInitialization_libPostAnalysisLib_Impl();
}
