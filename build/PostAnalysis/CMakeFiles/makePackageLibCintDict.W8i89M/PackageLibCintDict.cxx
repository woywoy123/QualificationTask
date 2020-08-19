// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PackageLibCintDict

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
#include "PostAnalysis/Functions.h"
#include "PostAnalysis/UnitClosures.h"
#include "PostAnalysis/Verification.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *Functions_Dictionary();
   static void Functions_TClassManip(TClass*);
   static void *new_Functions(void *p = 0);
   static void *newArray_Functions(Long_t size, void *p);
   static void delete_Functions(void *p);
   static void deleteArray_Functions(void *p);
   static void destruct_Functions(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Functions*)
   {
      ::Functions *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Functions));
      static ::ROOT::TGenericClassInfo 
         instance("Functions", "PostAnalysis/Functions.h", 36,
                  typeid(::Functions), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Functions_Dictionary, isa_proxy, 4,
                  sizeof(::Functions) );
      instance.SetNew(&new_Functions);
      instance.SetNewArray(&newArray_Functions);
      instance.SetDelete(&delete_Functions);
      instance.SetDeleteArray(&deleteArray_Functions);
      instance.SetDestructor(&destruct_Functions);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Functions*)
   {
      return GenerateInitInstanceLocal((::Functions*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Functions*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Functions_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Functions*)0x0)->GetClass();
      Functions_TClassManip(theClass);
   return theClass;
   }

   static void Functions_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_Functions(void *p) {
      return  p ? new(p) ::Functions : new ::Functions;
   }
   static void *newArray_Functions(Long_t nElements, void *p) {
      return p ? new(p) ::Functions[nElements] : new ::Functions[nElements];
   }
   // Wrapper around operator delete
   static void delete_Functions(void *p) {
      delete ((::Functions*)p);
   }
   static void deleteArray_Functions(void *p) {
      delete [] ((::Functions*)p);
   }
   static void destruct_Functions(void *p) {
      typedef ::Functions current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Functions

namespace {
  void TriggerDictionaryInitialization_libPackageLib_Impl() {
    static const char* headers[] = {
"PostAnalysis/Functions.h",
"PostAnalysis/UnitClosures.h",
"PostAnalysis/Verification.h",
0
    };
    static const char* includePaths[] = {
"/home/tnom6927/CTIDE/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/CTIDE/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/CTIDE/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/CTIDE/QualificationTask/PostAnalysis/PostAnalysis",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/include",
"/home/tnom6927/CTIDE/QualificationTask/build/PostAnalysis/CMakeFiles/makePackageLibCintDict.W8i89M/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libPackageLib dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$PostAnalysis/Functions.h")))  Functions;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libPackageLib dictionary payload"

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
#include "PostAnalysis/Functions.h"
#include "PostAnalysis/UnitClosures.h"
#include "PostAnalysis/Verification.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"Functions", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libPackageLib",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libPackageLib_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libPackageLib_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libPackageLib() {
  TriggerDictionaryInitialization_libPackageLib_Impl();
}
